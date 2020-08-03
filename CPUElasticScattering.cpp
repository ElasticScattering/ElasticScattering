#include <windows.h>

#include <random>
#include <limits>
#include <vector>
#include <assert.h>

#include "ElasticScattering.h"

void CPUElasticScattering::Init(SimulationParameters p_sp)
{
    sp = p_sp;
    if (sp.angular_speed == 0) sp.particle_max_lifetime = sp.tau;
    else  {
        double bound_time = GetBoundTime(true, false);
        sp.particle_max_lifetime = MIN(sp.tau, bound_time);
    }
    std::cout << "Particle max lifetime: " << sp.particle_max_lifetime << std::endl;

    // Initialize arrays.
    impurities.clear();
    impurities.resize(sp.impurity_count); 
    
    std::cout << "Impurity region: " << -sp.particle_speed * sp.tau << ", " << sp.region_size + sp.particle_speed * sp.tau << std::endl;
    std::uniform_real_distribution<double> unif(-sp.particle_speed * sp.tau, sp.region_size + sp.particle_speed * sp.tau);
    std::random_device r;
    std::default_random_engine re(0);
    
    for (int i = 0; i < sp.impurity_count; i++)
        impurities[i] = { unif(re), unif(re) };

    lifetimes.clear();
    lifetimes.resize(sp.particle_count, 0);
}

void CPUElasticScattering::Compute()
{
    double total_time;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    std::cout << "Simulating elastic scattering on the CPU..." << std::endl;

    QueryPerformanceCounter(&beginClock);
    for (int j = 0; j < sp.particle_row_count; j++)
        for (int i = 0; i < sp.particle_row_count; i++)
        {
            cl_double2 pos;
            pos.x = sp.region_size * (double(i) / sp.particle_row_count);
            pos.y = sp.region_size * (double(j) / sp.particle_row_count);

            double lifetime = sp.particle_max_lifetime;

            cl_double2 vel = { sp.particle_speed * cos(sp.phi), sp.particle_speed * sin(sp.phi) };

            double res = (sp.angular_speed == 0) ? ComputeA(pos, vel, sp) : ComputeB(pos, vel, sp);
            lifetimes[j * sp.particle_row_count + i] = res;
        }
    QueryPerformanceCounter(&endClock);

    double total = 0;
    for (int i = 0; i < sp.particle_count; i++) {
        total += lifetimes[i];
    }
    std::cout << "Average lifetime:     " << total / sp.particle_count << " s" << std::endl;

    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "CPU calculation time: " << total_time * 1000 << " ms" << std::endl;

    std::cout << "\n\nSorted results:" << std::endl;
    for (int i = 0; i < MIN(lifetimes.size() - 1, 200); i++)
        std::cout << lifetimes[i] << ", ";
    std::cout << "..." << std::endl;

    MakeTexture(sp);
}

double CPUElasticScattering::ComputeA(const cl_double2 pos, const cl_double2 vel, const SimulationParameters sp)
{
    const cl_double2 unit = { cos(sp.phi), sin(sp.phi) };
    double lifetime = sp.particle_max_lifetime;

    for (int k = 0; k < sp.impurity_count; k++)
    {
        const cl_double2 ip = impurities[k];
        const cl_double inner = (ip.x - pos.x) * unit.x + (ip.y - pos.y) * unit.y;
        const cl_double2 projected = { pos.x + inner * unit.x, pos.y + inner * unit.y };

        const double a = pow(projected.x - ip.x, 2) + pow(projected.y - ip.y, 2);
        if (a > sp.impurity_radius_sq) {
            continue; //@Speedup, if distance is greater than current min continue as well?
        }

        double L = sqrt(sp.impurity_radius_sq - a);

        cl_double2 time_taken;
        if (vel.x != 0) time_taken = { -((projected.x - L * unit.x) - pos.x) / vel.x, -((projected.x + L * unit.x) - pos.x) / vel.x };
        else            time_taken = { -((projected.y - L * unit.y) - pos.y) / vel.y, -((projected.y + L * unit.y) - pos.y) / vel.y };

        if ((time_taken.s0 * time_taken.s1) < 0)
            return 0;

        if (time_taken.s0 > 0 && time_taken.s0 < lifetime) {
            lifetime = time_taken.s0;
        }
        if (time_taken.s1 > 0 && time_taken.s1 < lifetime) {
            lifetime = time_taken.s1;
        }
    }

    return lifetime;
}

double CPUElasticScattering::ComputeB(const cl_double2 pos, const cl_double2 vel, const SimulationParameters sp)
{
    const bool clockwise = true;

    double vf = sqrt(vel.x * vel.x + vel.y * vel.y);
    double radius = vf / sp.angular_speed;
    auto center = GetCyclotronOrbit(pos, vel, radius, vf, clockwise);

    double lifetime = sp.particle_max_lifetime;
    for (int k = 0; k < sp.impurity_count; k++)
    {
        const cl_double2 ip = impurities[k];

        if (sp.impurity_radius_sq > pow(pos.x - ip.x, 2) + pow(pos.y - ip.y, 2))
        {
            lifetime = 0;
            break;
        }

        if (CirclesCross(center, radius, ip, sp.impurity_radius))
        {
            double t = GetCrossTime(center, pos, ip, radius, sp.impurity_radius, sp.angular_speed, clockwise);

            assert(t >= 0);

            if (t < lifetime)
                lifetime = t;
        }
    }

    return lifetime;
}

inline double smod(double a, double b)
{
    return a - b * floor(a / b);
}

double CPUElasticScattering::GetBoundTime(const bool is_electron, const bool is_future) const
{
    double remaining = smod(sp.phi + sp.alpha, PI * 0.5);

    double dphi;
    if (!is_electron && is_future) dphi = remaining;
    else if (!is_electron)         dphi = 2 * sp.alpha - remaining;
    else if (is_future)            dphi = 2 * sp.alpha - remaining;
    else                           dphi = remaining;

    return dphi / sp.angular_speed;
}

cl_double2 CPUElasticScattering::GetCyclotronOrbit(const cl_double2 p, const cl_double2 velocity, const double radius, const double vf, const bool is_electron) const
{
    cl_double2 shift = { radius * velocity.y / vf, -radius * velocity.x / vf };

    cl_double2 center;
    if (is_electron) center = { p.x - shift.x, p.y - shift.y };
    else             center = { p.x + shift.x, p.y + shift.y };

    return center;
}

bool CPUElasticScattering::CirclesCross(const cl_double2 p1, const double r1, const cl_double2 p2, const double r2) const
{
    double dist_squared = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
    if (dist_squared >= pow(r1 + r2, 2)) return false;
    if (dist_squared <= pow(r1 - r2, 2)) return false;

    return true;
}

std::pair<cl_double2, cl_double2> CPUElasticScattering::GetCrossPoints(const cl_double2 p1, const double r1, const cl_double2 p2, const double r2) const
{
    const double dist_squared = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
    const double dist = sqrt(dist_squared);
    const double xs = (dist_squared + r1 * r1 - r2 * r2) / (2.0 * dist);
    const double ys = sqrt(r1 * r1 - xs * xs);

    cl_double2 u = { (p2.x - p1.x) / dist, (p2.y - p1.y) / dist };

    std::pair<cl_double2, cl_double2> points;
    points.first = {
        p1.x + u.x * xs + u.y * ys,
        p1.y + u.y * xs + -u.x * ys
    };

    points.second = {
        p1.x + u.x * xs + -u.y * ys,
        p1.y + u.y * xs + u.x * ys
    };

    return points;
}

double CPUElasticScattering::GetPhi(const cl_double2 pos, const cl_double2 center, const double radius) const
{
    double p = (pos.x - center.x) / radius;
    assert(abs(p) < 1.0001);
    p = max(min(p, 1), -1);
    double phi = acos(p);

    if (pos.y < center.y)
        phi = PI2 - phi;

    return phi;
}

double CPUElasticScattering::GetCrossAngle(const double p, const double q, const bool clockwise) const
{
    double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

double CPUElasticScattering::GetCrossTime(const cl_double2 center, const cl_double2 pos, const cl_double2 ip, const double r, const double ir, const double w, const double clockwise) const
{
    const auto cross_points = GetCrossPoints(center, r, ip, ir);

    const double phi0 = GetPhi(pos, center, r);
    const double phi1 = GetPhi(cross_points.first, center, r);
    const double phi2 = GetPhi(cross_points.second, center, r);

    const double t1 = GetCrossAngle(phi0, phi1, clockwise) / w;
    const double t2 = GetCrossAngle(phi0, phi2, clockwise) / w;
    return min(t1, t2);
}

void CPUElasticScattering::MakeTexture(const SimulationParameters sp)
{
    double itau = 1.0 / sp.particle_max_lifetime;
    pixels.clear();
    pixels.resize(sp.particle_count * 3L);
    size_t j = 0;
    for (int i = 0; i < sp.particle_count; i++)
    {
        float k = float(lifetimes[i] * itau);
        if (k == 0) {
            pixels[j] = 1.0f;
            pixels[j + 1L] = 0.0f;
            pixels[j + 2L] = 0.0f;
        }
        else {
            pixels[j] = k;
            pixels[j + 1L] = k;
            pixels[j + 2L] = k;
        }
        j += 3;
    }
}
