__constant double PI = 3.141592653589793238463;
__constant double PI2 = 6.283185307179586;

#define ROW_SIZE 1000
#define GLINTEROP


double smod(double a, double b)
{
    return a - b * floor(a / b);
}

inline double GetBoundTime(const double phi, const double alpha, const double w, const bool is_electron, const bool is_future)
{
    double remaining = smod(phi + alpha, PI * 0.5);

    double dphi;
    if (!is_electron && is_future) dphi = remaining;
    else if (!is_electron)         dphi = 2.0 * alpha - remaining;
    else if (is_future)            dphi = 2.0 * alpha - remaining;
    else                           dphi = remaining;

    return dphi / w;
}

double2 GetCyclotronOrbit(double2 p, double2 velocity, double radius, double vf, bool is_electron)
{
    double2 shift = { radius * velocity.y / vf, -radius * velocity.x / vf };

    double2 center;
    if (is_electron)
    {
        center.x = p.x - shift.x;
        center.y = p.y - shift.y;
        }
    else { 
        center.x = p.x + shift.x;
        center.y = p.y + shift.y;
    }

    return center;
}

bool CirclesCross(double2 p1, double r1, double2 p2, double r2)
{
    double dist_squared = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
    if (dist_squared >= pow(r1 + r2, 2)) return false;
    if (dist_squared <= pow(r1 - r2, 2)) return false;

    return true;
}

double4 GetCrossPoints(double2 p1, double r1, double2 p2, double r2)
{
    double dist_squared = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
    double dist = sqrt(dist_squared);
    double xs = (dist_squared + r1 * r1 - r2 * r2) / (2.0 * dist);
    double ys = sqrt(r1 * r1 - xs * xs);

    double2 u = { (p2.x - p1.x) / dist, (p2.y - p1.y) / dist };

    double4 points = {
        p1.x + u.x * xs + u.y * ys,
        p1.y + u.y * xs + -u.x * ys,

        p1.x + u.x * xs + -u.y * ys,
        p1.y + u.y * xs + u.x * ys
    };

    return points;
}

double GetPhi(double2 pos, double2 center, double radius)
{
    double p = (pos.x - center.x) / radius;
    if (p > 1.0) p = 1.0;
    if (p < -1.0) p = -1.0;
    double phi = acos(p);

    if (pos.y < center.y)
        phi = PI2 - phi;

    return phi;
}

double GetCrossAngle(double p, double q, bool clockwise)
{
    double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

double GetFirstCrossTime(double2 center, double2 pos, double2 ip, double r, double ir, double w, double clockwise)
{
    double4 cross_points = GetCrossPoints(center, r, ip, ir);

    double2 p1 = {cross_points.x, cross_points.y};
    double2 p2 = {cross_points.z, cross_points.w};

    double phi0 = GetPhi(pos, center, r);
    double phi1 = GetPhi(p1, center, r);
    double phi2 = GetPhi(p2, center, r);

    double t1 = GetCrossAngle(phi0, phi1, clockwise) / w;
    double t2 = GetCrossAngle(phi0, phi2, clockwise) / w;
    return min(t1, t2);
}

double lifetime0(double max_lifetime, double2 pos, double2 vel, double angular_speed, bool clockwise, int impurity_count, double imp_radius,  __global double2 *imps)
{
    double lifetime = max_lifetime;

    double vf = sqrt(vel.x * vel.x + vel.y * vel.y);
    double radius = vf / angular_speed;
    double2 center = GetCyclotronOrbit(pos, vel, radius, vf, clockwise);

    double impurity_radius_sq = imp_radius * imp_radius;

    for (int i = 0; i < impurity_count; i++) {
        double2 imp_pos = imps[i];

        double2 d = { pos.x - imp_pos.x, pos.y - imp_pos.y };
        if (impurity_radius_sq > d.x*d.x + d.y*d.y)
        {
            lifetime = 0;
            break;
        }

        if (CirclesCross(center, radius, imp_pos, imp_radius))
        {
            double t = GetFirstCrossTime(center, pos, imp_pos, radius, imp_radius, angular_speed, clockwise);
            if (t < lifetime)
                lifetime = t;
		}
    }

    return lifetime;
}

__kernel void lifetime(double region_size,
                       double speed,
                       double imp_radius,
                       double tau,
                       double alpha,
                       double phi,
                       double angular_speed,
                       int impurity_count,
                       __global double2 *imps,
#ifdef GLINTEROP
                       __write_only image2d_t screen)
#else
                       __global double *lifetimes) 
#endif
{
    bool clockwise = false;
    int x = get_global_id(0);
    int y = get_global_id(1);
    
    double2 pos = {region_size * x / (ROW_SIZE-1), region_size * y / (ROW_SIZE-1)};

    double2 unit = { cos(phi), sin(phi) };
    double2 vel = { speed * unit.x, speed * unit.y };

    double bound_time = GetBoundTime(phi, alpha, angular_speed, clockwise, false);
    double max_lifetime =  min(tau, bound_time);
    double lifetime = lifetime0(max_lifetime, pos, vel, angular_speed, clockwise, impurity_count, imp_radius, imps);

#ifdef GLINTEROP
    float k = (float)(lifetime / max_lifetime);
    write_imagef(screen, (int2)(x, y), (float4)(k,k,k,1.0f));
#else
    lifetimes[y * ROW_SIZE + x] = lifetime;
#endif
}