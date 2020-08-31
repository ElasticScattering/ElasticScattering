#ifndef CL_COMMON_H
#define CL_COMMON_H

#ifndef DEVICE_PROGRAM
    #include "math.h"
    #include "escl/constants.h"
#else
    #include "src/escl/constants.h"
#endif



//Remove 1 from row_size to have an inclusive range, another because the kernel work dimension is even, but the integral requires uneven dimensions.
#define DECLARE_POS double2 pos = (double2)(x, y) * (sp->region_size / (row_size - 2));

#ifdef DEVICE_PROGRAM
    #define BUFFER_ARGS __global SimulationParameters* sp, __global double2* impurities
#else
    #define BUFFER_ARGS SimulationParameters* sp, double2* impurities

    struct v4 {
        double x, y, z, w;
    };

    struct v2 {
        double x, y;

        v2(double p_x = 0, double p_y = 0)
            : x(p_x), y(p_y)
        {
        }

        v2& operator=(const v2& a)
        {
            x = a.x;
            y = a.y;
            return *this;
        }

        v2 operator+(const v2& a) const
        {
            return v2(a.x + x, a.y + y);
        }
        
        v2 operator-(const v2& a) const
        {
            return v2(x - a.x, y - a.y);
        }

        v2 operator/(double a) const
        {
            return v2(x / a, y / a);
        }

        v2 operator*(double a) const
        {
            return v2(x * a, y * a);
        }
    };

    #define double2 v2
    #define double4 v4

    inline double dot(double2 a, double2 b) { return a.x * b.x + b.y * b.y; };
#endif
    //inline bool IsEdge(int i, int j, int dim) { return (i == 0) || (i == dim) || (j == 0) || (j == dim); }

    inline bool IsEdge(int i, int j, int dim) { return (i == 0) || (i == (dim - 2)) || (j == 0) || (j == (dim - 2)); }
    inline bool IsPadding(int i, int j, int dim) { return (i == (dim-1)) || (j == (dim-1)); }
    inline double GetWeight(int i, int j, int dim) {
        double w = IsPadding(i, j, dim) ? 0.0 : 1.0;
        if (!IsEdge(i, j, dim))
        {
            w  = ((i % 2) == 0) ? 2.0 : 4.0;
            w *= ((j % 2) == 0) ? 2.0 : 4.0;
        }

        return w;
    }
    inline bool IsSigma(int m) { return (m == MODE_SIGMA_XX || m == MODE_SIGMA_XY); }

typedef struct
{
    int mode;
    int dim; //-                           
    int particle_count; //-
    int impurity_count;
    int integrand_steps; //~
    int clockwise;
    unsigned int impurity_seed;

    double region_size;
    double region_extends;
    double particle_speed;      // v
    double particle_mass;       // m,-
    double impurity_radius;     // r
    double impurity_radius_sq;  // r^2,~
    double tau;

    double alpha;
    double phi;
    double magnetic_field;      // B
    double angular_speed;       // w
} SimulationParameters;


inline double smod(double a, double b)
{
    return a - b * floor(a / b);
}

inline double GetBoundTime(double phi, double alpha, double w, bool is_electron, bool is_future)
{
    double remaining = smod(phi + alpha, PI * 0.5);

    double dphi;
    if (!is_electron && is_future) dphi = remaining;
    else if (!is_electron)         dphi = 2.0 * alpha - remaining;
    else if (is_future)            dphi = 2.0 * alpha - remaining;
    else                           dphi = remaining;

    return dphi / w;
}

inline double2 GetCyclotronOrbit(double2 p, double2 velocity, double radius, double vf, bool is_electron)
{
    double2 shift = (double2)(velocity.y, -velocity.x) * radius / vf;

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

inline bool CirclesCross(double2 p1, double r1, double2 p2, double r2)
{
    double2 q = p1 - p2;

    double dist_squared = q.x * q.x + q.y * q.y;
    if (dist_squared >= (r1 + r2) * (r1 + r2)) return false;
    if (dist_squared <= (r1 - r2) * (r1 - r2)) return false;

    return true;
}

inline double4 GetCrossPoints(double2 p1, double r1, double2 p2, double r2)
{
    double2 q = p1 - p2;

    double dist_squared = dot(q, q);
    double dist = sqrt(dist_squared);
    double xs = (dist_squared + r1 * r1 - r2 * r2) / (2.0 * dist);
    double ys = sqrt(r1 * r1 - xs * xs);

    double2 u = (p2 - p1) / dist;

    double4 points = {
        p1.x + u.x * xs + u.y * ys,
        p1.y + u.y * xs + -u.x * ys,

        p1.x + u.x * xs + -u.y * ys,
        p1.y + u.y * xs + u.x * ys
    };

    return points;
}

inline double GetPhi(double2 pos, double2 center, double radius)
{
    double p = (pos.x - center.x) / radius;
    if (p > 1.0) p = 1.0;
    if (p < -1.0) p = -1.0;
    double phi = acos(p);

    if (pos.y < center.y)
        phi = PI2 - phi;

    return phi;
}

inline double GetCrossAngle(double p, double q, bool clockwise)
{
    double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

inline double GetFirstCrossTime(double2 center, double2 pos, double2 ip, double r, double ir, double w, double clockwise)
{
    double4 cross_points = GetCrossPoints(center, r, ip, ir);

    double2 p1 = { cross_points.x, cross_points.y };
    double2 p2 = { cross_points.z, cross_points.w };

    double phi0 = GetPhi(pos, center, r);
    double phi1 = GetPhi(p1, center, r);
    double phi2 = GetPhi(p2, center, r);

    double t1 = GetCrossAngle(phi0, phi1, clockwise) / w;
    double t2 = GetCrossAngle(phi0, phi2, clockwise) / w;
    return min(t1, t2);
}


inline double lifetime0(double2 pos, BUFFER_ARGS)
{
    double2 unit = { cos(sp->phi), sin(sp->phi) };
    double2 vel = unit * sp->particle_speed;

    double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double lifetime = sp->tau;

    for (int i = 0; i < sp->impurity_count; i++) {
        double2 imp_pos = impurities[i];
        double inner = (imp_pos.x - pos.x) * unit.x + (imp_pos.y - pos.y) * unit.y;
        double2 projected = pos + unit * inner;

        double2 d = projected - imp_pos;
        double diff = impurity_radius_sq - dot(d, d);
        if (diff < 0.0) {
            continue;
        }

        double L = sqrt(diff);

        double2 time_taken;
        if (vel.x != 0) {
            time_taken.x = -((projected.x - L * unit.x) - pos.x) / vel.x;
            time_taken.y = -((projected.x + L * unit.x) - pos.x) / vel.x;
        }
        else {
            time_taken.x = -((projected.y - L * unit.y) - pos.y) / vel.y;
            time_taken.y = -((projected.y + L * unit.y) - pos.y) / vel.y;
        }

        if ((time_taken.x * time_taken.y) < 0) {
            lifetime = 0;
            break;
        }

        if (time_taken.x > 0 && time_taken.x < lifetime) {
            lifetime = time_taken.x;
        }
        if (time_taken.y > 0 && time_taken.y < lifetime) {
            lifetime = time_taken.y;
        }
    }

    return lifetime;
}

inline double lifetimeB(double max_lifetime, double2 pos, bool clockwise, BUFFER_ARGS)
{
    double orbit_radius = sp->particle_speed / sp->angular_speed;
    double2 vel = (cos(sp->phi), sin(sp->phi)) * sp->particle_speed;
    double2 center = GetCyclotronOrbit(pos, vel, orbit_radius, sp->particle_speed, clockwise);

    double lifetime = max_lifetime;

    for (int i = 0; i < sp->impurity_count; i++) {
        double2 imp_pos = impurities[i];

        double2 d = pos - imp_pos;

        if (CirclesCross(center, orbit_radius, imp_pos, sp->impurity_radius))
        {
            if (sp->impurity_radius_sq > dot(d, d))
            {
                lifetime = 0;
            }

            double t = GetFirstCrossTime(center, pos, imp_pos, orbit_radius, sp->impurity_radius, sp->angular_speed, clockwise);

            if (t < lifetime)
                lifetime = t;
        }
    }

    return lifetime;
}

inline double single_lifetime(double2 pos, BUFFER_ARGS) {
    if (sp->angular_speed != 0) {
        bool clockwise = (sp->clockwise == 1);
        double bound_time = GetBoundTime(sp->phi, sp->alpha, sp->angular_speed, clockwise, false);
 
        return lifetimeB(min(sp->tau, bound_time), pos, clockwise, sp, impurities);
    }
    else {
        return lifetime0(pos, sp, impurities);
    }
}

inline double phi_lifetime(double2 pos, BUFFER_ARGS)
{
    bool clockwise = (sp->clockwise == 1);
    if (clockwise) sp->angular_speed *= -1;

    double angle_area = sp->alpha * 2.0;
    double step_size = angle_area / (sp->integrand_steps - 1);
    
    double integral = 0;
    for (int j = 0; j < 4; j++)
    {
        double start = -sp->alpha + j * (PI * 0.5);
        double total = 0.0;

        bool is_even = true;
        for (int i = 0; i < sp->integrand_steps; i++)
        {
            sp->phi = start + i * step_size;

            double result;
            if (sp->angular_speed != 0) {
                double bound_time = GetBoundTime(sp->phi, sp->alpha, sp->angular_speed, clockwise, false);
                result = lifetimeB(min(sp->tau, bound_time), pos, clockwise, sp, impurities);
            }
            else {
                result = lifetime0(pos, sp, impurities);
            }

            if (sp->mode == MODE_SIGMA_XX || sp->mode == MODE_SIGMA_XY)
            {
                double z = exp(-result / sp->tau);

                double r = cos(sp->phi) - cos(sp->phi + sp->angular_speed * result) * z;
                r += sp->angular_speed * sp->tau * sin(sp->phi + sp->angular_speed * result) * z;
                r -= sp->angular_speed * sp->tau * sin(sp->phi);
                r *= sp->tau;

                double v;
                if (sp->mode == MODE_SIGMA_XX) v = cos(sp->phi);
                else                           v = sin(sp->phi);

                result = r * v;
            }

            bool edge_item = i == 0 || i == (sp->integrand_steps - 1);
            is_even = (i % 2) == 0;
            double w = 1.0;

            if (!edge_item) {
                w = is_even ? 2.0 : 4.0;
            }

            total += result * w;
        }

        integral += total * angle_area / ((sp->integrand_steps - 1) * 3.0);
    }

    return integral;
}

#endif // CL_COMMON_H
