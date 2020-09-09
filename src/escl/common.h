#ifndef CL_COMMON_H
#define CL_COMMON_H

#ifdef DEVICE_PROGRAM
    #include "src/escl/util.h"

    #define BUFFER_ARGS __global SimulationParameters* sp, __global double2* impurities
    #define MIN(p_a, p_b) min((p_a), (p_b))
#else
    #include "math.h"
    #include "escl/constants.h"
    #define BUFFER_ARGS SimulationParameters *sp, std::vector<v2> &impurities
    
    #define MIN(p_a, p_b) ((p_a) < (p_b)) ? (p_a) : (p_b)

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
    
inline bool IsEdge(int i, int dim) {
    return (i == 0) || (i == (dim - 2)); 
}

inline bool IsPadding(int i, int dim) {
    return (i == (dim-1)); 
}

inline double GetWeight1D(int i, int dim) {
    double w = IsPadding(i, dim) ? 0.0 : 1.0;
    if (!IsEdge(i, dim))
    {
        w *= ((i % 2) == 0) ? 2.0 : 4.0;
    }

    return w;
}

inline double GetWeight2D(int i, int j, int dim) {
    double w = (IsPadding(i, dim) || IsPadding(j, dim)) ? 0.0 : 1.0;
    if (!IsEdge(i, dim))
    {
        w *= ((i % 2) == 0) ? 2.0 : 4.0;
    }
    if (!IsEdge(j, dim))
    {
        w *= ((j % 2) == 0) ? 2.0 : 4.0;
    }

    return w;
}

inline bool IsSigma(int m) {
    return (m == MODE_SIGMA_XX || m == MODE_SIGMA_XY);
}

typedef struct
{
    int mode;
    int dim; //-
    int impurity_count;
    int integrand_steps;
    int clockwise;
    unsigned int impurity_seed;

    double region_size;
    double region_extends;
    double particle_speed;      // v
    double impurity_radius;     // r
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
    double2 shift = { velocity.y, -velocity.x };
    shift = shift * radius / vf;

    //return is_electron ? (p - shift) :  (p + shift);
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
    double xs = (dist_squared + r1*r1 - r2*r2) / (2.0 * dist);
    double ys = sqrt(r1*r1 - xs*xs);

    double2 u = (p2 - p1) / dist;

    double4 points = {
        p1.x + u.x * xs +  u.y  * ys,
        p1.y + u.y * xs + -u.x  * ys,

        p1.x + u.x * xs + -u.y  * ys,
        p1.y + u.y * xs +  u.x  * ys
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
    return MIN(t1, t2);
}

inline double lifetimeB(double max_lifetime, double2 pos, double phi, bool clockwise, BUFFER_ARGS)
{
    double orbit_radius = sp->particle_speed / sp->angular_speed;
    double2 vel = { cos(phi), sin(phi) };
    vel = vel * sp->particle_speed;
    double2 center = GetCyclotronOrbit(pos, vel, orbit_radius, sp->particle_speed, clockwise);

    double lifetime = max_lifetime;

    double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    for (int i = 0; i < sp->impurity_count; i++) {
        double2 impurity = impurities[i];

        double2 d = pos - impurity;
        if (impurity_radius_sq > dot(d, d))
        {
            lifetime = 0;
            break;
        }

        if (CirclesCross(center, orbit_radius, impurity, sp->impurity_radius))
        {
            double t = GetFirstCrossTime(center, pos, impurity, orbit_radius, sp->impurity_radius, sp->angular_speed, clockwise);

            if (t < lifetime)
                lifetime = t;
        }
    }

    return lifetime;
}

inline double lifetime0(double2 pos, double phi, BUFFER_ARGS)
{
    double2 unit = { cos(phi), sin(phi) };
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
        if (fabs(vel.x) > (fabs(vel.y) * 1e-9)) {
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

////////////////////////////////////////////////////
// Main functions, for single or integrated lifetime
////////////////////////////////////////////////////

inline double single_lifetime(double2 pos, double phi, BUFFER_ARGS) {
    if (sp->angular_speed != 0) {
        bool clockwise = (sp->clockwise == 1);
        double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, clockwise, false);
 
        return lifetimeB(MIN(sp->tau, bound_time), pos, phi, clockwise, sp, impurities);
    }
    else {
        return lifetime0(pos, phi, sp, impurities);
    }
}

inline double phi_lifetime(double2 pos, BUFFER_ARGS)
{
    bool clockwise = (sp->clockwise == 1);
    //if (clockwise) sp->angular_speed *= -1;

    double angle_area = sp->alpha * 2.0;
    double step_size = angle_area / (sp->integrand_steps - 1);
    
    double integral = 0;
    for (int j = 0; j < 4; j++)
    {
        double start = -sp->alpha + j * (PI * 0.5);
        double total = 0.0;

        for (int i = 0; i < sp->integrand_steps; i++)
        {
            double phi = start + i * step_size;

            double result = single_lifetime(pos, phi, sp, impurities);

            if (sp->mode == MODE_SIGMA_XX || sp->mode == MODE_SIGMA_XY)
            {
                double z = exp(-result / sp->tau);

                double r = cos(phi) - cos(phi + sp->angular_speed * result) * z;
                
                //todo abs(angular_speed)
                r += sp->angular_speed * sp->tau * sin(phi + sp->angular_speed * result) * z;
                r -= sp->angular_speed * sp->tau * sin(phi);

                double v;
                if (sp->mode == MODE_SIGMA_XX) v = cos(phi);
                else                           v = sin(phi);

                result = r * v;
            }
      
            double w = 1.0;
            if (!(i == 0 || i == (sp->integrand_steps-1)))
            {
                w = ((i % 2) == 0) ? 2.0 : 4.0;
            }

            total += result * w;
        }

        integral += total * angle_area / ((sp->integrand_steps - 1) * 3.0);
    }

    return integral;
}

inline float GetColor(double v, double scale, int mode) {
    if      (mode == MODE_DIR_LIFETIME) v /= scale;
    else if (mode == MODE_PHI_LIFETIME) v /= (scale * 6.0);
    else if (mode == MODE_SIGMA_XX)     v /= (scale * 3.0);
    else if (mode == MODE_SIGMA_XY)     v /= scale;

    return (float)v;
}

#ifdef DEVICE_PROGRAM
__kernel void add_integral_weights_2d(__global double* A)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    int i = y * row_size + x;
    A[i] *= GetWeight2D(x, y, row_size);
}

__kernel void to_texture(__global double* lifetimes, int mode, double scale, __write_only image2d_t screen)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    float k = GetColor((float)(lifetimes[y * row_size + x]), scale, mode);
    float4 c;
    if (mode != MODE_SIGMA_XY) {
        c = (float4)(k, k, k, 1.0f);
    }
    else {
        if (k < 0.0) c = (float4)(0, 0, -k, 1.0f);
        else 		 c = (float4)(k, 0, 0, 1.0f);
    }

    write_imagef(screen, (int2)(x, y), c);
}

#endif

#endif // CL_COMMON_H
