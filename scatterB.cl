__constant double PI = 3.141592653589793238463;
__constant double PI2 = 6.283185307179586;

#define ROW_SIZE 1000

double smod(double a, double b)
{
    return a - b * floor(a / b);
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

__kernel void scatterB(double region_size,
                       double max_lifetime,
                       double speed,
                       double mass,
                       double imp_radius,
                       double tau,
                       double alpha,
                       double phi,
                       double magnetic_field,
                       double angular_speed,
                       int impurity_count,
                       __global double2 *imps,
                       __global double *lifetimes) 
{
    bool clockwise = true;
    int idx = get_global_id(0);
    
    int x = idx % ROW_SIZE;
    int y = idx / ROW_SIZE;
    double2 pos = {region_size * x / ROW_SIZE, region_size * y / ROW_SIZE};

    double2 unit = { cos(phi), sin(phi) };
    double2 vel = { speed * unit.x, speed * unit.y };

    double impurity_radius_sq = imp_radius * imp_radius;

    double vf = sqrt(vel.x * vel.x + vel.y * vel.y);
    double radius = vf / angular_speed;
    double2 center = GetCyclotronOrbit(pos, vel, radius, vf, clockwise);

    double lifetime = max_lifetime;

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

    lifetimes[idx] = lifetime;
}