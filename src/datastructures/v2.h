#pragma once

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
