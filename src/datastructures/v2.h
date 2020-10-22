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

struct v2i {
    int x, y;

    v2i(int p_x = 0, int p_y = 0)
        : x(p_x), y(p_y)
    {
    }

    v2i& operator=(const v2i& a)
    {
        x = a.x;
        y = a.y;
        return *this;
    }

    friend bool operator==(const v2i& a, const v2i& b)
    {
        return (a.x == b.x && a.y == b.y);
    }

    friend bool operator!=(const v2i& a, const v2i& b)
    {
        return !(a == b);
    }


    v2i operator+(const v2i& a) const
    {
        return v2i(a.x + x, a.y + y);
    }

    v2i operator-(const v2i& a) const
    {
        return v2i(x - a.x, y - a.y);
    }

    v2i operator/(int a) const
    {
        return v2i(x / a, y / a);
    }

    v2i operator*(int a) const
    {
        return v2i(x * a, y * a);
    }
};


#define double2 v2
#define double4 v4
#define int2 v2i

inline double dot(double2 a, double2 b) { return a.x * b.x + b.y * b.y; };
