/* -------------------------------------------------------------------------
    This code is part of ElasticScattering.

    Copyright(C) 2022 Elastic Scattering developers

    This program is free software : you can redistribute it and /or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

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

    bool operator==(const v2& a) const
    {
        return a.x == x && a.y == y;
    }

    bool operator<(const v2& a)
    {
        return (x < a.x) || ((!(a.x < x)) && (y < a.y));
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

    friend bool operator<(const v2i& a, const v2i& b)
    {
        return (a.x < b.x) || ((!(b.x < a.x)) && (a.y < b.y));
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
