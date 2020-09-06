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

__kernel void many_sqrt(__global double *a, int n, __global double *b) {
	int id = get_global_id(0);
	for (int i = 0; i < n; i++) {
		double w = fabs(sin((double)(i)));
	    b[id] += a[id] * sqrt(w);
	}
}

__kernel void test_double2(__global double *a) {
    int id = get_global_id(0);
    double2 x = (double2)(id, id) * 1.5;

	a[id] = x.x + x.y;
}