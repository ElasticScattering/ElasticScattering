#pragma once

const double PI   = 3.141592653589793238463;
const double PI2  = PI * 2.0;
const double M0   = 9.109e-31;
const double E    = 1.602e-19;
const double HBAR = 1.055e-34;
const double C    = 1.15e-9;

#define MIN(m_a, m_b) (((m_a) < (m_b)) ? (m_a) : (m_b))

typedef struct {
    double x, y;
} v2;
