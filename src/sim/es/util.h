/* -------------------------------------------------------------------------
	This code is part of ElasticScattering.

	Copyright(C) 2022 Stijn Hinlopen

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

#include "src/sim/es/shared_macros.h"

ESCL_INLINE double SimpsonWeight(const int i, const int dim) {
    const double main_multiplier = (i % 2 == 0) ? 2.0 : 4.0;
    const bool is_edge = (i == 0) || i == (dim - 1);
    return is_edge ? 1.0 : main_multiplier;
}

inline double SimpsonWeight2D(int i, int j, int dim) {
    return SimpsonWeight(i, dim) * SimpsonWeight(j, dim);
}

ESCL_INLINE double GetSigma(double lt, double phi, double tau, double w)
{
    double z = exp(-lt / tau);

    double r = cos(phi) - cos(phi + w * lt) * z;

    r += w * tau * sin(phi + w * lt) * z;
    r -= w * tau * sin(phi);

    double o = w * tau * (sin(phi + w * lt) * z - sin(phi));

    return r;
}
