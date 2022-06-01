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

#include <iostream>

#define CL_FAIL_CONDITION(err, msg)                                                                                                                                                          \
    if (err != CL_SUCCESS) {                                                                                                                                                                    \
        std::cout << "\x1B[31mCL Error [" << CLErrorString(err) << "]\x1B[0m: " << msg << "\n\tat: " << "L" << __LINE__ << ": " << __FILE__ << "\n\tin function: " << __FUNCTION__ << std::endl;   \
        exit(0);                                                                                                                                                                                \
    }

#define FAIL_CONDITION(cond, msg)                                                                                                                            \
	if (cond) {                                                                                                                                                 \
		std::cout << "\x1B[30mError\x1B[0m: " << msg << "\n\tat: " << "L" << __LINE__ << ": " << __FILE__ << "\n\tin function: " << __FUNCTION__ << std::endl;  \
		exit(0);                                                                                                                                                \
	}

#define CFG_EXIT_CONDITION(cond, msg)                                     \
	if (cond) {                                                           \
		std::cout << "\x1B[95mConfig Error\x1B[0m: " << msg << std::endl; \
		exit(0);                                                          \
	}

