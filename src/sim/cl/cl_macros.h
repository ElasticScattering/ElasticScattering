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

#define GET_BASE_INDEX(m_i, m_j) ((m_j * get_global_size(0) * ss->particles_per_position) + (m_i * ss->particles_per_position))

#define GET_INDEX(m_i, m_j, m_v) ((m_j * get_global_size(0)*get_global_size(2)) + (m_i * get_global_size(2)) + m_v)

#define GET_POSITION_INDEX(m_i, m_j) ((m_j * get_global_size(0)) + (m_i))

#define GET_PHI(m_v) (ps->phi_start + floor(v / (double)ss->particles_per_quadrant) * HALF_PI + (double)(v % ss->particles_per_quadrant) * ps->phi_step_size)
