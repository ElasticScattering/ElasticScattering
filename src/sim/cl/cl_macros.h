#define GET_BASE_INDEX(m_i, m_j) ((m_j * get_global_size(0) * ss->particles_per_position) + (m_i * ss->particles_per_position))

#define GET_INDEX(m_i, m_j, m_v) ((m_j * get_global_size(0)*get_global_size(2)) + (m_i * get_global_size(2)) + m_v)

#define GET_POSITION_INDEX(m_i, m_j) ((m_j * get_global_size(0)) + (m_i))

#define GET_PHI(m_v) (ps->phi_start + floor(v / (double)ss->particles_per_quadrant) * HALF_PI + (double)(v % ss->particles_per_quadrant) * ps->phi_step_size)
