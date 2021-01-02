//#define GET_INDEX(m_i, m_j, m_q, m_p) ((m_j * ss->values_per_row) + (m_i * ss->values_per_particle) + (m_q * ss->values_per_quadrant) + m_p)
#define GET_INDEX(m_i, m_j, m_v) ((m_j * ss->values_per_row) + (m_i * ss->values_per_particle) + m_v)

#define GET_PHI(m_q, m_p) (ps->phi_start + m_q * HALF_PI + m_p * ps->phi_step_size)
