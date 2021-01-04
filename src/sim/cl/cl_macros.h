//#define GET_INDEX(m_i, m_j, m_q, m_p) ((m_j * ss->particles_per_row) + (m_i * ss->values_per_particle) + (m_q * ss->values_per_quadrant) + m_p)
#define GET_INDEX(m_i, m_j, m_v) ((m_j * ss->particles_per_row) + (m_i * ss->particles_per_position) + m_v)
#define GET_POSITION_INDEX(m_i, m_j) ((m_j * ss->particles_per_row) + (m_i))
//#define GET_PARTICLE_LT_INDEX(m_i, m_j) ((m_j * ss->particles_per_row) + (m_i * ss->values_per_particle))

//#define GET_PHI(m_q, m_p) (ps->phi_start + m_q * HALF_PI + m_p * ps->phi_step_size)

#define GET_PHI(m_v) (ps->phi_start + floor(v / (double)ss->particles_per_quadrant) * HALF_PI + (double)(v % ss->particles_per_quadrant) * ps->phi_step_size)
//double tau = (ps.is_coherent) ? ss.coherent_tau : HBAR / (KB * temperature);
