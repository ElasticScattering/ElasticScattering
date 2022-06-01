# /* -------------------------------------------------------------------------
#     This code is part of ElasticScattering.
#     Copyright(C) 2022 Elastic Scattering developers
#     This program is free software : you can redistribute it and /or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with this program.If not, see < http://www.gnu.org/licenses/>.
#    ------------------------------------------------------------------------ */

from numba import njit
import numpy as np

HBAR = 1.054571628e-34
KB = 1.38E-23

@njit()
def exp(phi, centre, width):

    phi_red = min(abs(phi - centre),
                  abs(phi - centre - 2 * np.pi),
                  abs(phi - centre + 2 * np.pi))
    phi_eff = phi_red / width
    if phi_eff < -50:
        return 0
    return np.exp(-phi_eff**2 / 2)


@njit()
def get_tau(phi, T, tau_hot):

    width = 0.1
    f = exp(phi, 0, width)
    f += exp(phi, np.pi / 2, width)
    f += exp(phi, np.pi, width)
    f += exp(phi, 3 * np.pi / 2, width)
    tau0 = HBAR / (KB * T)
    tau = tau0 * (1 - (1 - tau_hot / tau0) * f)
    return tau


@njit()
def lorentzian(pphi, kf, mass, T, tau_hot):
    """ Isotropic with effective mass.
    Use this for testing against Drude.
    """

    vf = HBAR * kf / mass

    nr = len(pphi)
    kkx = np.zeros(nr)
    kky = np.zeros(nr)
    vvx = np.zeros(nr)
    vvy = np.zeros(nr)
    ttau = np.zeros(nr)
    for i, phi in enumerate(pphi):
        vvx[i] = vf * np.cos(phi)
        vvy[i] = vf * np.sin(phi)
        kkx[i] = kf * np.cos(phi)
        kky[i] = kf * np.sin(phi)
        ttau[i] = get_tau(phi, T, tau_hot)

    return kkx, kky, vvx, vvy, ttau
