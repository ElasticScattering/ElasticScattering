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

# Remove the lorentzian altogether and go for a hard boundary figure.

import os
import matplotlib.pyplot as plt
from numba import njit
import numpy as np
import scipy.integrate as si

# Constants
M0 = 9.109e-31
M = 5 * M0  # effective mass
E = 1.602e-19
HBAR = 6.626e-34 / 2 / np.pi
KB = 1.38e-23

p = 0.20  # doping
c = 1e-9  # 11.541e-10, lattice
a = 5.45e-10 / np.sqrt(2) # lattice parameter
n = (1 + p) / (a * a * c)  # charge density
kf = np.sqrt(n * c * 2 * np.pi)  # corresponding kf cylindrical Fermi surface

HT_GOAL = 20
MAX_H = 150

print('Charge density:', n)


@njit()
def get_planckian_tau(T):
    return HBAR / (KB * T)


@njit()
def omegac(b, m):
    return E * b / m


@njit()
def vf(kf, m):
    return HBAR * kf / m


@njit()
def vx(phi, kf, m):
    return vf(kf, m) * np.cos(phi)


@njit()
def vy(phi, kf, m):
    return vf(kf, m) * np.sin(phi)


@njit()
def odds(t, tau):
    """ Odds to survive a time t """
    return np.exp(-t / tau)


@njit()
def const(kf, m, c):
    """ Constant for sigma when doing the integral. """
    jacobian = kf  # phi integral
    dfde = 1 / (HBAR * vf(kf, m))  # dk/dE really.
    theory = E**2 / (4 * np.pi**3)
    zdirec = 2 * np.pi / c
    units = 1e-8  # (Omega-m)^-1 --> (mu-Omega-cm)^-1
    return jacobian * dfde * theory * zdirec * units


@njit()
def integrandxx(t, phi, B, kf, m, tau):
    return vx(phi, kf, m) * vx(phi - omegac(B, m) * t, kf, m) * odds(t, tau)


@njit()
def integrandxy(t, phi, B, kf, m, tau):
    return vx(phi, kf, m) * vy(phi - omegac(B, m) * t, kf, m) * odds(t, tau)


@njit()
def integrandyy(t, phi, B, kf, m, tau):
    return vy(phi, kf, m) * vy(phi - omegac(B, m) * t, kf, m) * odds(t, tau)


def integrate_1fluid(h, kf, m, tau, is_xx, pen=0):

    def timebound(angle):
        return min(15 * tau, (angle - minphi + pen) / omegac(h, m))

    if is_xx:
        formula = integrandxx
        sgn = 1
    else:
        formula = integrandxy
        sgn = -1

    minphi = 0
    res1 = si.dblquad(formula, 0, np.pi / 2, 0, timebound,
                      args=((h, kf, m, tau)))
    minphi = np.pi / 2
    res2 = si.dblquad(formula, np.pi / 2, np.pi, 0, timebound,
                      args=((h, kf, m, tau)))
    minphi = np.pi
    res3 = si.dblquad(formula, np.pi, np.pi * 1.5, 0, timebound,
                      args=((h, kf, m, tau)))
    minphi = 3 * np.pi / 2
    res4 = si.dblquad(formula, np.pi * 1.5, np.pi * 2, 0, timebound,
                      args=((h, kf, m, tau)))

    return sgn * (res1[0] + res2[0] + res3[0] + res4[0]) * const(kf, m, c)


# Assert that at B~=0 Drude is recovered as a check.
base_xx = integrate_1fluid(0.00001, kf, M, 1e-12, True)
base_xy = integrate_1fluid(0.00001, kf, M, 1e-12, False)
drude_xx = n * E**2 * 1e-12 / 5 / M0 * 1e-8
assert(abs(base_xy) < base_xx * 1e-5)
assert(abs(base_xx - drude_xx) < drude_xx * 1e-5)

##################################################################
# Theory
# Application
##################################################################


def make_1fluid_graphs(axXX, axXY, axDXX, axNH):

    tt = [5, 10, 20, 30, 40, 65, 80, 110, 140, 180]
    colors = plt.get_cmap('inferno', len(tt) + 2).colors[:-2]

    for color, temp in zip(colors, tt):
        tau = get_planckian_tau(temp)
        hh = np.linspace(0.01, MAX_H, 200)[1:]
        sxx = np.array([integrate_1fluid(h, kf, M, tau, True) for h in hh])
        sxy = np.array([integrate_1fluid(h, kf, M, tau, False) for h in hh])
        rxx = sxx / (sxx**2 + sxy**2)
        rxy = sxy / (sxx**2 + sxy**2)

        drxx = (rxx[1:] - rxx[:-1]) / (hh[1:] - hh[:-1])
        dhh = 0.5 * (hh[1:] + hh[:-1])

        kw = {'color': color, 'lw': 8, 'label': f'{temp:.0f} K'}
        if axDXX is not None:
            axDXX.plot(dhh / temp, drxx, **kw)
        if axXX is not None:
            axXX.plot(hh, rxx, **kw)
        if axXY is not None:
            axXY.plot(hh / temp, rxy / hh, **kw)
        if axNH is not None:
            axNH.plot(hh / temp, hh / rxy / E * 1e8, **kw)


def quad(x, a, c):
    return a * np.sqrt(1 + (c * x / a)**2)


def dquad(x, a, c):
    return c**2 * x / a / np.sqrt(1 + (c * x / a)**2)


##################################################################
# Application
# Lay-out
##################################################################

def add_approx_to_graphs(axDXX):

    Hlinear_slope = 2 * np.pi / (4 + (np.pi - 2)**2) / n / E * 1e8
    RH_saturation = np.pi * (np.pi - 2) / (4 + (np.pi - 2)**2) / n / E * 1e8
    Hlinear_inset = 2 / (n * E * np.pi) * 1e8
    rxx1 = M / (n * E**2 * get_planckian_tau(0.1)) * 1e8
    rxx100 = M / (n * E**2 * get_planckian_tau(100)) * 1e8
    approx_color = [0.319809, 0.770914, 0.411152, 1.]

    axDXX.plot([0, MAX_H * 2], [Hlinear_slope] * 2,
               color=approx_color, lw=4, dashes=[6, 2], zorder=-1)
    axDXX.scatter(0, Hlinear_inset, color='black', s=150, zorder=10)
    cut = KB * 5 * M0 * np.pi / 2 / E / HBAR
    print('Field scale:', cut, 'K/T')
    axDXX.plot([cut, cut], [0, plt.ylim()[1] * 2],
               color=approx_color, lw=4, dashes=[6, 2], zorder=-1)

    # plt.figure('Rxx')
    # plt.plot([0, 50], [rxx100, rxx100 + 50 * Hlinear_inset],
    #          color=approx_color, lw=8)
    # plt.plot([0, max(hh)], [rxx1, rxx1 + Hlinear_slope * max(hh)],
    #          color=approx_color, lw=8)

    # plt.figure('Rxy')
    # plt.plot([0, MAX_H], [RH_saturation] * 2, color=approx_color, lw=8)

    # plt.figure('nH')
    # plt.plot([0, MAX_H], [n] * 2, color='green', lw=8)
    # plt.plot([0, MAX_H], [1 / RH_saturation / E * 1e8]
    #          * 2, color=approx_color, lw=2)


def paper_fig_layout(fig, ax_rxx, ax_deriv):

    ax_deriv.set_xlabel('$B/T$ (T/K)')
    ax_deriv.set_ylabel('$d\u03C1_{xx}/dB$ (\u03BC\u03A9cm/T)')
    ax_deriv.set_xlim(0, 10)
    ax_deriv.set_ylim(0, 0.11)

    ax_rxx.set_xlabel('$B$ (T)')
    ax_rxx.set_ylabel('$\u03C1_{xx}$ (\u03BC\u03A9cm)')
    ax_rxx.set_ylim(0, 65)
    ax_rxx.set_xlim(0, MAX_H)
    ax_rxx.legend(loc='upper left', fontsize=30,
                  bbox_to_anchor=[1, 0.5], frameon=False)
    ax_rxx.set_xticks([0, 25, 50, 75, 100, 125])

    fig.subplots_adjust(right=0.75, bottom=0.1,
                        top=0.95, left=0.15, hspace=0.3)

    approx_color = [0.319809, 0.770914, 0.411152, 1.]
    ax_rxx.annotate('b)', xy=(0.01, 0.94), xycoords='figure fraction')
    ax_rxx.annotate('c)', xy=(0.01, 0.46), xycoords='figure fraction')
    ax_deriv.annotate('$B^*/T$', xy=(0.52, 0.13),
                      xycoords='figure fraction', color=approx_color)
    ax_deriv.annotate('$B\u279D\u221e$', [0.1, 0.095], color=approx_color)


def savefig(fig):
    path = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(path, 'code_last_iteration.png')
    fig.savefig(path, dpi=300)

    print(f'The figure has been exported to {path}')


plt.rc('font', size=30)
fig, axs = plt.subplots(nrows=2, figsize=(12, 13))
add_approx_to_graphs(axs[1])
make_1fluid_graphs(axs[0], None, axs[1], None)
paper_fig_layout(fig, *axs)
savefig(fig)
