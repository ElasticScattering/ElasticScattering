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

# Re-used from graph for theorists. 16-11-2020
# Re-used from writing v3. 09-04-2021

import os
import matplotlib.pyplot as plt
import scipy.optimize as so
from common import *

# Define the incoherent fraction,
# Coherent tau
# Temperatures (incoherent tau) and fields
F_INC = 0.5
tau_coh = 3e-8
T0 = 0
hh = np.linspace(0.01, 100, 1000)
tt = np.array([5, 10, 20, 30, 40, 50, 65, 80, 110, 140, 180])
print(f'Equality temperature is {HBAR / KB / tau_coh:.1f}')

plt.rc('font', size=30)
colors = plt.get_cmap('viridis', len(tt)).colors
LW_THICK = 10
LW_THIN = 4
GREEN = [0.319809, 0.770914, 0.411152, 1.]
RED = [221 / 256, 57 / 256, 57 / 256]


def get_crossover():
    return np.pi / 2 * F_INC * KB * M / (E * HBAR)


def get_slope():
    a2 = np.pi / 2 * F_INC
    xi = (1 - np.cos(a2)) / a2
    eta = 1 - np.sin(a2) / a2
    ninc = F_INC * n
    ncoh = (1 - F_INC) * n
    return xi * ninc / E / (xi**2 * ninc**2 + (ncoh + ninc * eta)**2) * 1e8


def coherent(hh, tau_coh):

    sxx_coh_base = n * E**2 * tau_coh / M * (1 - F_INC) / 1e8
    wct_coh = E * hh / M * tau_coh
    crossover_ht = get_crossover()
    sxx_coh = sxx_coh_base / (1 + wct_coh**2)
    sxy_coh = sxx_coh * wct_coh
    print(f' Resistivity coherent alone B=0 is {1/sxx_coh_base:.2e} muOhmcm')
    print(f' Crossover is at H/T={crossover_ht:.2f} T/K')
    print(f' Tau equal at T={HBAR / KB / tau_coh:.1f} K')

    return sxx_coh, sxy_coh


def incoherent(hh, temp, T0):

    sxx_inc = np.zeros(len(hh))
    sxy_inc = np.zeros(len(hh))
    tau_inc = HBAR / KB / (T0 + temp)
    for i, h in enumerate(hh):
        sxx_inc[i] = sigma(h, None, tau_inc, 'xx', F_INC, 0, True)
        sxy_inc[i] = -sigma(h, None, tau_inc, 'xy', F_INC, 0, True)

    return sxx_inc, sxy_inc


class Transport():

    def __init__(self, hh):
        self.hh = hh

    sxx_coh = None
    sxx_inc = None
    sxy_coh = None
    sxy_inc = None

    @property
    def sxx(self):
        return self.sxx_coh + self.sxx_inc

    @property
    def sxy(self):
        return self.sxy_coh + self.sxy_inc

    @property
    def rxx(self):
        return self.sxx / (self.sxx**2 + self.sxy**2)

    @property
    def rxy(self):
        return self.sxy / (self.sxx**2 + self.sxy**2)

    @property
    def dhh(self):
        return 0.5 * (self.hh[1:] + self.hh[:-1])

    @property
    def drxxdh(self):
        return (self.rxx[1:] - self.rxx[:-1]) / (self.hh[1:] - self.hh[:-1])


def get_transport(hh, tau_coh, temp, T0):

    tr = Transport(hh)
    tr.sxx_coh, tr.sxy_coh = coherent(hh, tau_coh)
    tr.sxx_inc, tr.sxy_inc = incoherent(hh, temp, T0)
    return tr


def plotR(axR, axQ, tr, temp, color, quadfit=False):

    axR.plot(tr.hh, tr.rxx, lw=LW_THIN, label=f'{temp} K', color=color)
    axQ.plot(tr.dhh / temp, tr.drxxdh, lw=LW_THICK,
             label=f'{temp} K', color=color, zorder=50)

    if quadfit:
        xx = np.linspace(0, max(tr.hh), 1000)
        p, perr, yy = fit_dquadrature(tr.dhh, tr.drxxdh, temp, xeval=xx)
        print('Fit parameters quadrature:')
        print(f'  alpha: {abs(p[0]):.5f} +- {perr[0]:.5f}')
        print(f'  gamma: {abs(p[1]):.5f} +- {perr[0]:.5f}')
        print(f'  g/a  : {abs(p[1] / p[0]):.2f}')

        xxquad = np.linspace(0, max(xx) * 1.2, 1000)
        yyquad = dquadraturedh(xxquad, *p, temp)
        axQ.set_ylim(bottom=0, top=get_slope() * 1.2)
        axQ.plot(xxquad / temp, yyquad, lw=2, color=RED,
                 zorder=100, dashes=[6, 2])
        axQ.plot([p[1] / p[0]] * 2, axQ.get_ylim(),
                 color=RED, lw=LW_THIN, zorder=-1, dashes=[6, 2])
        axQ.plot(axQ.get_xlim(), [get_slope()] * 2,
                 color=GREEN, lw=LW_THIN, zorder=-2, dashes=[6, 2])
        axQ.plot([get_crossover()] * 2, axQ.get_ylim(),
                 color=GREEN, lw=LW_THIN, zorder=-2, dashes=[6, 2])


def plotKohler():
    pass


def plotHall(axH, tr, temp, color, case):

    rH = tr.rxy / hh
    nH = 1 / rH / E * 1e8
    axH.plot(hh, nH / 1e27, label=f'{temp} K', lw=LW_THIN, color=color)


def layout(figMain, axR, axQ, axH):
    """ Responsible for the axes labels, limits and general layout. """

    # Main
    figMain.subplots_adjust(right=0.77, top=0.98, bottom=0.05, left=0.15,
                            hspace=0.3)

    # Top panel resistivity
    axR.set_xlabel('$B$ (T)')
    axR.set_ylabel('$\u03C1_{xx}$ (\u03BC\u03A9cm)')
    axR.set_xlim(0, 100)
    axR.set_ylim(bottom=0)
    axR.set_xticks([0, 20, 40, 60, 80])

    # Middle panel quadrature
    axQ.set_ylabel('$d\u03C1_{xx}/dB$ (\u03BC\u03A9cm/T)')
    axQ.set_xlabel('$B/T$ (T/K)')
    axQ.set_xlim(0, 5)
    axQ.set_ylim(bottom=0)
    axQ.plot([], [], color=RED, dashes=[6, 2],
             lw=2, label='fit')
    axQ.legend(loc='upper left', bbox_to_anchor=[1, 1.2], frameon=False)

    # Bottom panel Hall
    axH.plot([min(hh), max(hh)], [n / 1e27] * 2, color=GREEN,
             lw=LW_THIN, dashes=[6, 2])
    axH.plot([min(hh), max(hh)], [n * (1 - F_INC) / 1e27] * 2,
             lw=LW_THIN, color=RED, dashes=[6, 2])
    axH.set_xlabel('$B$ (T)')
    axH.set_ylabel('$n_H$ (10$^{27}$ m$^{-3}$)')
    axH.set_ylim(bottom=0, top=n * 1.1 / 1e27)
    axH.set_xlim(0, 100)
    axH.set_yticks([0, 2.5, 5, 7.5, 10])

    axR.annotate('a)', xy=(0.01, 0.98), xycoords='figure fraction')
    axQ.annotate('b)', xy=(0.01, 0.64), xycoords='figure fraction')
    axH.annotate('c)', xy=(0.01, 0.31), xycoords='figure fraction')
    axQ.annotate('$B^*$', xy=[3.1, 0.005], color=GREEN)
    axH.annotate('$n_{tot}$', xy=[8, 8.5], color=GREEN)
    axH.annotate('$n_{coh}$', xy=[8, 3.2], color=RED)
    axQ.annotate('$B\u279D\u221e$', [0.1, 0.045], color=GREEN)


def guidance_lin_quad(axR, trmaxT):
    """ Make a dashed B-lin and B-quad line. """

    def func(x, a, b):
        return a + b * x**2

    wh = trmaxT.hh < 90
    p, _ = so.curve_fit(func, trmaxT.hh[wh], trmaxT.rxx[wh])
    print(f'Expected B^2 coeff {0.03833 / 180:.2e}, got {p[1]:.2e}')
    xx = np.linspace(0, 70, 1000)
    yy = func(xx, *p) - 0.07
    xx, yy = xx[yy > 0], yy[yy > 0]

    axR.plot(xx, yy, color='black', dashes=[6, 2], lw=LW_THIN)
    axR.annotate('$\sim$B$^2$', [50, 0.2])

    xx2 = np.linspace(30, 80, 1000)
    yy2 = xx2 * 0.04301 + 0.05
    axR.plot(xx2, yy2, color='black', dashes=[6, 2], lw=LW_THIN)
    axR.annotate('$\sim$B', [xx2[400], yy2[400] + 0.6])


def main(temperatures, hh, tau_coh, T0):

    case = f'tau_coh={tau_coh:.2e}, T0={T0:.0f}'
    Tfit = temperatures[2]
    figMain, (axR, axQ, axH) = plt.subplots(nrows=3, figsize=(12, 22))

    sxx_coh, sxy_coh = coherent(hh, tau_coh)
    for color, temperature in zip(colors, temperatures):
        print()
        tr = get_transport(hh, tau_coh, temperature, T0)
        plotR(axR, axQ, tr, temperature, color, quadfit=temperature == Tfit)
        plotHall(axH, tr, temperature, color, case)
    layout(figMain, axR, axQ, axH)

    tr = get_transport(hh, tau_coh, max(temperatures), T0)
    guidance_lin_quad(axR, tr)

    path = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join('code_last_iteration.png')
    figMain.savefig(path, dpi=300)


print(f'At crossover, wct_coh={E*get_crossover()*tau_coh / M:.2f}')
main(tt, hh, tau_coh, T0)
