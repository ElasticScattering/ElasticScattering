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

import matplotlib.pyplot as plt
from numba import njit
import numpy as np
import os
import scipy.optimize as so
import time

from drude import simple_drude
from lorentzian import lorentzian

###############
# Global
###############

E = 1.602176487e-19
HBAR = 1.054571628e-34
ME = 9.10938215e-31
KB = 1.38E-23

A = 3.777e-10  # From Cava1987, table 2/3
C = 13.226e-10 / 2

plt.rc('font', size=30)
LW = 5

compute_FS = None


# for paper plot
LW_DATA = 5
LW_GUIDE = 3
C_GUIDE = 'black'
GREEN = [0.319809, 0.770914, 0.411152, 1.]
RED = [221 / 256, 57 / 256, 57 / 256]
DASH_GUIDE = [6, 2]


# Compiling is the default, my benchmark show a 150x speedup
# Not compiling can be useful if errors occur or new code is added.
COMPILE = True
if not COMPILE:
    def njit():
        def wrapper(func):
            return func
        return wrapper

###############
# Part 1
###############


# If the backend is not numba then disable the marked njit.
# marked
@njit()
def precompute(nr, *fsargs):
    """ Load all angles, kx, ky, vx, vy, tau """

    assert(nr % 2 == 0)
    pphi = np.linspace(0, 2 * np.pi * (nr - 1) / nr, nr)
    r = compute_FS(pphi, *fsargs)
    return pphi, r[0], r[1], r[2], r[3], r[4]

# marked


@njit()
def stepanalysis(pphi, *fsargs):
    """ Return this sigma component for the given B values.
    Also returns an array matching pphi with omegac values at B=1
    """

    # Determine exps [survival odds to go from each point to the next at B=1]
    # and tsteps [time it takes from one point to the next at B=1]
    N = len(pphi)
    exps = np.zeros(N)
    tsteps = np.zeros(N)
    for i, phi in enumerate(pphi):
        if i < N - 1:
            next_phi = pphi[i + 1]
        else:
            next_phi = pphi[0] + 2 * np.pi  # should just be 2pi
        tsteps[i], exps[i] = get_w_and_dt(phi, next_phi, *fsargs)

    # Just to make sure this is reasonable. Usually 1-100 ps or so.
    circulation_time = np.sum(tsteps)
    assert(1e-20 < circulation_time < 1e-8)

    omegac = np.zeros(len(pphi))
    omegac[:-1] = (pphi[1:] - pphi[:-1]) / tsteps[:-1]
    omegac[-1] = (pphi[0] - pphi[-1] + 2 * np.pi) / tsteps[-1]
    return tsteps, exps, omegac


###############
# Part 2
###############


@njit()
def integrate_zero_field(pphi, vvx, vvy, kkf, ttau):
    """ For the case B=0 there is a special case as here
    the wc representation diverges.
    Note that low field there is no taylor implemented,
    I tried but it lead to discontinuities whereever I terminated it.
    """

    assert(len(pphi) % 2 == 0 and "has to be *even* length for Simpson")
    assert(len(vvx) == len(pphi))
    assert(len(vvy) == len(pphi))
    assert(len(kkf) == len(pphi))
    assert(len(ttau) == len(pphi))

    total_xx = 0
    for i, phi in enumerate(pphi):
        integrand = vvx[i]**2 * ttau[i] * \
            kkf[i] / np.sqrt(vvx[i]**2 + vvy[i]**2)
        if i == 0:
            factor = 2  # combined the first and last, since phi is cyclic
        elif i % 2:
            factor = 4
        else:
            factor = 2
        total_xx += factor * integrand

    assert(factor == 4)  # guarantee that we finish with 4 before we wrap back.
    dphi = pphi[1] - pphi[0]

    return total_xx * dphi / 3, 0


@njit()
def integrate_phi(pphi, vvx, vvy, kkf, ww, ddt, ttau, wwc):
    """ Perform the phi integral and outsource the internal integrals. """

    # Even because element 0 is also the last element.
    assert(len(pphi) % 2 == 0 and "has to be *even* length for Simpson")
    assert(len(vvx) == len(pphi))
    assert(len(vvy) == len(pphi))
    assert(len(kkf) == len(pphi))
    assert(len(ww) == len(pphi))
    assert(len(ddt) == len(pphi))

    # Note that usually there is 1424242424242...41
    # Phi is cyclic, the first element is the last element,
    # so we simply visit everything once and let it be 242424242...4
    total_xx = 0
    total_xy = 0
    for i, phi in enumerate(pphi):
        vxx, vxy = integrate_time(i, pphi, vvx, vvy, kkf, ww, ddt, ttau, wwc)
        if i == 0:
            factor = 2
        elif i % 2:  # odd indices
            factor = 4
        else:
            factor = 2
        total_xx += vxx * factor
        total_xy += vxy * factor

    assert(factor == 4)  # guarantee that we finish with 4 before we wrap back.
    dphi = pphi[1] - pphi[0]

    return total_xx * dphi / 3, total_xy * dphi / 3


@njit()
def integrate_time(index, pphi, vvx, vvy, kkf, ww, ddt, ttau, wwc):
    """ Get the time integral for the phi value at index. """

    v0x = vvx[index]
    v0y = -vvy[index]
    vf0 = np.sqrt(vvx[index]**2 + vvy[index]**2)
    if ww[index] < 1e-1:
        wct = wwc[index] * ttau[index]
        drude_locally = kkf[index] * vvx[index] * ttau[index] / vf0
        drude_locally /= 1 + wct**2
        return v0x * drude_locally, v0y * drude_locally * wct

    N = len(pphi)
    survivalodds = 1
    position = index

    # Run down the orbit. Terminate on 2 conditions:
    #   1) The particle is dead. Negligible survival odds.
    #   2) We come full circle.
    # Use trapezoidal integration.
    integral = 0
    while survivalodds > 1e-8 and (position != index or integral == 0):
        # Previous position
        leftbound = vvx[position] * survivalodds

        # New position
        # Note that the ww and ddt array are 1 less in length,
        # [position] there dictates what happens between [position and position+1]
        # in the longer arrays of vx, vy, kx, ky. Hence ask for them
        # *before* incrementing.
        survivalodds *= ww[position]
        time_step = ddt[position]
        position = (position + 1) % N

        rightbound = vvx[position] * survivalodds
        integral += 0.5 * (leftbound + rightbound) * time_step

    # Relevant for when full orbits are completed
    integral /= 1 - survivalodds

    wc = (pphi[(index + 1) % N] - pphi[index]) / ddt[index]
    real = -ttau[index] / (1 + wc**2 * ttau[index]**2)
    real *= wc * ttau[index] * np.sin(pphi[index]) - np.cos(pphi[index])
    real *= np.sqrt(vvx[index]**2 + vvy[index]**2)

    # Velocity at time 0, determines if we do sigma_xx or sigma_xy
    # Really sigma_yx so we add the minus from Onsager's relations
    # using sigma_xy = -sigma_yx here right away to get that over with.
    integralxx = integral * v0x * kkf[index] / vf0
    integralxy = integral * v0y * kkf[index] / vf0
    return integralxx, integralxy


# marked
@njit()
def get_w_and_dt(phi1, phi2, *fsargs):
    """ Obtain the time t and survival oods exp(integral(-dt/tau)) """

    assert(0 < phi2 - phi1 < 0.1)

    ssubphi = np.linspace(phi1, phi2, 15)
    kkx, kky, vvx, vvy, ttau = compute_FS(ssubphi, *fsargs)

    timer = 0  # integral dt = dk/kdot
    exponent = 0  # integral dk/(kdot * tau)
    for index in range(len(ssubphi) - 1):
        dk = np.sqrt((kkx[index] - kkx[index + 1])**2 +
                     (kky[index] - kky[index + 1])**2)

        wc1 = E / HBAR * np.sqrt(vvx[index]**2 + vvy[index]**2)  # B=1
        wct1 = wc1 * ttau[index]
        wc2 = E / HBAR * np.sqrt(vvx[index + 1]**2 + vvy[index + 1]**2)
        wct2 = wc2 * ttau[index + 1]
        dt = dk * 0.5 * (1 / wc1 + 1 / wc2)

        timer += dt
        exponent += dk * 0.5 * (1 / wct1 + 1 / wct2)
    return timer, np.exp(-exponent)


@njit()
def prefactors():
    return E**2 / (2 * np.pi**2 * C * HBAR)


def sigma(Bvals, pphi, kkx, kky, vvx, vvy, ttau, exps, tsteps, wwc):
    """ Compute sigma_xx or sigma_xy for all magnetic field values given. """

    kkf = np.sqrt(kkx**2 + kky**2)
    sxx = np.zeros(len(Bvals))
    sxy = np.zeros(len(Bvals))
    st = time.time()
    elapsed = 0
    for i, B in enumerate(Bvals):
        if B == 0:
            sxx[i], sxy[i] = integrate_zero_field(pphi, vvx, vvy, kkf, ttau)
        else:
            ww = exps ** (1 / B)
            ddt = tsteps / B

            sxx[i], sxy[i] = integrate_phi(pphi, vvx, vvy, kkf,
                                           ww, ddt, ttau, wwc)

        iter_time = time.time() - st - elapsed
        if iter_time == 0:
            print('0 iter time ???')
            continue
        elapsed = time.time() - st
        print(f'Finished B={B:.3f}. '
              f'Busy for {elapsed:.3f} s, '
              f'expect another {(len(Bvals) - i - 1) * iter_time:.2f} s')

    sxx *= prefactors()
    sxy *= prefactors()
    return sxx, sxy


###############
# Execution
###############

def quad(x, a, c):
    return a * np.sqrt(1 + (c * x / a)**2)


def dquad(x, a, c):
    return c**2 * x / a / np.sqrt(1 + (c * x / a)**2)


def deriv(xx, yy, distance):
    """ Create running derivative, symmetrical on each point.
    Shortens the array by distance on both ends. """

    dx = xx[distance:-distance]
    dy = np.zeros(len(dx))
    for i in range(distance, len(xx) - distance):
        x = xx[i - distance:i + distance + 1]
        y = yy[i - distance: i + distance + 1]
        p = np.polyfit(x, y, 1)
        dy[i - distance] = p[0]
    return dx, dy


def run_Drude_test():
    """ Testing ground. Basically a visual unit test. """

    global compute_FS
    compute_FS = simple_drude
    compute_FS.recompile()

    KF = 0.75e10
    n = KF**2 / (2 * np.pi * C)
    print(f'Isotropic FS with p={n * A**2 * C:.2f} and mass 5m0.')

    # The entire computation is these few lines
    Bvals = np.linspace(0, 100, 131)
    pphi, kkx, kky, vvx, vvy, ttau = precompute(3000, KF, 5 * ME)
    ttsteps, eexps, wwc = stepanalysis(pphi, KF, 5 * ME)
    sxx, sxy = sigma(Bvals, pphi, kkx, kky, vvx,
                     vvy, ttau, eexps, ttsteps, wwc)

    # Test that wc is isotropic and the right value
    assert(all(np.abs(wwc - E / 5 / ME) < 1e-6 * E / 5 / ME))

    # Test that sigmaxx and sigmaxy have the right values at all fields
    drude = n * E**2 * ttau[0] / (5 * ME)
    print(
        f'Situation is n={n:.2e} /m3 and sxxDrude={drude/1e8:.2e} (muOhmcm)^-1 at B=0')
    wct = E * Bvals * ttau[0] / (5 * ME)
    drude_xx = drude / (1 + wct**2)
    drude_xy = drude_xx * wct

    plt.figure('Situation')
    plt.plot(pphi, kkx / max(kkx), label='kx')
    plt.plot(pphi, kky / max(kky), label='ky')
    plt.plot(pphi, vvx / max(vvx) + 0.1, label='vx + 0.1')
    plt.plot(pphi, vvy / max(vvy) + 0.2, label='vy + 0.2')
    plt.plot(pphi, ttau / max(ttau) + 0.3, label='tau + 0.3')
    plt.plot(pphi, wwc / max(wwc) + 0.4, label='wc + 0.4')
    plt.xlabel('angle (rad)')
    plt.ylabel('normed')
    plt.legend()

    plt.figure('Conductivity')
    plt.plot(Bvals, sxx / 1e8, label='sxx', lw=4, color='red')
    plt.plot(Bvals, sxy / 1e8, label='sxy', lw=4, color='blue')
    plt.plot(Bvals, drude_xx / 1e8, label='drude xx',
             dashes=[6, 2], color='black')
    plt.plot(Bvals, drude_xy / 1e8, label='drude xy',
             dashes=[6, 2], color='white')
    plt.legend()
    plt.xlabel('Field (T)')
    plt.ylabel('Conductivity ($\mu\Omega cm^{-1}$)')

    plt.figure('MR')
    plt.plot(Bvals, sxx / (sxx**2 + sxy**2) * 1e8, label='rho_xx')
    plt.plot(Bvals, sxy / (sxx**2 + sxy**2) * 1e8, label='rho_xy')
    plt.xlabel('B (T)')
    plt.ylabel('rho ($\Omega m$)')

    plt.show()


def estimate_breakdown_field(pphi, ttau, wwc, phi_from, phi_to):
    """ Bbr = integral_phi  of m*/(E tau) over this hotspot """

    Bstar = 0
    for i, (wc, tau) in enumerate(zip(wwc, ttau)):
        if phi_from < pphi[i] < phi_to:
            if i % 2:
                factor = 4
            else:
                factor = 2

            Bstar += factor / (wc * tau)
    Bstar *= (pphi[1] - pphi[0]) / 3
    return Bstar


def estimate_Bstar(mass, T):
    """ When is wc*tau = pi/2 """
    return np.pi * KB * T * mass / (2 * HBAR * E)

def simple_lorentz(Bvals, T, tau_hot, *, only_breakdown=False):

    global compute_FS
    compute_FS = lorentzian
    compute_FS.recompile()

    # Parameters from lee hone at p=0.25 from Caitlin's code
    kf = 7e9
    m = 5 * ME
    print(f'n={kf**2 / 2 / np.pi / C:.2e}')

    print('Precomputing Fermi surface')
    pphi, kkx, kky, vvx, vvy, ttau = precompute(20000, kf, m, T, tau_hot)
    ttsteps, eexps, wwc = stepanalysis(pphi, kf, m, T, tau_hot)
    Bbr = estimate_breakdown_field(
        pphi, ttau, wwc, np.pi / 2 - 0.2, np.pi / 2 + 0.2)
    print(f'Breakdown field is {Bbr:.0f} T')
    print(f'Bstar is {estimate_Bstar(m, T):.1f} T')

    # print(f'B* = {quick_Bstar(m, )')
    if only_breakdown:
        return Bbr

    print('Computing sigma')
    sxx, sxy = sigma(Bvals, pphi, kkx, kky, vvx,
                     vvy, ttau, eexps, ttsteps, wwc)

    make_situation_plot(pphi, kkx, kky, vvx, vvy, ttau, wwc)
    return sxx, sxy


def make_situation_plot(pphi, kkx, kky, vvx, vvy, ttau, wwc):
    """ Show the anisotropy of various quantitites involved. """

    plt.figure('Scattering rate')
    plt.semilogy(pphi, ttau * 1e12)
    plt.xlabel('$\phi$ (rad)')
    plt.ylabel('$\u03C4_0$ (ps)')
    plt.ylim(bottom=1e-7)

    plt.figure('Situation')
    plt.plot(pphi, kkx / max(kkx), label='kx')
    plt.plot(pphi, kky / max(kky), label='ky')
    plt.plot(pphi, vvx / max(vvx), label='vx')
    plt.plot(pphi, vvy / max(vvy), label='vy')
    plt.plot(pphi, ttau / max(ttau), label='tau')
    plt.plot(pphi, wwc / max(wwc), label='wc')
    # mass = E * B / wc, B=1 so you get mass = E / wc
    # Normalising for anisotropy drops E as well for 1/wc.
    plt.plot(pphi, E / wwc / ME, label='cycl mass (M0)', lw=LW)
    plt.xlabel('angle (rad)')
    plt.ylabel('normed')
    plt.legend()

    plt.figure('FS')
    plt.plot(kkx, kky)
    plt.xlabel('k$_x$ (m^${-1}$)')
    plt.xlabel('k$_y$ (m^${-1}$)')
    plt.subplots_adjust(left=0.15, bottom=0.15)

    plt.figure('L orbit')
    plt.plot(vvx * ttau, vvy * ttau)
    plt.xlabel('l$_x$ (m)')
    plt.xlabel('l$_y$ (m)')
    plt.subplots_adjust(left=0.15, bottom=0.15)


def show_resistivity(Bvals, sxx, sxy):
    """ Show the results of the computation. """

    plt.figure('Conductivity')
    plt.plot(Bvals, sxx / 1e8, label='sxx', lw=LW)
    plt.plot(Bvals, sxy / 1e8, label='sxy', lw=LW)
    plt.legend()
    plt.xlabel('B (T)')
    plt.ylabel('$\sigma$ ($\mu\Omega cm^{-1}$)')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    plt.figure('Resistivity')
    rxx = sxx / (sxx**2 + sxy**2) * 1e8
    rxy = sxy / (sxx**2 + sxy**2) * 1e8
    plt.plot(Bvals, rxx, label='\u03C1$_{xx}$', lw=LW)
    plt.plot(Bvals, rxy, label='\u03C1$_{xy}$', lw=LW)
    plt.xlabel('B (T)')
    plt.ylabel('\u03C1 (\u03BC\u03A9cm)')
    plt.legend()
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    plt.figure('Quadrature')
    dB, drxx = deriv(Bvals, rxx, 1)
    p, pcov = so.curve_fit(dquad, dB, drxx)
    xx = np.linspace(0, max(dB) * 1.05, 10000)
    yy = dquad(xx, *p)
    plt.plot(dB, drxx, lw=LW)
    plt.plot(*deriv(Bvals, rxy, 1), label='rxy', lw=LW)
    plt.xlabel('H (T)')
    plt.ylabel('d\u03C1/dB (\u03BC\u03A9cm/T)')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.xlim(left=0)
    plt.ylim(bottom=0)


def get_filepath(identifier):
    """ Make the filepath to this particular ID in storage. """

    if not isinstance(identifier, str):
        identifier = f'{identifier:02d}'
    if not len(identifier) == 2:
        raise ValueError(
            'identifier has to be 2 digits or an int <100 and >=0')

    scriptpath = os.path.realpath(__file__)
    scriptdir = os.path.dirname(scriptpath)
    return os.path.join(scriptdir, f'{identifier} data.dat')


def store_in_file(BB, sxx, sxy):
    """ Put the results in a new file

    Useful for long computations because then you can
    quickly reload while adjusting the plots.
    """

    direc = os.path.dirname(get_filepath(0))
    existing = os.listdir(direc)
    newid = 0
    while any(f.startswith(f'{newid:02d}') for f in existing):
        newid += 1
    filepath = get_filepath(newid)

    string = 'B sxx sxy\n'
    for b, sx, sy in zip(BB, sxx, sxy):
        string += f'{b:.5f} {sx:.15f} {sy:.15f}\n'

    with open(filepath, 'w') as f:
        f.write(string)
    return newid


def load_from_file(identifier):
    """ Get the results from this old computation """

    filepath = get_filepath(identifier)

    B, sxx, sxy = [], [], []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('B'):
                break
        for line in f:
            if line:
                parts = line.split()
                B.append(float(parts[0]))
                sxx.append(float(parts[1]))
                sxy.append(float(parts[2]))

    print(f'Loaded {identifier} from file')
    return np.array(B), np.array(sxx), np.array(sxy)


################

def _plot_tau(ax1, T, tau_hot):
    """ Make a plot (phi, tau0) of this case. """

    global compute_FS
    compute_FS = lorentzian

    kf = 7e9
    m = 5 * ME
    pphi, kkx, kky, vvx, vvy, ttau = precompute(20000, kf, m, T, tau_hot)

    pphi = np.array(list(pphi) + [2 * np.pi])
    ttau = np.array(list(ttau) + [ttau[0]])

    ax1.semilogy(pphi, ttau * 1e12, color='black', lw=LW_GUIDE)
    ax1.set_ylim(1e-6, 2)
    ax1.set_xlim(0, 2 * np.pi)
    ax1.set_xlabel('$\phi$ (rad)')
    ax1.set_ylabel('$\u03C4_0$ (ps)')


def _plot_half_data(ax2, ax3, fileid, divT):
    """ Plot either breakdown OR near-infinite bound results.
    Made such that there is no desync.
    Does not do the tau, only the resistivity results.

    divT: if set, divide by temperature for the horizontal ax3
    """

    # The code is structured such that the data is generated in one run
    # using sigma() and stored in a file. The next run then loads
    # these results from file. This way, the plots were much quicker
    # to re-make and optimize for the paper while guaranteeing that
    # none of the settings changed.


    # Import
    T = 10
    Bvals, sxx, sxy = load_from_file(fileid)
    sxx /= 1e8
    sxy /= 1e8
    rxx = sxx / (sxx**2 + sxy**2)
    rxy = sxy / (sxx**2 + sxy**2)
    dBxx, drxx = deriv(Bvals, rxx, 1)
    dBxy, drxy = deriv(Bvals, rxy, 1)
    dBxx = np.array([0] + list(dBxx))
    drxx = np.array([0] + list(drxx))
    maxH = np.max(Bvals)
    maxR = max([np.max(rxx), np.max(rxy)])
    maxDR = 0.07  # shared between plots

    # Plot data
    ax2.plot(Bvals, rxx, color='tab:blue',
             label='$\u03C1_{xx}$', lw=LW_DATA)
    ax2.plot(Bvals, rxy, color='tab:orange',
             label='$\u03C1_{xy}$', lw=LW_DATA)
    xx = dBxx / T if divT else dBxx
    ax3.plot(xx, drxx, color='tab:blue',
             label='$\u03C1_{xx}$', lw=LW_DATA)
    xx = dBxy / T if divT else dBxy
    ax3.plot(xx, drxy, color='tab:orange',
             label='$\u03C1_{xy}$', lw=LW_DATA)

    # Plot guides to the eye
    Bstar = estimate_Bstar(5 * ME, 10)
    Bbr = simple_lorentz([], 10, 1e-16, only_breakdown=True)
    yy = [-1, maxR * 0.75] if divT else [-1, 2 * maxR]
    ax2.plot([Bstar, Bstar], yy, color=GREEN,
             lw=LW_GUIDE, dashes=DASH_GUIDE)
    xx = [Bstar / T] * 2 if divT else [Bstar] * 2
    ax3.plot(xx, [-1, 100], color=GREEN,
             lw=LW_GUIDE, dashes=DASH_GUIDE)
    ax2.plot([Bbr, Bbr], [-1, 100],
             lw=LW_GUIDE, dashes=DASH_GUIDE, color=RED)
    xx = [Bbr / T] * 2 if divT else [Bbr] * 2
    ax3.plot(xx, [-1, 100],
             lw=LW_GUIDE, dashes=DASH_GUIDE, color=RED)
    pos = [Bstar + maxH / 25, maxR * 1]
    pos[1] = maxR * 0.05 if divT else pos[1]
    ax2.annotate('$B^*$', pos, color=GREEN)
    ax2.annotate('$B_t$', [Bbr - maxH / 7, maxR * 0.1], color=RED)
    txt = '$B^*/T$' if divT else '$B^*$'
    pos = [Bstar + maxH / 25, maxDR * 0.05]
    pos[0] /= T if divT else 1
    ax3.annotate(txt, pos, color=GREEN)
    txt = '$B_t/T$' if divT else '$B_t$'
    pos = [Bbr - maxH / 6, maxDR * 0.5]
    pos[0] /= T if divT else 1
    ax3.annotate(txt, pos, color=RED)

    # Plot limits / labels
    ax2.set_xlim(0, maxH * 1.1)
    ax3.set_xlim(0, maxH * 1.1 / T if divT else maxH * 1.1)
    ax2.set_ylim(0, maxR * 1.1)
    ax3.set_ylim(0, maxDR)

    ax2.set_xlabel('$B$ (T)')
    ax3.set_xlabel('$B/T$ (T/K)' if divT else '$B$ (T)')
    ax2.set_ylabel('$\u03C1$ (\u03BC\u03A9cm)')
    ax3.set_ylabel('$d\u03C1/dB$ (\u03BC\u03A9cm/T)')

    xticks = np.arange(4)
    xticks *= 50 if maxH < 250 else 500
    ax2.set_xticks(xticks)
    ax3.set_xticks(xticks / T if divT else xticks)


def _plot_quadrature_fit(ax3, fileid, T):
    """ Fit the quadrature form and plot it. T=10K original curve. """

    def func(x, a, b):
        return np.sqrt((a * T)**2 + (b * x)**2)

    def dfunc(x, a, b):
        return b**2 * x / np.sqrt((a * T)**2 + (b * x)**2)

    T = 10
    Bvals, sxx, sxy = load_from_file(fileid)
    sxx /= 1e8
    sxy /= 1e8
    rxx = sxx / (sxx**2 + sxy**2)
    dB, drxx = deriv(Bvals, rxx, 1)

    p, _ = so.curve_fit(dfunc, dB, drxx)
    print(p)
    print(f'Quadrature slope {p[1]:.4f} muOhmcm/T '
          f'and scale {p[0] / p[1]:.4f} T/K (here {p[0] / p[1] * T:.2f} T)')
    xx = np.linspace(0, max(Bvals) * 1.1, 1000)
    yy = dfunc(xx, *p)

    ax3.plot(xx / T, yy, color=C_GUIDE, dashes=DASH_GUIDE, lw=LW_GUIDE)


def _plot_breakdown_fit(ax3, fileid):
    """ Do an exponential fit in high field and plot.
    Also print the exponential characteristic field scale.
    """

    Bvals, sxx, sxy = load_from_file(fileid)
    sxx /= 1e8
    sxy /= 1e8
    rxx = sxx / (sxx**2 + sxy**2)
    dB, drxx = deriv(Bvals, rxx, 1)

    def func(x, a, b, c):
        return a * np.exp(-x / b) + c

    wh = dB > 400
    p, _ = so.curve_fit(func, dB[wh], drxx[wh], p0=[0.1, 500, 0])
    print(f'Exponential fit breakdown: field scale {p[1]:.0f} T')
    xx = np.linspace(0, 5 * max(Bvals), 1000)
    yy = func(xx, *p)

    ax3.plot(xx, yy, color='black', lw=LW_GUIDE, dashes=DASH_GUIDE)


def _get_color_scaling(temperature):

    TT = [0.01, 1, 2, 5, 10, 20, 50, 100, 200]
    colors = plt.get_cmap('cool', len(TT))
    colors = [colors(i) for i in range(len(TT))]
    return colors[TT.index(temperature)]


def _plot_scaling(ax):
    """ Extra runs with different temperatures added to this plot. """

    # Manually register the runs
    # Also add to the temperatures in color_scaling!
    file_temp_map = {'20': 1, '21': 2, '22': 5, '18': 10,
                     '24': 20, '25': 50, '26': 100, '27': 200}

    for fid, temp in file_temp_map.items():
        Bvals, sxx, sxy = load_from_file(int(fid))

        # Eliminate the lowest field values where they are straight off.
        # This is due to the discretisation of the Fermi surface,
        # and in the low B limit the chance to hop becomes near 0.
        if temp == 200:
            wh = Bvals > Bvals[6]
            Bvals, sxx, sxy = Bvals[wh], sxx[wh], sxy[wh]
        if temp == 100:
            wh = Bvals > Bvals[2]
            Bvals, sxx, sxy = Bvals[wh], sxx[wh], sxy[wh]
        if temp == 50:
            wh = Bvals > Bvals[1]
            Bvals, sxx, sxy = Bvals[wh], sxx[wh], sxy[wh]

        rxx = sxx / (sxx**2 + sxy**2) * 1e8
        dB, drxx = deriv(Bvals, rxx, 1)
        dB = np.array([0] + list(dB))
        drxx = np.array([0] + list(drxx))
        c = _get_color_scaling(temp)
        ax.plot(dB / temp, drxx, lw=LW_DATA, color=c, label=f'{temp:.0f} K')
    # ax.legend()


def plots_paper():
    """ """

    fig, subplots = plt.subplots(ncols=2, nrows=3, figsize=(12, 19))

    _plot_tau(subplots[0][1], 10, 1e-16)
    _plot_tau(subplots[0][0], 10, 1e-18)
    _plot_half_data(subplots[1][0], subplots[2][0], 18, True)
    _plot_half_data(subplots[1][1], subplots[2][1], 19, False)
    _plot_scaling(subplots[2][0])
    _plot_breakdown_fit(subplots[2][1], 19)
    _plot_quadrature_fit(subplots[2][0], 18, 10)

    fig.subplots_adjust(right=0.98, left=0.15, bottom=0.05, top=0.98,
                        wspace=0.4, hspace=0.22)
    subplots[1][0].legend(loc='upper left', frameon=False)
    subplots[0][0].annotate('a)', xy=[0.02, 0.98], xycoords='figure fraction')
    subplots[0][0].annotate('d)', xy=[0.52, 0.98], xycoords='figure fraction')
    subplots[0][0].annotate('b)', xy=[0.02, 0.65], xycoords='figure fraction')
    subplots[0][0].annotate('e)', xy=[0.52, 0.65], xycoords='figure fraction')
    subplots[0][0].annotate('c)', xy=[0.02, 0.33], xycoords='figure fraction')
    subplots[0][0].annotate('f)', xy=[0.52, 0.33], xycoords='figure fraction')

    # Save
    here = os.path.dirname(__file__)
    path = os.path.join(here, 'lorentzian.png')
    fig.savefig(path, dpi=300)



# Mode 1: just breakdown
# sxx, sxy = simple_lorentz(Bvals, only_breakdown=True)


# Mode 2: computation
# This generates the files that are then inserted in mode 3 to show.
# Bvals = np.linspace(0, 150, 301)
# for T in [0.01]:
#     sxx, sxy = simple_lorentz(Bvals, T, 1e-18)
#     fileid = store_in_file(Bvals, sxx, sxy)
#     print(f'Stored the results in id {fileid} for T = {T} K')
#     Bvals, sxx, sxy = load_from_file(fileid)
#     show_resistivity(Bvals, sxx, sxy)
#     plt.show()

# Mode 3: paper results
# One run at 10 K and with 1e-16 hot spot, Bvals linspace(0 1500 3001)
# One run at 10 K and with 1e-18 hot spot, Bvals linspace(0 150 301)
plots_paper()

# Mode 4: Drude
run_Drude_test()
