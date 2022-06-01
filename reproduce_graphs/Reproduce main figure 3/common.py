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

# Containing the actual code to do Boltzmann theoryimport numpy as np
import numpy as np
import scipy.optimize as so
import scipy.integrate as si
import numba

M0 = 9.10938188e-31
E = 1.602176634e-19
HBAR = 1.05457148e-34
KB = 1.380649e-23

# SETINGS MATCH 1fluid!!!!
# Namely 1e-9 c-axis and p=0.2 filling for n=8e27
#
# For Tl2201:
#       5 m0 [see Culo paper for sources]
#       k00 = 0.728e10 [quadrature paper for sources]
#       a = 5.45 / c = 23.1 A [Peets2018 arxiv:0910.3457]
M = 5 * M0
K = 7125232706
C = 1e-9  # 23.1e-10 / 2
A = 5.45e-10 / np.sqrt(2)
n = K**2 / C / 2 / np.pi


@numba.njit
def omegac(b):
    return E * b / M


@numba.njit
def vf():
    return HBAR * K / M


@numba.njit
def vx(phi):
    return vf() * np.cos(phi)


@numba.njit
def vy(phi):
    return vf() * np.sin(phi)


@numba.njit
def odds(t, tau):
    """ Probability to continue surviving until time t """
    return np.exp(-t / tau)


@numba.njit
def const():
    """ Constant for sigma when doing the integral. """
    jacobian = K  # phi integral
    dfde = 1 / (HBAR * vf())  # dk/dE really.
    theory = E**2 / (4 * np.pi**3)
    zdirec = 2 * np.pi / C
    units = 1e-8  # (Omega-m)^-1 --> (mu-Omega-cm)^-1
    return jacobian * dfde * theory * zdirec * units


@numba.njit
def integrandxx(t, phi, B, tau):
    return vx(phi) * vx(phi - omegac(B) * t) * odds(t, tau)


@numba.njit
def integrandxy(t, phi, B, tau):
    return vx(phi) * vy(phi - omegac(B) * t) * odds(t, tau)


@numba.njit
def integrandyy(t, phi, B, tau):
    return vy(phi) * vy(phi - omegac(B) * t) * odds(t, tau)


def getintegrand(component):
    if component == 'xx':
        return integrandxx
    if component == 'xy' or component == 'yx':
        return integrandxy
    if component == 'yy':
        return integrandyy
    raise ValueError(f'Unknown conductivity component {component}.')

#################
# Boltzmann
# Quadrature
#################


def quadrature(h, a, g, r0, T):
    return r0 + np.sqrt((a * T)**2 + (g * h)**2)


def dquadraturedh(h, a, g, T):
    return g**2 * h / np.sqrt((a * T)**2 + (g * h)**2)


def fit_quadrature(hh, rho, temperature):
    """ Get the fit parameters and their errors. """
    p, pcov = so.curve_fit(quadrature, hh, rho,
                           p0=[1, 1, rho[0]],
                           args=(temperature,))
    perr = np.sqrt(np.trace(pcov))
    return p, perr


def derivative(hh, rho):
    """ Get a new hh and drho/dh """
    return (rho[1:] - rho[:-1]) / (hh[1:] - hh[:-1])


def fit_dquadrature(hh, drhodh, temperature, *, xeval=[]):
    """ Get the fit parameters. """

    def fitfunc(h, a, g):
        return dquadraturedh(h, a, g, temperature)

    p, pcov = so.curve_fit(fitfunc, hh, drhodh,
                           p0=[1, 1])
    perr = np.array([np.sqrt(row[i]) for i, row in enumerate(pcov)])

    if len(xeval) == 0:
        return p, perr
    yeval = fitfunc(xeval, *p)
    return p, perr, yeval


#################
# Quadrature
# 2-fluid
#################


def _timeupper(angle1, angle2,
               tauCoh, tauInc, h, penetration, isUnbounded):
    """ Return an object to pass to si.quad for
    the upper bound of the integral.

    Physics:
    Incoherent particles are bounded,
    coherent particles are unbounded.
    i.e. this is "isCoherent"
    """

    max_tau_multiple = 20
    assert(angle2 > angle1)

    # Inserted rather than np.infty for
    # numerical stability as well as speed.
    # Overall an infinite integral asks for
    # a truncation point, this does that.
    if isUnbounded:
        t_infty = max_tau_multiple * tauCoh
    else:
        t_infty = max_tau_multiple * tauInc

    ############

    if isUnbounded:
        return t_infty

    ############
    # Either scattering or the end of the zone
    def timeupper(angle):
        wc = omegac(h)
        if wc == 0:
            return t_infty

        # Penetration adds to dangle in absolute value
        if wc > 0:
            dangle = angle - angle1 + penetration
        else:
            dangle = angle - angle2 - penetration
        assert(dangle / wc > 0)
        return min(t_infty, dangle / wc)

    return timeupper

def integrate(angle1, angle2, isCoherent, component,
              tauCoh, tauInc, h, penetration):

    timeupper = _timeupper(angle1, angle2, tauCoh, tauInc,
                           h, penetration, isCoherent)
    integrand = getintegrand(component)
    res, _ = si.dblquad(integrand, angle1, angle2, 0, timeupper,
                        args=(h, tauCoh if isCoherent else tauInc))
    return res


def sigma(h, tauCoh, tauInc, component, fractionInc, penetration, onlyInc=False):
    """ Returns the total and incoherent sigma. """

    assert(0 < fractionInc < 1)

    # X is half the angle the incoherent sector spans accross, radians.
    # The incoherent sectors are at [N pi/2-X, N pi/2+X]
    x = np.pi / 4 * fractionInc
    assert(x > 0)

    kw = {'component': component,
          'h': h, 'penetration': penetration,
          'tauCoh': tauCoh, 'tauInc': tauInc}

    res1 = integrate(-x, x, False, **kw)
    if not onlyInc:
        res2 = integrate(x, np.pi / 2 - x, True, **kw)
    res3 = integrate(np.pi / 2 - x, np.pi / 2 + x, False, **kw)
    if not onlyInc:
        res4 = integrate(np.pi / 2 + x, np.pi - x, True, **kw)

    # Factor 2 for [pi-x, 2pi-x] interval of the integral
    # which is by symmetry identical to [-x, pi-x]
    inc = 2 * const() * (res1 + res3)
    if not onlyInc:
        coh = 2 * const() * (res2 + res4)
        return coh + inc, coh, inc
    else:
        return inc


#################
# 2-fluid
# tests
#################


# Basically coherent
_vTot, _vCoh, _vInc = sigma(0.01, 1.5e-12, 1.5e-13, 'xx', 0.0001, 0)
_drude = n * E**2 * 1.5e-12 / M / 1e8
assert(_vCoh / _vInc > 100)
assert(abs(_vTot - _drude) < abs(_drude) * 1e-3)

# Basically incoherent
_tau = HBAR / KB / 20
_vTot, _vCoh, _vInc = sigma(0.01, 10 * _tau, _tau, 'yy', 0.9999, 0)
_drude = n * E**2 * _tau / M / 1e8
assert(_vInc / _vCoh > 100)
assert(abs(_vTot - _drude) < abs(_drude) * 1e-3)

print('> Unit tests pasts')
print(f'Carrier density (both spins) is ~{n} m-3')
print(f'Equivalent to p={n * A**2 * C}')
print(f'Fermi velocity vf={HBAR * K / M:.2e} m/s')
