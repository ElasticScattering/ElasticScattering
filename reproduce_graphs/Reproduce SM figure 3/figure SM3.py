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

import os
import numpy as np
import matplotlib.pyplot as plt


def raw_deriv(X, rho):
    """ Create a derivative between neighbouring points. """

    dX = 0.5 * (X[1:] + X[:-1])
    drhodX = (rho[1:] - rho[:-1]) / (X[1:] - X[:-1])
    return dX, drhodX


wct = np.linspace(0.003, 2 * np.pi, 1000)
D = np.pi / (4 * wct)
MR = D * (2 * D * np.cosh(2 * D) - np.sinh(2 * D))
MR /= (2 * D ** 2 + 1) * np.cosh(2 * D) - 2 * D * np.sinh(2 * D) - 1
MR -= 1

dX, dY = raw_deriv(wct, MR)

plt.rc('font', size=30)
fig, (ax1, ax2) = plt.subplots(figsize=(12, 15), nrows=2)
fig.subplots_adjust(bottom=0.07, top=0.92, right=0.98, left=0.15, hspace=0.03)

ax1.plot(wct, MR, lw=3, color='black')
ax1.set_ylabel('$\Delta\u03C1/\u03C1_{xx,0}$')
ax1.set_xticklabels([])
ax1.set_xlim(0, np.pi * 3 / 4)
ax1.set_ylim(bottom=0)
ax1.plot([np.pi / 2] * 2, [-1, 1000], color='tab:green',
         lw=3, label='scale', dashes=[6, 2])
ax1.annotate('scale', [1.2, 0.05], color='tab:green')
ax1.annotate('$B\u279D\u221e$', [0.1, 0.29], color='tab:green')
ax1.plot([0, 1000], [1 / 3] * 2, color='tab:green', lw=3, dashes=[6, 2])
axY = ax1.twiny()
axY.set_xticks([0, np.pi / 4, np.pi / 2, np.pi * 3 / 4])
axY.set_xticklabels(['0', '1', '2', '3'])
axY.set_xlabel('1/D = 4$\gamma/\pi$')

ax2.plot(dX, dY, lw=3, color='black')
ax2.set_ylabel('$d(\Delta\u03C1/\u03C1_{xx,0})/d\omega_c\u03C4$ (1/rad)')
ax2.set_xlabel('$\omega_c\u03C4$ (rad)')
ax2.set_xlim(0, np.pi * 3 / 4)
ax2.set_ylim(bottom=0)
ax2.plot([np.pi / 2] * 2, [-1, 1000], color='tab:green',
         lw=3, label='scale', dashes=[6, 2])
ax2.annotate('scale', [1.2, 0.5], color='tab:green')
axY = ax2.twiny()
axY.set_xticks([0, np.pi / 4, np.pi / 2, np.pi * 3 / 4])
axY.set_xticklabels([])

here = os.path.dirname(__file__)
path = os.path.join(here, 'last_code_execution.png')
fig.savefig(path, dpi=300)
