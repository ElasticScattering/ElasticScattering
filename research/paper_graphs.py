# Whereas view_results is aimed at individual results
# as they are generated,
# this code makes the graphs for the paper from specific results.

import os
import numpy as np
import matplotlib.pyplot as plt
import view_result as vr


plt.rc('font', size=30)
fig, subplots = plt.subplots(nrows=3, ncols=2, figsize=(12, 26), num='paper')

# High RR
paths = vr.get_resfiles(vr.maindir, 42)
data = vr.download(paths)
temperatures = np.unique(data['T'])
vr.show_fieldsweeps(data, temperatures, ax=subplots[0][0])
p = vr.show_resistivity(data, temperatures)
vr.show_kohler(data, temperatures, ax=subplots[1][0])
vr.show_mod_kohler(data, temperatures, p, ax=subplots[2][0])
vr.show_sigma_drift(data, temperatures)

# Low RR
paths = vr.get_resfiles(vr.maindir, 43)
data = vr.download(paths)
temperatures = np.unique(data['T'])
vr.show_fieldsweeps(data, temperatures, ax=subplots[0][1])
p = vr.show_resistivity(data, temperatures)
vr.show_kohler(data, temperatures, ax=subplots[1][1])
vr.show_mod_kohler(data, temperatures, p, ax=subplots[2][1])
vr.show_sigma_drift(data, temperatures)

# Custom adjustments
subplots[0][0].set_ylim(0, 160)
subplots[0][1].set_ylim(0, 160)
fig.subplots_adjust(bottom=0.12, left=0.13, right=0.935, top=0.98,
                    wspace=0.45, hspace=0.17)
subplots[2][0].legend(bbox_to_anchor=[1.1, -0.3],
                      loc='center', ncol=4, frameon=False)
subplots[1][1].set_xticks([0, 5, 10, 15])
X = 0.01
Y = 0.02
plt.figure('paper')
plt.annotate('a)', xy=[X, 1 - Y], xycoords='figure fraction')
plt.annotate('b)', xy=[X, 0.695 - Y], xycoords='figure fraction')
plt.annotate('c)', xy=[X, 0.39 - Y], xycoords='figure fraction')
plt.annotate('d)', xy=[0.48 + X, 1 - Y], xycoords='figure fraction')
plt.annotate('e)', xy=[0.48 + X, 0.695 - Y], xycoords='figure fraction')
plt.annotate('f)', xy=[0.48 + X, 0.39 - Y], xycoords='figure fraction')

# Store
here = os.path.dirname(__file__)
path = os.path.join(here, 'disorder.png')
fig.savefig(path, dpi=300)
# plt.show()
