import os
import numpy as np
import matplotlib.pyplot as plt


# Change this to show either the low residual or high residual
# results of figure 2 of the supplementary material
high_residual = False


###################


plt.rc('font', size=50)
maindir = os.path.join(os.path.dirname(__file__))
if high_residual:
    maindir = os.path.join(maindir, 'high residual')
else:
    maindir = os.path.join(maindir, 'low residual')
assert(os.path.isdir(maindir))


def get_resfiles(maindir):
    """ Get the Tx files available. """

    paths = []
    for obj in os.listdir(maindir):
        if obj.startswith('T') and obj.endswith('.dat'):
            paths.append(os.path.join(maindir, obj))

    return paths

def read_T0_file(path):

    labels = []
    data = []
    temperature = -1
    with open(path) as f:
        for line in f:
            if 'Temperature' in line:
                temperature = float(line.split()[-1])
            if line.startswith('B') and 'sigma' in line:
                break
        else:
            raise ValueError(f'No labels found in file {path}')
        if temperature < 0:
            raise ValueError(f'No temperature in file {path}')
        labels = line.split()

        for line in f:
            if line:
                data.append([float(number) for number in line.split()])
                assert(len(data[-1]) == len(labels))
    return temperature, labels, np.array(data)

def download(paths):

    temps = []
    magns = []
    sxxcs = []
    sxxis = []
    sxycs = []
    sxyis = []
    for path in paths:
        temp, labels, data = read_T0_file(path)
        temps += [temp] * len(data)
        magns += [col[labels.index('B')] for col in data]
        sxxcs += [col[labels.index('sigma_xx_coh')] for col in data]
        sxycs += [col[labels.index('sigma_xy_coh')] for col in data]
        sxxis += [col[labels.index('sigma_xx_inc')] for col in data]
        sxyis += [col[labels.index('sigma_xy_inc')] for col in data]

    result = {}
    result['T'] = np.array(temps)
    result['B'] = np.array(magns)
    result['sxxc'] = np.array(sxxcs)
    result['sxyc'] = np.array(sxycs)
    result['sxxi'] = np.array(sxxis)
    result['sxyi'] = np.array(sxyis)
    result['sxx'] = result['sxxc'] + result['sxxi']
    result['sxy'] = result['sxyc'] + result['sxyi']
    result['rxx'] = result['sxx'] / (result['sxx']**2 + result['sxy']**2)
    result['rxxi'] = result['sxxi'] / (result['sxxi']**2 + result['sxyi']**2)
    result['rxy'] = result['sxy'] / (result['sxx']**2 + result['sxy']**2)
    result['dB'] = 0.5 * (result['B'][1:] - result['B'][:-1])
    result['drxx'] = (result['rxx'][1:] - result['rxx'][:-1]
                      ) / (result['B'][1:] - result['B'][:-1])
    result['drude'] = result['sxxc'] / (result['sxxc']**2 + result['sxyc']**2)
    return result


def _derivative(X, rho, distance):
    """ Create a derivative over distance points to either side. """

    N = len(X)
    dX = X[distance:-distance]
    dR = np.zeros(len(dX))
    for i in range(len(dX)):
        stop = i + 2 * distance + 1
        assert(stop <= N)
        xlocal = X[i:stop]
        resss = rho[i:stop]
        p = np.polyfit(xlocal, resss, 1)
        dR[i] = p[0]
    assert(dX[-1] != 0)
    return dX, dR


def derivative(data, temperature):
    """ Get dH and dRxx for this temperature """

    Bvals = data['B'][data['T'] == temperature]
    Rvals = data['rxx'][data['T'] == temperature]
    dB, dR = _derivative(Bvals, Rvals, 2)
    # dB = 0.5 * (Bvals[1:] + Bvals[:-1])
    # dR = (Rvals[1:] - Rvals[:-1]) / (Bvals[1:] - Bvals[:-1])
    return dB, dR


def get_color(T, temperatures):
    """ Get the color of this temperature given the bunch. """
    assert(T in temperatures)
    index = list(temperatures).index(T)
    colors = plt.get_cmap('cividis', len(temperatures) + 1).colors[:-1]
    return colors[index]

#####################


def show_fieldsweeps(data, temperatures, *, ax=None):
    """ For these temperatures, make a (B,rho) plot. """

    if ax is None:
        _, ax = plt.subplots(num='Resistivity')

    for T in temperatures:
        wh = data['T'] == T
        label = f'{T:.1f} K' if T < 10 else f'{T:.0f} K'
        c = get_color(T, temperatures)
        ax.plot(data['B'][wh], data['rxx'][wh],
                label=label, lw=6, color=c)
    ax.set_xlabel('$B$ (T)')
    ax.set_ylabel('$\u03C1_{xx}$ (\u03BC\u03A9cm)')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    # plt.legend(bbox_to_anchor=[1, 1], loc='upper left')
    plt.gcf().subplots_adjust(bottom=0.2, left=0.2)


def show_resistivity(data, temperatures, *, ax=None):
    """ Show the T-linearity in 1 T field. Returns the fit. """

    if ax is None:
        _, ax = plt.subplots(num='Zero field')

    resistivity = []
    for T in temperatures:
        wh = data['T'] == T
        ind = np.argmin(np.abs(data['B'][wh] - 1))
        resistivity.append(data['rxx'][wh][ind])
    resistivity = np.array(resistivity)

    ax.semilogx(temperatures, resistivity,
                ms=10, color='red', marker='o', lw=0)

    p = np.polyfit(temperatures, resistivity, deg=1)
    xx = np.logspace(np.log10(min(temperatures)) * 1.2,
                     np.log10(max(temperatures)) * 1.1, 1000)
    yy = p[0] * xx + p[1]
    ax.semilogx(xx, yy, color='black', lw=3)
    ax.set_xlabel('$T$ (K)')
    ax.set_ylabel('$\u03C1_{xx,0}$ (\u03BC\u03A9cm)')

    T0 = p[1] / p[0]
    print(f'The residual resitsivity is equivalent to T0={T0:.1f} K')
    return p


def show_kohler(data, temperatures, *, ax=None):
    """ Make the scaling (H/rho0, DeltaRho/rho0) plot """

    legend = ax is None
    if ax is None:
        _, ax = plt.subplots(num='Kohler')

    for T in temperatures:
        wh = data['T'] == T
        label = f'{T:.1f} K' if T < 10 else f'{T:.0f} K'
        ind = np.argmin(np.abs(data['B'][wh] - 1))
        rho0 = data['rxx'][wh][ind]
        c = get_color(T, temperatures)

        ax.plot(data['B'][wh] / rho0, (data['rxx'][wh] - rho0) / rho0,
                label=label, lw=6, color=c)
    ax.set_xlabel('$B/\u03C1_{xx,0}$ (T/\u03BC\u03A9cm)')
    ax.set_ylabel('$\Delta\u03C1/\u03C1_{xx,0}$')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    plt.gcf().subplots_adjust(right=0.8)
    if legend:
        ax.legend(bbox_to_anchor=[1, 1], loc='upper left')


def show_mod_kohler(data, temperatures, zero_field_fit, *, ax=None):
    """ Make the failing scaling of kohler, but rho0 -> rho0-rho_res """

    if ax is None:
        _, ax = plt.subplots(num='Mod Kohler')

    for T in temperatures:
        wh = data['T'] == T
        label = f'{T:.1f} K' if T < 10 else f'{T:.0f} K'
        ind = np.argmin(np.abs(data['B'][wh] - 1))
        rho0 = data['rxx'][wh][ind]
        delta = rho0 - zero_field_fit[1]
        c = get_color(T, temperatures)

        ax.plot(data['B'][wh] / delta, (data['rxx'][wh] - rho0) / delta,
                label=label, lw=6, color=c)
    ax.set_xlabel('$B/(\u03C1_{xx,0}-\u03C1_{xx,res})$ (T/\u03BC\u03A9cm)')
    ax.set_ylabel('$\Delta\u03C1/(\u03C1_{xx,0}-\u03C1_{xx,res})$')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)


def show_sigma_drift(data, temperatures, *, ax=None):
    """ Tau is not constant with field, indicated by changing sigma. """

    # It does not matter which temperature, since coherent has
    # tau=infty always and no temperature dependence.
    #
    # The coherent sector is tiny, so the conductivity is tiny,
    # but its field dependence gives an idea of the changing
    # tau since rho_xx_drude has zero MR and is equal to
    # rho(B) = rho(0) = m*/ne^e tau.
    norm = data['sxxc'][0]

    if ax is None:
        _, ax = plt.subplots(num='Resistivity drift')

    for T in temperatures:
        wh = data['T'] == T
        label = f'{T:.1f} K' if T < 10 else f'{T:.0f} K'
        c = get_color(T, temperatures)
        rhoc = data['sxxc'][wh] / (data['sxxc'][wh]**2 + data['sxyc'][wh]**2)
        ax.plot(data['B'][wh], rhoc * norm,
                label=label, lw=6, color=c)
    ax.set_xlabel('$B$ (T)')
    ax.set_ylabel('$\u03C1_{xx,c}$ (\u03BC\u03A9cm)')
    ax.legend(bbox_to_anchor=[1, 1], loc='upper left')
    plt.gcf().subplots_adjust(right=0.75)
    ax.set_xlim(left=0)
    ax.set_title('Should be constant in Drude theory')


if __name__ == '__main__':
    paths = get_resfiles(maindir)
    data = download(paths)
    temperatures = np.unique(data['T'])

    show_fieldsweeps(data, temperatures)
    p = show_resistivity(data, temperatures)
    show_kohler(data, temperatures)
    show_mod_kohler(data, temperatures, p)
    show_sigma_drift(data, temperatures)
    plt.show()
