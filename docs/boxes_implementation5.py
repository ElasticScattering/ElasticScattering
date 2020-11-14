# Do array evaluation of impurity intersects.
# Work towards sigma.

import time
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import namedtuple
from numba import njit

HBAR = 1.05457148e-34
KB = 1.380649e-23
E = 1.60217646e-19

region = 1e-6
extends = 3e-6
imp_rad = 1e-8
imp_nr = 10000
isEl = True
avg_box = 10
mass = 5 * 9.109e-31
vf = 1.68e5
alpha = 0.5
dimension = 15
nr_phi = 13

COUNTER = 0

##############################
# Existing functionality
##############################

def generate_impurities(number, low, high):
    """ Get (x,y) pairs for all the impurity locations. """

    np.random.seed(123875)
    impurities = np.zeros((number, 2))
    for i in range(number):
        impurities[i, 0] = np.random.uniform(low, high)
        impurities[i, 1] = np.random.uniform(low, high)

    return impurities


@njit()
def angle_position(angle_vel, isEl):
    """ This is a very easily overlooked and mixed up angle transformation.

    RULE: to prevent confusion about what angle is what, ALWAYS hold on
        to angle_vel and isEl and only derive angle_pos when you need
        it to perform maths.

    The incoming particle is [x, y, phi] where phi is the velocity angle.
    But for box motion we are not interested in velocity but position angle.

    If we have an electron, which moves ccw with time and cw to the past,
      then phi_pos = phi_vel - pi/2
    If we have a hole, which moves cw with time and ccw to the past,
      then phi_pos = phi_vel + pi/2
    """
    if isEl:
        return (angle_vel - np.pi / 2) % (2 * np.pi)
    else:
        return (angle_vel + np.pi / 2) % (2 * np.pi)


@njit()
def angle_difference(phi1, phi2, isEl):
    """ How many radians from phi1 towards phi2?

    Simple, but with 2pi modularity and clockwise, counter clockwise,
    it can be quite confusing.

    Assume motion to the past, so electrons are clockwise
    """
    if isEl:
        r = (phi1 - phi2) % (2 * np.pi)
    else:
        r = (phi2 - phi1) % (2 * np.pi)

    if r > 2 * np.pi - 1e-12:
        return 0
    return r


@njit()
def angle_in_range(phi1, phirange, isEl):
    """ Is phi1 in the interval given?

    This one is tricky because you have to consider to go
    cw or ccw and keep in mind the modulo.

    As usual, assume we are looking towards the past
    phirange is [from, to]

    NOTE: all input is assumed to be [0,2pi]
    """

    # Consider this to be full circle, not empty circle.
    if abs(phirange[0] - phirange[1]) < 1e-10:
        return True

    start, end = phirange

    # El cw so go *down*
    if isEl:
        low = end
        high = start
    # Hole then ccw so go *up*
    if not isEl:
        low = start
        high = end

    # Detect modulo triggers
    if high < low:
        high += 2 * np.pi

    if low < phi1 < high:
        return True
    if low < phi1 + 2 * np.pi < high:
        return True
    return False


def test_angle_in_range():
    """ Trigger all the if statements in angle_in_range. """

    # Direct inequality
    assert(angle_in_range(1, np.array([0.5, 1.5]), False))

    # Move 0.5 through 0 down to 1.5
    assert(not angle_in_range(1, np.array([0.5, 1.5]), True))
    assert(angle_in_range(2, np.array([0.5, 1.5]), True))
    assert(angle_in_range(0.1, np.array([0.5, 1.5]), True))

    # Go up from 1.5 through 2pi to 0.5
    assert(not angle_in_range(1, np.array([1.5, 0.5]), False))

    # Go down from 1.5 to 0.5, direct inequality
    assert(angle_in_range(1, np.array([1.5, 0.5]), True))

    # Full circle
    assert(angle_in_range(1, np.array([1, 1]), True))
    assert(angle_in_range(2, np.array([1, 1]), True))


test_angle_in_range()


@njit()
def angle_velocity(angle_pos, isEl):
    """ This is rare functionality, see angle_pos for the main function """
    if isEl:
        return (angle_pos + np.pi / 2) % (2 * np.pi)
    else:
        return (angle_pos - np.pi / 2) % (2 * np.pi)


def cyclotron_orbit(x, y, vx, vy, w, isEl):
    """ Compute the associated cyclotron radius and center coordinates
    of this starting position.

    isEl == isClockwise to future
    """

    if w == 0:
        raise ValueError('No cyclotron orbit without magnetic field.')
    elif w < 0:
        raise ValueError('wc has to be positive.')

    vf = np.sqrt(vx**2 + vy**2)
    radius = vf / w
    xshift = radius * vy / vf
    yshift = -radius * vx / vf

    if isEl:
        centerX = x - xshift
        centerY = y - yshift
    else:
        centerX = x + xshift
        centerY = y + yshift

    return radius, centerX, centerY


def plot_situation(axis, impurities, radius, region, extends):
    """ Plot the impurities and the region box. """

    axis.set_xlim(-extends, region + extends)
    axis.set_ylim(-extends, region + extends)

    for imp in impurities:
        circle = mpatches.Circle(imp, radius=radius)
        axis.add_patch(circle)


def plot_region(axis, region):
    axis.plot([0, region], [0, 0], lw=2, color='black')
    axis.plot([0, 0], [0, region], lw=2, color='black')
    axis.plot([0, region], [region, region], lw=2, color='black')
    axis.plot([region, region], [0, region], lw=2, color='black')


def plot_cyclotron(axis, particle, vf, wc, isEl):
    """ Plot the cyclotron orbit of this particle. """

    vx = vf * np.cos(particle[2])
    vy = vf * np.sin(particle[2])
    rad, xc, yc = cyclotron_orbit(particle[0], particle[1], vx, vy, wc, isEl)

    theta = np.linspace(0, 2 * np.pi, 500)
    xx = xc + rad * np.cos(theta)
    yy = yc + rad * np.sin(theta)
    axis.plot(xx, yy, color='red', lw=1)
    axis.scatter(particle[0], particle[1], s=500, color='black')

    # dx = rad / 3 * np.cos(particle[2])
    # dy = rad / 3 * np.sin(particle[2])
    # axis.arrow(particle[0], particle[1], dx, dy, color='red',
    #            width=1e-9, head_width=1e-9, head_length=7e-10, lw=3)


##############################
# Griddivision
##############################

# Note: After determining it, I take nr1D to be the defining characteristic,
#       the other choice is to take boxdim to be leading, but it has to exactly
#       match the region size at a multiple - and that coincidence seems more
#       natural by taking the integer than having perfectly divisible floats.

@njit()
def get_boxdim(region, extends, nr1D):
    return (region + 2 * extends) / nr1D


@njit()
def get_xind(x, region, extends, nr1D):
    return (x + extends) // get_boxdim(region, extends, nr1D)


@njit()
def get_yind(y, region, extends, nr1D):
    return (y + extends) // get_boxdim(region, extends, nr1D)


@njit()
def get_boxid(xind, yind, nr1D):
    assert(xind >= 0)
    assert(xind < nr1D)
    assert(yind >= 0)
    assert(yind < nr1D)
    return int(yind * nr1D + xind)


@njit()
def get_boxid_full(x, y, region, extends, nr1D):
    xind = get_xind(x, region, extends, nr1D)
    yind = get_yind(y, region, extends, nr1D)
    return get_boxid(xind, yind, nr1D)


@njit()
def get_xlow(xind, region, extends, nr1D):
    return -extends + (region + 2 * extends) * xind / nr1D


@njit()
def get_ylow(yind, region, extends, nr1D):
    return -extends + (region + 2 * extends) * yind / nr1D


@njit()
def get_xind_id(id, nr1D):
    assert(0 <= id < nr1D**2)
    return id % nr1D


@njit()
def get_yind_id(id, nr1D):
    assert(0 <= id < nr1D**2)
    return id // nr1D


def grid_divide(impurities, rad, region, extends, avg_per_box):
    """ From [n, 2] to [B, n, 2] where B is the box and n~avg """

    aim_nr1D = np.sqrt(len(impurities) / avg_per_box)
    nr1D = int(np.ceil(aim_nr1D))
    assert(get_boxdim(region, extends, nr1D) > rad)

    boxdim = get_boxdim(region, extends, nr1D)
    bins = [[] for _ in range(nr1D**2)]
    for impurity in impurities:

        # Add impurity to the box it lives in
        xind = get_xind(impurity[0], region, extends, nr1D)
        yind = get_yind(impurity[1], region, extends, nr1D)

        box_id = get_boxid(xind, yind, nr1D)
        bins[box_id].append(impurity)

        dx = (impurity[0] + extends) % boxdim
        dy = (impurity[1] + extends) % boxdim

        # And figure out if it also lives in a neighbour.
        # Horizontal
        # Box to left
        if dx < rad and xind > 0:
            bins[get_boxid(xind - 1, yind, nr1D)].append(impurity)

        # Box to bottom
        if dy < rad and yind > 0:
            bins[get_boxid(xind, yind - 1, nr1D)].append(impurity)

        # Box to right
        if dx > boxdim - rad and xind < nr1D - 1:
            bins[get_boxid(xind + 1, yind, nr1D)].append(impurity)

        # Box to top
        if dy > boxdim - rad and yind < nr1D - 1:
            bins[get_boxid(xind, yind + 1, nr1D)].append(impurity)

        # Diagonals
        # Topleft
        d_topleft2 = dx**2 + (boxdim - dy)**2
        if d_topleft2 < rad**2 / 2 and xind > 0 and yind < nr1D - 1:
            bins[get_boxid(xind - 1, yind + 1, nr1D)].append(impurity)

        d_topright2 = (boxdim - dx)**2 + (boxdim - dy)**2
        if d_topright2 < rad**2 / 2 and xind < nr1D - 1 and yind < nr1D - 1:
            bins[get_boxid(xind + 1, yind + 1, nr1D)].append(impurity)

        d_bottomleft2 = dx**2 + dy**2
        if d_bottomleft2 < rad**2 / 2 and xind > 0 and yind > 0:
            bins[get_boxid(xind - 1, yind - 1, nr1D)].append(impurity)

        d_bottomright2 = (boxdim - dx)**2 + dy**2
        if d_bottomright2 < rad**2 / 2 and xind < nr1D - 1 and yind > 0:
            bins[get_boxid(xind + 1, yind - 1, nr1D)].append(impurity)

    return bins, nr1D


def show_grid(axis, impurities3D, rad, region, extends):

    nr1D = int(np.sqrt(len(impurities3D)))
    boxdim = (region + 2 * extends) / nr1D

    all_colors = list(plt.get_cmap('tab20', 20).colors)
    colors = all_colors * (len(impurities3D) // 20) + \
        all_colors[:len(impurities3D) % 20]
    assert(len(colors) == len(impurities3D))
    for i, (c, box) in enumerate(zip(colors, impurities3D)):
        for imp in box:
            rad_now = rad  # * (1 - i / 10 / nr1D)
            circle = mpatches.Circle(
                imp, radius=rad_now, color=c, lw=0)  # , alpha=0.5)
            axis.add_patch(circle)

        x_id = i % nr1D
        y_id = i // nr1D
        xL = -extends + x_id * boxdim + boxdim * 0.01
        yL = -extends + y_id * boxdim + boxdim * 0.01
        L = boxdim * 0.98
        plt.plot([xL, xL], [yL, yL + L], color=c, lw=1)
        plt.plot([xL, xL + L], [yL, yL], color=c, lw=1)
        plt.plot([xL, xL + L],
                 [yL + L, yL + L], color=c, lw=1)
        plt.plot([xL + L, xL + L],
                 [yL, yL + L], color=c, lw=1)

    impurities = sum(impurities3D, [])
    impurities = sorted(impurities, key=lambda x: x[0] + x[1])
    impurities.append([1e99, 1e99])
    ind = 0
    while ind < len(impurities) - 1:
        imp = impurities[ind]
        multi = 1
        while abs(impurities[ind + 1][0] - imp[0]) < 1e-5 * (2 * extends + region) \
                and abs(impurities[ind + 1][1] - imp[1]) < 1e-5 * (2 * extends + region):
            multi += 1
            ind += 1

        # plt.annotate(f'{multi}', imp)
        ind += 1

##############################
# Gridmovement
##############################


@njit()
def boundary_intersects(x1, y1, x2, y2, xc, yc, rc, L):
    """ Find the intersections between the circle and the connection
    line between the two points. L is the distance between
    point 1 and point 2. See the document for the maths. """

    ux = (x2 - x1) / L
    uy = (y2 - y1) / L
    projection_distance = ux * (xc - x1) + uy * (yc - y1)
    x_proj = x1 + projection_distance * ux
    y_proj = y1 + projection_distance * uy

    Empty = np.ones(2) * 1e99
    if (x_proj - xc)**2 + (y_proj - yc)**2 >= rc**2:
        return [Empty, Empty]

    dist_proj_c_squared = (x_proj - xc)**2 + (y_proj - yc)**2
    dist_proj_intersect = np.sqrt(rc**2 - dist_proj_c_squared)
    xa = x_proj + dist_proj_intersect * ux
    ya = y_proj + dist_proj_intersect * uy
    xb = x_proj - dist_proj_intersect * ux
    yb = y_proj - dist_proj_intersect * uy

    if abs(ux) > abs(uy):
        intersect_a_valid = ((x1 < xa) and (x2 > xa)) or \
            ((x2 < xa) and (x1 > xa))
        intersect_b_valid = ((x1 < xb) and (x2 > xb)) or \
            ((x2 < xb) and (x1 > xb))
    else:
        intersect_a_valid = ((y1 < ya) and (y2 > ya)) or \
            ((y2 < ya) and (y1 > ya))
        intersect_b_valid = ((y1 < yb) and (y2 > yb)) or \
            ((y2 < yb) and (y1 > yb))

    if intersect_a_valid:
        intersect_a = np.array([xa, ya])
    else:
        intersect_a = Empty
    if intersect_b_valid:
        intersect_b = np.array([xb, yb])
    else:
        intersect_b = Empty
    return [intersect_a, intersect_b]


def test_intersect():
    """ Unit tests created in the word doc before implementation """

    r = boundary_intersects(0, -3, 0, 3, 0, 0, 2, 6)
    assert(r[0][0] == r[1][0] == 0)
    assert(r[0][1] == 2 or r[1][0] == 2)
    assert(r[0][1] == -2 or r[1][1] == -2)

    r = boundary_intersects(0, -3, 0, 3, 0.2, 0.2, 20, 6)
    assert(r[0][0] > 1e98)
    assert(r[1][0] > 1e98)
    assert(r[0][1] > 1e98)
    assert(r[1][1] > 1e98)

    r = boundary_intersects(1, 0, 1, 6, 0.7, 0, 0.5, 6)
    assert(r[0][0] > 1e98 or r[0][0] == 1)
    assert(r[0][1] > 1e98 or abs(r[0][1] - 0.4) < 1e-10)
    assert(r[1][0] > 1e98 or r[1][0] == 1)
    assert(r[1][1] > 1e98 or abs(r[1][1] - 0.4) < 1e-10)

    r = boundary_intersects(5, -50, 5, 60, 0, 0, 1, 110)
    assert(r[0][0] > 1e98)
    assert(r[1][0] > 1e98)
    assert(r[0][1] > 1e98)
    assert(r[1][1] > 1e98)

    r = boundary_intersects(0, 0, 0, 2, 1, 1, 1.1, 2)
    A = 1 - np.sqrt(1.1**2 - 1)
    B = 1 + np.sqrt(1.1**2 - 1)
    assert(abs(r[0][1] - A) < 1e-8 or abs(r[0][1] - B) < 1e-8)
    assert(abs(r[1][1] - A) < 1e-8 or abs(r[1][1] - B) < 1e-8)
    assert(r[0][0] == 0)
    assert(r[1][0] == 0)


test_intersect()


@njit()
def get_angle(x, y, xc, yc, rc, isEl):
    """ Return the VELOCITY angle of a point on a circle.
    Assumes the point is on the circle. """

    if abs(y - yc) > 1e-10 * rc:
        # within numerical precision y==yc
        angle = np.arcsin((y - yc) / rc)
        if x < xc:
            angle = np.pi - angle
        return angle_velocity(angle % (2 * np.pi), isEl)
    elif x > xc:
        return angle_velocity(0, isEl)
    else:
        return angle_velocity(np.pi, isEl)


@njit()
def box_intersects(xlow, ylow, L, xc, yc, rc, isEl):
    """ Find the intersects of this particle’s orbit with the given box. """

    ints = np.zeros((8, 3))

    ints[0, :2], ints[1, :2] = \
        boundary_intersects(xlow, ylow, xlow + L, ylow, xc, yc, rc, L)
    ints[2, :2], ints[3, :2] = \
        boundary_intersects(xlow, ylow, xlow, ylow + L, xc, yc, rc, L)
    ints[4, :2], ints[5, :2] = \
        boundary_intersects(xlow, ylow + L, xlow + L, ylow + L, xc, yc, rc, L)
    ints[6, :2], ints[7, :2] = \
        boundary_intersects(xlow + L, ylow, xlow + L, ylow + L, xc, yc, rc, L)

    ints = [i for i in ints if i[0] < 1e98]

    for i in range(len(ints)):
        ints[i][2] = get_angle(ints[i][0], ints[i][1], xc, yc, rc, isEl)
    return ints


def test_angle():
    """ Tests created before function implementation

    Then adjusted to reflect velocity angle.
    This now tests both get_angle as well as
    the conversion between position and velocity.
    These tests are an authority, they are correct.

    As a starting point to see what is happening,
    realise that vel_phi=0 (moving to positive x)
    is at the bottom of an electron orbit where
    pos_phi=3pi/2.
    On the other hand, for holes vel_phi=0 happens
    at the top or pos_phi=pi/2.
    """

    assert(get_angle(1, 0, 0, 0, 1, True) == np.pi / 2)
    assert(get_angle(-1, 0, 0, 0, 1, True) == 3 * np.pi / 2)
    assert(get_angle(0, 1, 0, 0, 1, True) == np.pi)
    assert(get_angle(0, -1, 0, 0, 1, True) == 0)
    assert(get_angle(5, 5, 6, 6, np.sqrt(2), True) == 7 * np.pi / 4)

    assert(get_angle(1, 0, 0, 0, 1, False) == 3 * np.pi / 2)
    assert(get_angle(-1, 0, 0, 0, 1, False) == np.pi / 2)
    assert(get_angle(0, 1, 0, 0, 1, False) == 0)
    assert(get_angle(0, -1, 0, 0, 1, False) == np.pi)
    assert(get_angle(5, 5, 6, 6, np.sqrt(2), False) == 3 * np.pi / 4)


test_angle()


def find_exit_intersect(intersects, xc, yc, rc, phi, start_intersect, isEl):

    if not len(intersects):
        angle = angle_position(phi, isEl)
        x0 = xc + rc * np.cos(angle)
        y0 = yc + rc * np.sin(angle)
        return [-1, [x0, y0, phi]]

    # This has to stay, you can get completely wrong results otherwise.
    for intersect in intersects:
        if abs(intersect[2] - phi) < 1e-11:
            raise ValueError('Intersect and particle (grid) align, change the grid.')

    closest = -1
    dphi_closest = angle_difference(start_intersect[2], phi, isEl)
    nextone = -2
    dphi_nextone = 1e99

    for ind, intersect in enumerate(intersects):
        dphi_now = angle_difference(start_intersect[2], intersect[2], isEl)
        if dphi_now < dphi_closest:
            nextone = closest
            dphi_nextone = dphi_closest
            closest = ind
            dphi_closest = dphi_now
        elif dphi_now < dphi_nextone:
            nextone = ind
            dphi_nextone = dphi_now

    assert(nextone > -2)
    if nextone >= 0:
        return [nextone, intersects[nextone]]
    else:  # -1
        angle = angle_position(phi, isEl)
        x0 = xc + rc * np.cos(angle)
        y0 = yc + rc * np.sin(angle)
        return [-1, [x0, y0, phi]]


def find_exit_intersect_old(intersects, xc, yc, rc, phi, start_intersect, isEl):

    raise NotImplementedError()

    # DEPRECATED:
    #   This is just staying around for reference.
    #   This version works perfectly, but it fails when the circles
    #   are large enough that angle differences are less than
    #   1e-7 radians. I could have reduced this number, but that
    #   is just waiting for a finer mesh to violate it anyways.

    # If we are at the starting point, we have to find one of the
    # intersects to move towards - we cannot finish at phi.
    if abs(start_intersect[2] - phi) < 1e-7:
        angle_step = 1e99


    # If we are somewhere on the trajectory,
    # then one of the termination points is guaranteed to be the phi
    # or the starting point. Hence set this as the initial value
    # and then loop and take an impurity intersect whenever it is earlier.
    else:
        angle_step = angle_difference(start_intersect[2], phi, isEl)

    next_index = -1
    for ind, intersect in enumerate(intersects):
        angle_diff = angle_difference(start_intersect[2], intersect[2], isEl)

        # Of course, the distance from the current point to the
        # current point is 0, just like the starting point above,
        # this has to be excluded. We want the closest *next in line* impurity.
        if angle_diff < angle_step and angle_diff > 1e-7:
            angle_step = angle_diff
            next_index = ind

    if next_index == -1:
        angle = angle_position(phi, isEl)
        x0 = xc + rc * np.cos(angle)
        y0 = yc + rc * np.sin(angle)
        return [-1, [x0, y0, phi]]
    else:
        return [next_index, intersects[next_index]]


def test_exit():

    phi = 0
    radians = [0.4, 1, 4.7, 2.3, np.pi, 6]
    points_on_circle = [[np.cos(r), np.sin(r), r] for r in radians]

    # From east to northeast
    start_intersect = [np.cos(phi), np.sin(phi), phi]
    result = find_exit_intersect(
        points_on_circle, 0, 0, 1, phi, start_intersect, False)

    assert(result[0] == 0)
    assert(result[1][0] == points_on_circle[0][0])
    assert(result[1][1] == points_on_circle[0][1])
    assert(result[1][2] == 0.4)

    # From east to southeast
    result = find_exit_intersect(
        points_on_circle, 0, 0, 1, phi, start_intersect, True)
    assert(result[0] == 5)
    assert(result[1][0] == points_on_circle[5][0])
    assert(result[1][1] == points_on_circle[5][1])
    assert(result[1][2] == 6)

    # Follow this particle around
    result = find_exit_intersect(
        points_on_circle, 0, 0, 1, phi, result[1], True)
    assert(result[0] == 2)

    result = find_exit_intersect(
        points_on_circle, 0, 0, 1, phi, result[1], True)
    assert(result[0] == 4)

    result = find_exit_intersect(
        points_on_circle, 0, 0, 1, phi, result[1], True)
    assert(result[0] == 3)

    result = find_exit_intersect(
        points_on_circle, 0, 0, 1, phi, result[1], True)
    assert(result[0] == 1)

    result = find_exit_intersect(
        points_on_circle, 0, 0, 1, phi, result[1], True)
    assert(result[0] == 0)

    result = find_exit_intersect(
        points_on_circle, 0, 0, 1, phi, result[1], True)
    assert(result[0] == -1)
    assert(result[1][2] == 0)

    # Nothing to be found, go full circle
    result = find_exit_intersect([], 0, 0, 1, phi, start_intersect, True)
    assert(result[0] == -1)
    assert(result[1][2] == 0)


test_exit()


@njit()
def find_shared_side(point, id, region, extends, nr1D):
    """ Given a point on the border of box id, get the box on the other side. """

    # First we need to check if we are going out of bounds.
    boxdim = get_boxdim(region, extends, nr1D)
    delta = boxdim * 1e-6
    if point[0] < -extends + delta:
        return -1
    if point[0] > extends + region - delta:
        return -1
    if point[1] < -extends + delta:
        return -1
    if point[1] > extends + region - delta:
        return -1

    # Now that is out of the way, it is both:
    # 1) Guaranteed that there is a pair
    # 2) Safe to get the boxID
    xind = get_xind(point[0], region, extends, nr1D)
    yind = get_yind(point[1], region, extends, nr1D)
    boxid1 = get_boxid(xind, yind, nr1D)

    # Find the other one
    dx = (extends + point[0]) % boxdim
    dy = (extends + point[1]) % boxdim
    if dx < delta:
        boxid2 = get_boxid(xind - 1, yind, nr1D)
    elif dy < delta:
        boxid2 = get_boxid(xind, yind - 1, nr1D)
    elif dx > boxdim - delta:
        boxid2 = get_boxid(xind + 1, yind, nr1D)
    elif dy > boxdim - delta:
        boxid2 = get_boxid(xind, yind + 1, nr1D)
    else:
        raise ValueError('Not at a boundary')

    if id != boxid1 and id != boxid2:
        raise ValueError('Given value was not on the boundary or at a 4-point.')
    if id == boxid1:
        return boxid2
    return boxid1


def next_box(boxID, entry_intersect, xc, yc, rc, phi, isEl, region, extends, nr1D):
    """ Find the ID of the next box and the intersection where we cross over.

    Returns 2 values.
    The first is the ID of the next box in line.
    The second is a 3-point (x, y, phi) indicating the intersection point
        with the next box. This 3-point is the entryIntersect
        for the next iteration.
        NOTICE: PHI is at all times the velocity angle, not the position angle.
            See also: angle_position

    In case there is no box to move towards,
        the returned ID will be -1 to signal a termination of the recursion.
        This happens when we stop at the starting point OR if we leave the field.
        If we do leave the field, we do NOT re-enter, which is fine only because
        this scenario should only happen when out of bounds is at such long times
        that the particle is completely irrelevant due to tau.

    The input arguments speak for themselves except for entryIntersect.
        entryIntersect is the starting point of the current piece of
        the path we are going to investigate for impurities.
        Usually, this is the point where we enter the box.
        The only exception is the first iteration,
        in which case it is [x, y, phi] of the particle’s starting position.
        Yes, this means phi exists twice in the arguments,
        but that is how we know this is the first iteration.

    For impurity intersects, the interval of validity on angle is
        [entryIntersect.phi, exitIntersect.phi]
        that holds for the first iteration as well as for the final
        one when nextID==-1. Remember, if -1 is returned you still #
        have to still search the final interval once more!
    """

    xind = get_xind_id(boxID, nr1D)
    xlow = get_xlow(xind, region, extends, nr1D)
    yind = get_yind_id(boxID, nr1D)
    ylow = get_ylow(yind, region, extends, nr1D)

    L = get_boxdim(region, extends, nr1D)
    intersects = box_intersects(xlow, ylow, L, xc, yc, rc, isEl)


    nextID, exit_intersect = find_exit_intersect(
        intersects, xc, yc, rc, phi, entry_intersect, isEl)

    if nextID == -1:
        return [nextID, exit_intersect]
    else:
        nextBoxID = find_shared_side(
            exit_intersect, boxID, region, extends, nr1D)
        return [nextBoxID, exit_intersect]


def boxmotion(imp3d, particle, vf, wc, isEl, region, extends, nr1D):
    vx = vf * np.cos(particle[2])
    vy = vf * np.sin(particle[2])
    rc, xc, yc = cyclotron_orbit(particle[0], particle[1], vx, vy, wc, isEl)

    boxes = [get_boxid_full(particle[0], particle[1], region, extends, nr1D)]
    points = [particle]
    while len(boxes) == 1 or boxes[-1] != -1:
        i, p = next_box(boxes[-1], points[-1],
                        xc, yc, rc, particle[2], isEl, region, extends, nr1D)
        boxes.append(i)
        points.append(p)

    return boxes, points


def show_boxmotion(axis, points, particle, vf, isEl):

    vx = vf * np.cos(particle[2])
    vy = vf * np.sin(particle[2])
    rc, xc, yc = cyclotron_orbit(particle[0], particle[1], vx, vy, wc, isEl)
    colors = plt.get_cmap('viridis', len(points) - 1).colors

    # Guarantee El is ccw and hole is cw,
    # which goes wrong at the point where you go through 0 back to 2pi
    # without this correction.
    for i, p in enumerate(points[1:], 1):
        if isEl and p[2] > points[i - 1][2]:
            points[:i] = [[p[0], p[1], p[2] + 2 * np.pi] for p in points[:i]]
        if not isEl and p[2] < points[i - 1][2]:
            points[:i] = [[p[0], p[1], p[2] - 2 * np.pi] for p in points[:i]]

    for p1, p2, color in zip(points[:-1], points[1:], colors):
        theta = np.linspace(p1[2], p2[2], 100)

        xx = xc + rc * np.cos(angle_position(theta, isEl))
        yy = yc + rc * np.sin(angle_position(theta, isEl))
        axis.plot(xx, yy, lw=2, color=color)
        axis.scatter(p1[0], p1[1], s=200, color=color)


def get_test_impurities():
    path = os.path.join(os.path.dirname(__file__), 'test_impurities.dat')
    imps = []
    with open(path) as f:
        line = f.readline()
        while 'x' not in line or 'y' not in line:
            line = f.readline()
        for line in f:
            if line:
                parts = line.split()
                imps.append([float(parts[0]), float(parts[1])])

    return np.array(imps)


def test_boxmotion():
    """ Trace out full orbits and make sure the motion is correct. """

    # For these settings, you get that:
    #   - B=2 leaves the region and does not fulfill the orbit, left side
    #   - B=4 leaves the region along the bottom/top side.
    #   - B=10 is conventional circle
    #   - B=42 retraces through the starting square halfway through
    #   - B=50 stays within a single square
    # This concludes all the regimes and boundary cases.
    # The tests also run past every significant if statement.
    region = extends = 1e-6
    imps = get_test_impurities()
    imp3d, nr1D = grid_divide(imps, 1e-8, 1e-6, 1e-6, 10)
    particle = np.array([7e-7, 2.98e-7, 0])
    m = 5 * 9.109e-31
    vf = 1.68e5

    # Crash into the negative x boundary for electrons
    wc = E * 2 / m
    boxids, points = boxmotion(
        imp3d, particle, vf, wc, True, region, extends, nr1D)
    assert(boxids[-1] == -1)
    assert(boxids[0] != boxids[-2])
    assert(get_xind_id(boxids[-2], nr1D) == 0)
    assert(get_yind_id(boxids[-2], nr1D) > get_yind_id(boxids[0], nr1D))
    assert(points[-1][0] != points[0][0])
    assert(points[-1][1] != points[0][1])

    # And the same happens for holes, but the y value is lower
    # rather than higher. (Deflects to bottomleft, electrons topleft)
    boxids2, points2 = boxmotion(
        imp3d, particle, vf, wc, False, region, extends, nr1D)
    assert(boxids2[-1] == -1)
    assert(boxids2[0] != boxids2[-2])
    assert(get_xind_id(boxids2[-2], nr1D) == 0)
    assert(get_yind_id(boxids2[-2], nr1D) < get_yind_id(boxids2[0], nr1D))
    assert(get_yind_id(boxids2[-2], nr1D) < get_yind_id(boxids[-2], nr1D))
    assert(points2[-1][0] == -extends)

    # Crash into the top side for electrons once the orbit is sharp enough
    # that the left boundary is not touched
    wc = E * 4 / m
    boxids, points = boxmotion(
        imp3d, particle, vf, wc, True, region, extends, nr1D)
    assert(boxids[-1] == -1)
    assert(boxids[0] != boxids[-2])
    assert(get_yind_id(boxids[-2], nr1D) == nr1D - 1)
    assert(get_xind_id(boxids[-2], nr1D) < get_xind_id(boxids[0], nr1D))
    assert(points[-1][1] == region + extends)

    # And the same happens for holes, but the y value is lower
    # rather than higher. (Deflects to bottomleft, electrons topleft)
    boxids2, points2 = boxmotion(
        imp3d, particle, vf, wc, False, region, extends, nr1D)
    assert(boxids2[-1] == -1)
    assert(boxids2[0] != boxids2[-2])
    assert(get_yind_id(boxids2[-2], nr1D) == 0)
    assert(get_xind_id(boxids2[-2], nr1D) < get_xind_id(boxids2[0], nr1D))
    assert(points2[-1][1] == -extends)

    # Create a normal orbit where each square is traversed once
    # except for the starting point when we terminate.
    wc = E * 10 / m
    boxids, points = boxmotion(
        imp3d, particle, vf, wc, True, region, extends, nr1D)
    assert(boxids[-1] == -1)
    assert(boxids[0] == boxids[-2])
    assert(abs(points[-1][0] - points[0][0]) < 1e-15)
    assert(abs(points[-1][1] - points[0][1]) < 1e-15)

    # And the same happens for holes, but the y value is lower
    # rather than higher. (Deflects to bottomleft, electrons topleft)
    boxids2, points2 = boxmotion(
        imp3d, particle, vf, wc, False, region, extends, nr1D)
    assert(boxids2[-1] == -1)
    assert(boxids2[0] == boxids2[-2])
    assert(len(boxids) == len(boxids2))

    # Trace through the starting square three times, namely
    # start, finish and cutting off a corner halfway through.
    # That corner should not terminate.
    wc = E * 42 / m
    boxids, points = boxmotion(
        imp3d, particle, vf, wc, True, region, extends, nr1D)
    assert(boxids[-1] == -1)
    assert(boxids[0] == boxids[-2])
    assert(boxids[0] in boxids[1:-2])

    # Stay within the box and never leave it
    wc = E * 50 / m
    boxids, points = boxmotion(
        imp3d, particle, vf, wc, True, region, extends, nr1D)
    assert(boxids[-1] == -1)
    assert(len(boxids) == 2)


test_boxmotion()

##############################
# Old verified intersection code
##############################

# Adapted from code2, 20-07-2020
# Expanded for is_coherent, but NOT for is_diag_regions.


@njit()
def is_crossed(xA, yA, rA, xB, yB, rB):
    """ Find out if and how these circles interact.

    Returns 0 if the circles are disjoint.
    Returns 1 if the first circle is fully inside the second.
    Returns 2 if the second circle is fully inside the first.
    Returns 3 if the circles touch
    Returns 4 if they cross.
    """

    distance2 = (xA - xB)**2 + (yA - yB)**2
    if distance2 > (rA + rB)**2:
        return 0

    if distance2 < (rA - rB)**2:
        if rA < rB:
            return 1
        else:
            return 2

    if distance2 == (rA + rB)**2 or distance2 == (rA - rB)**2:
        return 3

    return 4


@njit()
def crosspoints(xA, yA, rA, xB, yB, rB):
    """ Blindly assuming the circles cross (not touch),
    find the (x1,y1,x2,y2) where this happens.
    """

    dist = np.sqrt((xA - xB)**2 + (yA - yB)**2)
    xstar = (dist**2 + rA**2 - rB**2) / (2 * dist)
    assert(xstar > 0)  # non-crossing or touching (==) excluded.
    ystar = rA**2 - xstar**2
    assert(ystar > 0)
    ystar = np.sqrt(ystar)

    ux = (xB - xA) / dist
    uy = (yB - yA) / dist
    vx = uy
    vy = -ux

    x1 = xA + ux * xstar + vx * ystar
    x2 = xA + ux * xstar - vx * ystar
    y1 = yA + uy * xstar + vy * ystar
    y2 = yA + uy * xstar - vy * ystar

    return x1, y1, x2, y2


@njit()
def crossangle(phi0, phi1, isCw):
    """ Private. Find how many radians you have to move to get
    from phi0 to phi1.

    This takes care of mod 2pi as well as rotation direction.
    """
    if isCw:
        return (phi0 - phi1) % (2 * np.pi)
    else:
        return (phi1 - phi0) % (2 * np.pi)


@njit()
def starts_inside(x, y, impurity, radius):
    """ Test if this coordinate is inside the impurity or not. """
    distance2 = (x - impurity[0])**2 + (y - impurity[1])**2
    return distance2 < radius**2


def test_is_crossing():
    """ Strong test that claims what the situation is.

    I am not going to test it 100% rigorously, that is
    very hard to do, but I can test each case.
    """

    status = is_crossed(1, 2, 1, 5, 5, 1)
    assert(status == 0)
    status = is_crossed(1, 2, 1, 1, 3, 5)
    assert(status == 1)
    status = is_crossed(1, 2, 5, 1, 3, 1)
    assert(status == 2)
    status = is_crossed(0, 0, 1, 0, 0.5, 0.5)
    assert(status == 3)
    status = is_crossed(0, 0, 1, 0, 0.5, 0.6)
    assert(status == 4)

    # One extreme case that I can imagine is
    # near the edge of one of the inequalities.
    status = is_crossed(0, 0, 1000, 2, 0, 1001)
    assert(status == 4)


test_is_crossing()


def test_starts_inside():
    """ Starting point in/out of an impurity. """

    yesno = starts_inside(1, 1, np.array([0, 0]), 2)
    assert(yesno)

    yesno = starts_inside(0, 0, np.array([1, 1]), 2)
    assert(yesno)

    yesno = starts_inside(2, 1, np.array([0, 0]), 2)
    assert(not yesno)

    # Edge = outside
    # A definition, but at least guaranteed.
    yesno = starts_inside(1, 0, np.array([0, 0]), 1)
    assert(not yesno)


test_starts_inside()


def test_crosspoint():
    """ Given two circles intersect. Find where. """

    # Symmetric circles
    x1, y1, x2, y2 = crosspoints(-1, 0, 1.5, 1, 0, 1.5)
    assert(abs(x1) < 1e-7)
    assert(abs(x2) < 1e-7)
    assert(abs(abs(y1) - (1.5**2 - 1)**0.5) < 1e-7)
    assert(abs(abs(y2) - (1.5**2 - 1)**0.5) < 1e-7)
    assert(abs(y1 + y2) < 1e-7)

    # Somewhat asymmetric
    x1, y1, x2, y2 = crosspoints(0, 0, 1, 1, 1, 1.5)
    assert(abs(x1**2 + y1**2 - 1) < 1e-7)
    assert(abs((x1 - 1)**2 + (y1 - 1)**2 - 1.5**2) < 1e-7)
    assert(abs(x2**2 + y2**2 - 1) < 1e-7)
    assert(abs((x2 - 1)**2 + (y2 - 1)**2 - 1.5**2) < 1e-7)

    # Extreme asymmetry, just outside
    x1, y1, x2, y2 = crosspoints(100, 0, 100, -1, 0, 1.5)
    assert(abs((x1 - 100)**2 + y1**2 - 100**2) < 1e-7)
    assert(abs((x1 + 1)**2 + y1**2 - 1.5**2) < 1e-7)
    assert(abs((x2 - 100)**2 + y2**2 - 100**2) < 1e-7)
    assert(abs((x2 + 1)**2 + y2**2 - 1.5**2) < 1e-7)

    # Extreme asymmetry, just inside
    x1, y1, x2, y2 = crosspoints(100, -0.5, 100, 1, 0.5, 1.5)
    assert(abs((x1 - 100)**2 + (y1 + 0.5)**2 - 100**2) < 1e-7)
    assert(abs((x1 - 1)**2 + (y1 - 0.5)**2 - 1.5**2) < 1e-7)
    assert(abs((x2 - 100)**2 + (y2 + 0.5)**2 - 100**2) < 1e-7)
    assert(abs((x2 - 1)**2 + (y2 - 0.5)**2 - 1.5**2) < 1e-7)


test_crosspoint()


def test_crossangle():
    """ Make sure direction and modulo are correct. """

    angle = crossangle(6.1, 0.1, True)
    assert(abs(angle - 6) < 1e-5)

    angle = crossangle(6.1, 0.1, False)
    assert(abs(angle - 0.28) < 0.01)


test_crossangle()

##############################
# Intersects in a box
##############################


@njit()
def boundangle(phi, alpha, isEl):
    """ By moving to the past, find at what angle the first bound is found """

    assert(alpha <= np.pi / 4)
    multiple = (phi + alpha) // (np.pi / 2)
    bound1 = multiple * (np.pi / 2) - alpha
    bound2 = multiple * (np.pi / 2) + alpha

    dangle1 = angle_difference(phi, bound1, isEl)
    dangle2 = angle_difference(phi, bound2, isEl)

    assert(min(dangle1, dangle2) < 2 * alpha)
    assert(max(dangle1, dangle2) > np.pi / 3 or min(dangle1, dangle2) == 0)

    if dangle1 < dangle2:
        return bound1
    else:
        return bound2


@njit()
def traversal_time(angle1, angle2, isEl, wc):
    """ Find the time to go 1->2 """

    dangle = angle_difference(angle1, angle2, isEl)
    dangle = dangle % (2 * np.pi)
    return dangle / abs(wc)


@njit()
def lifetime_one(particle, imp, imp_rad, xc, yc, rc, wc, isEl, anglerange):
    """ Find in B!=0 if there is ever an intersect.
    If yes, return intersect time. Else return 1.
    """

    lifetime = 1
    if is_crossed(xc, yc, rc, imp[0], imp[1], imp_rad) == 4:
        x1, y1, x2, y2 = crosspoints(xc, yc, rc,
                                     imp[0], imp[1], imp_rad)

        angle1 = get_angle(x1, y1, xc, yc, rc, isEl)
        angle2 = get_angle(x2, y2, xc, yc, rc, isEl)

        if angle_in_range(angle1, anglerange, isEl):
            t = traversal_time(particle[2], angle1, isEl, wc)
            lifetime = min(lifetime, t)

        if angle_in_range(angle2, anglerange, isEl):
            t = traversal_time(particle[2], angle2, isEl, wc)
            lifetime = min(lifetime, t)

    return lifetime


def box_intersect(particle, xc, yc, rc, wc, isEl, imp3D, boxID, imp_rad, anglerange):
    """ Find if there is an intersect in this crossing of boxID.
    If so, return the lowest crossing time. If not, return 1. """

    global COUNTER
    T = 1
    for impurity in imp3D[boxID]:
        COUNTER+=1
        t = lifetime_one(particle, impurity, imp_rad,
                         xc, yc, rc, wc, isEl, anglerange)
        T = min(T, t)
    return T


def lifetime(particle, xc, yc, rc, wc, region, extends, nr1D,
             imp3d, imp_rad, alpha, isEl, iscoh):
    """ For one particle at one angle, compute the lifetime.
    particle: 3-array of [x, y, velocity angle]
    """

    # First do the starting box to see if we start inside an impurity
    boxID = get_boxid_full(particle[0], particle[1], region, extends, nr1D)
    for impurity in imp3d[boxID]:
        if starts_inside(particle[0], particle[1], impurity, imp_rad):
            return 0

    # then determine how far we have to go for bound intersect.
    if iscoh:
        bound_distance = 1e99
        tbound = 1e99
    else:
        bound_angle = boundangle(particle[2], alpha, isEl)
        bound_distance = angle_difference(particle[2], bound_angle, isEl)
        tbound = traversal_time(particle[2], bound_angle, isEl, wc)

    # Then loop through the boxes in chronologic order / following the orbit
    entrypoint = particle
    nowID = boxID
    while nowID != -1:
        # Find out the part of the trajectory in this box ...
        nextID, exitpoint = next_box(nowID, entrypoint,
                                     xc, yc, rc, particle[2],
                                     isEl, region, extends, nr1D)

        # ... to find which angles are valid, for in case we
        # traverse this box multiple times we only take the current
        # segment.
        # The last time this is done, is the first time ID==-1,
        # which means we did not overlook the last interval.
        anglerange = np.array([entrypoint[2], exitpoint[2]])
        t = box_intersect(particle, xc, yc, rc, wc, isEl, imp3d,
                          nowID, imp_rad, anglerange)

        # Because it is chronologic, we can terminate in the first
        # box we find an intersection. 1 is no intersection.
        if t < 1:
            break

        # If we go beyond the bound, then we can stop early.
        nowID = nextID
        entrypoint = exitpoint
        if bound_distance < angle_difference(particle[2], entrypoint[2], isEl):
            break

    t = min(t, tbound)
    return t


def lifetime_full(particle, vf, wc,
                  region, extends, nr1D, imp3d, imp_rad,
                  alpha, isEl, iscoh):

    # I put this here separate, because the suggestion was to do
    # all angles within a single thread. The cyclotron orbit can be
    # computed once for all angles. It is not a big deal overall.
    #
    # We'll likely do 13*4 incoherent values in one sweep
    # Do that for all of it and do the integral to get
    # sigma_inc_xx and sigma_inc_xy
    #
    # Then 13*4 coherent values in one sweep.
    # Do that for all of it and do the integral to get
    # sigma_coh_xx and sigma_inc_xy
    vx = vf * np.cos(particle[2])
    vy = vf * np.sin(particle[2])
    rc, xc, yc = cyclotron_orbit(particle[0], particle[1], vx, vy, wc, isEl)
    return lifetime(particle, xc, yc, rc, wc, region, extends, nr1D,
                    imp3d, imp_rad, alpha, isEl, iscoh)


def plot_lifetime(ax, particle, vf, wc, time, isEl):
    """ Given the result of how long this guy has lived,
    plot its trajectory thus far. """

    vx = vf * np.cos(particle[2])
    vy = vf * np.sin(particle[2])
    rc, xc, yc = cyclotron_orbit(particle[0], particle[1], vx, vy, wc, isEl)
    sgn = -1 if isEl else 1

    if time < 1:
        phi_scatter = particle[2] + sgn * wc * time
        color = 'black'
    else:
        phi_scatter = particle[2] + 2 * np.pi
        color = 'red'

    theta = np.linspace(particle[2], phi_scatter, 100)
    xx = xc + rc * np.cos(theta + sgn * np.pi / 2)
    yy = yc + rc * np.sin(theta + sgn * np.pi / 2)

    ax.plot(xx, yy, lw=5, color=color)


def test_lifetime():
    """ These are as hard tests as I can make them.

    The main elements added by scattering compared to just
    motion through boxes, is the early termination,
    no intersect, staying within a single box,
     """

    region, extends = 1e-6, 1e-6
    particle = np.array([7e-7, 2.99e-7, 0])
    vf = 1.68e5
    wc = E * 10 / (5 * 9.109e-31)
    imps = get_test_impurities()
    imp3d, nr1D = grid_divide(imps, 1e-8, region, extends, 10)
    alpha = np.pi / 4
    imp_rad = 1e-8

    # Intersect in the second box
    t = lifetime_full(particle, vf, wc,
                      region, extends, nr1D, imp3d, imp_rad,
                      alpha, True, True)
    radians = wc * t
    assert(0.88 < radians < 0.89)

    # Intersect in the starting box
    t = lifetime_full(particle, vf, wc,
                      region, extends, nr1D, imp3d, imp_rad,
                      alpha, False, True)
    radians = wc * t
    assert(0.078 < radians < 0.079)

    # Boundary limited, even if the boundary is in the
    # same box as the first impurity intersect
    t = lifetime_full(particle, vf, wc,
                      region, extends, nr1D, imp3d, imp_rad,
                      alpha, True, False)
    radians = wc * t
    assert(abs(radians - np.pi / 4) < 1e-7)

    # For coherent, most of a circle through many boxes
    particle = np.array([7e-7, 4.3e-7, 0])
    t = lifetime_full(particle, vf, wc,
                      region, extends, nr1D, imp3d, imp_rad,
                      alpha, True, True)
    radians = wc * t
    assert(3.82 < radians < 3.83)

    # No intersect, simulated at high field as this is hard
    particle = np.array([6.5e-7, 1.35e-6, 0])
    wc = E * 50 / (5 * 9.109e-31)
    t = lifetime_full(particle, vf, wc,
                      region, extends, nr1D, imp3d, imp_rad,
                      alpha, False, True)
    assert(t == 1)

    # Make sure that within a single box does work out,
    # which traces through a different codeline in
    # angle_in_range than other kinds of orbits
    #
    # This one also tests 2 impurities that are close together
    # within the same box and both crossed,
    # but it is the first one that counts.
    particle = np.array([6.15e-7, 1.35e-6, 0])
    t = lifetime_full(particle, vf, wc,
                      region, extends, nr1D, imp3d, imp_rad,
                      alpha, False, True)
    radians = wc * t
    assert(2.21 < radians < 2.215)



    # CORNER CASE
    # Start on a horizontal edge.
    # Python will not assign you to the above or below ID based on floating
    # point precision, instead it will assign you to the higher index
    # through integer division.
    #
    # Take then a hole, which moves clockwise to the past.
    # The initial nextbox then cannot get stuck on the particle itself,
    # and the box to the south intersects AT the particle itself and
    # thus also cannot be chosen. Subsequently you are doomed to choose
    # the backtracking edge at about 5 radians.
    #
    # Then you look for intersects, which will be ignored since you
    # are supposed to reach 5 radians by taking the long way around
    # and you are looking for intersects the short way around.
    #
    # Then you go to the next iteration and you will find that the
    # closest point to move to is both the boxintersect and starting
    # position. The starting point is considered first and
    # equality ignored, hence you go to the starting point and terminate.
    #
    # It is a complicated story, but the bottom line is that NO particle
    # should EVER start on a gridline exactly. If there is floating point
    # precision, that is fine, but without it (like in Python) you get
    # stuck by indecision.
    wc = E * 8 / (5 * 9.109e-31)
    particle = np.array([1.00000000e-06, 5.0000000e-07, 1.57079633e+00])
    try:
        t = lifetime_full(particle, vf, wc,
                         region, extends, nr1D, imp3d, imp_rad * 2,
                          alpha, False, True)
    except ValueError as e:
        assert('grid' in str(e))
    else:
        assert(False and "Particles on gridlines go undetected.")

    wc = E * 8 / (5 * 9.109e-31)
    particle = np.array([1.00000000e-06, 5.0000100e-07, 1.57079633e+00])
    t = lifetime_full(particle, vf, wc,
                      region, extends, nr1D, imp3d, imp_rad * 2,
                      alpha, False, True)
    assert(t < 1)

    wc = E * 8 / (5 * 9.109e-31)
    particle = np.array([1.00000000e-06, 4.9999900e-07, 1.57079633e+00])
    t = lifetime_full(particle, vf, wc,
                      region, extends, nr1D, imp3d, imp_rad * 2,
                      alpha, False, True)
    assert(t < 1)


test_lifetime()

##############################
# Gridevaluation
##############################


@njit()
def get_xcoordinate(index, dimension, region):
    """ For x or y since it is square.

    These are the particle positions.
    Slightly offset
    """
    off = region / dimension * 1e-2
    assert(0 <= index < dimension)
    return index / (dimension - 1) * region + off


@njit()
def get_ycoordinate(index, dimension, region):
    """ For x or y since it is square.

    These are the particle positions.
    Slightly offset
    """
    off = region / dimension * 5e-3
    assert(0 <= index < dimension)
    return index / (dimension - 1) * region + off


@njit()
def get_phi(index, nr_phi, alpha, quadrant, iscoh):
    """ Currently ONLY for horizontal regions.

    quadrant skips by angles 90 forward starting at
    -alpha for iscoh and at alpha for incoherent.
    """
    assert(0 <= index < nr_phi)
    assert(0 <= alpha <= np.pi / 4)

    if iscoh:
        lower_bound = alpha + quadrant * np.pi / 2
        res = lower_bound + index / (nr_phi - 1) * (np.pi / 2 - 2 * alpha)
    else:
        lower_bound = -alpha + quadrant * np.pi / 2
        res = lower_bound + index / (nr_phi - 1) * 2 * alpha
    return res % (2 * np.pi)


def grid_lifetime(dimension, region, extends,
                  imp3d, nr1D, imp_rad,
                  wc, isEl, alpha, nr_phi, iscoh):

    Exists = False

    particle = np.zeros(3)
    result = np.zeros((dimension, dimension, 4 * nr_phi))
    for xind in range(dimension):
        particle[0] = get_xcoordinate(xind, dimension, region)
        for yind in range(dimension):
            particle[1] = get_ycoordinate(yind, dimension, region)
            for pind in range(nr_phi):
                for quadrant in range(4):
                    particle[2] = get_phi(pind, nr_phi, alpha, quadrant, iscoh)
                    vx = vf * np.cos(particle[2])
                    vy = vf * np.sin(particle[2])
                    rc, xc, yc = cyclotron_orbit(particle[0], particle[1],
                                                 vx, vy, wc, isEl)
                    pq_ind = pind + quadrant * nr_phi
                    try:
                        result[xind, yind, pq_ind] = \
                            lifetime(particle, xc, yc, rc, wc, region, extends, nr1D,
                                     imp3d, imp_rad, alpha, isEl, iscoh)
                    except:
                        # DEBUG code only
                        fig, ax = plt.subplots()
                        print(particle)
                        show_grid(ax, imp3d, imp_rad, region, extends)
                        plot_cyclotron(ax, particle, vf, wc, isEl)
                        plt.show()
                        raise

                    if result[xind, yind, pq_ind] == 1:
                        if not Exists:
                            fig, ax = plt.subplots()
                            show_grid(ax, imp3d, imp_rad, region, extends)
                            Exists = True
                        plot_cyclotron(ax, particle, vf, wc, isEl)
    if Exists:
        plt.show()
    return result


@njit()
def simpson_weight(index, length):
    """ Get the weight (1, 2, 4) at the given index position. """

    assert(length % 2)

    if index == 0 or index == length - 1:
        return 1
    if index % 2:
        return 4
    return 2


@njit()
def eval_lifetime(lifetimes, dimension, nr_phi, region, extends, alpha, iscoh):
    """ Integrate lifetimes over space and angle.

    This DOES divide out the total space and angle.
    """

    integral = 0
    weight = 1
    for xind in range(dimension):
        weightX = simpson_weight(xind, dimension)
        for yind in range(dimension):
            weightY = simpson_weight(yind, dimension)
            for pind in range(nr_phi):
                weightP = simpson_weight(pind, dimension)
                for quadrant in range(4):
                    pq_ind = pind + quadrant * nr_phi
                    weight = weightX * weightY * weightP
                    integral += lifetimes[xind, yind, pq_ind] * weight

    phirange = np.pi / 2 - 2 * alpha if iscoh else 2 * alpha
    integral *= phirange / (6**3 * 4 * (nr_phi - 1) * (dimension - 1)**2)
    average = integral / phirange
    return average


@njit()
def eval_sigma(lifetimes, tau, isEl, wc, mass, vf, caxis,
               dimension, nr_phi,
               alpha, iscoh):

    if isEl:
        wc *= -1

    integralxx = 0
    integralxy = 0
    weight = 1
    for xind in range(dimension):
        weightX = simpson_weight(xind, dimension)
        for yind in range(dimension):
            weightY = simpson_weight(yind, dimension)
            for pind in range(nr_phi):
                weightP = simpson_weight(pind, dimension)
                for quadrant in range(4):
                    pq_ind = pind + quadrant * nr_phi
                    phi = get_phi(pind, nr_phi, alpha, quadrant, iscoh)
                    weight = weightX * weightY * weightP

                    # See the 20/07/2020 introduction
                    # I moved one factor tau to the end, outside the loop.
                    t = lifetimes[xind, yind, pq_ind]
                    exp = np.exp(-t / tau)
                    integrand = np.cos(phi)
                    integrand -= np.cos(phi + wc * t) * exp
                    integrand += wc * tau * np.sin(phi + wc * t) * exp
                    integrand -= wc * tau * np.sin(phi)

                    integralxx += integrand * weight * np.cos(phi)
                    integralxy += integrand * weight * np.sin(phi)

    # These factors are stepsize as given by "interval / (nr - 1)"
    # and for the x/y integrals we divide by the interval (=region) to get the
    # average rather than integral, leading to "1/(dimension-1)"
    phirange = np.pi / 2 - 2 * alpha if iscoh else 2 * alpha
    integralxx *= phirange / (4 * (nr_phi - 1) * (dimension - 1)**2)
    integralxy *= phirange / (4 * (nr_phi - 1) * (dimension - 1)**2)

    kf = mass * vf / HBAR
    factor = E**2 * kf**2 * tau / (6**3 * 2 * np.pi**2 * mass * caxis)
    factor /= 1 + (wc * tau)**2
    sxx = integralxx * factor
    sxy = integralxy * factor
    return sxx, sxy


def computation(Bvals, Tvals, taucoh, vf, caxis,
                dimension, region, extends,
                nr_imp, avg_box, imp_rad,
                mass, isEl, alpha, nr_phi):
    """ Compute for all T and B the average lifetime, sxx, sxy
    split for coherent and incoherent.

    Return an array with:
        1st dimension is T
        2nd dimension is B
        3rd dimension is time_coh, time_inc, sxx_c, sxx_i, sxy_c, sxy_i
    """

    global COUNTER
    imps = generate_impurities(nr_imp, -extends, region + extends)
    imp3d, nr1D = grid_divide(imps, imp_rad, region, extends, avg_box)

    fig, ax = plt.subplots()
    show_grid(ax, imp3d, imp_rad, region, extends)

    result = np.zeros((len(Tvals), len(Bvals), 6))
    for bind, b in enumerate(Bvals):
        wc = E * b / mass
        st = time.time()
        COUNTER = 0
        incoherent = grid_lifetime(dimension, region, extends,
                                   imp3d, nr1D, imp_rad,
                                   wc, isEl, alpha, nr_phi, False)
        coherent = grid_lifetime(dimension, region, extends,
                                 imp3d, nr1D, imp_rad,
                                 wc, isEl, alpha, nr_phi, True)

        timer = time.time() - st
        efficiency = COUNTER / (dimension**2 * len(imps) * nr_phi * 4 * 2) * 100
        print(f'Computation for B={b:.2f} lasted {timer:.2f} s, {COUNTER / timer:.0f} '
              f'intersects/s, {efficiency:.2f} % of full intersect')
        st = time.time()
        for tind, t in enumerate(Tvals):
            result[tind, bind, 0] = eval_lifetime(coherent, dimension,
                                                  nr_phi, region, extends,
                                                  alpha, True)
            result[tind, bind, 1] = eval_lifetime(incoherent, dimension,
                                                  nr_phi, region, extends,
                                                  alpha, False)

            tauinc = HBAR / (KB * t)
            sxxc, sxyc = eval_sigma(coherent, taucoh, isEl, wc, mass, vf, caxis,
                                    dimension, nr_phi, alpha, True)
            sxxi, sxyi = eval_sigma(incoherent, tauinc, isEl, wc, mass, vf, caxis,
                                    dimension, nr_phi, alpha, False)
            result[tind, bind, 2] = sxxc
            result[tind, bind, 3] = sxxi
            result[tind, bind, 4] = sxyc
            result[tind, bind, 5] = sxyi
        print(f'Generalised over temperature in {time.time() - st:.2f} s')

    return result

##############################
# Main
##############################


Bvals = np.linspace(0.1, 10, 10)
Tvals = [0.1, 4, 10]
tcoh = 1e-11
caxis = 11.5e-10

r = computation(Bvals, Tvals, tcoh, vf, caxis,
                dimension, region, extends,
                imp_nr, avg_box, imp_rad,
                mass, isEl, alpha, nr_phi)

for i, T in enumerate(Tvals):
    sxx = r[i, :, 2] + r[i, :, 3]
    sxy = r[i, :, 4] + r[i, :, 5]

    kf = mass * vf / HBAR
    n = kf**2 / (2 * np.pi * caxis)
    drude = r[i, :, 0] * n * E**2 / mass
    sxxc = r[i, :, 2]
    print(r[i,:,0])

    rxx = sxx / (sxx**2 + sxy**2)
    rxy = sxy / (sxx**2 + sxy**2)

    plt.figure('sigma')
    plt.plot(Bvals, sxxc, label=f'T={T:.2f} K')
    plt.plot(Bvals, drude, label=f'T={T:.2f} K Drude')
    plt.legend()

    plt.figure('rho')
    plt.plot(Bvals, rxx, label=f'T={T:.2f} K')

plt.xlabel('B (T)')
plt.ylabel('\u03C1$_{xx}$ (T)')
plt.show()
