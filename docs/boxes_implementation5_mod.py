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

##############################
# Griddivision
##############################

# Note: After determining it, I take nr1D to be the defining characteristic,
#       the other choice is to take boxdim to be leading, but it has to exactly
#       match the region size at a multiple - and that coincidence seems more
#       natural by taking the integer than having perfectly divisible floats.

 
def get_boxdim(region, extends, nr1D):
    return (region + 2 * extends) / nr1D


 
def get_xind(x, region, extends, nr1D):
    return (x + extends) // get_boxdim(region, extends, nr1D)


 
def get_yind(y, region, extends, nr1D):
    return (y + extends) // get_boxdim(region, extends, nr1D)


 
def get_boxid(xind, yind, nr1D):
    assert(xind >= 0)
    assert(xind < nr1D)
    assert(yind >= 0)
    assert(yind < nr1D)
    return int(yind * nr1D + xind)


 
def get_boxid_full(x, y, region, extends, nr1D):
    xind = get_xind(x, region, extends, nr1D)
    yind = get_yind(y, region, extends, nr1D)
    return get_boxid(xind, yind, nr1D)


 
def get_xlow(xind, region, extends, nr1D):
    return -extends + (region + 2 * extends) * xind / nr1D 


 
def get_ylow(yind, region, extends, nr1D):
    return -extends + (region + 2 * extends) * yind / nr1D


 
def get_xind_id(id, nr1D):
    assert(0 <= id < nr1D**2)
    return id % nr1D


 
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

##############################
# Gridmovement
##############################


 
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


##############################
# Old verified intersection code
##############################

# Adapted from code2, 20-07-2020
# Expanded for is_coherent, but NOT for is_diag_regions.


 
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

 
def crossangle(phi0, phi1, isCw):
    """ Private. Find how many radians you have to move to get
    from phi0 to phi1.

    This takes care of mod 2pi as well as rotation direction.
    """
    if isCw:
        return (phi0 - phi1) % (2 * np.pi)
    else:
        return (phi1 - phi0) % (2 * np.pi)

 
def starts_inside(x, y, impurity, radius):
    """ Test if this coordinate is inside the impurity or not. """
    distance2 = (x - impurity[0])**2 + (y - impurity[1])**2
    return distance2 < radius**2

##############################
# Intersects in a box
##############################

 
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


def traversal_time(angle1, angle2, isEl, wc):
    """ Find the time to go 1->2 """

    dangle = angle_difference(angle1, angle2, isEl)
    dangle = dangle % (2 * np.pi)
    return dangle / abs(wc)


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

    T = 1
    for impurity in imp3D[boxID]:
        t = lifetime_one(particle, impurity, imp_rad,
                         xc, yc, rc, wc, isEl, anglerange)
        T = min(T, t)
    return T


def lifetime(particle, xc, yc, rc, wc, region, extends, nr1D, imp3d, imp_rad, alpha, isEl, iscoh):
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


def lifetime_full(particle, vf, wc, region, extends, nr1D, imp3d, imp_rad, alpha, isEl, iscoh):

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


##############################
# Gridevaluation
##############################


 
def get_xcoordinate(index, dimension, region):
    """ For x or y since it is square.

    These are the particle positions.
    Slightly offset
    """
    off = region / dimension * 1e-2
    assert(0 <= index < dimension)
    return index / (dimension - 1) * region + off


 
def get_ycoordinate(index, dimension, region):
    """ For x or y since it is square.

    These are the particle positions.
    Slightly offset
    """
    off = region / dimension * 5e-3
    assert(0 <= index < dimension)
    return index / (dimension - 1) * region + off


 
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


def grid_lifetime(dimension, region, extends, imp3d, nr1D, imp_rad, wc, isEl, alpha, nr_phi, iscoh):

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


 
def simpson_weight(index, length):
    """ Get the weight (1, 2, 4) at the given index position. """

    assert(length % 2)

    if index == 0 or index == length - 1:
        return 1
    if index % 2:
        return 4
    return 2


 
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
                weightP = simpson_weight(pind, dimension) #@Fout, nr_steps.
                for quadrant in range(4):
                    pq_ind = pind + quadrant * nr_phi
                    weight = weightX * weightY * weightP
                    integral += lifetimes[xind, yind, pq_ind] * weight

    phirange = np.pi / 2 - 2 * alpha if iscoh else 2 * alpha
    integral *= phirange / (3**3 * 4 * (nr_phi - 1) * (dimension - 1)**2) #@Vraag, 6^3?
    average = integral / phirange
    return average


 
def eval_sigma(lifetimes, tau, isEl, wc, mass, vf, caxis, dimension, nr_phi, alpha, iscoh):

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
    integralxx *= phirange / (4 * (nr_phi - 1) * (dimension - 1)**2) #@Vraag, volgorde?
    integralxy *= phirange / (4 * (nr_phi - 1) * (dimension - 1)**2)

    kf = mass * vf / HBAR
    factor = E**2 * kf**2 * tau / (6**3 * 2 * np.pi**2 * mass * caxis) #@Vraag, 6^3?
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

    imps = generate_impurities(nr_imp, -extends, region + extends)
    imp3d, nr1D = grid_divide(imps, imp_rad, region, extends, avg_box)

    result = np.zeros((len(Tvals), len(Bvals), 6))
    for bind, b in enumerate(Bvals):
        wc = E * b / mass

        incoherent = grid_lifetime(dimension, region, extends,
                                   imp3d, nr1D, imp_rad,
                                   wc, isEl, alpha, nr_phi, False)
        coherent = grid_lifetime(dimension, region, extends,
                                 imp3d, nr1D, imp_rad,
                                 wc, isEl, alpha, nr_phi, True)

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

