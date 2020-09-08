"""Identification of the outflow of warm conveyor belts
"""

from math import radians, cos, sin, asin, sqrt

import numpy as np
import matplotlib as mpl

import iris.plot as iplt
from iris.analysis import cartography


def get_points(dtheta, pv):
    # Extract the contour surrounding the outflow region
    criteria = np.logical_and(pv.data < 2, dtheta.data > 0)
    criteria = pv.copy(data=criteria.astype(int))
    cs = iplt.contour(criteria, [0.5], colors="g")
    contours = cs.allsegs[0]
    closed_loop = get_longest_closed_contour(contours)
    path = mpl.path.Path(closed_loop)

    # Create an array containing all the grid points within the outflow region
    lon, lat = cartography.get_xy_grids(dtheta)
    points = np.transpose(np.array([lon.flatten(), lat.flatten()]))
    points = points[np.where(path.contains_points(points))]

    return closed_loop, points


def get_longest_closed_contour(contours, threshold=100):
    """Find the tropopause contour that wraps around the globe

    Args:
        contours (list):
            A list of contours. Each individual contour is an array of
            coordinates for each point along the contour of shape (N, 2).

    """
    # Array of lengths (in km) of the contour boundaries
    lengths = np.array([contour_length(x) for x in contours])

    # Set the length of any non-closed contours to zero
    closed = np.array([is_closed_contour(x, threshold) for x in contours])
    lengths *= closed.astype(int)

    # Extract the longest closed contour
    imax = lengths.argmax()
    long_contour = contours[imax]

    if lengths[imax] == 0:
        raise ValueError("No closed contours found")

    return long_contour


def contour_length(points):
    """Contour length in kilometres

    Args:
        points: Nx2 array of longitude and latitude points around a contour (degrees)

    Returns:
        float: The total length of the contour in kilometres
    """
    conlen = haversine(points[-1], points[0])
    for n in range(len(points) - 1):
        conlen += haversine(points[n], points[n+1])

    return conlen


def is_closed_contour(contour_section, threshold=100):
    """Checks that a contour is closed

    Checks that the final point along a contour is sufficiently close to the
    initial point on a countour to determine if it is closed

    Args:
        contour_section (np.Array):
            An array of coordinates for each point along the contour of shape
            (N, 2).

    Returns:
        True: If contour is closed

        False: If contour is open
    """
    return haversine(contour_section[0], contour_section[-1]) < threshold


def increase_circuit_resolution(points, resolution):
    """Add points around the circuit loop
    """
    # Loop over all points around the circuit.
    n = 0
    while n < len(points):
        # Allow the loop to connect the final point with the first point
        # i.e np1 = 0
        np1 = (n+1) % len(points)

        # Calculate distances from current point to next point and
        # Insert a point if the next point is further away than the required
        # resolution
        distance = haversine(points[n], points[np1])
        if distance > resolution:
            # Put the next point along the same line
            dlon = points[np1, 0] - points[n, 0]
            dlat = points[np1, 1] - points[n, 1]

            new_lon = points[n, 0] + (resolution / distance) * dlon
            new_lat = points[n, 1] + (resolution / distance) * dlat

            points = np.insert(points, np1, [new_lon, new_lat], axis=0)

        # Always jump to the next point. This will either be the inserted point
        # or the next point along that is within the required resolution
        n += 1

    return points


def haversine(x1, x2):
    """ Calculate the great circle distance between two points on the earth
    (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [x1[0], x1[1], x2[0], x2[1]])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of earth in kilometers
    return c * r


