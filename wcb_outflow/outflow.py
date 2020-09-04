"""Identification of the outflow of warm conveyor belts
"""

from math import cos
import numpy as np
import matplotlib as mpl
from scipy.ndimage import filters
import iris.plot as iplt
from irise import convert, grid, plot
from irise.constants import r
from irise.diagnostics.rossby_waves import tropopause_contour


def get_points(dtheta, pv, filter_sizes=dict(dtheta=30, pv=30)):
    # Smooth the theta_adv and pv data first
    dtheta.data = filters.median_filter(dtheta.data, size=filter_sizes["dtheta"])
    pv.data = filters.median_filter(pv.data, size=filter_sizes["pv"])

    # Extract the contour surrounding the outflow region
    criteria = np.logical_and(pv.data < 2, dtheta.data > 0)
    pv.data = criteria.astype(int)
    cs = iplt.contour(pv, [0.5])
    contours = tropopause_contour.get_contour_verts(cs)
    closed_loop = tropopause_contour.get_tropopause_contour(contours[0])
    path = mpl.path.Path(closed_loop)

    # Create an array containing all the grid points within the outflow region
    lon, lat = grid.get_xy_grids(dtheta)
    points = np.transpose(np.array([lon.flatten(), lat.flatten()]))
    points = points[np.where(path.contains_points(points))]

    return closed_loop, points


def increase_circuit_resolution(points, resolution):
    """Add points around the circuit loop
    """
    # Loop over all points around the circuit.
    n = 0
    while n < len(points) - 1:
        # Calculate distances from current point to next point
        lat = (points[n, 1] + points[n + 1, 0]) / 2
        dlon = points[n + 1, 0] - points[n, 0]
        dlat = points[n + 1, 1] - points[n, 1]
        dx = r * cos(np.deg2rad(lat)) * np.deg2rad(dlon)
        dy = r * np.deg2rad(dlat)
        distance = (dx ** 2 + dy ** 2) ** 0.5
        # Insert a point if the next point is further away than the required
        # resolution
        if distance > resolution:
            # Put the next point along the same line
            frac = resolution / distance
            new_lon = points[n, 0] + frac * dlon
            new_lat = points[n, 1] + frac * dlat
            points = np.insert(points, n + 1, [new_lon, new_lat], 0)

        # Always jump to the next point. This will either be the inserted point
        # or the next point along that is within the required resolution
        n += 1

    return points