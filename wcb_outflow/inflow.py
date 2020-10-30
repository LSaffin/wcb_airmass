"""Identification of WCB inflow region
"""
import numpy as np
from scipy.ndimage import filters

import iris
import iris.plot as iplt

from pylagranto import trajectory

from wcb_outflow import outflow


def from_3d_trajectories(case):
    # Select trajectories that do not leave the domain
    tr3d = trajectory.load(case.data_path / "3d_trajectories.pkl")
    tr3d = tr3d.select("air_potential_temperature", ">", 0)

    # Select the trajectories starting on designated outflow levels
    tr3d = tr3d.select(
        "air_potential_temperature", ">", case.outflow_theta[0] - 1,
        time=[tr3d.relative_times[0]]
    )

    tr3d = tr3d.select(
        "air_potential_temperature", "<", case.outflow_theta[-1] + 1,
        time=[tr3d.relative_times[0]]
    )

    # Select 3d trajectories below the outflow level at the inflow time
    # (i.e trajectories that ascend into the outflow region)
    tr3d = tr3d.select(
        "air_potential_temperature", "<", case.outflow_theta[0],
        time=[tr3d.relative_times[-1]]
    )

    # Load an example cube (this is just used to get the xy grid
    cube = iris.load_cube(
        case.filename_theta(case.start_time),
        iris.Constraint(name="ertel_potential_vorticity", time=case.start_time)
    )
    x, y = cube.coord(axis="x").points, cube.coord(axis="y").points
    bins = (_bin_edges(x), _bin_edges(y))

    # Produce a contour around all the points in the inflow
    # Use the same median filter as for the outflows
    h, xed, yed = np.histogram2d(tr3d.x[:, -1], tr3d.y[:, -1], bins=bins)
    h_sm = filters.median_filter(h.transpose(), size=30)
    h_sm = cube[0].copy(data=h_sm)
    cs = iplt.contour(h_sm, [0.5])

    closed_loop = outflow.get_longest_closed_contour(cs.allsegs[0])
    np.save(str(case.data_path / "inflow_boundaries.npy"), closed_loop)


def _bin_edges(x):
    xed = np.zeros(len(x) + 1)
    dx = np.diff(x).mean()
    xed[:-1] = x - dx
    xed[-1] = x[-1] + dx

    return xed
