import numpy as np
import matplotlib as mpl
from matplotlib.path import Path
import cartopy.crs as ccrs

import iris.plot as iplt


projection = ccrs.NearsidePerspective(central_longitude=-20, central_latitude=50)


def background_map(ax, pv):
    ax.coastlines()
    ax.gridlines()

    # With the projection the axis limits don't set well, so shrink it down as much as
    # possible
    x, y = pv.coord(axis="x"), pv.coord(axis="y")
    crs = x.coord_system.as_cartopy_crs()

    vertices = [[x.points[0], yp] for yp in y.points] + \
               [[xp, y.points[-1]] for xp in x.points] + \
               [[x.points[-1], yp] for yp in y.points[::-1]] + \
               [[xp, y.points[0]] for xp in x.points[::-1]]
    path = Path(np.array(vertices))
    ax.set_boundary(path, transform=crs)

    # Add a contour for PV=2 and also shade PV>2
    iplt.contour(pv, [2], colors='k', zorder=30)
    pv.data = (pv.data > 2).astype(float)
    iplt.contourf(pv, [0.9, 1.1], colors="k", alpha=0.25, zorder=20)


def set_map_rcparams():
    mpl.rcParams["figure.subplot.left"] = 0.0
    mpl.rcParams["figure.subplot.right"] = 1.0
    mpl.rcParams["figure.subplot.wspace"] = 0.0


def set_plot_rcparams():
    mpl.rcParams["figure.subplot.right"] = 0.95
    mpl.rcParams["figure.subplot.bottom"] = 0.05
    mpl.rcParams["figure.subplot.top"] = 0.95
