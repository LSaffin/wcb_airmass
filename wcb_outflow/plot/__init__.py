import iris.plot as iplt


def background_map(ax, pv):
    ax.coastlines()
    ax.gridlines()

    # Add a contour for PV=2 and also shade PV>2
    iplt.contour(pv, [2], colors='k')
    pv.data = (pv.data > 2).astype(float)
    iplt.contourf(pv, [0.9, 1.1], colors="k", alpha=0.25)

    # With the projection the axis limits don't set well, so shrink it down as much as
    # possible
    x, y = pv.coord(axis="x"), pv.coord(axis="y")
    crs = x.coord_system.as_cartopy_crs()
    ax.set_extent([x.points.min(), x.points.max(), y.points.min(), y.points.max()], crs=crs)
