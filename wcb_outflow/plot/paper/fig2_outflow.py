"""
> python -m wcb_outflow.plot.paper.fig2_outflow

Maps showing the identified outflow regions with a background of the fields used to
identify it (dtheta + PV)
"""

import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
import cartopy.crs as ccrs

from pylagranto import trajectory

from ... import case_studies


def main():
    projection = ccrs.Mollweide()
    fig = plt.figure(figsize=[16, 10])

    for n, case_name in enumerate(case_studies):
        case = case_studies[case_name]
        theta_level = case.outflow_theta[1]

        # Load the forecast at the outflow time
        cubes = iris.load(
            case.filename_theta(case.outflow_time),
            iris.Constraint(time=case.outflow_time)
        )
        pv = cubes.extract_strict(
            iris.Constraint(
                "ertel_potential_vorticity",
                air_potential_temperature=theta_level)
        )
        dtheta = cubes.extract_strict(
            iris.Constraint(
                "total_minus_adv_only_theta",
                air_potential_temperature=theta_level)
        )

        # Load the isentropic trajectories on the isentropic level of the forecast data
        tr = trajectory.load(case.data_path / "isentropic_trajectories.pkl")
        tr = tr.select("air_potential_temperature", "==", theta_level,
                       time=[tr.relative_times[0]])

        ax = plt.subplot(2, 2, n+1, projection=projection)
        im = make_plot(ax, dtheta, pv, tr)

        ax.set_title("{}: {}K".format(case.name, int(theta_level)))

    plt.subplots_adjust(bottom=0.2)
    cax = plt.axes([0.1, 0.1, 0.8, 0.05])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label(r"$\theta - \theta_{adv}$ (K)")

    plt.show()

    return


def make_plot(ax, dtheta, pv, tr):
    im = iplt.pcolormesh(dtheta, vmin=-30, vmax=30, cmap="seismic")
    ax.coastlines()
    ax.gridlines()

    # Add a contour for PV=2 and also shade PV>2
    iplt.contour(pv, [2], colors='k')
    pv.data = (pv.data > 2).astype(float)
    iplt.contourf(pv, [0.9, 1.1], colors="k", alpha=0.25)

    # Need to use the cube's coordinate reference system to properly add regular
    # coordinates on top
    crs = dtheta.coord(axis="x").coord_system.as_cartopy_crs()
    plt.plot(tr.x[:, 0] - 360, tr.y[:, 0], transform=crs, lw=3, color='cyan')

    # With the projection the axis limits don't set well, so shrink it down as much as
    # possible
    x, y = pv.coord(axis="x"), pv.coord(axis="y")
    ax.set_extent([x.points.min(), x.points.max(), y.points.min(), y.points.max()],
                  crs=crs)

    return im


if __name__ == '__main__':
    main()
