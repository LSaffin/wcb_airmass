"""
> python -m wcb_outflow.plot.paper.fig2_outflow

Maps showing the identified outflow regions with a background of the fields used to
identify it (dtheta + PV)
"""

import matplotlib.pyplot as plt
import iris
import iris.plot as iplt

from pylagranto import trajectory

from wcb_outflow.plot import background_map, projection, set_map_rcparams
from wcb_outflow import case_studies


def main():
    set_map_rcparams()
    fig = plt.figure(figsize=[8, 9])

    for n, case_name in enumerate(case_studies):
        case = case_studies[case_name]
        theta_level = case.outflow_theta[1]

        # Load the forecast at the outflow time
        cubes = iris.load(
            case.filename_theta(case.outflow_time),
            iris.Constraint(time=case.outflow_time)
        )
        pv = cubes.extract_cube(
            iris.Constraint(
                "ertel_potential_vorticity",
                air_potential_temperature=theta_level)
        )
        dtheta = cubes.extract_cube(
            iris.Constraint(
                "total_minus_adv_only_theta",
                air_potential_temperature=theta_level)
        )

        # Load the isentropic trajectories on the isentropic level of the forecast data
        tr = trajectory.load(case.data_path / "isentropic_trajectories_backward.pkl")
        tr = tr.select("air_potential_temperature", "==", theta_level,
                       time=[tr.relative_times[0]])

        ax = plt.subplot(2, 2, n+1, projection=projection)
        im = make_plot(ax, dtheta, pv, tr)

        ax.set_title("{}: {} K \n {} ($t_0 + {}$h)".format(
            case.name,
            int(theta_level),
            case.start_time.strftime("%d/%m %H UTC"),
            int(case.outflow_lead_time.total_seconds() // 3600),
        ))

    plt.subplots_adjust(bottom=0.2)
    cax = plt.axes([0.1, 0.1, 0.8, 0.05])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label(r"$\theta - \theta_{adv}$ (K)")

    plt.show()

    return


def make_plot(ax, dtheta, pv, tr):
    background_map(ax, pv)

    im = iplt.pcolormesh(dtheta, vmin=-30, vmax=30, cmap="seismic", zorder=1)

    # Need to use the cube's coordinate reference system to properly add regular
    # coordinates on top
    crs = dtheta.coord(axis="x").coord_system.as_cartopy_crs()
    plt.plot(tr.x[:, 0] - 360, tr.y[:, 0], transform=crs, lw=3, color='cyan', zorder=40)

    return im


if __name__ == '__main__':
    main()
