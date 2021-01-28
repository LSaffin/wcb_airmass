"""
> python -m wcb_outflow.plot.paper.fig5_inflow

Maps showing identified inflow region using both Lagrangian and isentropic trajectories
"""
import numpy as np
import matplotlib.pyplot as plt

import iris

from pylagranto import trajectory

from wcb_outflow import case_studies, outflow
from wcb_outflow.plot import background_map, projection


def main():
    plt.figure(figsize=[8, 5])
    for m, case_name in enumerate(case_studies):
        case = case_studies[case_name]

        pv = iris.load_cube(
            case.filename_theta(case.start_time),
            iris.Constraint(
                name="ertel_potential_vorticity",
                air_potential_temperature=case.outflow_theta[1],
                time=case.start_time)
        )
        crs = pv.coord(axis="x").coord_system.as_cartopy_crs()
        tr = trajectory.load(case.data_path / "isentropic_trajectories_from_volume.pkl")

        # Select trajectories for theta level that do not leave the domain
        tr = tr.select("air_potential_temperature", ">=", case.outflow_theta[0])
        tr = tr.select("air_potential_temperature", "<=", case.outflow_theta[-1])

        inflow_points = np.load(str(case.data_path / "inflow_boundaries.npy"))

        ax = plt.subplot(2, 2, m+1, projection=projection)
        background_map(ax, pv)

        try:
            outflow.contour_around_points(
                tr.x[:, -1], tr.y[:, -1],
                pv,
                colors="cyan",
                linewidths=3
            )
        except ValueError:
            # IOP3 doesn't return a closed contour but we're only using the function to
            # produce the plot so just ignore the error
            pass

        ax.plot(
            inflow_points[:, 0] - 360, inflow_points[:, 1],
            color="magenta",
            lw=3,
            alpha=0.5,
            transform=crs
        )

        ax.set_title("{}".format(case.name))

    plt.show()


def _bin_edges(x):
    xed = np.zeros(len(x) + 1)
    dx = np.diff(x).mean()
    xed[:-1] = x - dx
    xed[-1] = x[-1] + dx

    return xed


if __name__ == '__main__':
    main()
