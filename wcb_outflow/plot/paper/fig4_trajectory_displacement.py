"""
> python -m wcb_outflow.plot.paper.fig4_trajectory_displacement

Joint probability distributions of the horizontal and vertical (dtheta) displacement
between Lagrangian and isentropic trajectories initialised from the same starting points
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from pylagranto import trajectory

from ... import case_studies
from ...outflow import haversine


# Adjustable plot parameters
norm = LogNorm(vmin=1, vmax=350)
bins = 51
dx_range = [-2500, 2500]
dtheta_range = [-40, 10]
time_index = -1


def main():
    fig, axes = plt.subplots(4, 3, sharex="all", sharey="all", figsize=(16, 10))
    for m, case_name in enumerate(case_studies):
        case = case_studies[case_name]
        tr = trajectory.load(case.data_path / "isentropic_trajectories_from_volume.pkl")
        tr3d = trajectory.load(case.data_path / "3d_trajectories.pkl")

        # Filter trajectories that leave the domain (in either set of trajectories)
        idx = np.where(
            np.logical_and(
                tr.x[:, time_index] != -1000,
                tr3d.x[:, time_index] != -1000,
            )
        )[0]

        print("{} trajectories out of bounds for {}".format(
            len(tr) - len(idx), case.name)
        )

        tr = tr[idx]
        tr3d = tr3d[idx]

        for n, theta_level in enumerate(case.outflow_theta):
            idx = np.where(tr.data[:, 0, 2] == theta_level)[0]

            # Make sure the initial positions match
            assert (tr[idx].x[:, 0] == tr3d[idx].x[:, 0]).all()
            assert (tr[idx].y[:, 0] == tr3d[idx].y[:, 0]).all()

            # Calculate horizonatal and vertical (theta) displacements
            displacement = haversine([tr[idx].x, tr[idx].y], [tr3d[idx].x, tr3d[idx].y])
            theta = tr3d["air_potential_temperature"][idx] - theta_level

            # Use the x (longitude) difference to plot the sign of the displacement
            # (purely for visualisation)
            dx = tr3d[idx].x - tr[idx].x

            h, xed, yed, im = axes[m, n].hist2d(
                displacement[:, time_index]*np.sign(dx)[:, time_index],
                theta[:, time_index],
                bins=bins,
                range=[dx_range, dtheta_range],
                norm=norm,
            )

            axes[m, n].set_title("{}: {}K".format(case.name, theta_level))

            # Print out numbers for adjusting plot parameters
            print(case.name, theta_level)
            print("Max trajectories: ", h.max())
            print("Max dx: ", displacement.max())
            print("Max dtheta: ", theta.max())

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Number of Trajectories")

    axes[-1, 1].set_xlabel(r"$|\Delta \mathbf{x}|$ (km)")
    fig.text(0.05, 0.5, r"$\Delta \theta$ (K)", rotation="vertical")
    plt.savefig("fig4_trajectory_displacement.png")


if __name__ == '__main__':
    main()
