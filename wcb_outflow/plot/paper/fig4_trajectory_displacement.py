"""
> python -m wcb_outflow.plot.paper.fig4_trajectory_displacement

Joint probability distributions of the horizontal and vertical (dtheta) displacement
between Lagrangian and isentropic trajectories initialised from the same starting points
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from pylagranto import trajectory

from wcb_outflow import case_studies
from wcb_outflow.outflow import haversine


# Adjustable plot parameters
norm = LogNorm(vmin=1, vmax=500)
bins = 51
dx_range = [-3000, 3000]
dtheta_range = [-40, 10]
time_index = -1


def main():
    fig, axes = plt.subplots(4, 1, sharex="all", sharey="all", figsize=(8, 10))
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

        # Select trajectories that start in the outflow
        idx = np.where(np.logical_and(
            tr.data[:, 0, 2] <= case.outflow_theta[-1],
            tr.data[:, 0, 2] >= case.outflow_theta[0],
        ))[0]

        tr = tr[idx]
        tr3d = tr3d[idx]

        # Make sure the initial positions match
        assert (tr.x[:, 0] == tr3d.x[:, 0]).all()
        assert (tr.y[:, 0] == tr3d.y[:, 0]).all()

        # Calculate horizonatal and vertical (dtheta) displacements
        displacement = haversine([tr.x, tr.y], [tr3d.x, tr3d.y])[:, time_index]
        dtheta = tr3d["air_potential_temperature"][:, time_index] - \
                 tr3d["air_potential_temperature"][:, 0]

        # Use the x (longitude) difference to plot the sign of the displacement
        # (purely for visualisation)
        dx = tr3d.x[:, time_index] - tr.x[:, time_index]

        h, xed, yed, im = axes[m].hist2d(
            displacement*np.sign(dx),
            dtheta,
            bins=bins,
            range=[dx_range, dtheta_range],
            norm=norm,
        )

        axes[m].set_title(case.name)

        # Print out numbers for adjusting plot parameters
        print(case.name)
        print("Max trajectories: ", h.max())
        print("Max dx: ", displacement.max())
        print("Max dtheta: ", dtheta.max())

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Number of Trajectories")

    axes[-1].set_xlabel(r"$|\Delta \mathbf{x}|$ (km)")
    fig.text(0.05, 0.5, r"$\Delta \theta$ (K)", rotation="vertical", va="center")
    plt.show()


if __name__ == '__main__':
    main()
