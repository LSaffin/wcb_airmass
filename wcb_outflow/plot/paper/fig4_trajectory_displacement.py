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
from wcb_outflow.plot import set_plot_rcparams


# Adjustable plot parameters
norm = LogNorm(vmin=1, vmax=500)
bins = 51
dx_range = [-3000, 3000]
dtheta_range = [-40, 10]
time_index = -1


def main():
    set_plot_rcparams()
    fig, axes = plt.subplots(2, 2, sharex="all", sharey="all", figsize=(8, 6))
    for m, case_name in enumerate(case_studies):
        ax = axes[m // 2, m % 2]
        case = case_studies[case_name]
        tr = trajectory.load(case.data_path / "isentropic_trajectories_from_volume.pkl")
        tr3d = trajectory.load(case.data_path / "3d_trajectories.pkl")

        # Select trajectories that start in the outflow
        idx = np.where(np.logical_and(
            tr.data[:, 0, 2] <= case.outflow_theta[-1],
            tr.data[:, 0, 2] >= case.outflow_theta[0],
        ))[0]

        tr = tr[idx]
        tr3d = tr3d[idx]

        # Filter trajectories that leave the domain (in either set of trajectories)
        idx = np.where(
            np.logical_and(
                tr.x[:, time_index] != -1000,
                tr3d.x[:, time_index] != -1000,
            )
        )[0]

        out_of_bounds = "{:.1f}% trajectories out of bounds".format(
            (len(tr) - len(idx))/len(tr) * 100
        )
        ax.text(0, 0, out_of_bounds, transform=ax.transAxes)

        tr = tr[idx]
        tr3d = tr3d[idx]

        # Calculate horizonatal and vertical (dtheta) displacements
        displacement = haversine([tr.x, tr.y], [tr3d.x, tr3d.y])[:, time_index]
        dtheta = tr3d["air_potential_temperature"][:, time_index] - \
                 tr3d["air_potential_temperature"][:, 0]

        # Use the x (longitude) difference to plot the sign of the displacement
        # (purely for visualisation)
        dx = tr3d.x[:, time_index] - tr.x[:, time_index]

        h, xed, yed, im = ax.hist2d(
            displacement*np.sign(dx),
            dtheta,
            bins=bins,
            range=[dx_range, dtheta_range],
            norm=norm,
        )

        ax.axvline(color="k", lw=1)
        ax.axhline(color="k", lw=1)
        ax.set_title(case.name + " (-{:d}h from outflow)".format(
            int((case.outflow_time - case.start_time).total_seconds() // 3600)
        ))

        # Print out numbers for adjusting plot parameters
        print(case.name)
        print("Max trajectories: ", h.max())
        print("Max dx: ", displacement.max())
        print("Max dtheta: ", dtheta.max())

    plt.subplots_adjust(bottom=0.25)
    cax = plt.axes([0.1, 0.1, 0.8, 0.05])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("Number of Trajectories")

    fig.text(0.5, 0.175, r"$|\Delta \mathbf{x}|$ (km)", ha="center")
    fig.text(0.05, 0.55, r"$\Delta \theta$ (K)", rotation="vertical", va="center")
    plt.show()


if __name__ == '__main__':
    main()
