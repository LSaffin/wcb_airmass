"""
Replicate the figure 2 from Methven (2015). The schematic showing the isentropic volumes
of warm conveyor belt inflow and outflow with the trajectories connecting them
"""

import numpy as np
import matplotlib as mpl
import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from pylagranto import trajectory

from wcb_outflow import case_studies


def main():
    mpl.rcParams["figure.subplot.left"] = 0.0
    mpl.rcParams["figure.subplot.right"] = 1.0
    mpl.rcParams["figure.subplot.bottom"] = 0.2
    mpl.rcParams["figure.subplot.top"] = 1.0

    case = case_studies["IOP3"]

    cmap = plt.cm.inferno

    outflow = trajectory.load(case.data_path / "isentropic_trajectories_backward.pkl")
    outflow = outflow.select("air_potential_temperature", "==", case.outflow_theta[1], time=[outflow.relative_times[0]])
    inflow = trajectory.load(case.data_path / "isentropic_trajectories_from_inflow_forward.pkl")
    inflow = inflow.select("air_potential_temperature", "==", 300)

    tr_3d = trajectory.load(case.data_path / "3d_trajectories.pkl")

    tr_3d = tr_3d.select("air_potential_temperature", "<", case.outflow_theta[0], time=[tr_3d.relative_times[-1]])
    tr_3d = tr_3d.select("air_potential_temperature", ">", case.outflow_theta[1] - 1, time=[tr_3d.relative_times[0]])
    tr_3d = tr_3d.select("air_potential_temperature", "<", case.outflow_theta[1] + 1, time=[tr_3d.relative_times[0]])
    tr_3d = tr_3d.select("air_potential_temperature", ">", 0)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.gca(projection='3d')

    plot_trajectory_selection(ax, tr_3d.x-360, tr_3d.y, tr_3d["air_potential_temperature"], cmap((tr_3d["air_potential_temperature"] - 300)/25), spacing=100)

    ax.plot(outflow.x[:, 0]-360, outflow.y[:, 0], outflow.z[:, 0], color="c", lw=3)
    idx = np.where(outflow.x[:, -1] != -1000)[0]
    ax.plot(outflow.x[idx, -1]-360, outflow.y[idx, -1], outflow.z[idx, -1], color="c", lw=3)

    ax.plot(inflow.x[:, 0]-360, inflow.y[:, 0], inflow.z[:, 0], color="m", lw=3)
    ax.plot(inflow.x[:, -1]-360, inflow.y[:, -1], inflow.z[:, -1], color="m", lw=3)

    plot_trajectory_selection(
        ax,
        outflow.x[idx, :]-360,
        outflow.y[idx, :],
        outflow["air_potential_temperature"][idx, :],
        "c",
        spacing=10,
        alpha=0.25,
    )

    plot_trajectory_selection(
        ax,
        inflow.x-360,
        inflow.y,
        inflow["air_potential_temperature"],
        "m",
        spacing=10,
        alpha=0.25,
    )

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_zlabel(r"$\theta$ (K)")

    fs = 14
    ax.text(-40, 60, 325, r"$M_{out} (t=t_{out})$", fontsize=fs)
    ax.text(-75, 35, 325, r"$M_{out} (t=t_0)$", fontsize=fs)
    ax.text(-70, 35, 300, r"$M_{in} (t=t_0)$", fontsize=fs)
    ax.text(5, 55, 300, r"$M_{in} (t=t_{out})$", fontsize=fs)

    #norm = colors.Normalize(vmin=5, vmax=10)
    cax = fig.add_axes([0.1, 0.1, 0.8, 0.05])
    cbar = plt.colorbar(matplotlib.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(vmin=300, vmax=325),
        cmap=cmap,
    ), cax=cax, orientation="horizontal")

    cbar.set_label(r"$\theta$ (Lagrangian Trajectories)")

    plt.show()


def plot_trajectory_selection(ax, x, y, z, color, spacing=100, **kwargs):
    for n in range(1, len(x), spacing):
        for i in range(len(x[n]-1)):
            if type(color) == str:
                ax.plot(x[n, i:i + 2], y[n, i:i + 2], z[n, i:i + 2], color=color, **kwargs)
            else:
                ax.plot(x[n, i:i+2], y[n, i:i+2], z[n, i:i+2], color=color[n, i], **kwargs)


if __name__ == '__main__':
    main()
