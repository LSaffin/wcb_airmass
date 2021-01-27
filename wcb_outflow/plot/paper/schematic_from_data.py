"""
Replicate the figure 2 from Methven (2015). The schematic showing the isentropic volumes
of warm conveyor belt inflow and outflow with the trajectories connecting them
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from pylagranto import trajectory

from wcb_outflow import case_studies


def main():
    case = case_studies["IOP3"]

    outflow = trajectory.load(case.data_path / "isentropic_trajectories_backward.pkl")
    outflow = outflow.select("air_potential_temperature", "==", case.outflow_theta[1], time=[outflow.relative_times[0]])
    inflow = trajectory.load(case.data_path / "isentropic_trajectories_from_inflow_forward.pkl")
    inflow = inflow.select("air_potential_temperature", "==", 300)

    tr_3d = trajectory.load(case.data_path / "3d_trajectories.pkl")

    tr_3d = tr_3d.select("air_potential_temperature", "<", case.outflow_theta[0], time=[tr_3d.relative_times[-1]])
    tr_3d = tr_3d.select("air_potential_temperature", ">", case.outflow_theta[1] - 1, time=[tr_3d.relative_times[0]])
    tr_3d = tr_3d.select("air_potential_temperature", "<", case.outflow_theta[1] + 1, time=[tr_3d.relative_times[0]])
    tr_3d = tr_3d.select("air_potential_temperature", ">", 0)

    fig = plt.figure(figsize=(8,5))
    ax = fig.gca(projection='3d')

    plot_trajectory_selection(ax, tr_3d.x, tr_3d.y, tr_3d["air_potential_temperature"], plt.cm.cubehelix(tr_3d.z/10000), spacing=100)

    ax.plot(outflow.x[:, 0], outflow.y[:, 0], outflow.z[:, 0], color="C0", lw=3)
    idx = np.where(outflow.x[:, -1] != -1000)[0]
    ax.plot(outflow.x[idx, -1], outflow.y[idx, -1], outflow.z[idx, -1], color="C0", lw=3)

    ax.plot(inflow.x[:, 0], inflow.y[:, 0], inflow.z[:, 0], color="C1", lw=3)
    ax.plot(inflow.x[:, -1], inflow.y[:, -1], inflow.z[:, -1], color="C1", lw=3)

    plot_trajectory_selection(
        ax,
        outflow.x[idx, :],
        outflow.y[idx, :],
        outflow["air_potential_temperature"][idx, :],
        "C0",
        spacing=10,
        alpha=0.25,
    )

    plot_trajectory_selection(
        ax,
        inflow.x,
        inflow.y,
        inflow["air_potential_temperature"],
        "C1",
        spacing=10,
        alpha=0.25,
    )

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
