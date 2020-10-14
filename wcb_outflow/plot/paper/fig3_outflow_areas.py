"""
> python -m wcb_outflow.plot.paper.fig3_outflow_areas

The area of the outflow area on each isentropic level
"""

import numpy as np
import matplotlib.pyplot as plt

from pylagranto import trajectory

from wcb_outflow import case_studies


def main():
    theta_levels = np.arange(300, 345, 5)

    fig, axes = plt.subplots(2, 2, sharex="all", sharey="all", figsize=(8, 5))
    for n, case_name in enumerate(case_studies):
        ax = axes[n // 2, n % 2]
        case = case_studies[case_name]

        tr = trajectory.load(case.data_path / "3d_trajectories.pkl")
        panel(ax, tr, theta_levels)
        ax.set_title(case_name)

        ax.axhspan(case.outflow_theta[0]-2.5, case.outflow_theta[-1]+2.5, color="r", alpha=0.25)

    fig.text(0.5, 0.01, "Outflow Area (number of gridpoints)", ha="center")
    fig.text(0.05, 0.5, "Isentropic Level (K)", rotation="vertical", va="center")
    plt.show()


def panel(ax, tr, theta_levels):
    sizes = []

    for theta_level in theta_levels:
        tr_theta = tr.select("air_potential_temperature", ">", theta_level - 1,
                             [tr.relative_times[0]])
        tr_theta = tr_theta.select("air_potential_temperature", "<", theta_level + 1,
                                   [tr.relative_times[0]])

        sizes.append(len(tr_theta))

    ax.plot(sizes, theta_levels, 'kx')


if __name__ == '__main__':
    main()
