import numpy as np
import matplotlib.pyplot as plt

from pylagranto import trajectory

from wcb_outflow import case_studies
from wcb_outflow.plot import set_plot_rcparams


def main():
    #set_plot_rcparams()
    fig, axes = plt.subplots(2, 2, sharex="all", sharey="all", figsize=(8, 5))
    for n, case in enumerate(case_studies):
        plt.axes(axes[n // 2, n % 2])
        make_plot(case_studies[case])

    fig.text(0.5, 0.01, "Forecast Lead Time (hours)", ha="center")
    fig.text(0.05, 0.5, r"$\theta$ (K)", rotation="vertical", va="center")
    plt.legend()

    plt.show()


def make_plot(case):
    tr3d = trajectory.load(case.data_path / "3d_trajectories.pkl")

    # Trajectories that stay in the domain at all times
    tr3d = tr3d.select("air_potential_temperature", ">", 0)

    # Trajectories that start in the outflow
    tr3d = tr3d.select(
        "air_potential_temperature", ">=", case.outflow_theta[0]-1,
        time=[tr3d.relative_times[0]]
    )
    tr3d = tr3d.select(
        "air_potential_temperature", "<=", case.outflow_theta[-1]+1,
        time=[tr3d.relative_times[0]]
    )

    # Trajectories that end below the lowest outflow level
    tr3d = tr3d.select(
        "air_potential_temperature", "<", case.outflow_theta[0],
        time=[tr3d.relative_times[-1]]
    )
    times = [(t - case.start_time).total_seconds() // 3600 for t in tr3d.times]
    theta = tr3d["air_potential_temperature"]

    xMed = np.median(theta, axis=0)
    xMean = theta.mean(axis=0)
    x100 = np.percentile(theta, 100, axis=0)
    x95 = np.percentile(theta, 95, axis=0)
    x75 = np.percentile(theta, 75, axis=0)
    x25 = np.percentile(theta, 25, axis=0)
    x05 = np.percentile(theta, 5, axis=0)
    x00 = np.percentile(theta, 0, axis=0)

    # Make the plot
    plt.fill_between(times, x00, x100, color="b", alpha=0.1, label="0-100%")
    plt.fill_between(times, x05, x95, color="b", alpha=0.25, label="5-95%")
    plt.fill_between(times, x25, x75, color="grey", alpha=0.75, label="25-75%")
    plt.plot(times, xMean, color="k", lw=3, label="Mean")
    plt.plot(times, xMed, color="C1", lw=3, alpha=0.5, label="Median")
    plt.title(case.name)
    plt.grid(True)


if __name__ == '__main__':
    main()
