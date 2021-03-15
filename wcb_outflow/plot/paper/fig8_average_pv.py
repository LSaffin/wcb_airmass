"""
Plot the average PV in the inflow and outflow regions at the inflow and outflow times
"""


from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

import iris

from wcb_outflow import case_studies
from wcb_outflow.plot.paper.fig7_bulk_fraction import time_zero
from wcb_outflow.plot import set_plot_rcparams


# Lowest isentropic layer where the inflow doesn't intersect the ground at any time
theta_min = dict(
    IOP3=300,
    IOP5=305,
    IOP6=290,
    IOP7=290,
)


def main():
    set_plot_rcparams()
    fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex="all", sharey="row")
    for n, case in enumerate(case_studies):
        ax = plt.axes(axes[n // 2, n % 2])
        make_figure(case_studies[case])
        ax.set_title(case)

    _, xmax = axes[0, 0].get_xlim()
    axes[0, 0].set_xlim(0, xmax)

    legend_elements = [
        Line2D([0], [0], color="k", marker="o", linestyle="", label=r"$t=t_{out}$"),
        Line2D([0], [0], color="k", marker="s", linestyle="", label=r"$t=t_{out} - \Delta t$"),
        Line2D([0], [0], color="k", marker="D", linestyle="", label=r"$t=t_0$")
    ]

    # Create the figure
    axes[0, 0].legend(handles=legend_elements)

    fig.text(0.5, 0.01, "Average PV (PVU)", ha="center")
    fig.text(0.05, 0.5, r"$\theta$ (K)", va="center", rotation="vertical")

    plt.show()


def make_figure(case):
    pv_outflow = load_PV(str(case.data_path / "circulation.nc"))

    for theta in case.outflow_theta:
        cube_outflow_time = pv_outflow.extract(iris.Constraint(
            air_potential_temperature=theta,
            time=case.outflow_time,
        ))

        # Inflow time for the outflow volume is actually just the earlies time at which
        # the data still looks reasonable
        cube_inflow_time = pv_outflow.extract(iris.Constraint(
            air_potential_temperature=theta,
            time=case.start_time + time_zero[case.name][theta] * case.timestep,
        ))

        plt.plot(
            [cube_inflow_time.data, cube_outflow_time.data],
            [theta, theta],
            "-k"
        )

        plt.plot((cube_inflow_time.data + cube_outflow_time.data) / 2, theta, "k<")
        plt.plot(cube_inflow_time.data, theta, "cs")
        plt.plot(cube_outflow_time.data, theta, "co")

    # Dividing line between inflow and outflow
    plt.axhline(case.outflow_theta[0] - 2.5, color="k", linewidth=1)

    pv_inflow = load_PV(str(case.data_path / "circulation_inflow.nc"))

    for theta in range(theta_min[case.name], case.outflow_theta[0], 5):
        cube_outflow_time = pv_inflow.extract(iris.Constraint(
            air_potential_temperature=theta,
            time=case.outflow_time,
        ))

        cube_inflow_time = pv_inflow.extract(iris.Constraint(
            air_potential_temperature=theta,
            time=case.start_time,
        ))

        plt.plot(
            [cube_inflow_time.data, cube_outflow_time.data],
            [theta, theta],
            "-k"
        )

        plt.plot((cube_inflow_time.data + cube_outflow_time.data) / 2, theta, "k>")
        plt.plot(cube_inflow_time.data, theta, "mD")
        plt.plot(cube_outflow_time.data, theta, "mo")

    return


def load_PV(filename):
    mass = iris.load_cube(filename, "mass")
    circulation = iris.load_cube(filename, "mass_integrated_circulation")

    return (circulation / mass) * 1e6


if __name__ == "__main__":
    main()
