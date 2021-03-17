"""
> python -m wcb_outflow.plot.paper.fig6_circulation

Circulation and mass vs time for the identified outflow regions in each case study
"""

import matplotlib.pyplot as plt

import iris

from wcb_outflow import case_studies, calc
from wcb_outflow.plot import set_plot_rcparams


def main():
    set_plot_rcparams()
    fig, axes = plt.subplots(4, 2, sharex="all", figsize=(8, 10))

    for n, case in enumerate(case_studies):
        i = (n // 2) * 2
        j = n % 2

        ax1 = axes[i, j]
        ax2 = axes[i + 1, j]
        ax1.set_title(case)

        print(case)
        case_study = case_studies[case]
        cubes = iris.load(str(case_study.data_path / "circulation_outflow.nc"))

        for m, theta_level in enumerate(case_study.outflow_theta):
            cubes_theta = cubes.extract(iris.Constraint(
                air_potential_temperature=theta_level)
            )

            circ = cubes_theta.extract_cube("circulation")
            dt, outflow_time = calc.timediff(circ, case_study)

            color = "C{}".format(m)
            ax1.plot(dt, circ.data, color=color, label="{}K".format(theta_level))

            circ_m = cubes_theta.extract_cube("mass_integrated_circulation")
            dt, outflow_time = calc.timediff(circ_m, case_study)
            ax1.plot(dt, circ_m.data, color=color, linestyle="--", alpha=0.5)

            mass = cubes_theta.extract_cube("mass")
            ax2.plot(dt, mass.data, color=color, linestyle="--")

        _, ymax = ax1.get_ylim()
        ax1.set_ylim([0, ymax])
        _, ymax = ax2.get_ylim()
        ax2.set_ylim([0, ymax])

        ax1.axvline(outflow_time, color="k")
        ax2.axvline(outflow_time, color="k")

        ax1.legend(loc="lower right")

    axes[0, 0].set_xlim(0, 72)
    axes[0, 0].set_xticks(range(0, 73, 12))

    axes[0, 0].set_ylabel(r"Circulation (m$^2$s$^{-1}$)")
    axes[1, 0].set_ylabel(r"Mass (Kg)")
    axes[2, 0].set_ylabel(r"Circulation (m$^2$s$^{-1}$)")
    axes[3, 0].set_ylabel(r"Mass (Kg)")

    fig.text(0.5, 0.01, "Forecast Lead Time (hours)", ha="center")

    plt.show()


if __name__ == '__main__':
    main()
