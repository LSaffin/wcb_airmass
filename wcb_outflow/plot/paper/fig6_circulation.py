"""
> python -m wcb_outflow.plot.paper.fig6_circulation

Circulation and mass vs time for the identified outflow regions in each case study
"""

import matplotlib.pyplot as plt

import iris

from wcb_outflow import case_studies, calc


def main():
    fig, axes = plt.subplots(2, 4, sharex="all", figsize=(16, 5))

    for n, case in enumerate(case_studies):
        print(case)
        case_study = case_studies[case]
        cubes = iris.load(str(case_study.data_path / "circulation.nc"))

        axes[0, n].set_title(case)

        for m, theta_level in enumerate(case_study.outflow_theta):
            cubes_theta = cubes.extract(iris.Constraint(
                air_potential_temperature=theta_level)
            )

            circ = cubes_theta.extract_strict("circulation")
            dt, outflow_time = calc.timediff(circ, case_study)

            color = "C{}".format(m)
            ax1 = plt.axes(axes[0, n])
            ax1.plot(dt, circ.data, color=color, label="{}K".format(theta_level))

            circ_m = cubes_theta.extract_strict("mass_integrated_circulation")
            dt, outflow_time = calc.timediff(circ_m, case_study)
            ax1.plot(dt, circ_m.data, color=color, linestyle="--", alpha=0.5)

            ax2 = plt.axes(axes[1, n])
            mass = cubes_theta.extract_strict("mass")
            ax2.plot(dt, mass.data, color=color, linestyle="--")

        _, ymax = axes[0, n].get_ylim()
        axes[0, n].set_ylim([0, ymax])
        _, ymax = axes[1, n].get_ylim()
        axes[1, n].set_ylim([0, ymax])

        axes[0, n].set_xlim(0, 90)
        axes[0, n].set_xticks(range(0, 90, 12))
        axes[1, n].set_xlim(0, 90)
        axes[1, n].set_xticks(range(0, 90, 12))

        axes[0, n].axvline(outflow_time, color="k")
        axes[1, n].axvline(outflow_time, color="k")

        axes[0, n].legend(loc="lower right")

    axes[0, 0].set_ylabel(r"Circulation (m$^2$s$^{-1}$)")
    axes[1, 0].set_ylabel(r"Mass (Kg)")

    fig.text(0.5, 0.01, "Forecast Lead Time (hours)", ha="center")

    plt.show()


if __name__ == '__main__':
    main()
