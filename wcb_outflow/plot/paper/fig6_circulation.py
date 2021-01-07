"""
> python -m wcb_outflow.plot.paper.fig6_circulation

Circulation and mass vs time for the identified outflow regions in each case study
"""

import numpy as np
import matplotlib.pyplot as plt

import iris

from wcb_outflow import case_studies


def main():
    fig, axes = plt.subplots(4, 2, sharex="all", figsize=(8, 10))

    for n, case in enumerate(case_studies):
        print(case)
        case_study = case_studies[case]
        cubes = iris.load(str(case_study.data_path / "circulation.nc"))

        for m, theta_level in enumerate(case_study.outflow_theta):
            print(theta_level)
            cubes_theta = cubes.extract(iris.Constraint(
                air_potential_temperature=theta_level)
            )

            circ = cubes_theta.extract_strict("circulation")
            dt = circ.coord("time")
            dt.convert_units("Hours Since {}".format(case_study.start_time.strftime("%Y-%m-%d %H:%M:%S")))

            plt.axes(axes[n, 0])
            plt.plot(dt.points, circ.data, color="C{}".format(m))

            circ_m = cubes_theta.extract_strict("mass_integrated_circulation")
            dt = circ_m.coord("time")
            dt.convert_units("Hours Since {}".format(case_study.start_time.strftime("%Y-%m-%d %H:%M:%S")))
            plt.plot(dt.points, circ_m.data, color="C{}".format(m), linestyle="--", alpha=0.5)

            plt.axes(axes[n, 1])
            mass = cubes_theta.extract_strict("mass")
            print(np.diff(mass.data))
            plt.plot(dt.points, mass.data, color="C{}".format(m), linestyle="--", label="{}K".format(theta_level))

        _, ymax = axes[n, 0].get_ylim()
        axes[n, 0].set_ylim([0, ymax])
        _, ymax = axes[n, 1].get_ylim()
        axes[n, 1].set_ylim([0, ymax])
        axes[n, 1].yaxis.set_label_position("right")
        axes[n, 1].yaxis.tick_right()

        axes[n, 0].set_xlim(0, 90)
        axes[n, 0].set_xticks(range(0, 90, 12))
        axes[n, 1].set_xlim(0, 90)
        axes[n, 1].set_xticks(range(0, 90, 12))

        axes[n, 0].axvline(case_study.outflow_lead_time.total_seconds() / 3600, color="k")
        axes[n, 1].axvline(case_study.outflow_lead_time.total_seconds() / 3600, color="k")

        axes[n, 1].legend(loc="lower right")
        axes[n, 0].set_ylabel(case_study.name)

    fig.text(0.05, 0.5, r"Circulation (m$^2$s$^{-1}$)", rotation="vertical", va="center")
    fig.text(0.95, 0.5, r"Mass (Kg)", rotation="vertical", va="center")

    fig.text(0.5, 0.05, "Forecast Lead Time (hours)", ha="center")

    plt.show()


if __name__ == '__main__':
    main()
