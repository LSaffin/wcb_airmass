import matplotlib.pyplot as plt

import iris
import iris.plot as iplt

from wcb_outflow import case_studies


def main():
    fig, axes = plt.subplots(3, 4, figsize=(16, 10), sharex="col", sharey="all")

    for i, case in enumerate(case_studies):
        case_study = case_studies[case]

        # Load trajectories
        circulation = iris.load(str(case_study.data_path / "circulation_outflow.nc"))

        # Repeat for all theta levels
        for j, theta_level in enumerate(case_study.outflow_theta):
            plt.axes(axes[j, i])

            cs = iris.Constraint(air_potential_temperature=theta_level)
            subcubes = circulation.extract(cs)

            for cube in subcubes:
                if "circ" in cube.name():
                    iplt.plot(cube, label=cube.name())

            axes[j, i].axhline(color="k")
            axes[j, i].axvline(case_study.outflow_time, color="k")
            axes[j, i].set_title("{} ({} K)".format(case_study.name, theta_level))

    fig.autofmt_xdate()
    axes[0, 0].legend()
    plt.show()

    return


if __name__ == '__main__':
    main()
