import matplotlib.pyplot as plt

import iris

from wcb_outflow import case_studies, calc


def main(diagnostics, ylims=(-1, 1), colours=['b', 'r', 'm', 'g']):
    ldg = len(diagnostics) // 2
    fig, axes = plt.subplots(2, ldg, figsize=(ldg*3 + 1, 10), sharex="all")

    for m, name in enumerate(diagnostics):
        ax = axes[m // ldg, m % ldg]
        panel(ax, name, ylims, colours)

    # plt.suptitle(r'$\frac{\Delta X(t)}{X(0)}$', fontsize = 30)
    axes[0, 0].legend()
    plt.show()


def panel(ax, name, ylims, colours):
    for n, case_name in enumerate(case_studies):
        case_study = case_studies[case_name]
        colour = colours[n]
        cubes = iris.load(str(case_study.data_path / "circulation.nc"))

        # pick the middle theta level
        theta_level = case_study.outflow_theta[1]

        cubes_theta = cubes.extract(iris.Constraint(
            air_potential_temperature=theta_level
        ))

        if name == "mass area anomaly ratio":
            diff = calc.ratio(cubes_theta)
            ylims = (0, 3)
            x_line = 1

        else:
            x_line = 0

            cube_theta = calc.get_cube_by_name(cubes_theta, name)

            # Fractional change from start to end, relative to initial value
            diff = calc.diff(cube_theta)

        # time elapsed from initial time
        timediff, out_time = calc.timediff(diff, case_study)

        ax.plot(timediff, diff.data, zorder=2, color=colour,
                label=case_name + ' ' + str(theta_level) + 'K')

        # plot outflow times
        ax.axvline(out_time, linestyle='-.', color=colour)

    ax.set_ylim(ylims)
    ax.axhline(x_line, color="k")
    ax.set_title(name)


if __name__ == '__main__':
    main(['mass', 'area', 'density', 'mass area anomaly ratio'])
    # main(['mass', 'area', 'density', 'mass area anomaly ratio',
    #       'PV', 'vorticity', 'circulation'])
