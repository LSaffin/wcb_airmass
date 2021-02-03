"""
Timeseries of fractional differences in bulk quantities as laid out in eqn.25
But just current time quantity instead of notional adiabatic state
"""

import matplotlib.pyplot as plt

import iris

from wcb_outflow import case_studies, calc

# list of diagnostics to plot
diagnostics_all = [
    "mass",
    "area",
    "density",
    "mass area anomaly ratio",
    "circulation",
    "relative_circulation",
    "planetary_circulation",
    "mass_integrated_circulation",
    "PV",
    "relative_PV",
    "planetary_PV",
    "mass_integrated_PV",
    "vorticity",
    "relative_vorticity",
    "planetary_vorticity",
    "mass_integrated_vorticity"
]


time_zero = dict(
    IOP3={320: 1, 325: 2, 330: 4},
    IOP5={325: 1, 330: 1, 335: 4},
    IOP6={310: 3, 315: 3, 320: 2},
    IOP7={310: 0, 315: 0, 320: 1}
)


def main(diagnostics):
    ldg = len(diagnostics)

    fig, axes = plt.subplots(ldg, 4, figsize=(16, 10), sharex="all", sharey="row")

    for n, case_name in enumerate(case_studies):
        case_study = case_studies[case_name]
        cubes = iris.load(str(case_study.data_path / "circulation.nc"))

        for m, name in enumerate(diagnostics):
            ax = plt.axes(axes[m, n])
            if name == "mass" or name == "area":
                ylims = [0, 1]
            else:
                ylims = [-1, 1]
            panel(ax, cubes, name, case_study, ylims)

            if m == 0:
                plt.title(case_name)
                plt.legend()
            if n == 0:
                plt.ylabel(name)

    axes[0, 0].set_xlim(0, 90)
    axes[0, 0].set_xticks(range(0, 90, 12))
    # plt.suptitle(r'$\frac{\Delta X(t)}{X(0)}$', fontsize = 30)
    plt.show()


def panel(ax, cubes, name, case_study, ylims):
    cube = calc.get_cube_by_name(cubes, name)
    dt, outflow_time = calc.timediff(cube, case_study)

    for theta_level in case_study.outflow_theta:
        cube_theta = cube.extract(iris.Constraint(
            air_potential_temperature=theta_level
        ))

        idx_0 = time_zero[case_study.name][theta_level]

        # Fractional change from start to end, relative to initial value
        diff = calc.diff(cube_theta, idx_0)

        ax.plot(dt, diff.data, label=str(theta_level) + 'K')

    # plot outflow time
    ax.axvline(outflow_time, color="k")

    # plot zero line
    ax.axhline(color="k", lw=1)
    ax.set_ylim(ylims)


if __name__ == '__main__':
    main(['mass', 'area', 'PV', 'vorticity', 'density'])
