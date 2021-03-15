"""
Show the outflow region as function of time when tracked backwards with the isentropic
trajectories. It is useful to see how much the contour breaks down over time

Usage:
    outflow_tracking.py <case> <theta_level>
    outflow_tracking.py  (-h | --help)

Arguments:
    <case>         Name of case study
    <theta_level>  Which outflow theta level to show (K)

Options:
    -h --help      Show help
"""
import docopt

import matplotlib.pyplot as plt

import iris

from pylagranto import trajectory

from wcb_outflow import case_studies
from wcb_outflow.plot import background_map, projection


def main():
    arguments = docopt.docopt(__doc__)

    case = case_studies[arguments["<case>"]]
    theta_level = float(arguments["<theta_level>"])

    fig = plt.figure(figsize=[8, 5])

    ntimes = (case.outflow_time - case.start_time) // case.timestep
    times = [case.start_time + n*case.timestep for n in range(ntimes+1)]

    tr = trajectory.load(case.data_path / "isentropic_trajectories_backward.pkl")
    tr_theta = tr.select(
        "air_potential_temperature", "==", theta_level,
        time=[tr.relative_times[0]]
    )

    tr3d = trajectory.load(case.data_path / "3d_trajectories.pkl")
    tr3d = tr3d.select(
        "air_potential_temperature", ">", case.outflow_theta[0] - 1,
        time=[tr.relative_times[0]]
    )

    tr3d = tr3d.select(
        "air_potential_temperature", "<", case.outflow_theta[-1] + 1,
        time=[tr.relative_times[0]]
    )

    for n, time in enumerate(times):
        print(time)
        pv = iris.load_cube(
            case.filename_theta(time=time),
            iris.Constraint(
                name="ertel_potential_vorticity",
                time=time,
                air_potential_temperature=theta_level,
            )
        )
        crs = pv.coord(axis="x").coord_system.as_cartopy_crs()
        ax = plt.subplot(ntimes // 2, 3, n+1, projection=projection)
        background_map(ax, pv)

        idx = tr.times.index(time)
        plt.plot(tr_theta.x[:, idx] - 360, tr_theta.y[:, idx], transform=crs)
        tr3d_inflow = tr3d.select(
            "air_potential_temperature", "<", case.outflow_theta[0] - 5,
            time=[tr.relative_times[idx]]
        )

        plt.scatter(tr3d_inflow.x[:, idx] - 360, tr3d_inflow.y[:, idx], c=tr3d_inflow.z[:, idx]/10000, transform=crs, marker=".")
        ax.set_title(time)
    fig.suptitle("{}: {:.0f}K".format(case.name, theta_level))
    plt.show()


if __name__ == '__main__':
    main()
