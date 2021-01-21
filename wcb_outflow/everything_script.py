"""
Rerunning through everything in Jake's everything_script
"""

import warnings

import numpy as np

from wcb_outflow import case_studies, outflow, trajectories, inflow, circulation


def main():
    # Define initial outflow volumes
    theta_levels = np.arange(300, 350, 5)
    outflow.outflow_th(
        case_studies["IOP3"],
        theta_levels,
        lon_bounds=(335, 400),
        lat_bounds=(40, 75),
        filtersize_t=25,
        filtersize_p=1,
        save=True,
    )

    outflow.outflow_th(
        case_studies["IOP5"],
        theta_levels,
        lon_bounds=(300, 375),
        lat_bounds=(20, 65),
        filtersize_t=30,
        filtersize_p=30,
    )
    outflow.outflow_th(
        case_studies["IOP6"],
        theta_levels[1:],
        lon_bounds=(310, 365),
        lat_bounds=(40, 75),
        filtersize_t=30,
        filtersize_p=30,
    )
    outflow.outflow_th(
        case_studies["IOP7"],
        theta_levels,
        lon_bounds=(280, 350),
        lat_bounds=(30, 65),
        filtersize_t=30,
        filtersize_p=5,
    )

    # Trajectory calculations
    for case_name in case_studies:
        trajectories.preprocess_winds(case_studies[case_name])

        # Isentropic circuit trajectories
        # Backward
        trajectories.isentropic_trajectories(case_studies[case_name])
        # and forward
        trajectories.isentropic_trajectories(case_studies[case_name], fbflag=1)

        # 3d back trajectories
        trajectories.lagrangian_trajectories(case_studies[case_name])

        # Isentropic trajectories from the same points as 3d trajectories
        trajectories.isentropic_trajectories_from_volume(case_studies[case_name])

    # Inflow regions
    inflow.from_3d_trajectories(case_studies["IOP3"])
    inflow.from_3d_trajectories(case_studies["IOP5"])
    inflow.from_3d_trajectories(case_studies["IOP6"])
    inflow.from_3d_trajectories(case_studies["IOP7"])

    outflow.at_inflow_time(case_studies["IOP3"])
    outflow.at_inflow_time(case_studies["IOP5"])
    outflow.at_inflow_time(case_studies["IOP6"])
    outflow.at_inflow_time(case_studies["IOP7"])

    # Circulation
    warnings.filterwarnings("ignore", category=UserWarning)
    circulation.calc_circulation(case_studies["IOP3"])
    circulation.calc_circulation(case_studies["IOP5"])
    circulation.calc_circulation(case_studies["IOP6"])
    circulation.calc_circulation(case_studies["IOP7"])


if __name__ == "__main__":
    main()
