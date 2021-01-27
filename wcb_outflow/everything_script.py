"""
Rerunning through everything in Jake's everything_script
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np

from wcb_outflow import case_studies, outflow, trajectories, inflow, circulation


def main():
    warnings.filterwarnings("ignore", category=UserWarning)

    # Define initial outflow volumes
    outflow.outflow_th(
        case_studies["IOP3"],
        theta_levels = np.arange(320, 335, 5),
        lon_bounds=(335, 400),
        lat_bounds=(30, 75),
        filtersize_t=25,
        filtersize_p=5,
        save=True,
    )

    outflow.outflow_th(
        case_studies["IOP5"],
        theta_levels=np.arange(315, 345, 5),
        lon_bounds=(320, 375),
        lat_bounds=(30, 65),
        filtersize_t=25,
        filtersize_p=5,
        save=True,
    )

    outflow.outflow_th(
        case_studies["IOP6"],
        theta_levels=np.arange(305, 330, 5),
        lon_bounds=(310, 365),
        lat_bounds=(30, 75),
        filtersize_t=25,
        filtersize_p=5,
        save=True,
    )

    outflow.outflow_th(
        case_studies["IOP7"],
        theta_levels=np.arange(310, 350, 5),
        lon_bounds=(280, 360),
        lat_bounds=(30, 65),
        filtersize_t=25,
        filtersize_p=5,
        save=True,
    )

    plt.show()

    # Trajectory calculations
    for case_name in case_studies:
        case = case_studies[case_name]
        trajectories.preprocess_winds(case)

        # Isentropic circuit trajectories
        # Backward
        trajectories.isentropic_trajectories(case, fbflag=-1)
        # and forward
        trajectories.isentropic_trajectories(case, fbflag=1)

        # 3d back trajectories
        trajectories.lagrangian_trajectories(case)

        # Isentropic trajectories from the same points as 3d trajectories
        trajectories.isentropic_trajectories_from_volume(case)

        # Circulation
        circulation.calc_circulation(case)

        # Inflow regions
        inflow.from_3d_trajectories(case)

        # Isentropic trajectories from the inflow
        trajectories.isentropic_trajectories(case, fbflag=1, region="inflow")
        circulation.calc_circulation(case, region="inflow")


if __name__ == "__main__":
    main()
