"""
Rerunning through everything in Jake's everything_script
"""

import warnings

import numpy as np

from . import case_studies, outflow, trajectories, circulation


def main():
    # Define initial outflow volumes
    theta_levels = np.arange(300, 340, 5)
    outflow.outflow_th(case_studies["IOP3"], theta_levels=theta_levels)
    outflow.outflow_th(case_studies["IOP5"], theta_levels=theta_levels)
    outflow.outflow_th(case_studies["IOP6"], theta_levels=theta_levels[1:])
    outflow.outflow_th(case_studies["IOP7"], theta_levels=theta_levels)

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

    # Circulation
    warnings.filterwarnings("ignore", category=UserWarning)
    circulation.calc_circulation(case_studies["IOP3"])
    circulation.calc_circulation(case_studies["IOP5"])
    circulation.calc_circulation(case_studies["IOP6"])
    circulation.calc_circulation(case_studies["IOP7"])


if __name__ == "__main__":
    main()
