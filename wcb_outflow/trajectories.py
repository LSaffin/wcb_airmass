import numpy as np

from pylagranto import caltra

from . import case_studies, data_path


def isentropic_trajectories(case, theta_levels):
    # Set up the mapping of times to files:
    mapping = case_studies[case].time_to_filename_mapping()

    # Get the outflow points for the selected theta levels
    trainp = np.load(str(case_studies[case].data_path / "outflow_boundaries.npy"))
    trainp = np.vstack(
        trainp[np.where(trainp[:, 2] == theta)] for theta in theta_levels
    )

    levels = ("air_potential_temperature", theta_levels)

    # Calculate the trajectories
    traout = caltra.caltra(trainp, mapping, nsubs=12, fbflag=-1, levels=levels)

    # Save the trajectories
    traout.to_pickle(str(data_path / case / "isentropic_trajectories.pkl"))


def lagrangian_trajectories(case, theta_levels):
    # Set up the mapping of times to files:
    mapping = case_studies[case].time_to_filename_mapping()

    # Get the outflow points for the selected theta levels
    trainp = np.load(str(case_studies[case].data_path / "outflow_volume.npy"))
    trainp = np.vstack(
        trainp[np.where(trainp[:, 2] == theta)] for theta in theta_levels
    )

    # Calculate the trajectories
    traout = caltra.caltra(trainp, mapping, nsubs=12, fbflag=-1)

    # Save the trajectories
    traout.to_pickle(str(data_path / case / "3d_trajectories.pkl"))
