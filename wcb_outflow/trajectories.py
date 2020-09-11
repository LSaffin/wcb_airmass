import numpy as np

import iris

from irise import convert, interpolate
from pylagranto import caltra

from . import data_path


def preprocess_winds(case):
    """Interpolate wind fields to a common grid ready for trajectory calculations
    """
    # Don't load pressure on model theta levels
    pressure_cs = iris.Constraint(
        cube_func=lambda x: str(x.attributes["STASH"]) != "m01s00i407")

    mapping = case.time_to_filename_mapping()

    for time in mapping:
        print(time)
        cubes = iris.load(mapping[time], iris.Constraint(time=time) & pressure_cs)

        w = convert.calc("upward_air_velocity", cubes)
        theta = convert.calc("air_potential_temperature", cubes)

        newcubes = iris.cube.CubeList([w, theta])
        for variable in ["x_wind", "y_wind"]:
            cube = convert.calc(variable, cubes)
            cube = interpolate.remap_3d(cube, w)
            newcubes.append(cube)

        iris.save(newcubes, case.filename_winds(time))


def isentropic_trajectories(case):
    # Get the outflow points for the selected theta levels
    trainp = np.load(str(case.data_path / "outflow_boundaries.npy"))
    trainp = np.vstack(
        trainp[np.where(trainp[:, 2] == theta)] for theta in case.theta_levels
    )

    levels = ("air_potential_temperature", case.theta_levels)

    # Calculate the trajectories
    traout = caltra.caltra(
        trainp,
        case.time_to_filename_winds_mapping(),
        nsubs=12,
        fbflag=-1,
        levels=levels,
        tracers=["x_wind", "y_wind", "upward_air_velocity", "altitude"]
    )

    # Save the trajectories
    traout.save(str(data_path / case / "isentropic_trajectories.pkl"))


def lagrangian_trajectories(case):
    # Get the outflow points for the selected theta levels
    trainp = np.load(str(case.data_path / "outflow_volume.npy"))

    # Calculate the trajectories
    traout = caltra.caltra(
        trainp,
        case.time_to_filename_winds_mapping(),
        nsubs=12,
        fbflag=-1,
        tracers=["air_potential_temperature"]
    )

    # Save the trajectories
    traout.save(str(data_path / case / "3d_trajectories.pkl"))
