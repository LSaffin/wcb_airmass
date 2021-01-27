import numpy as np

import iris

from irise import convert, interpolate
from pylagranto import caltra


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


def isentropic_trajectories(case, fbflag=-1, region="outflow"):
    # Get the outflow points for the selected theta levels
    trainp = np.load(str(case.data_path / "{}_boundaries.npy".format(region)))

    if region == "outflow":
        theta_levels = case.outflow_theta
        trainp = np.vstack(
            trainp[np.where(trainp[:, 2] == theta)] for theta in theta_levels
        )

    elif region == "inflow":
        theta_levels = list(range(280, case.outflow_theta[-1]+1, 5))

        trainp2 = np.zeros([trainp.shape[0], trainp.shape[1]+1])
        trainp2[:, :2] = trainp
        trainp2 = np.vstack([trainp2 for theta in theta_levels])

        npoints = trainp.shape[0]
        for n, theta in enumerate(theta_levels):
            trainp2[npoints*n:npoints*(n+1), 2] = theta
        trainp = trainp2

    levels = ("air_potential_temperature", theta_levels)

    if fbflag == 1:
        if region == "outflow":
            fname = "forward"
            mapping = case.time_to_filename_winds_mapping(
                start=case.outflow_time, end=case.forecast_end_time
            )
        elif region == "inflow":
            fname = "from_inflow_forward"
            mapping = case.time_to_filename_winds_mapping()
    elif fbflag == -1:
        fname = "backward"
        mapping = case.time_to_filename_winds_mapping()

    # Calculate the trajectories
    traout = caltra.caltra(
        trainp,
        mapping,
        nsubs=12,
        fbflag=fbflag,
        levels=levels,
        tracers=["x_wind", "y_wind", "upward_air_velocity", "altitude"]
    )

    # Save the trajectories
    traout.save(str(case.data_path / "isentropic_trajectories_{}.pkl".format(fname)))


def lagrangian_trajectories(case):
    # Get the outflow points for the selected theta levels
    trainp = np.load(str(case.data_path / "outflow_volume.npy"))[:, 0:3]

    # Calculate the trajectories
    traout = caltra.caltra(
        trainp,
        case.time_to_filename_winds_mapping(),
        nsubs=12,
        fbflag=-1,
        tracers=["air_potential_temperature"]
    )

    # Save the trajectories
    traout.save(str(case.data_path / "3d_trajectories.pkl"))


def isentropic_trajectories_from_volume(case):
    # Get the outflow points for the selected theta levels
    trainp = np.load(str(case.data_path / "outflow_volume.npy"))

    trainp_theta = np.zeros([trainp.shape[0], 3])
    trainp_theta[:, 0:2] = trainp[:, 0:2]
    trainp_theta[:, 2] = trainp[:, 3]

    levels = ("air_potential_temperature", case.outflow_theta)

    # Calculate the trajectories
    traout = caltra.caltra(
        trainp_theta,
        case.time_to_filename_winds_mapping(),
        nsubs=12,
        fbflag=-1,
        levels=levels,
        tracers=["x_wind", "y_wind", "upward_air_velocity", "altitude"]
    )

    # Save the trajectories
    traout.save(str(case.data_path / "isentropic_trajectories_from_volume.pkl"))
