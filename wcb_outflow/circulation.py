"""Calculate integrals following an isentropic circuit

Calculate length and circulation around the isentropic circuit.
Calculate area/volume/mass/pv substance of the enclosed area of the isentropic
circuit.
"""

import datetime

import numpy as np
from matplotlib.path import Path

import iris
from iris.analysis import SUM, cartography

from irise import convert
from irise.constants import Omega

from pylagranto import trajectory

a = 6378100

base_time = datetime.datetime(1970, 1, 1)
base_time_units = "hours since 1970-01-01 00:00:00"


def calc_circulation(case, dtheta=1, region="outflow"):
    if region == "outflow":
        tr = trajectory.load(case.data_path / "isentropic_trajectories_backward.pkl") + \
             trajectory.load(case.data_path / "isentropic_trajectories_forward.pkl")
    elif region == "inflow":
        tr = trajectory.load(case.data_path / "isentropic_trajectories_from_inflow_forward.pkl")

    ntra, nt, ndim = tr.data.shape
    theta = np.unique(tr["air_potential_temperature"])

    time = [(t - base_time).total_seconds() // 3600 for t in tr.times]
    time = iris.coords.DimCoord(
        points=time,
        standard_name="time",
        units=base_time_units,
    )

    theta_levels = iris.coords.DimCoord(
        points=theta,
        standard_name="air_potential_temperature",
        units="K",
    )

    circ = iris.cube.Cube(
        data=np.empty([len(theta), nt]),
        long_name="circulation",
        units="m2 s-1",
        dim_coords_and_dims=[(theta_levels, 0), (time, 1)]
    )

    circ_r = circ.copy()
    circ_r.rename("relative_circulation")

    circ_p = circ.copy()
    circ_p.rename("planetary_circulation")

    results = iris.cube.CubeList([circ, circ_r, circ_p])

    # Repeat for all theta levels
    for j, theta_level in enumerate(theta):
        tr_theta = tr.select("air_potential_temperature", "==", theta_level)

        results_intermediate = iris.cube.CubeList()

        for n in range(nt):
            lon = tr_theta.x[:-2, n]
            lat = tr_theta.y[:-2, n]
            alt = tr_theta['altitude'][:-2, n]
            u = tr_theta['x_wind'][:-2, n]
            v = tr_theta['y_wind'][:-2, n]
            w = tr_theta['upward_air_velocity'][:-2, n]

            # Integrals are invalid once trajectories leave the domain but we don't
            # want to stop the script so put NaNs in the output instead
            if (alt == -1000).any():
                circ_r.data[j, n], circ_p.data[j, n], circ.data[j, n] = \
                    np.nan, np.nan, np.nan

            else:
                circ_r.data[j, n], circ_p.data[j, n], circ.data[j, n] = \
                    circuit_integrals(u, v, w, lon, lat, alt)

            # Calculate enclosed area integrals
            try:
                cubes = iris.load(
                    case.filename_theta(tr_theta.times[n], tracer_files=["c"]),
                    iris.Constraint(time=tr_theta.times[n])
                )
                print("{}K: {}".format(theta_level, tr_theta.times[n]))
            except OSError:
                print(str(tr_theta.times[n]) + " not available")
                break

            # Remove duplicate altitudes
            z = cubes.extract("altitude")
            cubes.remove(z[0])

            dlambda = np.deg2rad(np.diff(z[1].coord("longitude").points).mean())
            dphi = np.deg2rad(np.diff(z[1].coord("latitude").points).mean())

            integrals = mass_integrals(cubes, lon, lat, theta_level, dtheta, dlambda, dphi)

            for icube in integrals:
                # Set integrals to NaN if trajectories have left the domain
                if (alt == -1000).any():
                    icube.data = np.nan
                results_intermediate.append(icube)
        for cube in results_intermediate.merge():
            results.append(cube)

    iris.save(results.merge(), str(case.data_path / "circulation_{}.nc".format(region)))
    return


def circuit_integrals(u, v, w, lon, lat, z):
    """

    Args:
        u (numpy.ndarray):
        v (numpy.ndarray):
        w (numpy.ndarray):
        lon (numpy.ndarray):
        lat (numpy.ndarray):
        z (numpy.ndarray):

    Returns:

    """
    # Convert to radians
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    # u_planetary = \Omega r cos(lat)
    u_planetary = Omega.data * (a + z) * np.cos(lat)

    # Integrate dl around the circuit of trajectories
    dx, dy, dz = [], [], []
    for n in range(len(u)):
        # Allow a complete loop (back to zero at the end)
        np1 = (n + 1) % len(u)

        # dx = r cos(lat) dlon
        dx.append((a + z[n]) * np.cos(lat[n]) * 0.5 * (lon[np1] - lon[n - 1]))

        # dy = r dlat
        dy.append((a + z[n]) * 0.5 * (lat[np1] - lat[n - 1]))

        # dz is independent of grid rotation
        dz.append(0.5 * (z[np1] - z[n - 1]))

    dx = np.array(dx)
    dy = np.array(dy)
    dz = np.array(dz)

    # Circulation and separate components
    relative_circulation = np.sum(u*dx + v*dy + w*dz)
    planetary_circulation = np.sum(u_planetary * dx)
    circulation = relative_circulation + planetary_circulation

    return relative_circulation, planetary_circulation, circulation


def mass_integrals(cubes, lon, lat, theta_level, dtheta, dlambda, dphi):
    """

    Args:
        cubes (iris.cube.CubeList):
        lon (np.Array): Circuit longitudes (degrees)
        lat (np.Array): Circuit latitudes (degrees)
        theta_level: Isentropic level of the circuit
        dtheta (int): Isentrope spacing used to calculate volume integrals
        dlambda (float): Longitudinal grid spacing in radians
        dphi (float): Latitudinal grid spacing in radians

    Returns:
        iris.cube.Cube: A tuple containing four `iris.cube.Cube`s. Area, volume, mass
            and circulation calculated from the integrals.

    """
    # Get height at [theta - dtheta/2, theta, theta + dtheta/2]
    levels = ('air_potential_temperature',
              [theta_level - dtheta / 2.0, theta_level,
               theta_level + dtheta / 2.0])
    zth = convert.calc('altitude', cubes, levels=levels)

    # Extract model output grid
    glon, glat = cartography.get_xy_grids(zth)
    gridpoints = np.array([glon.flatten(), glat.flatten()]).transpose()

    # Mask all points that are not contained in the circuit
    points = np.array([lon, lat]).transpose()
    pth = Path(points)
    mask = np.logical_not(pth.contains_points(gridpoints).reshape(glat.shape))

    # Area = r**2 cos(phi) dlambda dphi
    area = (a + zth[1]) ** 2 * np.cos(np.deg2rad(glat)) * dlambda * dphi
    area.units = 'm^2'
    area.rename('area')
    total_area = integrate(area, mask)

    # Volume = area * dz
    volume = area * (zth[2] - zth[0])
    # theta coordinate disappears here
    volume.add_aux_coord(area.coord("air_potential_temperature"))
    volume.rename('volume')
    total_volume = integrate(volume, mask)

    # Mass = density * volume
    levels = ('air_potential_temperature', [theta_level])
    density = convert.calc('air_density', cubes, levels=levels)[0]
    mass = density * volume
    mass.rename('mass')
    total_mass = integrate(mass, mask)

    # Circulation = \int rho.pv.dv / dtheta
    pv = convert.calc('ertel_potential_vorticity', cubes, levels=levels)[0]
    pv.convert_units('m^2 K s-1 kg-1')
    pv_substance = pv * mass
    circulation = integrate(pv_substance, mask) / dtheta

    circulation.rename('mass_integrated_circulation')

    return total_area, total_volume, total_mass, circulation


def integrate(cube, mask):
    cube = cube.copy()
    cube.data = np.ma.masked_where(mask, cube.data)
    result = cube.collapsed(['longitude', 'latitude'], SUM)

    return result
