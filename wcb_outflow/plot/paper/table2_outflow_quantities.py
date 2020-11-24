"""Generate a table showing percentage changes in quantities in the outflow

Creates a table that can be pasted in to a LaTeX document
"""

import datetime

import numpy as np
import iris

from wcb_outflow import case_studies


# List of variables to include in the table
diagnostics = [
    "circulation",
    "mass_integrated_circulation",
    "mass",
    "area",
    "vorticity"
]


def main():
    # Create the table header
    print(" & ".join(["Case"] + ["Isentropic Level (K)"] + diagnostics))
    for case in case_studies:
        case_study = case_studies[case]
        cubes = iris.load(str(case_study.data_path / "circulation.nc"))

        dt = case_study.outflow_time - datetime.datetime(1970, 1, 1)
        idx_outflow = (dt.total_seconds() // 3600 == cubes[2].coord("time").points).argmax()

        for theta_level in case_study.outflow_theta:
            # Prepend the case study name and theta level to each line of results
            results = [case_study.name, "{}".format(theta_level)]
            for name in diagnostics:
                cubes_theta = cubes.extract(iris.Constraint(
                    air_potential_temperature=theta_level
                ))
                if name == "vorticity":
                    cube = get_vorticity(cubes_theta)
                else:
                    cube = cubes_theta.extract_strict(name)

                # Percentage change from start to end, relative to initial value
                diff = percentage_increase(cube, idx_outflow)
                try:
                    # Round to 2 decimal places to make the table neater
                    results.append(round(diff.data[()], 2))
                except TypeError:
                    # If any NaNs are present, round doesn't work because the
                    # calculation of derived variables (i.e. vorticity) creates a masked
                    # array
                    results.append("--")

            print(" & ".join([str(x) for x in results]) + r" \\")


def get_vorticity(cubes):
    circulation = cubes.extract_strict("mass_integrated_circulation")
    mass = cubes.extract_strict("area")

    return circulation / mass


def percentage_increase(cube, idx_outflow):
    # Use the earliest non-NaN values
    idx_inflow = np.isfinite(cube.data).argmax()

    return ((cube[idx_outflow] - cube[idx_inflow]) / cube[idx_inflow]) * 100


if __name__ == '__main__':
    main()

