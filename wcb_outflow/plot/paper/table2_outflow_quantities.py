"""Generate a table showing percentage changes in quantities in the outflow

Creates a table that can be pasted in to a LaTeX document
"""

import datetime

import iris

from wcb_outflow import case_studies, calc
from wcb_outflow.plot.paper.fig7_bulk_fraction import time_zero


# List of variables to include in the table
diagnostics = [
    "circulation",
    "mass_integrated_circulation",
    "mass",
    "area",
    "mass_integrated_vorticity",
    "planetary_vorticity",
    "relative_vorticity",
]


def main():
    # Create the table header
    print(" & ".join(
        ["Case"] + [r"$\Delta$ t"] + ["Isentropic Level (K)"] + diagnostics
    ))
    for case in case_studies:
        case_study = case_studies[case]
        cubes = iris.load(str(case_study.data_path / "circulation_outflow.nc"))

        dt = case_study.outflow_time - datetime.datetime(1970, 1, 1)
        idx_outflow = (
                dt.total_seconds() // 3600 == cubes[2].coord("time").points
        ).argmax()

        results = [case_study.name]

        for theta_level in case_study.outflow_theta:
            cubes_theta = cubes.extract(iris.Constraint(
                air_potential_temperature=theta_level
            ))

            # Prepend the time difference and theta level to each line of results
            idx_inflow = time_zero[case][theta_level]
            results.append("{:0d}h".format((idx_outflow - idx_inflow)*6))
            results.append("{}".format(theta_level))

            for name in diagnostics:
                cube = calc.get_cube_by_name(cubes_theta, name)

                if "vorticity" in name:
                    if name == "mass_integrated_vorticity":
                        vort_0 = cube[idx_inflow]

                    if name == "vorticity" or name == "mass_integrated_vorticity":
                        diff = calc.diff(cube, idx_0=idx_inflow) * 100
                    else:
                        diff = ((cube - cube[idx_inflow]) / vort_0) * 100

                else:
                    # Percentage change from start to end, relative to initial value
                    diff = calc.diff(cube, idx_0=idx_inflow) * 100

                try:
                    # Round to 2 decimal places to make the table neater
                    results.append(round(diff.data[idx_outflow], 1))
                except TypeError:
                    # If any NaNs are present, round doesn't work because the
                    # calculation of derived variables (i.e. vorticity) creates a masked
                    # array
                    results.append("--")

            print(" & ".join([str(x) for x in results]) + r" \\")
            results = [""]


if __name__ == '__main__':
    main()
