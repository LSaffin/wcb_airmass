"""Generate a table showing percentage changes in quantities in the outflow

Creates a table that can be pasted in to a LaTeX document
"""

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
        cubes = iris.load(
            str(case_study.data_path / "circulation.nc"),
            iris.Constraint(time=[case_study.inflow_time, case_study.outflow_time])
        )

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
                diff = ((cube[-1] - cube[0]) / cube[0]) * 100
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


if __name__ == '__main__':
    main()

