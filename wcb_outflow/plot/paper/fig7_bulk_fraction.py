"""
Timeseries of fractional differences in bulk quantities as laid out in eqn.25
But just current time quantity instead of notional adiabatic state
"""

import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
import numpy as np

from wcb_outflow import case_studies

# list of diagnostics to plot
diagnostics7 = [
    "mass",
    "area",
    "density"
]

def main(ylims = [-.3, .8]):
    
    fig, axes = plt.subplots(4, 3, figsize=(12, 16))
    
    for n, case_name in enumerate(case_studies):
        case_study = case_studies[case_name]
        cubes = iris.load(
            #'/home/users/bn826011/WCB_Leo_Git_clone/wcb_airmass/wcb_outflow/data/' 
            #+ case_name + '/circulation.nc'
            str(case_study.data_path / "circulation.nc")
            #iris.Constraint(time=[case_study.inflow_time, case_study.forecast_end_time])
        )
        
        for m, name in enumerate(diagnostics7):
            
            ax = axes[n, m]
            
            panel(ax, cubes, name, case_study, ylims)
            
            if m == 0:
                plt.ylabel(case_name, fontsize = 16)
            if n == 0:
                plt.title(name, fontsize = 16)
                
            
    #plt.suptitle(r'$\frac{\Delta X(t)}{X(0)}$', fontsize = 30)
    #plt.savefig('/home/users/bn826011/figure7_WCB_bulk_fractions')
            
            
def panel(ax, cubes, name, case_study, ylims):
    
    # don't know how to do ax.iplt or equivalent
    plt.sca(ax)

    if name == "density":
        cube = get_density_1Klayer(cubes)
    else:
        cube = cubes.extract_strict(name)
        
    for theta_level in case_study.outflow_theta:
        
        cube_theta = cube.extract(iris.Constraint(
            air_potential_temperature=theta_level
        ))
        # define earliest non-zero value
        cube_initial = cube_theta[~np.isnan(cube_theta.data)][0]
        
        # Fractional change from start to end, relative to current value
        #diff = ((cube_theta - cube_initial) / cube_theta)
        
        # Fractional change from start to end, relative to initial value
        diff = ((cube_theta - cube_initial) / cube_initial)
        
        iplt.plot(diff, zorder = 2, label = str(theta_level) + 'K')
        
    # plot outflow time
    ax.plot([case_study.outflow_time, case_study.outflow_time],
            ylims, 'k', zorder = 3)
    current_xlim = ax.get_xlim()
    # plot zero line
    ax.plot(current_xlim, [0, 0], 'k', zorder = 1)
    ax.set_ylim(ylims)
    plt.xticks(rotation = 25)
        
    if name == "density":
        plt.legend()
            
            
def get_density_1Klayer(cubes):
    mass = cubes.extract_strict("mass")
    area = cubes.extract_strict("area")
    
    return mass / area


if __name__ == '__main__':
    main()