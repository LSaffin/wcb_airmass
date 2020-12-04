"""
Timeseries of fractional differences in bulk quantities as laid out in eqn.25
But just current time quantity instead of notional adiabatic state
"""

import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
import numpy as np
import datetime
from wcb_outflow import case_studies

# list of diagnostics to plot
diagnostics_all = [
    "mass",
    "area",
    "density",
    "mass area anomaly ratio",
    "circulation",
    "relative_circulation",
    "planetary_circulation",
    "mass_integrated_circulation",
    "PV",
    "relative_PV",
    "planetary_PV",
    "mass_integrated_PV",
    "vorticity",
    "relative_vorticity",
    "planetary_vorticity",
    "mass_integrated_vorticity"
]

def main(diagnostics, ylims = [-1, 1]):
    
    ldg = len(diagnostics)
    
    fig, axes = plt.subplots(4, ldg, figsize=(ldg*2 + 4, 16))
    
    for n, case_name in enumerate(case_studies):
        case_study = case_studies[case_name]
        cubes = iris.load(
            #'/home/users/bn826011/WCB_Leo_Git_clone/wcb_airmass/wcb_outflow/data/' 
            #+ case_name + '/circulation.nc'
            str(case_study.data_path / "circulation.nc")
            #iris.Constraint(time=[case_study.inflow_time, case_study.forecast_end_time])
        )
        
        for m, name in enumerate(diagnostics):
            
            ax = axes[n, m]
            
            panel(ax, cubes, name, case_study, ylims)
            
            if m == 0:
                plt.ylabel(case_name, fontsize = 16)
            if n == 0:
                plt.title(name, fontsize = 16)
                
            
    #plt.suptitle(r'$\frac{\Delta X(t)}{X(0)}$', fontsize = 30)
    plt.savefig('/home/users/bn826011/figure7_WCB_bulk_fractions_new')
            
            
def panel(ax, cubes, name, case_study, ylims):
    
    # don't know how to do ax.iplt or equivalent
    plt.sca(ax)

    if name == "density":
        cube = get_density_1Klayer(cubes)
    elif "vorticity" in name:
        cube = get_vorticity(cubes, name[:-9])
    elif "PV" in name:
        cube = get_PV(cubes, name[:-2])
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
        
        
def main_alt(diagnostics, ylims = (-1, 1), colours = ['b', 'r', 'm', 'g']):
    
    ldg = len(diagnostics) // 2
    
    fig, axes = plt.subplots(2, ldg, figsize=(ldg*3 + 1, 10))
        
    for m, name in enumerate(diagnostics):
            
        ax = axes[m // ldg, m % ldg]
            
        panel_alt(ax, name, ylims, colours)
                
            
    #plt.suptitle(r'$\frac{\Delta X(t)}{X(0)}$', fontsize = 30)
    plt.savefig('/home/users/bn826011/figure7_WCB_bulk_fractions_fewerpanels.png')
            
            
def panel_alt(ax, name, ylims, colours):
    
    # don't know how to do ax.iplt or equivalent
    #plt.sca(ax)
    
    for n, case_name in enumerate(case_studies):
        case_study = case_studies[case_name]
        colour = colours[n]
        cubes = iris.load(
            #'/home/users/bn826011/WCB_Leo_Git_clone/wcb_airmass/wcb_outflow/data/' 
            #+ case_name + '/circulation.nc'
            str(case_study.data_path / "circulation.nc")
            #iris.Constraint(time=[case_study.inflow_time, case_study.forecast_end_time])
        )
        
        # pick the middle theta level
        theta_level = case_study.outflow_theta[1]
        
        cubes_theta = cubes.extract(iris.Constraint(
            air_potential_temperature=theta_level
        ))

        if name == "mass area anomaly ratio":
            
           cube_theta = get_ratio(cubes_theta)
           diff = cube_theta
           ylims = (0, 3)
           x_line= [1, 1]
              
        else:
            
            x_line = [0, 0]
            
            if name == "density":
                cube_theta = get_density_1Klayer(cubes_theta)
                #ylims = (-.2, .4)
            elif "vorticity" in name:
                cube_theta = get_vorticity(cubes_theta, name[:-9])
            elif "PV" in name:
                cube_theta = get_PV(cubes_theta, name[:-2])
            else:
                cube_theta = cubes_theta.extract_strict(name)
    
            # Fractional change from start to end, relative to initial value
            diff = get_diff(cube_theta)
            
        # time elapsed from initial time
        timediff, out_time = get_timediff(cube_theta, case_study)
            
        ax.plot(timediff, diff.data, zorder = 2, color = colour,
                label = case_name + ' ' + str(theta_level) + 'K')        
        
        # plot outflow times
        ax.plot([out_time, out_time],
                ylims, '-.', zorder = 3, color = colour)
        
    current_xlim = ax.get_xlim()
    # plot zero line
    ax.plot(current_xlim, x_line, 'k', zorder = 1)
    ax.set_ylim(ylims)
    #plt.xticks(rotation = 25)
    ax.set_title(name)
        
    #if name == "density":
    plt.legend()
            
            
def get_density_1Klayer(cubes):
    mass = cubes.extract_strict("mass")
    area = cubes.extract_strict("area")
    
    return mass / area

def get_vorticity(cubes, circulation_name = ""):
    area = cubes.extract_strict("area")
    circulation = cubes.extract_strict(circulation_name + "circulation")
    
    area, circulation = equate_times(area.copy(), circulation.copy())
    
    vort = circulation/area
    
    return vort

def get_PV(cubes, circulation_name = ""):
    mass = cubes.extract_strict("mass")
    circulation = cubes.extract_strict(circulation_name + "circulation")
    mass, circulation = equate_times(mass.copy(), circulation.copy())
    
    pv = circulation/mass
    
    return pv

def get_diff(cube_theta):
    # define earliest non-zero value
    if not np.isnan(cube_theta.data).all():
        cube_initial = cube_theta[~np.isnan(cube_theta.data)][0]
    else:
        cube_initial = cube_theta[0]
    
    # Fractional change from start to end, relative to current value
    #diff = ((cube_theta - cube_initial) / cube_theta)
    
    # Fractional change from start to end, relative to initial value
    diff = ((cube_theta - cube_initial) / cube_initial)
    
    return diff

def get_ratio(cubes):
    mass = cubes.extract_strict("mass")
    area = cubes.extract_strict("area")
    alpha = get_diff(mass)
    epsilon = get_diff(area)
    return alpha/epsilon

def get_timediff(cube, case_study):
    
    times = cube.coord('time').points
    
    if not np.isnan(cube.data).all():
        initial_time = times[~np.isnan(cube.data)][0]
    else:
        initial_time = times[0]
    
    out_time = case_study.outflow_time
    oftd = out_time - datetime.datetime(1970, 1, 1)
    outime_hrs = oftd.days*24 + oftd.seconds/3600
    out_time_rel = outime_hrs - initial_time
    
    return (times - initial_time)/24, out_time_rel/24


def equate_times(cube, reference_cube):
    
    ref_dates = [cell.point for cell in reference_cube.coord('time').cells()]
    ref_cst = iris.Constraint(time = lambda cell: cell in ref_dates)
    cube = cube.extract(ref_cst)
    cube_dates = [cell.point for cell in cube.coord('time').cells()]
    eqtim_cst = iris.Constraint(time = lambda cell: cell in cube_dates)
    reference_cube = reference_cube.extract(eqtim_cst)
    
    new_cube = reference_cube.copy()
    
    new_cube.rename(cube.name)
    new_cube.units = cube.units
    new_cube.data= cube.data
    
    return new_cube, reference_cube
    
    

if __name__ == '__main__':
    main(['PV', 'vorticity', 'density'])
    #main(['mass', 'area', 'PV', 'vorticity', 'density'])
    #main_alt(['mass', 'area', 'density', 'mass area anomaly ratio'])
    main_alt(['mass', 'area', 'density', 'mass area anomaly ratio', 'PV', 'vorticity', 'circulation'])
