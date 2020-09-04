# -*- coding: utf-8 -*-
"""
A copy of stuff from circuit_integral, to produce plots of circulation split
"""

import numpy as np
import iris
from matplotlib.path import Path
from scripts.circuit_integral import integrate
from mymodule import convert

import matplotlib.pyplot as plt
import iris.plot as iplt

from matplotlib.dates import DayLocator, HourLocator, DateFormatter


a = 6371229 # JB changed from 6378100

def circuit_integral_split(u, u_plan, v, w, lon, lat, z):
    """
    Args:
        u: Zonal wind
        v: Meridional Wind
        w: Vertical wind
        lon: Longitude
        lat: Latitude
        z: Altitude
    Returns:
    """
    circ_u, circ_u_plan, circ_v, circ_w = (0, 0, 0, 0)
    # Elements 0, -1 and -2 are identical
    for n in range(len(u)-1):
        # u.r.cos(phi).dlambda
        dx = (a+z[n]) * np.cos(lat[n]) * (lon[n+1] - lon[n-1])/2
        circ_u += u[n] * dx
        
        circ_u_plan += u_plan[n] * dx

        # v.r.dphi
        dy = (a + z[n]) * (lat[n+1] - lat[n-1])/2
        circ_v += v[n] * dy

        # w.dz
        dz = (z[n+1] - z[n-1])/2
        circ_w += w[n]*dz

    circulation_rel = circ_u + circ_v + circ_w
    circulation_plan = circ_u_plan
    circulation = circulation_rel + circulation_plan

    return [circulation, circulation_rel, circulation_plan]


def mass_integrated_circulation_split(cubes, x, y, glat, gridpoints, theta_level, dtheta, dlambda, dphi, timeh, basetime_str):
    """
    Args:
        cubes (iris.cube.CubeList):
        x (np.Array): Circuit longitudes
        y (np.Array): Circuit latitudes
        glat (np.array): Grid latitudes
        gridpoints(np.Array): Array of gridpoint longitudes and latitudes of
            shape (ngp, 2)
        theta_level:
        dtheta (int): Isentrope spacing used to calculate volume integrals
        JB addition::
        dlambda (float): Longitudinal spacing of grid boxes (degrees)
        dphi (float): Latitudinal spacing of grid boxes (degrees)
    Returns:
        CubeList of 5 different mass integrated circulations
        # Risk taken, 'timeh' changed to 'time' at end
    """
    # Include points within circuit boundary
    points = np.array([x, y]).transpose()
    pth = Path(points)

    # Mask all points that are not contained in the circuit
    mask = np.logical_not(pth.contains_points(gridpoints).reshape(glat.shape))

    # Area = r**2 cos(phi) dlambda dphi
    levels = ('air_potential_temperature',
              [theta_level - dtheta / 2.0, theta_level,
               theta_level + dtheta / 2.0])
    zth = convert.calc('altitude', cubes, levels=levels)
    # JB remove area = (a + zth[1]) ** 2 * np.cos(np.deg2rad(glat)) * np.deg2rad(0.11) ** 2 #where does 0.11 come from?
    area = (a + zth[1]) ** 2 * np.cos(np.deg2rad(glat)) * np.deg2rad(dlambda) * np.deg2rad(dphi)
    area.units = 'm^2'
    area.rename('area')

    total_area = integrate(area, mask) ##

    # Volume = area * dz
    volume = area * (zth[2] - zth[0])
    volume.rename('volume')

    # Mass = density * volume
    levels = ('air_potential_temperature', [theta_level])
    density = convert.calc('air_density', cubes, levels=levels)[0]
    mass = density * volume
    mass.rename('mass')

    # Circulation = \int rho.pv.dv / dtheta
    pv = convert.calc('ertel_potential_vorticity', cubes, levels=levels)[0]
    pv.convert_units('m^2 K s-1 kg-1')
    pv_substance = pv * mass
    
    circulation = integrate(pv_substance, mask) / dtheta

    circulation.rename('mass_integrated_circulation')
    
    vty = pv_substance / area ##
    
    meanvort_1 = integrate(vty, mask) / dtheta ##
    
    meanvort_2 = circulation / total_area ##
    
    meanvort_1.rename('mass_integrated_circulation'+'darea_1') ##
    meanvort_2.rename('mass_integrated_circulation'+'darea_2') ##
    
    # and now for more, separated
    
    tlev_cstrnt = iris.Constraint(air_potential_temperature = theta_level)
    time_cstrnt = iris.Constraint(time = timeh)
    # constraint for loading pv devived quanitities
    # taken risk by changing 'timeh' to 'time' in above line at end
    
    pvnames = ['ertel_potential_vorticity', 'total_minus_adv_only_PV', 'dPV_tot', 'adv_only_PV']
    # names of pv files to extract
    
    pvlist = iris.load('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/' + basetime_str +
     '/prodm_op_gl-mn_' + basetime_str + '_c*_thsfcs_5K.nc', pvnames)
     
#    print pvlist
    
#    print pvlist[0].coord('time')
    
#    print timeh
     
    pvlist = pvlist.extract(tlev_cstrnt & time_cstrnt)
    
#    print pvlist
     
    circulist = iris.cube.CubeList([])
    
    circulist.append(circulation)
    circulist.append(meanvort_1) ##
    circulist.append(meanvort_2) ##
    
    circnames = ['total_circulation', 'tot-adv_circulation', 'diab_circulation', 'adv_only_circulation']

    for cn, cube in enumerate(pvlist):
        cube.convert_units('m^2 K s-1 kg-1')
        c_sub = cube * mass
        
        circln = integrate(c_sub, mask) / dtheta
        
        vty = c_sub / area ##

        meanvort_1 = integrate(vty, mask) / dtheta        ## 
        
        meanvort_2 = circln / total_area ##
        
        circln.rename(circnames[cn])
        
        meanvort_1.rename(circnames[cn]+'darea_1') ##
        meanvort_2.rename(circnames[cn]+'darea_2') ##
    
        circulist.append(circln)
        
        circulist.append(meanvort_1) ##
        circulist.append(meanvort_2) ##
        
#        print 'circulist'
#        
#        print circulist
        
    return circulist, total_area ##
    
    
def plot_timeseries_split(case = 0, levs = 3, tls = [320, 325, 310, 310], divarea = False, folder = 'IOP3/T42', strapp = ''):
    "Based on leo's plot_timeseries"
    
    save_dir = '/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/'
    
    #tls = [320, 325, 310, 310]
    # pre-defined starting theta levels for three cases, subject to change
    
    #colours = [[[.8, .8, .8], [.4, .4, .4], [0, 0, 0]], [[.5, .5, 1], [0, 0, 1], [0, 0, .3]], [[1, .5, .5], [1, 0, 0], [.3, 0, 0]]]
    colours = ['k', 'b', 'r', 'b']
    # list of greys, blues, reds of different shades 
    lines = ['-', '-', '-', '--']
    
    plt.figure(figsize=(10, 5*levs)) # JB impose size to make time labels readable
    
    for i in xrange(levs):
        # three theta levels par case, also subject to change
        theta = tls[case] + 5*i
        # theta levels are spaced by 5
        cubes = iris.load(save_dir + folder + '/circulation/circulations_' + str(theta) + 'K_' + strapp + '_split.nc')
        
        circnames = ['mass_integrated_circulation', 'tot-adv_circulation', 'adv_only_circulation', 'diab_circulation']
        #circnames = ['total_circulation', 'tot-adv_circulation', 'diab_circulation', 'adv_only_circulation']
    
        m = 0
        plt.subplot(levs, 2, i+1) # JB for this specific case with 6 plots
        for name in circnames:
            # then it's a mass integrated
        
            if divarea == True:
                name = name + 'darea_2' # or 'darea_1'
            cube = cubes.extract(name)[0]
            if cube == False:
                pass
            else:
#                cm = int(np.floor(m/2))
#                dm = int(2*(m/2 - cm))
#                print m
                iplt.plot(cube, label=name, color = colours[m], linestyle = lines[m], linewidth = 2)
                m += 1

        pg = plt.gca()
        fmt = DateFormatter('\n%m/%d')                          # <--- change this bit to get different formatting of dates/times
        fmt2 = DateFormatter('%H')
        majorLocator = DayLocator(interval=1)
        minorLocator = HourLocator(range(0, 24, 6))
        pg.xaxis.set_major_formatter(fmt)
        pg.xaxis.set_minor_formatter(fmt2)
        pg.xaxis.set_minor_locator(minorLocator)
        pg.xaxis.set_major_locator(majorLocator)
        if i == 0: 
            plt.legend(ncol=3, loc=3, bbox_to_anchor=(-.03, 1.07, .6, .102))
        # plt.savefig(plotdir + 'circulation_' + theta + 'K.png')
        #plt.title(folder + '_' + strapp)

        m = 0
        plt.subplot(levs, 2, i + 3) # JB for this specific case with 6 plots
        for name in ['circulation', 'relative_circulation', 'planetary_circulation']:  
            # then it's a line integral
        
            if divarea == True:
                name = name + 'darea'
            cube = cubes.extract(name)[0]
            if cube == False:
                pass
            else:
                #iplt.plot(cube, label = name + '_' + str(theta) + 'K', color = colours[m], linewidth = 2)
#                print m
                iplt.plot(cube, label = name, color = colours[m], linewidth = 2)
                m += 1
            
        pg = plt.gca()
        fmt = DateFormatter('\n%m/%d')                          # <--- change this bit to get different formatting of dates/times
        fmt2 = DateFormatter('%H')
        majorLocator = DayLocator(interval=1)
        minorLocator = HourLocator(range(0, 24, 6))
        pg.xaxis.set_major_formatter(fmt)
        pg.xaxis.set_minor_formatter(fmt2)
        pg.xaxis.set_minor_locator(minorLocator)
        pg.xaxis.set_major_locator(majorLocator)
        if i == levs - 1:
            plt.legend(ncol=3, loc=1, bbox_to_anchor=(0.02, -.22, 1., .102))
        
    if divarea == True:
        strapp = 'darea' + strapp
        
    plt.savefig('IOP3_T42_circsplit_divarea_2levs.png')
    #plt.savefig(save_dir + folder + '/integral_split_' + strapp + '.png')
    plt.show()

    return    
    

    
    
    