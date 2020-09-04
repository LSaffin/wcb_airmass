# -*- coding: utf-8 -*-
"""
Created on Wed May 30 15:32:52 2018

@author: bn826011
"""

from scripts.circuit_integral import circuit_integral_rotated, mass_integrals, plot_timeseries
#from scripts.circuit_integral import circuit_integrals, plot_timeseries, integrate, get_geographic_coordinates
from circ_split_2019rad import circuit_integral_split, mass_integrated_circulation_split, plot_timeseries_split
#from circ_split import circuit_integral_split, mass_integrated_circulation_split, plot_timeseries_split
import iris
from mymodule import convert, grid
import numpy as np
from mymodule.constants import omega
from lagranto.trajectory import load
import datetime
from mymodule import forecast

import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt

from matplotlib.path import Path
import tropopause_contour as trop
from outflow_area_2019 import len_con, increase_nodes
from scipy.ndimage import filters

from matplotlib.dates import DayLocator, HourLocator, DateFormatter

iris.FUTURE.cell_datetime_objects=True

a = 6371229 # JB changed from 6378100

#save_dir = '/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/'
save_dir= '/storage/silver/scenario/bn826011/WCB_outflow/Final/'
IOP = [3, 5, 6]

def main_cc(k = 0, filtersize = 30, theta_level = 325, dtheta = 1, split = False, folder = 'IOP3/T42', strapp = ''):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    #datadir = '/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    datadir = '/storage/silver/NCAS-Weather/ben/nawdex/mi-ar482//{}/'.format(basetime_str)
    
    TrB = load(save_dir + folder + '/{}_TrajectoryEnsemble_backward'.format(basetime_str) + strapp)
    TrF = load(save_dir + folder + '/{}_TrajectoryEnsemble_forward'.format(basetime_str) + strapp)
    
    start_time = [-1*TrB.relative_times[-1]]
    # time at which trajectories were initialised
    
    trajectories = TrB.__add__(TrF)
#    trajectories = load(save_dir + folder + '/{}_TrajectoryEnsemble_forward'.format(basetime_str) + strapp)
#    trajectories = TrB
#    start_time = [trajectories.relative_times[0]]
    
    mapping = {}
    for t in trajectories.times:
        leadtimehr = np.int((t - basetime[k]).total_seconds()) / 3600
        fn = 'prodm_op_gl-mn_{0}_d{1:03d}_thgrid.pp'.\
            format(basetime_str, 12 * (leadtimehr / 12))
        mapping[t] = datadir + fn
        
    fcast = forecast.Forecast(basetime[k], mapping)
    
    if split == True:
        
        strapp = strapp + '_splitLW'#'_split'
        
    
    calc_circulation(trajectories, fcast, theta_level, dtheta, start_time, filtersize, k, split, folder, strapp, basetime_str)
    #shadow(trajectories, fcast, theta_level, dtheta, start_time, 5, folder, strapp)
    # handy to tack on here to the setup
    

def calc_circulation(trajectories, fcast, theta_level, dtheta, start_time, filtersize, case, split, folder, strapp, basetime_str):
    "Copied directly from l.saffin circuit_integral.py"  
    
    # Select an individual theta level
    trajectories = trajectories.select(
        'air_potential_temperature', '==', theta_level, time = start_time)
    # JB this allows trajectories which leave the domain to be selected also
    # JB this is dealt with later (start_time addition)
    print len(trajectories)
    levels = ('air_potential_temperature', [theta_level])

    results = iris.cube.CubeList()
    for n, time in enumerate(trajectories.times):
        cubelist = fcast.set_time(time)
        #now for some reason I'm having real problems extracting the right time, so
#        u = iris.unit.Unit('hours since 1970-01-01 00:00:00', calendar=iris.unit.CALENDAR_STANDARD)
#        timeh = u.date2num(time)
        #this puts the time in the same format as in the cubelist
        cubes = cubelist.extract(iris.Constraint(time = time))
        # NOTE: THIS USED TO ONLY WORK WITH timeh, NOW ONLY WORKS WITH time ???
        # JB needed because each file contains two times
        #print cubelist
        print n
        if n == 0:
            # Load grid parameters
            example_cube = convert.calc('upward_air_velocity', cubes,
                                        levels=levels)

            # Create a 1d array of points for determining which gridpoints are
            # contained in the trajectory circuit when performing volume
            # integrals
            glon, glat = grid.get_xy_grids(example_cube)
            gridpoints = np.array([glon.flatten(), glat.flatten()]).transpose()
            cs = example_cube.coord_system()
            
            # JB to do integrals need dlambda, dphi, not always 0.11
            
            dlon = np.diff(glon)
            if np.diff(dlon).any():
                print 'longitudinal spacing not uniform'
                # there exists a non-zero difference between longitude spacings
            else:
                dlambda = dlon[0][0]
                #as they are all the same
                
            dlat = np.diff(glat.transpose())
            if np.diff(dlat).any():
                print 'latitudinal spacing not uniform'
                # there exists a non-zero difference between latitude spacings
            else:
                dphi = dlat[0][0]
                #as they are all the same
                

        # Load trajectory positions -(n+2) because the trajectories are
        # backwards in time. +2 to skip the analysis which is not in the
        # forecast object (i.e. n=0 corresponds to idx=-2 in the trajectories)

        #JB I'm making them forwards in time so no worries, load n
        x = trajectories.x[:, n]
        y = trajectories.y[:, n]
        z = trajectories['altitude'][:, n]
        u = trajectories['x_wind'][:, n]
        v = trajectories['y_wind'][:, n]
        w = trajectories['upward_air_velocity'][:, n]
        
        # Integrals are invalid once trajectories leave the domain but we don't
        # want to stop the script so just print out the number of trajectories
        # that have left the domain
        
        leftflag = (trajectories['air_pressure'][:, n] < 0).astype(int)
        leftcount = np.count_nonzero(leftflag)
        print leftcount
        
        # Calculate enclosed area integrals
        if split == False:
            integrals = mass_integrals(cubes, x, y, glat, gridpoints,
                                   theta_level, dtheta, dlambda, dphi)
        else:
            integrals, area = mass_integrated_circulation_split(cubes, x, y, glat, gridpoints,
                                   theta_level, dtheta, dlambda, dphi, time, basetime_str)
        # function options mass_integrals or mass_integrated_circulation_split
                                   # RISK TAKEN timeh CHANGED TO time
                                   
        #print area
        #print area.data
                                   

        for icube in integrals:
            # Set integrals to zero if trajectories have left the domain
            if leftcount > 0:
                icube.data = float('nan') # was 0.
            results.append(icube)

        # Convert to global coordinates in radians
        #u, v, lon, lat = get_geographic_coordinates(u, v, x, y, cs)
        # JB I don't think that is necessary as the data is not on a rotated grid
        lon = x
        lat = y
        # JB but I'm not sure if that is right

        # Unrotated coordinates in radians
        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)

        # Calculate the velocity due to Earth's rotation
        u_abs = omega.data * (a+z) * np.cos(lat)
        #u += u_abs JB remove this line, move into circ_int_rot

        # Integrate around the circuit
        if leftcount > 0:
            circulation = float('nan') # was 0
        else:
            if split == False:
                circulations = [circuit_integral_rotated(u, v, w, lon, lat, z)]
            else:
                #circulations = circuit_integral_split(u, u_abs, v, w, lon, lat, z)
                circulations = []
            # JB edit to cope with multiple output from above function 
            # function options circuit_integral_rotated, or circuit_integral_split
            
            cnames = ['circulation', 'relative_circulation', 'planetary_circulation']
            
            for ccn, circulation in enumerate(circulations):
            
                ccube = icube.copy(data=circulation)
                ccube.rename(cnames[ccn])
                ccube.units = 'm^2 K s-1 kg-1'
                results.append(ccube)
                
                if split == True:
                
                    acube = icube.copy(data=circulation/area.data)
                    acube.rename(cnames[ccn]+'darea')
                    acube.units = 'K s-1 kg-1'
                    results.append(acube)
                

#    print results

    iris.save(results.merge(),
               save_dir + folder + '/circulation/circulations_' + str(theta_level) + 'K_' + strapp + '.nc')
    #          'outflow/T42_mfs' + str(filtersize) + '/circulations_' + str(theta_level) + 'K_case' + str(case) + '.nc')
    # save on glusterfs now

    return
    
    
def load_from_files(theta, case):
    cubes = iris.load(save_dir + 'IOP' + str(IOP[case]) + '/circulations_' + str(theta) + '.nc')
    #cubes = iris.load('outflow/T42_mfs' + str(filtersize) + '/circulations_' + str(theta) + 'K_case' + str(case) + '.nc')
    #plot_timeseries(cubes, str(theta), save_dir + 'IOP' + str(IOP[case]) + '/')
    plot_timeseries(cubes, str(theta), 'circ_plots/IOP' + str(IOP[case]) + '/')
    #plot_timeseries(cubes, str(theta), 'outflow/T42_mfs' + str(filtersize) + '/Case_' + str(case) + '_timeseries/')
    return cubes
    
    
def load_all_files():
    list_of_cubes = [[[], [], []], [[], [], []], [[], [], []]]
    tls = [320, 325, 310]
    for k in [0, 1, 2]:
        for j, tps in enumerate([0, 5, 10]):
            cubes = load_from_files(tls[k]+tps, k)
            list_of_cubes[k][j].append(cubes)
                
    return list_of_cubes
    
    
    
def plot_timeseries_levs(case = 0, levs = 3, tls = [320, 325, 310, 310], folder = 'IOP3/T42', strapp = ''):
    "Based on leo's plot_timeseries"
    
    #tls = [320, 325, 310, 310]
    # pre-defined starting theta levels for three cases, subject to change
    
    if levs == 3:
        
        colours = [['c', 'b', [0, 0, .3]], [[0, 1, 0], 'g', [0, .3, 0]]]
        # list of blues and greens of different shades 
        
    else:
        
        colours = [[[0, 0, 1]], [[0, 1, 0]]]
        
        for lev in xrange(levs - 1):
            
            blu = [0, 0, 1 - (float(lev)+1)/levs]
            grn = [0, 1 - (float(lev)+1)/levs, 0]
            
            colours[0].append(blu)
            colours[1].append(grn)
    
    plt.figure(figsize=(10, 8)) # JB impose size to make time labels readable
    
    for i in xrange(levs):
        # three theta levels par case, also subject to change
        theta = tls[case] + 5*i
        # theta levels are spaced by 5
        cubes = iris.load(save_dir + folder + '/circulation/circulations_' + str(theta) + 'K_' + strapp + '.nc')
    
        n = 1
        m = 0
        plt.subplot(2, 2, n) # JB for this specific case with 4 plots
        for cube in cubes:
            if 'circulation' in cube.name(): # the '_' before can be added to remove the line integral
                iplt.plot(cube, label=cube.name(), color = colours[m][i], linewidth = 2.5)
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

        plt.legend(ncol=2, loc=3, bbox_to_anchor=(0., 1.07, .6, .102))
        # plt.savefig(plotdir + 'circulation_' + theta + 'K.png')
        #plt.title(folder + '_' + strapp)
        plt.title('Circulation')

        for cube in cubes:        
            if 'circulation' not in cube.name():
                n += 1
                plt.subplot(2, 2, n) # JB for this specific case with 4 plots
                qplt.plot(cube/cube.data[2], color = colours[0][i], label = str(theta) + 'K', linewidth = 2.5)
                # plt.savefig(plotdir + cube.name() + '_' + theta + 'K.png')
                
                pg = plt.gca()
                fmt = DateFormatter('\n%m/%d')                          # <--- change this bit to get different formatting of dates/times
                fmt2 = DateFormatter('%H')
                majorLocator = DayLocator(interval=1)
                minorLocator = HourLocator(range(0, 24, 6))
                pg.xaxis.set_major_formatter(fmt)
                pg.xaxis.set_minor_formatter(fmt2)
                pg.xaxis.set_minor_locator(minorLocator)
                pg.xaxis.set_major_locator(majorLocator)
                plt.xlabel('')
                
        plt.legend(loc='best')
        
    if levs == 1:
        strapp = str(int(tls[case])) + strapp

    #plt.savefig(save_dir + folder + '/circ_integral_ts_' + strapp + '.png')
    plt.savefig('IOP3_T42_integrals.png')
    plt.show()

    return
    
def percent_change(k = 0, theta = 320, folder = 'IOP3/T42', strapp = '',
    t1 = datetime.datetime(2016, 9, 23, 00), t2 = datetime.datetime(2016, 9, 24, 06)):
    "percentage change in timeseries between inflow and outflow times"
    
    ts = [[datetime.datetime(2016, 9, 23, 00), datetime.datetime(2016, 9, 24, 06)],
          [datetime.datetime(2016, 9, 26, 18), datetime.datetime(2016, 9, 28, 00)],
          [datetime.datetime(2016, 10, 1, 00), datetime.datetime(2016, 10, 2, 06)],
          [datetime.datetime(2016, 10, 3, 12), datetime.datetime(2016, 10, 4, 12)]]
    
    cubes = iris.load(save_dir + folder + '/circulation/circulations_' + str(theta) + 'K_' + strapp + '.nc')
    
    chdic = {}
    
    for cube in cubes:
        inf = cube.extract(iris.Constraint(time = ts[k][0])).data
        otf = cube.extract(iris.Constraint(time = ts[k][1])).data
        pchange = (1 - inf/otf)*100 
        chdic.update({cube.name():pchange})
        
    return chdic

def print_percent_change():
    
    print 'IOP3'
    for th in [320, 325]:
        print th
        p = percent_change(theta = th)
        print('{:03.2f}'.format(round(p['area'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['volume'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['mass'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['mass_integrated_circulation'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['circulation'], 2)) + ' \\\ ')
    print 'IOP5'
    for th in [325, 330, 335]:
        print th
        p = percent_change(k = 1, theta = th, folder = 'IOP5/T36')
        print('{:03.2f}'.format(round(p['area'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['volume'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['mass'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['mass_integrated_circulation'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['circulation'], 2)) + ' \\\ ')
    print 'IOP6'
    for th in [310, 315, 320]:
        print th
        p= percent_change(k = 2, theta = th, folder = 'IOP6', strapp = '0')
        print('{:03.2f}'.format(round(p['area'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['volume'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['mass'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['mass_integrated_circulation'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['circulation'], 2)) + ' \\\ ')
    print 'IOP7'
    for th in [310, 315, 320]:
        print th
        p = percent_change(k = 3, theta = th, folder = 'IOP7')
        print('{:03.2f}'.format(round(p['area'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['volume'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['mass'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['mass_integrated_circulation'], 2)) + ' & '
        + '{:03.2f}'.format(round(p['circulation'], 2)) + ' \\\ ')
        
        
def mass_change(k = 0, theta = 320, tht = 0, folder = 'IOP3/T42', strapp = ''):
    "percentage change in timeseries between inflow and outflow times"
    
    ts = [[datetime.datetime(2016, 9, 23, 00), datetime.datetime(2016, 9, 24, 06)],
          [datetime.datetime(2016, 9, 26, 18), datetime.datetime(2016, 9, 28, 00)],
          [datetime.datetime(2016, 10, 1, 00), datetime.datetime(2016, 10, 2, 06)],
          [datetime.datetime(2016, 10, 3, 12), datetime.datetime(2016, 10, 4, 12)]]
    
    cubes = iris.load(save_dir + folder + '/circulation/circulations_' + str(theta) + 'K_' + strapp + '.nc')
    
    rmass = np.load('easy_to_find/' + folder + strapp + '_ridge_mass.np.npy')[tht]  
    
    mass = cubes.extract(iris.Constraint(name = 'mass'))[0]
    
    minf = mass.extract(iris.Constraint(time = ts[k][0])).data
    motf = mass.extract(iris.Constraint(time = ts[k][1])).data
    mdiab = motf - minf
    madv = rmass - mdiab
    
    print('{:03.2f}'.format(round(rmass/1e14, 2)) + ' & '
        + '{:03.2f}'.format(round(minf/1e14, 2)) + ' & '
        + '{:03.2f}'.format(round(motf/1e14, 2)) + ' & '
        + '{:03.2f}'.format(round(mdiab/1e14, 2)) + ' & '
        + '{:03.2f}'.format(round(madv/1e14, 2)) + ' & '
        + '{:03.2f}'.format(round(mdiab*100/rmass, 2)) + ' & '
        + '{:03.2f}'.format(round(madv*100/rmass, 2)) + ' \\\ ')
    
def print_mass_change():
    
    print 'IOP3'
    for th in [320, 325, 330]:
        print th
        mass_change(k = 0, theta = th, tht = (th-320)/5)
    print 'IOP5'
    for th in [325, 330, 335]:
        print th
        mass_change(k = 1, theta = th, tht = (th-325)/5, folder = 'IOP5/T36')
    print 'IOP6'
    for th in [310, 315, 320]:
        print th
        mass_change(k = 2, theta = th, tht = (th-310)/5, folder = 'IOP6', strapp = '0')
    print 'IOP7'
    for th in [310, 315, 320]:
        print th
        mass_change(k = 3, theta = th, tht = (th-310)/5, folder = 'IOP7')
    
    
        
        
        
        
 # slightly out of place tack-ons       
        
def shadow(trajectorys, fcast, theta_levels, dtheta, start_time, chosen_time_index, folder, strapp):
    """
    Shitty hack to return sum of shadows
    Chosen_time_index 
    trajectorys deliberately misspelt
    definitely move somewhere else when done
    """
    masks = []
    for theta_level in theta_levels:
    
        # Select an individual theta level
        trajectories = trajectorys.select(
            'air_potential_temperature', '==', theta_level, time = start_time)
        # JB this allows trajectories which leave the domain to be selected also
        # JB this is dealt with later (start_time addition)
        print len(trajectories)
        levels = ('air_potential_temperature', theta_levels)

        time = trajectories.times[chosen_time_index]
        n = chosen_time_index
        cubelist = fcast.set_time(time)
        #now for some reason I'm having real problems extracting the right time, so
#        u = iris.unit.Unit('hours since 1970-01-01 00:00:00', calendar=iris.unit.CALENDAR_STANDARD)
#        timeh = u.date2num(time)
        #this puts the time in the same format as in the cubelist
        cubes = cubelist.extract(iris.Constraint(time = time))
        # NOTE: THIS USED TO ONLY WORK WITH timeh, NOW ONLY WORKS WITH time ???
        # JB needed because each file contains two times
        #print cubelist
        
        if 1:
        
            # Load grid parameters
            example_cube = convert.calc('upward_air_velocity', cubes,
                                        levels=levels)

            # Create a 1d array of points for determining which gridpoints are
            # contained in the trajectory circuit when performing volume
            # integrals
            glon, glat = grid.get_xy_grids(example_cube)
            gridpoints = np.array([glon.flatten(), glat.flatten()]).transpose()
            cds = example_cube.coord_system()
            
            # JB to do integrals need dlambda, dphi, not always 0.11
            
            dlon = np.diff(glon)
            if np.diff(dlon).any():
                print 'longitudinal spacing not uniform'
                # there exists a non-zero difference between longitude spacings
            else:
                dlambda = dlon[0][0]
                #as they are all the same
                
            dlat = np.diff(glat.transpose())
            if np.diff(dlat).any():
                print 'latitudinal spacing not uniform'
                # there exists a non-zero difference between latitude spacings
            else:
                dphi = dlat[0][0]
                #as they are all the same
                

        # Load trajectory positions -(n+2) because the trajectories are
        # backwards in time. +2 to skip the analysis which is not in the
        # forecast object (i.e. n=0 corresponds to idx=-2 in the trajectories)

        #JB I'm making them forwards in time so no worries, load n
        x = trajectories.x[:, n]
        y = trajectories.y[:, n]
        z = trajectories['altitude'][:, n]
        u = trajectories['x_wind'][:, n]
        v = trajectories['y_wind'][:, n]
        w = trajectories['upward_air_velocity'][:, n]
    
        # Include points within circuit boundary
        points = np.array([x, y]).transpose()
        pth = Path(points)

        # Mask all points that are not contained in the circuit
        mask = np.logical_not(pth.contains_points(gridpoints).reshape(glat.shape)) 
    
        masks.append(mask)
        
    sumask = np.array(masks[0], dtype = int)
    for mask in masks[1:]:
        sumask = sumask + np.array(mask, dtype = int)
    value = len(masks) - 0.5  
    
    #sumooth = filters.median_filter(sumask, size = 20)
    
    cs = plt.contour(glon, glat, sumask, [value])
    
    contours = trop.get_contour_verts(cs)
    # returns array of arrays of vertices for zero contours
    
    #ncon = np.size(contours[0])
    ncon = len(contours[0])
    #number of contours
    
    lencon = np.zeros(ncon)
    # empty array of lengths of contorus (in lat/long space)
    
    print ncon
    for j in xrange(ncon):
        print j
        
        lencon[j] = len_con(contours[0][j])
    
    imax = lencon.argmax()
    # index of longest contour
    lcontour = contours[0][imax]
    # longest contour
    # don't worry about closed-ness
    
    points = increase_nodes(lcontour, resolution = .25)
    # increase number of points on the contour such that they have minimum spacing resolution    
    
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/outflow_shadow_contour.npy', points)
    
    
def bound_inflow_3D(k = 0, folder = 'IOP3/T42', in_time = 5, min_lev = 320, strapp = '', gmod = False):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    Tr3 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_new'.format(basetime_str) + strapp)
    
    Tr3_inflow = Tr3.select('air_potential_temperature', '<', min_lev - 2.5, time = [datetime.timedelta(hours = -6*in_time)])
    
    inf_area = plt.hist2d(Tr3_inflow.data[:, in_time, 0], Tr3_inflow.data[:, in_time, 1], bins = [360, 90], range = [[0, 360], [0, 90]])
    plt.show()
    if gmod:
        ia = filters.gaussian_filter(inf_area[0], gmod)
    else:
        ia = inf_area[0]
    cs = plt.contour(np.transpose(ia), [.5])
    
    contours = trop.get_contour_verts(cs)
    ncon = len(contours[0])
    #number of contours

    lencon = np.zeros(ncon)
    # empty array of lengths of contorus (in lat/long space)

    for j in xrange(ncon):

       lencon[j] = len_con(contours[0][j])

    imax = lencon.argmax()
    # index of longest contour
    lcontour = contours[0][imax]
    # longest contour
    # don't worry about closed-ness

    points = increase_nodes(lcontour, resolution = .25)
    
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow_bound_3D_less_than_' + str(min_lev) + 'K_gmod.npy', points)
    
    
def bound_outflow_3D(k = 0, folder = 'IOP3/T42', in_time = 5, min_lev = 320, strapp = '', gmod = False):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    Tr3 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_new'.format(basetime_str) + strapp)
    
    Tr3_inflow = Tr3.select('air_potential_temperature', '>', min_lev - 2.5, time = [datetime.timedelta(hours = -6*in_time)])
    
    inf_area = plt.hist2d(Tr3_inflow.data[:, in_time, 0], Tr3_inflow.data[:, in_time, 1], bins = [360, 90], range = [[0, 360], [0, 90]])
    plt.show()
    if gmod:
        ia = filters.gaussian_filter(inf_area[0], gmod)
    else:
        ia = inf_area[0]
    cs = plt.contour(np.transpose(ia), [.5])
    
    contours = trop.get_contour_verts(cs)
    ncon = len(contours[0])
    #number of contours

    lencon = np.zeros(ncon)
    # empty array of lengths of contorus (in lat/long space)

    for j in xrange(ncon):
        print ncon
        print j

        lencon[j] = len_con(contours[0][j])

    imax = lencon.argmax()
    # index of longest contour
    lcontour = contours[0][imax]
    # longest contour
    # don't worry about closed-ness

    points = increase_nodes(lcontour, resolution = .25)
    
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/outflow_bound_3D_less_than_' + str(min_lev) + 'K_gmod.npy', points)
   
   
   
def compare_inflow_bounds(folder = 'IOP3/T42', min_theta = 320):
    
    shadow_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/'+ folder + '/shadow_contour.npy')
    
    smooth_shadow_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/'+ folder + '/smooth_shadow_contour.npy')

    bound_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow_bound_3D_less_than_' + str(min_theta) + 'K.npy')
    
    big_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow_bound_3D_less_than_' + str(min_theta) + 'K_gmod.npy')
    
    bound_outf = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/outflow_bound_3D_less_than_' + str(min_theta) + 'K.npy')

    plt.plot(shadow_cont[:, 0], shadow_cont[:, 1], label = 'sum of shadows of isentropic V2 at ts')
    plt.plot(smooth_shadow_cont[:, 0], smooth_shadow_cont[:, 1], linestyle = ':', linewidth = 4, label = 'smoothed sum of shadows of isentropic V2 at ts')
    plt.plot(bound_outf[:, 0], bound_outf[:, 1],  linestyle = '--', label = 'bound around 3D trajectories at ts within V2 (should == shadow)')
    plt.plot(bound_cont[:, 0], bound_cont[:, 1], linewidth = 1.5, label = 'bound around 3D trajectories at ts below V2')
    plt.plot(big_cont[:, 0], big_cont[:, 1], linewidth = 2, label = 'bound around 3D trajectories at ts below V2, slightly larger')
    
    plt.title(folder + '  comparison of inflow lateral bounds')
    plt.legend(bbox_to_anchor=(1.2, -.1))
    
    plt.savefig('/home/users/bn826011/NAWDEX/From N Drive/2019_figs/' + folder[:4] + '_bounds_comparison2.jpg')
    plt.show()