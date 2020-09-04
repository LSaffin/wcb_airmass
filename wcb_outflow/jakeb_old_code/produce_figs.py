# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 14:37:18 2018

@author: bn826011
"""

import numpy as np
import iris
import matplotlib
import matplotlib.pyplot as plt
import iris.plot as iplt
from scipy.ndimage import filters
import tropopause_contour as trop
import datetime
import operator

from lagranto import caltra

from lagranto.trajectory import Trajectory, TrajectoryEnsemble, load
import mymodule.plot as modplt

import iris.quickplot as qplt

from mymodule import forecast, convert, grid
from July3 import find_ridge
from scripts.circuit_integral import mass_integrals

import random
from matplotlib.dates import DayLocator, HourLocator, DateFormatter

PVU = 2 # location of the dynamical tropopause


def first_fig(k=0):

    datestr = ['0922_12', '0922_12', '0922_12', '1002_00', '1003_12']

    thetastart = [8, 9, 8, 6, 6]#8, 9, 6

    hourst = [24, 42, 60]
    
    llevs = 3    

    V = np.linspace(-10, 40, 11)

    filtersize = 30

    fig = plt.figure(figsize=(10, 5.5))
    
    cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])
    #for colorbar

    cldt = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_b*_thsfcs_5K.nc', 'total_minus_adv_only_theta') #  '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
    dtcube = cldt.concatenate_cube()
    # create cube of delta theta
    clpv = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
    pvcube = clpv.concatenate_cube()
    # create cube of delta theta
    
    times = len(hourst)
    
    for i in xrange(times):

        dtcub = dtcube[hourst[i]/3, thetastart[k]:thetastart[k]+llevs]
        pvcub = pvcube[hourst[i]/3, thetastart[k]:thetastart[k]+llevs]
        
        for j in xrange(llevs):
            
            dt = dtcub.copy()[j]
            pv = pvcub.copy()[j]
            
            
            
            #plt.subplot(fsz, 3, 3*j+i+1)#7
            plt.subplot(times, llevs, llevs*i + j + 1)
            
            ax = iplt.contourf(dt, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
            
            pv.data = filters.median_filter(pv.data, size = filtersize)
            dt.data = filters.median_filter(dt.data, size = filtersize)            
            
            iplt.contour(pv, [PVU], colors = ['m'])
            iplt.contour(dt, [0], colors = ['g'])
            
            criteria = np.logical_and(pv.data < 2, dt.data > 0)
            # from Leo's code
        
            pv.data = criteria.astype(int)        
        
            cs = iplt.contour(pv, [0.5], colors = ['k'])
            # additional requirement of pv < 2 on contour
            
            if i == 0:
                
                potemp = dtcube.coord('air_potential_temperature')
                
                plt.title(str(potemp.points[thetastart[k]+j]) + ' K')
                #label columns with isentropic levels
                
            if j == 0:
                
                reftime = datetime.datetime(1970, 1, 1)
                
                hourssince = dtcub.coord('time').points[0]
                
                timedatenow = reftime + datetime.timedelta(hours = hourssince)
                
                plt.text(-90, 70, timedatenow.strftime("%d/%m_%HUTC"), rotation = 'vertical')
                #label rows with dates and times
                
    cbar = fig.colorbar(ax, cax, orientation = 'horizontal')
    cbar.set_label('Delta theta on isentropic surfaces and dynamical tropopause')
            
    plt.savefig('IOP3_brief.png')
    
    
    
def second_fig(k = 0, levs = 3, filtersize = 30, field = 'total_minus_adv_only_theta', scatter3D = False, lim = 40, indiv = False, folder = 'IOP3/T42', strapp = ''):
    "Options for field:"
    "total_minus_adv_only_theta (default)"
    "ertel_potential_vorticity"
    "pv_substance" # recommend limits ~ 1e11
    "isentropic_density" # recommend limits ~ 1000
    "dPV_: cld, LW, mass, iau, bl, tot, mic, SW, gwd, adv, phl, sol, conv"
    "adv_only_PV"
    "total_minus_adv_only_PV"
    "Or, if 3Dscatter == True:"
    "air_potential_temperature"
    "air_pressure"
    "specific_humidity"
    "mass_fraction_of_cloud_ice_in_air"
    "mass_fraction_of_cloud_liquid_water_in_air"
    "lim is +/- limits of contours coloured distinctly for field, default 40"
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    #TrB = load('outflow/T42_mfs' + str(filtersize) + '/{}_TrajectoryEnsemble_backward'.format(basetime_str))
    #TrF = load('outflow/T42_mfs' + str(filtersize) + '/{}_TrajectoryEnsemble_forward'.format(basetime_str))
    TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/{}_TrajectoryEnsemble_backward'.format(basetime_str) + strapp)
    TrF = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/{}_TrajectoryEnsemble_forward'.format(basetime_str) + strapp)
    
    TrEn = TrB.__add__(TrF)
    
    tsteps = len(TrEn.times)
    
    theta_0 = TrEn.data[0, 0, 2]
    
    V = np.linspace(-lim, lim, 35)
    
    clpv = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/' + basetime_str +
    '/prodm_op_gl-mn_' + basetime_str + '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
    pvcube = clpv.concatenate_cube()
    # create cube of PV
    
    if field == 'total_minus_adv_only_theta':
    
        cldt = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/' + basetime_str +
        '/prodm_op_gl-mn_' + basetime_str + '_b*_thsfcs_5K.nc', 'total_minus_adv_only_theta')
        cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
        fieldcube = cldt.concatenate_cube()
        # create cube of delta theta
        
    elif field == 'ertel_potential_vorticity':
        
        fieldcube = pvcube
        
    elif field in ['pv_substance', 'isentropic_density']:
        
        IOPs = [3, 5, 6, 7]
        
        fieldcube = iris.load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP' + 
        str(IOPs[k]) + '/' +  basetime_str + '_isentropic_density.nc', field)[0]
        
    elif scatter3D == False:
        #assume user has input the name of a diabatic PV field
        
        cldpv = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/' + basetime_str +
        '/prodm_op_gl-mn_' + basetime_str + '_c*_thsfcs_5K.nc', field)
        cldpv[-1] = iris.util.new_axis(cldpv[-1], 'time')
        fieldcube = cldpv.concatenate_cube()
        # and for PV
        
    if scatter3D != False:
        
        Tr3 = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble'.format(basetime_str) + strapp)
        
        field_no = Tr3.names.index(field)     
        
        strapp = 'scatter3D_' + strapp
        
        
    if indiv == False:    
        plt.figure(figsize = (levs*5, 40))
    
    for tht in [1]:
        
        theta = theta_0 + 5*tht
        
        TrEntheta = TrEn.select('air_potential_temperature', '==', theta, time = [-1*TrB.relative_times[-1]])
        # This gives the levels at the designated start time
        
        if scatter3D != False:
            
            Tt = Tr3.select('air_potential_temperature', '>=', theta - 2.5, time = [datetime.timedelta(hours = 0)])
            Tt = Tt.select('air_potential_temperature', '<', theta + 2.5, time = [datetime.timedelta(hours = 0)])
            Tt = Tt.select(field, '!=', -1000)
            # points starting from theta level for which field is defined
                    
        for t in [12]:
            
            if indiv == False:
                plt.subplot(tsteps, levs, t*levs + tht + 1)
            else:
                pass
                #plt.figure(figsize = (5, 4))
                # for plotting individual frames
                
            pvort = pvcube.extract(iris.Constraint(time = TrEn.times[t], air_potential_temperature = theta))
            
            if scatter3D == False:
            
                c_field = fieldcube.extract(iris.Constraint(time = TrEn.times[t], air_potential_temperature = theta))            
            
                iplt.contourf(c_field, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-lim, lim))
                # plot contours of field, blue & red
            
            iplt.contour(pvort, [2], colors = ['g'])
            # plot 2 PVU contour, green
            plt.gca().coastlines(color = [.6, .4, 0])
            # plot coastlines, dark yellow
            
            if scatter3D == True:
                
                if t < len(Tt.times):
                    
                    if field == 'air_potential_temperature':
                        
                        norm = matplotlib.colors.Normalize(theta-lim, theta+lim)
                        
                    else:
                        
                        norm = matplotlib.colors.Normalize(-lim, lim)
                
                    plt.scatter(Tt.data[:, -1-t, 0]-360, Tt.data[:, -1-t, 1], c = Tt.data[:, -1-t, field_no], marker = '.', 
                            cmap = 'inferno', norm = norm, edgecolors = 'face')
            
            
            TrEnindomain = TrEn.select('air_potential_temperature', '==', theta, time = [datetime.timedelta(hours=6*t)])
            
            Leftdomain = TrEntheta.__len__() - TrEnindomain.__len__()
            
            if np.shape(TrEntheta) == (0,):
                
                pass
            
            elif Leftdomain != 0:
                #If any of the trajectories have left the domain
            
                plt.plot(TrEnindomain.data[:, t, 0]-360, TrEnindomain.data[:, t, 1], color = [.4, .4, .4], marker = '.', ms = .8, mec = 'k', mfc = 'k')
                
                plt.title(str(Leftdomain) + ' points have left the domain')
            
            else:
            
                plt.plot(TrEntheta.data[:, t, 0]-360, TrEntheta.data[:, t, 1], color = [.4, .4, .4], marker = '.', ms = .8, mec = 'k', mfc = 'k')
                # plot outflow area, black
            
            if indiv == True:
                pass
#                plt.savefig('IOP3_T42_325K_t12.png')
            #for plotting individual frames            
            
#            if t == 0:
#                
#                if Leftdomain != 0:
#                
#                    plt.title(str(theta) + ' K \n ' + str(Leftdomain) + ' points have left the domain')
#                    #label columns with isentropic levels
#                    
#                else:
#                    
#                    plt.title(str(theta) + ' K')
#                
#                
#            if tht == 0:
#                
#                timedatenow = TrEn.times[t]
#                
#                plt.text(-90, 70, timedatenow.strftime("%d/%m %HUTC"), rotation = 'vertical')
#                
#                #label rows with dates and times
                            
    #plt.savefig('outflow/T42_mfs' + str(filtersize) + '/' + basetime_str + 'fulltraj_' + field + '.jpg')
    if indiv == False:
        
        plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + basetime_str + 'fulltraj_' + field + strapp + '.png')

    
def prod2fig():
    
    plt.figure(figsize = (9, 4), dpi = 400)
    plt.subplot(1, 2, 1)
    second_fig(indiv = True)
    plt.subplot(1, 2, 2)
    second_fig(indiv = True, folder = 'IOP3/T60')
    plt.savefig('IOP3_325K_t12.png', dpi = 400)
    
def single_plume(k = 0, levs = 3, theta_0 = [320, 325, 310, 310], trajs = 50, folder = 'IOP3/T42', strapp = '', add = ''):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble'.format(basetime_str) + strapp)
    TrF = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_fw'.format(basetime_str) + strapp)
    
    Tr3 = TrB.__add__(TrF) 
#    Tr3 = TrB
    
    #deftime = len(TrB.times) - 1
    
    plt.figure(figsize = (5, 3))
    
    Tall = Tr3.select('air_potential_temperature', '!=', -1000)
    #only those which remain in domain
    
    for i in xrange(trajs):
        i = random.randint(0, len(Tall)-1)
        plt.plot(Tall.times, Tall.data[i, :, 5])#4   
        
#    Tmean = np.mean(Tall.data[:, :, 5], axis = 0)
#    plt.plot(Tall.times, Tmean, color = 'k', linewidth = 3)
        
    #plt.ylabel('theta, K')
    plt.ylabel('specific_humidity')
    plt.xlabel('time')
    
    pg = plt.gca()
    fmt = DateFormatter('\n%m/%d')                  
    fmt2 = DateFormatter('%H')
    majorLocator = DayLocator(interval=1)
    minorLocator = HourLocator(range(0, 24, 6))
    pg.xaxis.set_major_formatter(fmt)
    pg.xaxis.set_minor_formatter(fmt2)
    pg.xaxis.set_minor_locator(minorLocator)
    pg.xaxis.set_major_locator(majorLocator)
    #pg.invert_yaxis()
        
    plt.savefig('IOP3_T60' + strapp + '50_random_3D-trajectories_humidity' + add + '.png')
    plt.show()
    
def single_hist(k = 0, levs = 3, theta_0 = [320, 325, 310, 310], t = 0, folder = 'IOP3/T42', strapp = ''):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    Tr3 = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble'.format(basetime_str) + strapp)
    
    bins = np.linspace(-200, 650, 18)
    #bins = np.linspace(-10, 40, 26)
    
    plt.figure(figsize = (5, 3))
    
    Tall = Tr3.select('air_potential_temperature', '!=', -1000)
    #only those which remain in domain
    
    plt.hist((Tall.data[:, -(t+1), 4] - Tall.data[:, 0, 4])/100, bins = bins)
    #plt.title(folder + '_' + strapp)
    #plt.title('Started at all surfaces')
    #plt.xlabel('theta, K')
    plt.xlabel('pressure, hPa')
    plt.ylabel('number of trajectories')
    
    
    plt.savefig('IOP3_T60' + strapp + '_change_pressure_hist_t=' + str(t) + '.png')
    plt.show()
    
    
    
def plot_circ_mass(case = 0, levs = 3, tls = [320, 325, 310, 310], folder = 'IOP3/T42', strapp = ''):
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
    
    plt.figure(figsize=(10, 4)) # JB impose size to make time labels readable
    
    for i in xrange(levs):
        # three theta levels par case, also subject to change
        theta = tls[case] + 5*i
        # theta levels are spaced by 5
        cubes = iris.load(save_dir + folder + '/circulation/circulations_' + str(theta) + 'K_' + strapp + '.nc')
    
        n = 1
        m = 0
        plt.subplot(1, 2, n) # JB for this specific case with 4 plots
        for cube in cubes:
            if 'circulation' in cube.name():
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
                if 'mass' in cube.name():
                    n += 1
                    plt.subplot(1, 2, n) # JB for this specific case with 4 plots
                    qplt.plot(cube, color = colours[1][i], label = str(theta) + 'K', linewidth = 2.5)
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
    #plt.savefig('IOP6_' + strapp + '_circ_mass_2levs.png')
    plt.savefig('IOP3_legend')
    plt.show()

    return    
    
def out_ridge_mass(k = 0, levs = 3, thetastart = [8, 9, 6, 6], hoursafterinit = [42, 36, 42, 24], filtersize = 5, dtheta = 1, folder = 'IOP3/T42', strapp = ''):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    datadir = '/export/cloud/NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    
    #eqlat = [0, 0, 0, 0, 0, 0, 55, 50, 47, 43, 38, 36] # equivalent latitude indexed on theta level, for November
    eqlat = [0, 0, 0, 0, 0, 0, 55, 51, 48, 45, 42, 40] # a guess at equivalent latitude for September, because who knows
    
    TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/{}_TrajectoryEnsemble_backward'.format(basetime_str) + strapp)
    
    theta_0 = TrB.data[0, 0, 2]
    
    ridge_mass = np.zeros(levs)
    
    for tht in xrange(levs):
        
        theta = theta_0 + 5*tht
        
        TrEntheta = TrB.select('air_potential_temperature', '==', theta, time = [datetime.timedelta(hours = 0)])
        # This gives the levels at the designated start time
        
        mapping = {}
        for t in TrEntheta.times:
            leadtimehr = np.int((t - basetime[k]).total_seconds()) / 3600
            fn = 'prodm_op_gl-mn_{0}_d{1:03d}_thgrid.pp'.\
                format(basetime_str, 12 * (leadtimehr / 12))
            mapping[t] = datadir + fn
        
        fcast = forecast.Forecast(basetime[k], mapping)
        
        time = TrB.times[0]
        t = 0
            
        cubelist = fcast.set_time(time)
        cubes = cubelist.extract(iris.Constraint(time = time))
        levels = ('air_potential_temperature', [theta])
                   
        # Load grid parameters
        example_cube = convert.calc('upward_air_velocity', cubes,
                                        levels=levels)

        # Create a 1d array of points for determining which gridpoints are
        # contained in the trajectory circuit when performing volume
        # integrals
        glon, glat = grid.get_xy_grids(example_cube)
        gridpoints = np.array([glon.flatten(), glat.flatten()]).transpose()
            
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

        x = TrEntheta.x[:, t]
        y = TrEntheta.y[:, t]
                
        thtst = [0, 0, 0, 0]
                
        thtst[k] += tht + thetastart[k]
                
        ridge_area, pv = find_ridge(k = k, thlevs = 1, thetastart = thtst, hoursafterinit = hoursafterinit, filtersize = filtersize, eqlat = eqlat[thetastart[k]+tht], savet = False, folder = folder, strapp = strapp)
                
        xr = np.array(ridge_area[:, 0])
        yr = np.array(ridge_area[:, 1])
                
        ridge_mass[tht] = mass_integrals(cubes, xr, yr, glat, gridpoints, theta, dtheta, dlambda, dphi)[2].data
        
        iplt.contour(pv, [2], colours = ['m'])
        plt.plot(xr-360, yr, color = [.4, .4, .4], marker = '.', mec = 'k', mfc = 'k')
        plt.plot(x-360, y, color = [.4, .4, 1], marker = '.', mec = 'b', mfc = 'b')
        plt.gca().coastlines(color = [.6, .4, 0])
        
        plt.savefig('easy_to_find/' + folder + '_ridge_' + str(theta) + 'K_Nov.png')
        plt.show()
        
    np.save('easy_to_find/' + folder + strapp + '_ridge_mass.np', ridge_mass)
        
    return ridge_mass
    

#print 'IOP3/T42'
#print out_ridge_mass()
#print 'IOP5/T36'
#print out_ridge_mass(k = 1, folder = 'IOP5/T36')
#print 'IOP3/T60'
#print out_ridge_mass(hoursafterinit = [60], filtersize = 10, folder = 'IOP3/T60')


def pv_ascent(k = 0, levs = 3, theta_0 = [320, 325, 310, 310], trajs = 50, asc = 10, voltimes = [[12, 42], [6, 36], [12, 42], [0, 24]], folder = 'IOP3/T42', strapp = '', add = ''):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble'.format(basetime_str) + strapp)
#    TrF = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_fw'.format(basetime_str) + strapp)
    
#    Tr3 = TrB.__add__(TrF) 
    Tr3 = TrB
    
    
    #times = [datetime.timedelta(hours = voltimes[k][1]), datetime.timedelta(hours = voltimes[k][0])]
    times = [datetime.timedelta(hours = 0), datetime.timedelta(hours = (voltimes[k][0]-voltimes[k][1]))]
    
    #deftime = len(TrB.times) - 1
    
    plt.figure(figsize = (5, 3))
    
    Tall = Tr3.select('air_potential_temperature', '!=', -1000)
    #only those which remain in domain
    
    Tasc = Tall.select('air_potential_temperature', '<', asc, time = times)
    #only those which experience 10K of heating
    
    for i in xrange(trajs):
        i = random.randint(0, len(Tasc)-1)
        plt.plot(Tasc.times, Tasc.data[i, :, 8])#4   
        
    Tmean = np.mean(Tasc.data[:, :, 8], axis = 0)
    plt.plot(Tasc.times, Tmean, color = 'k', linewidth = 3)
        
    #plt.ylabel('theta, K')
    plt.ylabel('PV')
    plt.xlabel('time')
    
    pg = plt.gca()
    fmt = DateFormatter('\n%m/%d')                  
    fmt2 = DateFormatter('%H')
    majorLocator = DayLocator(interval=1)
    minorLocator = HourLocator(range(0, 24, 6))
    pg.xaxis.set_major_formatter(fmt)
    pg.xaxis.set_minor_formatter(fmt2)
    pg.xaxis.set_minor_locator(minorLocator)
    pg.xaxis.set_major_locator(majorLocator)
    #pg.invert_yaxis()
        
    plt.savefig('easy_to_find/' + folder + strapp + '50_random_3D-trajectories_PV_asc' + str(asc) + add + '.png')
    #plt.savefig('easy_to_find/' + folder + strapp + 'PV_asc' + str(asc) + add + '.png')
    plt.show()
            