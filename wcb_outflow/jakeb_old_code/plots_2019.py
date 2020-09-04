# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:29:49 2019

@author: bn826011
"""
from __future__ import division

import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
from matplotlib.colors import SymLogNorm, LogNorm
import numpy as np
import datetime

from mpl_toolkits.mplot3d import Axes3D
import matplotlib

from lagranto.trajectory import load

from lagranto.trajectory import Trajectory, TrajectoryEnsemble, load
from outflow_area_2019 import try_caltra_rev
from circ_int_test import main_cc

def area_by_theta(folder = 'IOP3/contours_2019'):
    # produce plots of area by theta level of contours identified at a specified start time

    plt.figure(figsize = (8, 8))
    
    save_dir = '/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/circulation/'
    
    filename = lambda theta: 'circulations_' + str(theta) + 'K_.nc'
    
    mass = iris.Constraint(name = 'mass')
    volume = iris.Constraint(name = 'volume')
    area = iris.Constraint(name = 'area')
    
    for theta in xrange(300, 340, 5):
        
        cubelist = iris.load(save_dir + filename(theta))
        
        plt.subplot(2, 2, 2)
        plt.scatter(cubelist.extract(mass)[0][0].data, theta)
        plt.xlabel('mass')
        plt.ylabel('theta')
        plt.subplot(2, 2, 3)
        plt.scatter(cubelist.extract(area)[0][0].data, theta)
        plt.xlabel('area')
        plt.ylabel('theta')
        plt.subplot(2, 2, 4)
        plt.scatter(cubelist.extract(volume)[0][0].data, theta)
        plt.xlabel('volume')
        plt.ylabel('theta')
        
    plt.savefig(folder + '_outflow_area_v_theta_at_tf.jpg')
        
    plt.show()
    
    
    
def trajectory_2D_3D_displacement_difference(k = 0, folder = 'IOP3/T42', strapp = ''):
    # figure 4
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    t3 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_new'.format(basetime_str) + strapp)
    t2 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_2DTrajectoryEnsemble'.format(basetime_str) + strapp)
    pts = len(t3.data[:, 0, 0])
    traj23diffs = np.zeros([pts, 3, 4])
    
    for i, t in enumerate([0, 5, -1]):
        traj23diffs[:, i, 0] = t2.data[:, t, 2] - t3.data[:, t, 3]
        traj23diffs[:, i, 1] = t2.data[:, t, 0] - t3.data[:, t, 0]
        traj23diffs[:, i, 2] = t2.data[:, t, 1] - t3.data[:, t, 1]
        traj23diffs[:, i, 3] = (traj23diffs[:, i, 1]**2 + traj23diffs[:, i, 2]**2)**0.5
        
    for i, t in enumerate([0, 5, -1]):
        plt.figure(figsize = (10, 8))
        plt.suptitle('Differences at ' + str(t3.times[t]))
        for j, name in enumerate([['longitude_difference', -40, 20], ['latitude_difference', -10, 10], ['total_displacement', -5, 40]]):
            plt.subplot(2, 2, j+1)
            plt.hist2d(traj23diffs[:, i, j+1], traj23diffs[:, i, 0], bins = 50, range = [[name[1], name[2]], [-10, 35]], cmap = 'binary', norm = LogNorm())
            plt.colorbar()#ticks = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000])
            plt.ylabel('theta_difference')
            plt.xlabel(name[0])
        plt.subplot(2, 2, 4)
        plt.hist2d(traj23diffs[:, i, 1], traj23diffs[:, i, 2], bins = 50, range = [[-40, 20], [-10, 10]], cmap = 'binary', norm = LogNorm())
        plt.colorbar()#ticks = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000])
        plt.xlabel('longitude_displacement')
        plt.ylabel('latitude_displacement')
        plt.savefig('/home/users/bn826011/NAWDEX/From N Drive/2019_figs/' + folder[:4] + '/traj23diffs_' + str(t) + 'x6_hours_before_new.jpg')
        plt.show()
    
    return t2, t3, traj23diffs

#for strapp in ['bound3DV1', 'bound3DV2', 'shadow', 'shadow_smooth']:
#    inflow_quantities(strapp)    
    
def inflow_quantities(strapp, IOP = 5, thmin = 285, thmax = 325):
    # figure 8
    
    cinf = iris.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/IOP' + str(IOP) +
                     '/contours_2019/circulation/circulations_' + str(thmin) + 'K_' + strapp + '.nc')
                     
    plt.figure(figsize = (10, 8))
                     
    for theta in range(thmin + 5, thmax, 5):
        cube = iris.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/IOP' + str(IOP) + '/contours_2019/circulation/circulations_' + str(theta) + 'K_' + strapp + '.nc')
        cinf.extend(cube)
    cinfm = cinf.merge()
    
    for i, name in enumerate(['area', 'mass']):
        
        cube = cinfm.extract(iris.Constraint(name = name))[-1][:, 0]
        # extract only start time
        
        plt.subplot(2, 2, i+1)
        
        plt.scatter(cube.data, cube.coord('air_potential_temperature').points)
        plt.xlabel(name)
        plt.ylabel('theta')
        
    print cinfm
        
    #cube1 = cinfm.extract(iris.Constraint(name = 'circulation'))[0][:, 0]
    cube2 = cinfm.extract(iris.Constraint(name = 'mass_integrated_circulation'))[-1][:, 0]
    
    plt.subplot(2, 2, 3)
    
    #plt.scatter(cube1.data, cube1.coord('air_potential_temperature').points, label = 'circuit_integrated', color = 'r', edgecolor = 'face')
    plt.scatter(cube2.data, cube2.coord('air_potential_temperature').points, label = 'mass_integrated')
    plt.xlabel('circulation')
    plt.ylabel('theta')
    
    plt.subplot(2, 2, 4)

    plt.scatter(cube2.data/cube.data, cube2.coord('air_potential_temperature').points)
    
    plt.xlim(0, 1e-6)
    
    plt.xlabel('circulation/mass')
    plt.ylabel('theta')
    
    plt.suptitle('inflow_' + strapp)
    #plt.legend()
    plt.savefig('/home/users/bn826011/NAWDEX/From N Drive/2019_figs/IOP' + str(IOP) + '/inflow_integrals_' + strapp + '.jpg')
    plt.show()
    
    return cinfm
    
    
        
def plot_timeseries_new(case = 0, levs = 3, tls = [320, 325, 310, 310], folder = 'IOP3/T42', strapp = '', ts = 2):
    "Based on leo's plot_timeseries"
    # figure 7
    
    save_dir = '/storage/silver/scenario/bn826011/WCB_outflow/Final/'
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
    
    plt.figure(figsize=(14, 8)) # JB impose size to make time labels readable
    
    for i in xrange(levs):
        # three theta levels par case, also subject to change
        theta = tls[case] + 5*i
        # theta levels are spaced by 5
        cubes = iris.load(save_dir + folder + '/circulation/circulations_' + str(theta) + 'K_' + strapp + '.nc')
        
        mass = cubes.extract(iris.Constraint(name = 'mass'))[0]
        area = cubes.extract(iris.Constraint(name = 'area'))[0]
        circ = cubes.extract(iris.Constraint(name = 'mass_integrated_circulation'))[0]
        
        time = mass.coord('time').points - 409608 # for IOP3 this gives time from ts
        
        dq =  1 - mass.data[ts]/mass.data
        dz =  1 - area.data[ts]/area.data
        dr = (mass.data*area.data[ts])/(mass.data[ts]*area.data) - 1
        dc = circ.data/mass.data
        
        plt.subplot(2, 3, 1)
        plt.plot(time, dq, color = colours[0][i], linewidth = 2.5)
        plt.title("-q2'/qs")
        #plt.title(folder + '  Circulation/Mass')
        plt.ylabel('1 - Ms/M2D')
        plt.subplot(2, 3, 2)
        plt.plot(time,dz,  color = colours[0][i], linewidth = 2.5)
        plt.title("-zeta2'/zetas")
        plt.ylabel('1 - As/A2D')
        plt.subplot(2, 3, 3)
        plt.plot(time, dr, color = colours[0][i], linewidth = 2.5)
        plt.title("r2'/rs")
        plt.ylabel('1 - M2D/Ms As/A2D')
        
#        plt.subplot(2, 3, 4)
#        plt.plot(time, dc, color = colours[0][i], linewidth = 2.5, label = str(theta) + 'K')
#        plt.title("circulation_div_by_mass")
#        plt.ylabel('C2D/M2D')
#        plt.xlabel('time from ts')
        
        alpha = dq/(1 - dq)
        epsilon = dz/(1 - dz)
        rder = (1 + alpha)/(1 + epsilon) - 1
        
        plt.subplot(2, 3, 4)
        plt.plot(time, alpha, color = colours[0][i], linewidth = 2.5)
        plt.title("alpha")
        plt.ylabel('dq/(1-dq)')
        plt.xlabel('time from ts')
        plt.subplot(2, 3, 5)
        plt.plot(time, epsilon,  color = colours[0][i], linewidth = 2.5)
        plt.title("epsilon")
        plt.ylabel('dz/(1-dz)')
        plt.xlabel('time from ts')
        plt.subplot(2, 3, 6)
        plt.plot(time, alpha/epsilon, color = colours[0][i], linewidth = 2.5)
        plt.title("alpha/epsilon")
        plt.ylabel('alpha/epsilon')
        plt.xlabel('time from ts')
        plt.ylim(-1, 3)
        
        
    plt.legend()
        
    #plt.savefig('/home/users/bn826011/NAWDEX/From N Drive/2019_figs/' + folder[:4] + '/_alpha_epsilon.png')
    plt.show()
    
    
def movement_compare(t2, t3, traj23diffs):
    # maps of a select few pairs of 2d and 3d trajectories from the same point to illustrate the positions within the 
    # outflow volume from which there are the largest and least difference in inflow position
    
    pvts = iris.load('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/20160922_12/prodm_op_gl-mn_20160922_12_c012_thsfcs_5K.nc', 'ertel_potential_vorticity')[0][0][9]
    # pv at ts on 325K
    pvtf = iris.load('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/20160922_12/prodm_op_gl-mn_20160922_12_c042_thsfcs_5K.nc', 'ertel_potential_vorticity')[0][0][9]
    # pv at tf on 325K
    
    
    
    move_mid = np.nonzero((traj23diffs[:, 2, 0] > 25)*(traj23diffs[:, 2, 0] < 100)*(traj23diffs[:, 2, 3] > 18)*(traj23diffs[:, 2, 3] < 22))[0][::16]
    # 8 trajectories move between 18 - 22 degrees
    move_most = np.nonzero((traj23diffs[:, 2, 0] > 25)*(traj23diffs[:, 2, 0] < 100)*(traj23diffs[:, 2, 3] > 25))[0][::40]#29.35))[0]
    # 8 trajectories move > 29.35 degrees
    move_least = np.nonzero((traj23diffs[:, 2, 0] > 25)*(traj23diffs[:, 2, 0] < 100)*(traj23diffs[:, 2, 3] < 15))[0][::11]#10))[0]
    # 8 trajectories move < 10 degrees
    
    for move in [[move_most, 'most'], [move_mid, 'mid'], [move_least, 'least']]:
        
        plt.figure(figsize = (13, 7))
        
        iplt.contour(pvts, [2], colors = [[.5, .5, .5]], linewidths = [.5])
        iplt.contour(pvtf, [2], colors = ['k'])
        
        for i in move[0]:
            plt.plot(t3.data[i, :, 0]-360, t3.data[i, :, 1], color = 'r', linewidth = 3)
        for i in move[0]:
            plt.plot(t2.data[i, :, 0]-360, t2.data[i, :, 1], color = 'c', linestyle = '--', linewidth = 2)
            
        plt.title('Difference in path taken by isentropic v full trajectories')
        plt.savefig('/home/users/bn826011/NAWDEX/From N Drive/2019_figs/IOP3/Path_differences_'+move[1]+'.png')
        plt.show()
        
        
def this_would_do_integrals():
    # a function I wrote to "conveniently" bulk calcualte circulation integrals on inflow volumes
    
    #for kfs in [[1, 'IOP5/T36', '', 325, 5, '0926_18'], [2, 'IOP6', '0', 310, 5, '1001_00'], [3, 'IOP7', '', 310, 4, '1003_12']]:
    for kfs in [[0, 'IOP3/T42', '', 320, 5, '0923_00'], [2, 'IOP6', '0', 310, 5, '1001_00'], [3, 'IOP7', '', 310, 4, '1003_12']]:
        folder = kfs[1][:4] + '/contours_2019'
        sp2 = ['wider_gauss']
        for i, strapp in enumerate(['inflow_bound_3D_less_than_' + str(kfs[3]) + 'K_gmod']):
            contour = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + kfs[1] + '/' +  strapp + '.npy')
            for thlev in range(285, kfs[3], 5):
                points = np.concatenate((contour, thlev*np.ones([len(contour), 1])), axis = 1)
                if thlev == 285:
                    inflow_volume = points
                else:
                    inflow_volume = np.concatenate((inflow_volume, points))

            np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/' + kfs[5] + sp2[i] + '.npy', inflow_volume)
            
            try_caltra_rev(kfs[0], levs = int((kfs[3] - 285)/5), hoursafterinit = [12, 6, 12, 0], hoursuntil = [18, 12, 18, 6], folder = folder, strapp = sp2[i], theta_max = kfs[3] - 5)
     
            for theta_level in range(285, kfs[3], 5):
                main_cc(kfs[0], theta_level = theta_level, folder = folder, strapp = sp2[i])
                
                
                
def plot_ridge_build_IOP6_3D(indiv = False, azm = -60, elv = 40):
    # produce a 3D animation of the buliding of scandinavian blocking ridge, split into 4 seperate air masses in different colours ascending at different time because it looks cool
    
    basetime = datetime.datetime(2016, 9, 30, 12)
    basetime_str = basetime.strftime('%Y%m%d_%H')
    
    TrEnList = []

    for k in xrange(4):#5
        
        #TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP6/{}_TrajectoryEnsemble_backward'.format(basetime_str) + str(k))
        #TrF = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP6/{}_TrajectoryEnsemble_forward'.format(basetime_str) + str(k))
        TrB = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/IOP6/inflow/{}_3DTrajectoryEnsemble'.format(basetime_str) + str(k))
        TrF = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/IOP6/inflow/{}_3DTrajectoryEnsemble_fw'.format(basetime_str) + str(k))
    
        TrEn = TrB.__add__(TrF)
        
        TrEnList.append(TrEn)
    
    tsteps = len(TrEnList[0].times)
    
    theta_0 = TrEnList[0].data[0, 0, 2]
    
    V = np.linspace(0, 40, 9)
    
#    cldt = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/' + basetime_str +
#    '/prodm_op_gl-mn_' + basetime_str + '_b*_thsfcs_5K.nc', 'total_minus_adv_only_theta')
#    cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
#    dtcube = cldt.concatenate_cube()
#    # create cube of delta theta
#    clpv = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/' + basetime_str +
#    '/prodm_op_gl-mn_' + basetime_str + '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
#    clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
#    pvcube = clpv.concatenate_cube()
    # and for PV
    
#    fig = plt.figure()
#    ax = fig.gca(projection='3d')
    
    scatarray = []
        
    for t in xrange(tsteps):
        print t
        scatarray.append([])
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
            
            #dtheta = dtcube.extract(iris.Constraint(time = TrEnList[0].times[t], air_potential_temperature = theta))
            
#        if t == tsteps-1:
#            
#            trop_theta = dtcube.extract(iris.Constraint(time = TrEnList[0].times[t], ertel_potential_vorticity = 2))
#                
#            lon = pvort.coord('longitude').points
#            lat = pvort.coord('latitude').points
#            x, y = np.meshgrid(lon, lat)
#            z = pvort.data
#                
#            ax.contour(x, y, z)
            
        #colour = [t/tsteps, 1-t/tsteps, .5*(1-t/tsteps)]
        #colour = [[.5*(t+1)/tsteps, 0, .1*(t+1)/tsteps], [0, 0, (t+1)/tsteps], [.9*(t+1)/tsteps, .5*(t+1)/tsteps, 0], [0, (t+1)/tsteps, 0]]
        colour = ['Oranges', 'Blues', 'Greys', 'Greens']   
           
        hoursuntilstart = [0, 42, 66, 60, 132]#[0, 0, 0, 0, 0]
        # times at which I want each contour to appear
        
        norm = matplotlib.colors.Normalize(-20, 40)
            
        for k in xrange(4):#5
                
            tidx = t - hoursuntilstart[k]/6
            #adjusted time from start of respective Trajectory Ensemble
                
            if tidx >= 0:
            
                TrEnindomain = TrEnList[k].select('air_potential_temperature', '>', 275, time = [datetime.timedelta(hours=t*6)])
                
                if TrEnindomain.__len__() != 0:
            
                    scatarray[t].append(ax.scatter(TrEnindomain.data[:, t, 0]-360, 
TrEnindomain.data[:, t, 1], TrEnindomain.data[:, t, 2], c = TrEnindomain.data[:, t, 3]-TrEnindomain.data[:, hoursuntilstart[k]/6, 3], 
marker = '.', edgecolors = 'face', alpha = .5, cmap = colour[k], norm = norm))

        ax.set_xlim([-100, 60])
        ax.set_ylim([10, 90])
        ax.set_zlim([-2000, 14000])
        ax.view_init(azim = azm, elev = elv)
        plt.savefig('/home/users/bn826011/NAWDEX/IOP6_series/all_3D_' + str('{0:02d}'.format(t)) + '_azm' + str(azm) + '.jpg')

        plt.show()
    
    return scatarray
    
                    
                    
            
        #plt.title(str(theta) + ' K')
                
#        timedatenow = TrEnList[0].times[t]
#                
#        plt.text(-90, 70, timedatenow.strftime("%d/%m %HUTC"), rotation = 'vertical')
#                
#        #label rows with dates and times
                
                
                
def PV_change_hist(k = 0, levs = 3, theta_0 = [320, 325, 310, 310], t = 0, folder = 'IOP3/T42', strapp = '', var_index = 5):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    #Tr3 = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble'.format(basetime_str) + strapp)
    Tr3 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_dPV'.format(basetime_str) + strapp)   
    
    bins = np.linspace(-5, 5, 21)
    
    plt.figure(figsize = (13, 9))
    
    Tall = Tr3.select('air_potential_temperature', '!=', -1000)
    #only those which remain in domain
    
    plt.subplot(int(levs/2)+1, 2, 1)
    curve = plt.hist(Tall.data[:, 0, var_index] - Tall.data[:, -(t+1), var_index], bins = bins)
    #plt.title(folder + '_' + strapp)
    plt.title('Started at all surfaces')
    plt.xlabel('Change in PV')
    
    for i in xrange(levs):
        
        theta = theta_0[k] + i*5
        
        Tt = Tall.select('air_potential_temperature', '>=', theta - 2.5, time = [datetime.timedelta(hours = 0)])
        Tt = Tt.select('air_potential_temperature', '<', theta + 2.5, time = [datetime.timedelta(hours = 0)])
        # those which remain in domain that start on desired theta surface
        
        plt.subplot(int(levs/2)+1, 2, 2+i)
        plt.hist(Tt.data[:, 0, var_index] - Tt.data[:, -(t+1), var_index], bins = bins)
        plt.title('Started at ' + str(theta) + 'K surface')
        plt.xlabel('Change in PV')
        
    #plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/PV_change_histograms_t=' + str(t) + strapp + '.png')
    #plt.savefig('IOP3_T42_2x2_hist_t=' + str(t) + '.png')
    #plt.savefig('/home/users/bn826011/NAWDEX/From N Drive/2019_figs/' + folder[:4] + 'PV_change_histograms_t=' + str(t) + strapp + '.png')
    plt.show()
    
    return curve
