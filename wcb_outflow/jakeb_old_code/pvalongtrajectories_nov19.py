#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 15:09:11 2019

@author: bn826011
"""
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from lagranto.trajectory import Trajectory, TrajectoryEnsemble, load
import numpy as np

def PV_along_trajectories(folder = 'IOP3/T42', time_string = '20160922_12', name = '',
                          theta_min = 24, no_of_trajs = 10, plotnotmean = True):
    
    TrEn = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/'+
                 folder + '/inflow/' + time_string + '_3DTrajectoryEnsemble_new' + name)
    # load arbitrary set of 3D trajectories
    
    times = TrEn.times
    # get array containing pertinent times
    
    clpv = iris.load('/storage/silver/NCAS-Weather/ben/nawdex/mi-ar482/' + time_string +
                     '/prodm_op_gl-mn_' + time_string + '_c*_thsfcs.nc', 'ertel_potential_vorticity')
    clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
    pvcube = clpv.concatenate_cube()
    # load full 3D PV fields for corresponding case
    # could restrict this to pertinent times to save processing time
    
    cldt = iris.load('/storage/silver/NCAS-Weather/ben/nawdex/mi-ar482/' + time_string +
                     '/prodm_op_gl-mn_' + time_string + '_b*_thsfcs.nc', 'total_minus_adv_only_theta') #  '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
    dtcube = cldt.concatenate_cube()
    # same for diabatic heating proxy
    
    delta_lat = np.mean(np.diff(pvcube.coord('latitude').points[:10]))/2
    # spacing of latitude grid
    delta_lon = np.mean(np.diff(pvcube.coord('longitude').points[:10]))/2
    # spacing of longitude grid
     
    trajectory_bin = []
    for traj in TrEn:
        if abs(traj.data[0, 3] - traj.data[-1, 3]) > theta_min and min(traj.data[:, 3]) > 300:
            trajectory_bin.append(traj)
    # make a list of trajectories which ascend the most
    # NOTE: the data I have for some reason only goes down to 300K - possible drawback
            
    n = int(max(np.floor(len(trajectory_bin)/no_of_trajs), 1))
    # interval of selection based on desired number of trajectories    
    
    for figno, trajex in enumerate(trajectory_bin[::n]):
        
        lat = trajex.data[:, 1]
        lon = trajex.data[:, 0]
        theta = trajex.data[:, 3]
        
        pvs = []
        dts = []
        
        for i in range(len(times)):
            lat_constraint = iris.Constraint(latitude=lambda cell: lat[i] - delta_lat < cell < lat[i] + delta_lat)
            lon_constraint = iris.Constraint(longitude=lambda cell: lon[i] - delta_lon < cell < lon[i] + delta_lon)
            time_constraint = iris.Constraint(time = times[i])
            pvs.append(pvcube.extract(lat_constraint & lon_constraint & time_constraint))
            dts.append(dtcube.extract(lat_constraint & lon_constraint & time_constraint))
            
        ### hack fix for points not being found    
        ncl = []
        tcl = []
        try:
            for cube in pvs:
                if cube.ndim == 1:
                    ncl.append(cube)
                elif cube.ndim == 2:
                    ncl.append(cube[:, 0])
                else:
                    ncl.append(cube[:, 0, 0])
            ### hack fix for points not being found 
            for cube in dts:
                if cube.ndim == 1:
                    tcl.append(cube)
                elif cube.ndim == 2:
                    tcl.append(cube[:, 0])
                else:
                    tcl.append(cube[:, 0, 0])
            ### hack fix for points not being found 
                
            pvtrajcubes = iris.cube.CubeList(ncl)
            dttrajcubes = iris.cube.CubeList(tcl)
            
            pvmerge = pvtrajcubes.merge_cube()
            dtmerge = dttrajcubes.merge_cube()
            
            if plotnotmean:
            
                plt.figure(figsize = (12, 12))
                plt.subplot(2, 2, 1)
                qplt.contourf(pvmerge, np.linspace(-3, 3, 25), cmap = 'RdBu_r')
                plt.plot(times, theta)
                plt.subplot(2, 2, 2)
                qplt.contourf(dtmerge, np.linspace(-25, 25, 26), cmap = 'RdBu_r')
                plt.plot(times, theta)
                plt.subplot(2, 1, 2)
                qplt.contour(pvcube[14, 10], [2])
                plt.gca().coastlines()
                plt.plot(lon-360, lat, linewidth = 3)
                plt.savefig('PV_dtheta_trajectory_crosssection_' + str(figno) + '_' + time_string + '.png')
                plt.show()
        except AttributeError as e:
            print e
            
        else:
            
            if figno == 0:               
                pvarray = np.array([pvmerge.data])
                dtarray = np.array([dtmerge.data])
                thetarray = np.array([theta])
                # for the first profile, initialise a numpy array
            else:
                pvarray = np.append(pvarray, [pvmerge.data], axis=0)
                dtarray = np.append(dtarray, [dtmerge.data], axis=0)
                thetarray = np.append(thetarray, [theta], axis = 0)
                
    if not plotnotmean:
        
        lts = len(times)
        
        pvmean = np.mean(pvarray, axis = 0)
        dtmean = np.mean(dtarray, axis = 0)
        thetamean = np.mean(thetarray, axis = 0)
        # create mean fields along trajectories
        
        ytheta = np.repeat([np.linspace(300, 340, 17)], lts, axis = 0)
        xtime = np.repeat([np.linspace(0, (lts-1)*6, lts)], 17, axis = 0).T
        # create arrays for axes
        
        plt.figure(figsize = (12, 8))
        plt.subplot(1, 2, 1)
        plt.contourf(xtime, ytheta, pvmean, np.linspace(-3, 3, 25), cmap = 'RdBu_r')
        plt.plot(np.linspace((lts-1)*6, 0, lts), thetamean)
        plt.title('Average PV along trajectory for > 20K ascent')
        plt.xlabel('time from start, hours')
        plt.ylabel('theta, kelvin')
        plt.subplot(1, 2, 2)
        plt.contourf(xtime, ytheta, dtmean, np.linspace(-25, 25, 26), cmap = 'RdBu_r')
        plt.plot(np.linspace((lts-1)*6, 0, lts), thetamean)
        plt.title('Average diabatic heating')
        plt.xlabel('time, hours')
        plt.savefig('PV_dtheta_trajectory_crosssection_mean_' + time_string + '.png')
        plt.show()
                
                
            
    
            
            


    
