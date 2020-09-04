# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 17:46:57 2018

@author: bn826011

file to run trajectories for IOP67 case, in a sensible way
"""

import numpy as np
import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import datetime

from lagranto import caltra

from lagranto.trajectory import load

def try_caltra_ridge(k = 0):
    
    # Create a mapping from datetime objects to filenames
    basetime = datetime.datetime(2016, 9, 30, 12)
    basetime_str = basetime.strftime('%Y%m%d_%H')
    datadir = '/export/cloud/NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    timestep = datetime.timedelta(hours=6)
    
    hoursafterinit = [42, 72, 96, 96]
    #time contour defined after 30/12
    hoursuntilstart = [0, 42, 66, 0]
    #time that feature appears, s.t. contour doesn't collapse to stupid nothing
    
    times = [basetime + timestep * i for i in range(hoursuntilstart[k]/6, hoursafterinit[k]/6 + 1)]
    print times[-1]
    #creates mapping up to and including selected outflow time    
    
    mapping = {}
    for t in times:
        leadtimehr = np.int((t - basetime).total_seconds()) / 3600
        fn = 'prodm_op_gl-mn_{0}_d{1:03d}_thgrid.pp'.\
            format(basetime_str, 12 * (leadtimehr / 12))
        mapping[t] = datadir + fn
        
    if k != 3:
        trainp_th = np.load('/glusterfs/msc/users_2018/bn826011/NAWDEX/IOP67ridge/' + times[-1].strftime("%m%d_%H") + '.npy')
        #trajectory input on theta levels
    else:
        trainp_th = np.load('/glusterfs/msc/users_2018/bn826011/NAWDEX/IOP67ridge/' + times[-1].strftime("%m%d_%H") + '_ridge2.npy')
    
    levels = ('air_potential_temperature', [trainp_th[0, 2] + 5*i for i in range(3)])
    # three isentropic levels including that of the first data point
    
    print levels
    
    tracers = ['altitude', 'x_wind', 'y_wind', 'upward_air_velocity', 'air_pressure']
    #need these for circulation integrals
    
    traout = caltra.caltra(trainp_th, mapping, fbflag=-1, nsubs = 12, tracers = tracers, levels=levels)
    # 12 steps between files = 30 mins apart 
        

    traout.save('/glusterfs/msc/users_2018/bn826011/NAWDEX/IOP67ridge/{}_TrajectoryEnsemble_backward'.format(basetime_str) + str(k+1))
    
    return traout
    

def try_caltra_ridge_rev(k = 0):
    
    # Create a mapping from datetime objects to filenames
    basetime = datetime.datetime(2016, 9, 30, 12)
    basetime_str = basetime.strftime('%Y%m%d_%H')
    datadir = '/export/cloud/NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    timestep = datetime.timedelta(hours=6)
    
    hoursafterinit = [42, 72, 96, 96]
    hoursuntil = [156, 156, 156, 156]
    
    times = [basetime + timestep * i for i in range(hoursafterinit[k]/6, hoursuntil[k]/6 + 1)]
    print times[0]
    #creates mapping from selected outflow time until selected end time 
    
    mapping = {}
    for t in times:
        leadtimehr = np.int((t - basetime).total_seconds()) / 3600
        fn = 'prodm_op_gl-mn_{0}_d{1:03d}_thgrid.pp'.\
            format(basetime_str, 12 * (leadtimehr / 12))
        mapping[t] = datadir + fn
        
    if k != 3:
        trainp_th = np.load('/glusterfs/msc/users_2018/bn826011/NAWDEX/IOP67ridge/' + times[0].strftime("%m%d_%H") + '.npy')
        #trajectory input on theta levels
    else:
        trainp_th = np.load('/glusterfs/msc/users_2018/bn826011/NAWDEX/IOP67ridge/' + times[0].strftime("%m%d_%H") + '_ridge2.npy')
    
    levels = ('air_potential_temperature', [trainp_th[0, 2] + 5*i for i in range(3)])
    # three isentropic levels including that of the first data point
    
    print levels
    
    tracers = ['altitude', 'x_wind', 'y_wind', 'upward_air_velocity', 'air_pressure']
    #need these for circulation integrals
    
    traout = caltra.caltra(trainp_th, mapping, fbflag=1, nsubs = 12, tracers = tracers, levels=levels)
    # 12 steps between files = 30 mins apart 
    

    traout.save('/glusterfs/msc/users_2018/bn826011/NAWDEX/IOP67ridge/{}_TrajectoryEnsemble_forward'.format(basetime_str) + str(k+1))
    
    return traout
    
    
def plot_ridge_build(indiv = False):
    
    basetime = datetime.datetime(2016, 9, 30, 12)
    basetime_str = basetime.strftime('%Y%m%d_%H')
    
    TrEnList = []

    for k in xrange(5):#4
        
        TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP6/{}_TrajectoryEnsemble_backward'.format(basetime_str) + str(k))
        TrF = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP6/{}_TrajectoryEnsemble_forward'.format(basetime_str) + str(k))
    
        TrEn = TrB.__add__(TrF)
        
        TrEnList.append(TrEn)
    
    tsteps = len(TrEnList[0].times)
    
    theta_0 = TrEnList[0].data[0, 0, 2]
    
    V = np.linspace(0, 40, 9)
    
    cldt = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/' + basetime_str +
    '/prodm_op_gl-mn_' + basetime_str + '_b*_thsfcs_5K.nc', 'total_minus_adv_only_theta')
    cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
    dtcube = cldt.concatenate_cube()
    # create cube of delta theta
    clpv = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/' + basetime_str +
    '/prodm_op_gl-mn_' + basetime_str + '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
    pvcube = clpv.concatenate_cube()
    # and for PV
    
    if indiv == False:
        plt.figure(figsize = (15, 80))
    
    for tht in xrange(3):
        
        theta = theta_0 + 5*tht
        
        for t in xrange(tsteps):
            
            if indiv == False:
                plt.subplot(tsteps, 3, t*3 + tht + 1)
            else:
                plt.figure(figsize = (10, 10))
                # for plotting individual frames
            
            dtheta = dtcube.extract(iris.Constraint(time = TrEnList[0].times[t], air_potential_temperature = theta))
            
            pvort = pvcube.extract(iris.Constraint(time = TrEnList[0].times[t], air_potential_temperature = theta))
            
            iplt.contourf(dtheta, V, cmap = 'binary')
            # plot contours of delta theta, gray
            plt.gca().coastlines(color = [.6, .4, 0])
            # plot coastlines dark yellow?
            iplt.contour(pvort, [2], colors = ['k'], linewidths = 1.7)
            # plot 2 PVU contour, black
            
            colour = [[.5, 0, .1], 'b', [.9, .5, 0], 'g', [1, 0, 1]]
            
            hoursuntilstart = [0, 0, 0, 0, 0]#[0, 42, 66, 0, 132]
            # times at which I want each contour to appear
            
            for k in xrange(5):#4
                
                tidx = t - hoursuntilstart[k]/6
                #adjusted time from start of respective Trajectory Ensemble
                
                if tidx >= 0:
            
                    TrEnindomain = TrEnList[k].select('air_potential_temperature', '==', theta, time = [datetime.timedelta(hours=t*6)])
                
                    if TrEnindomain.__len__() != 0:
            
                        plt.scatter(TrEnindomain.data[:, t, 0]-360, TrEnindomain.data[:, t, 1], s = 2, c = colour[k], marker = '.', edgecolor = colour[k])
            
            if indiv == True:
                plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP6/busy_plts/' + basetime_str + '_' + str(int(theta)) + 'K_t=' + str('{0:02d}'.format(t)) + '.png')
                # for plotting individual frames
            
            if t == 0:
                    
                plt.title(str(theta) + ' K')
                
            if tht == 0:
                
                timedatenow = TrEnList[0].times[t]
                
                plt.text(-90, 70, timedatenow.strftime("%d/%m %HUTC"), rotation = 'vertical')
                
                #label rows with dates and times
                
    if indiv == False:  
        pass                      
#        plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP6/' + basetime_str + 'fulltraj_all.png')