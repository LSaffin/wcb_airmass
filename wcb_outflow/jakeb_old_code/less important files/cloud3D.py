# -*- coding: utf-8 -*-
"""
Created on Sat Jul 28 17:30:45 2018

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

def plot_cloud(k = 0, ts = [[2, 4, 7], [1, 3, 6], [2, 4, 7], [0, 2, 4], [2, 6, 10]], levs = 3, filtersize = 30, field = 'air_potential_temperature', scatter3D = True, lim = 30, indiv = False, folder = 'IOP3/T42', strapp = ''):
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
    
    tsteps = 3
    
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
        plt.figure(figsize = (levs*5, tsteps*2.5), dpi = 600)
        #plt.figure(figsize = (12, 8), dpi = 300)
    
    for tht in xrange(levs):
        
        theta = theta_0 + 5*tht
        
        TrEntheta = TrEn.select('air_potential_temperature', '==', theta, time = [-1*TrB.relative_times[-1]])
        # This gives the levels at the designated start time
        
        if scatter3D != False:
            
            Tt = Tr3.select('air_potential_temperature', '>=', theta - 2.5, time = [datetime.timedelta(hours = 0)])
            Tt = Tt.select('air_potential_temperature', '<', theta + 2.5, time = [datetime.timedelta(hours = 0)])
            Tt = Tt.select(field, '!=', -1000)
            # points starting from theta level for which field is defined
                    
        for ti in xrange(tsteps):
            
            t = ts[k][ti]
            
            if indiv == False:
                plt.subplot(tsteps, levs, ti*levs + tht + 1)
            else:
                plt.figure(figsize = (10, 10))
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
                        
                        norm = matplotlib.colors.Normalize(theta-lim, theta+10)
                        
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
                plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/indv_plts/' + basetime_str + '_' + field + strapp + str(int(theta)) + 'K_t=' + str('{0:02d}'.format(t)) + '.png')
            #for plotting individual frames            
            
            if t == 0:
                
                if Leftdomain != 0:
                
                    plt.title(str(theta) + ' K \n ' + str(Leftdomain) + ' points have left the domain')
                    #label columns with isentropic levels
                    
                else:
                    
                    plt.title(str(theta) + ' K')
                
                
            if tht == 0:
                
                timedatenow = TrEn.times[t]
                
                plt.text(-90, 70, timedatenow.strftime("%d/%m %HUTC"), rotation = 'vertical')
                
                #label rows with dates and times
                            
    #plt.savefig('outflow/T42_mfs' + str(filtersize) + '/' + basetime_str + 'fulltraj_' + field + '.jpg')
    if indiv == False:
        
        #plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + basetime_str + 'fulltraj_' + field + strapp + '.png')
        plt.savefig('easy_to_find/' + folder + strapp + '_3Dcloud.png', dpi = 600)