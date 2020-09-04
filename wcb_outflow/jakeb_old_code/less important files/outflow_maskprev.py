# -*- coding: utf-8 -*-
"""
Created on Mon May 14 14:33:53 2018

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
from matplotlib.path import Path

from lagranto import caltra

from lagranto.trajectory import Trajectory, TrajectoryEnsemble, load

from mymodule import grid

import matplotlib.patches as patches

PVU = 2 # location of the dynamical tropopause


def compare_filter_size_AB(k=0):

    datestr = ['0922_12', '0926_12', '0930_12', '0930_12', '0930_12']

    thetastart = [6, 9, 6, 6, 6]#8, 9

    hoursafterinit = [60, 42, 42, 72, 150]#24, 42
    

    V = np.linspace(-10, 40, 11)

    filtersize = [40]#[1, 10, 20, 30, 40, 50, 60]
    
    fsz = len(filtersize)

    fig = plt.figure(figsize=(15, 10*fsz + 1))#15

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
    
    
    #################################################
    TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP6/2016{}_TrajectoryEnsemble_backward3'.format(datestr[k]))
    TrF = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP6/2016{}_TrajectoryEnsemble_forward3'.format(datestr[k]))
    
    TrEn = TrB.__add__(TrF)
    
    start_time = [-1*TrB.relative_times[-1]]
    ###############################################
    

    dtcub = dtcube[hoursafterinit[k]/3, thetastart[k]:thetastart[k]+5]
    pvcub = pvcube[hoursafterinit[k]/3, thetastart[k]:thetastart[k]+5]
    
    for i in xrange(3):
        
        for j in xrange(fsz):#7
            
            dt = dtcub.copy()[i]
            pv = pvcub.copy()[i]
            
            #pv.data = filters.median_filter(pv.data, size = filtersize[j])
            #dt.data = filters.median_filter(dt.data, size = filtersize[j])
            
            ##########################################
            trajectories = TrEn.select(
                'air_potential_temperature', '==', 280 + 5*(thetastart[k]+i), time = start_time)

            glon, glat = grid.get_xy_grids(dt)
            gridpoints = np.array([glon.flatten(), glat.flatten()]).transpose()
    
            x = trajectories.x[:, hoursafterinit[k]/6]
            y = trajectories.y[:, hoursafterinit[k]/6]
    
            plt.subplot(5*fsz, 3, 12*j+i+1)
            plt.scatter(x, y, c = 'k', marker = '.')
    
            # Include points within circuit boundary
            points = np.array([x, y]).transpose()
            pth = Path(points)
            
            print pth
            
            ax = fig.add_subplot(5*fsz, 3, 12*j+i+4)
            patch = patches.PathPatch(pth)
            ax.add_patch(patch)

            # Mask all points that are contained in the circuit
            mask = pth.contains_points(gridpoints).reshape(glat.shape)
            
            plt.subplot(5*fsz, 3, 12*j+i+7)
            plt.contourf(mask)
            
            plt.subplot(5*fsz, 3, 12*j+i+10)
            iplt.contourf(dt, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
            plt.plot(trajectories.data[:, hoursafterinit[k]/6, 0]-360, trajectories.data[:, hoursafterinit[k]/6, 1], color = 'k', marker = '.', ms = .8, mec = 'k', mfc = 'k')
    
            dt = dt.copy()
            #dt.data = np.ma.masked_where(mask, dt.data)
            dt.data = dt.data - 100*mask
            #pv = pv.copy()
            #pv.data = np.ma.masked_where(mask, pv.data)
            #########################################
            # Apply mask to the data
            # Have made decision to apply mask before and after filtering
            # and instead of traditional mask, trying to take 100 off affected areas instead
            
            pv.data = filters.median_filter(pv.data, size = filtersize[j])
            dt.data = filters.median_filter(dt.data, size = filtersize[j])
            
            dt = dt.copy()
            dt.data = dt.data - 100*mask
            
            plt.subplot(5*fsz, 3, 6*j+i+13)#7
            
            iplt.contourf(dt, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
            iplt.contour(pv, [PVU], colors = ['m'])
            iplt.contour(dt, [0], colors = ['g'])
            
            criteria = np.logical_and(pv.data < 2, dt.data > 0)
            # from Leo's code
        
            pv.data = criteria.astype(int)        
        
            cs = iplt.contour(pv, [0.5], colors = ['k'])
            # additional requirement of pv < 2 on contour
            
            if j == 0:
                
                potemp = dtcube.coord('air_potential_temperature')
                
                plt.title(str(potemp.points[thetastart[k]+i]) + ' K')
                #label columns with isentropic levels
                
            if i == 0:
                
                plt.text(-90, 70, 'filter size = ' + str(filtersize[j]), rotation = 'vertical')
                #label rows with filter size
                
    reftime = datetime.datetime(1970, 1, 1)
                
    hourssince = dtcub.coord('time').points[0]
                
    timedatenow = reftime + datetime.timedelta(hours = hourssince)
            
    plt.savefig('2016' + timedatenow.strftime("%d%m_%H") + '_delta_theta_' + str(PVU) + 'PVU_outflow__dtfilter_5.jpg')
    
            
def increase_nodes(points, resolution):
    "Taken from l.saffin old outflow_region.py, to avoid stupid __init__ files"
    "This has since been updated to account for different weightings at different latitudes"
    "I should also update this"
    
    n = 0
    while n < len(points) - 1:
        x1, y1 = points[n]
        x2, y2 = points[n + 1]
        dr = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        if dr > resolution:
            theta = np.arctan2(x2 - x1, y2 - y1)
            xp = x1 + resolution * np.sin(theta)
            yp = y1 + resolution * np.cos(theta)

            points = np.insert(points, n + 1, [xp, yp], axis=0)
        n += 1

    return points
    

def min_spac(points):
    "return minimum l2 spacing between two points along contour"
    
    dr = np.zeros(len(points)-1)
    
    for n in xrange(len(points)-1):
        x1, y1 = points[n]
        x2, y2 = points[n + 1]
        dr[n] = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
    return min(dr)
    
def len_con(points):
    "return sum of l2 spacings between all pairs of adjacent points along contour"
    "Probably needs changing to account for different weightings at different latitudes"
    
    conlen = 0
    
    for n in xrange(len(points)-1):
        x1, y1 = points[n]
        x2, y2 = points[n + 1]
        conlen += np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
    return conlen
            
    
def outflow_AB(k=0, thlevs = 3, thetastart = [6, 9, 6, 6, 6, 6], hoursafterinit = [60, 42, 42, 72, 96, 108], 
               filtersize = 30, resolution = 0.05, masks = 1, folder = 'IOP3/ABsplit', strapp = ''):
    "Returns numpy array of indices on contour around WCB outflow"
    "k is an index for case, of 0, 1, 2, 3, 4, 5"
    "thlevs is the number of isentropic levels considered"
    "thetastart is an array of four initial isentropic levels considered"
    "isentropic level = 280 + thetastart*5"
    "filtersize is size for median filter"
    "resolution is minimum L2 distance between two points on returned contour"
       

    datestr = ['0922_12', '0926_12', '0930_12', '0930_12', '0930_12', '0930_12', '0930_12']

    #hoursafterinit = [60, 42, 42, 72, 96, 96]# 42 for all 21, 18
    # start time in hours after first time in dataset
    #filtersize = [30, 30, 30, 30, 20, 35], the 20 & 35 to deal with
    #two separate outflows incredibly close together
    
    cldt = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_b*_thsfcs_5K.nc', 'total_minus_adv_only_theta')
    cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
    dtcube = cldt.concatenate_cube()
    # create cube of delta theta
    
        
    clpv = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
    pvcube = clpv.concatenate_cube()
    # and for PV
    
    #################################################
    TrEn = []
    start_time = []
    
    if k == 6:        
        rg = [3]        
    else:        
        rg = xrange(masks)

    for TE in rg:
    
        TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/2016{}_TrajectoryEnsemble_backward'.format(datestr[k]) + str(TE))
        TrF = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/2016{}_TrajectoryEnsemble_forward'.format(datestr[k]) + str(TE))
    
        TrEn.append(TrB.__add__(TrF))
    
        start_time.append([-1*TrB.relative_times[-1]])
    ###############################################
    
    # for labelling    
    potemp = dtcube.coord('air_potential_temperature')
    
    reftime = datetime.datetime(1970, 1, 1)
               
    hourssince = dtcube.coord('time').points[hoursafterinit[k]/3]
                
    timedatenow = reftime + datetime.timedelta(hours = hourssince)
    # for labelling
            
    plt.figure(figsize = (10, 14))
    
    
    V = np.linspace(-10, 40, 11)    
    
    outflow_volume = np.array([])
    
    for i in xrange(thlevs):
        
        cutoff = 250
        endoff = -1
        
        filtersize_t, filtersize_p = filtersize, filtersize
        
        if k == 3:
            cutoff = 50
            endoff = -50
        
        if k == 4: #3rd outflow into IOP6
            cutoff = 200
            filtersize_p = 20 #Pv division
            filtersize_t = 20 #narrow feature anyway
        
        if k == 5: #outflow into IOP7, bit of a weird one
            cutoff = 0
            endoff = -150
            filtersize_p = 5 #to maintain very fine PV structure dividing the ridges
            
        if k == 6: #why am I doing another outflow into IOP7?!!
            cutoff = 150
            
            
        dtcub = dtcube[hoursafterinit[k]/3, thetastart[k]+i, :, cutoff:endoff].copy()
        #restrict to outflow time and selected theta level
        #the cutoff artificially truncates the domain s.t. only the region we want is found
        #removes West: can't see new developing storms
        
        pvcub = pvcube[hoursafterinit[k]/3, thetastart[k]+i, :, cutoff:endoff].copy()
        #restrict to outflow time and selected theta level
        #the 150 artificially truncates the domain s.t. only the region we want is found
        
        
        ##########################################
        mask = []        
        
        for TE in xrange(masks):        
        
            trajectories = TrEn[TE].select(
                'air_potential_temperature', '==', 280 + 5*(thetastart[k]+i), time = start_time[TE])
#                'air_potential_temperature', '==', 280 + 5*(thetastart[k]+i), time = [datetime.timedelta(hours=hoursafterinit[k])])
#                option for if points which have left the domain are a problem, but it probably won't help                
            glon, glat = grid.get_xy_grids(dtcub)
            gridpoints = np.array([glon.flatten(), glat.flatten()]).transpose()
    
            x = trajectories.x[:, hoursafterinit[k]/6]
            y = trajectories.y[:, hoursafterinit[k]/6]
    
            # Include points within circuit boundary
            points = np.array([x, y]).transpose()
            pth = Path(points)

            # Mask all points that are contained in the circuit
            mask.append(pth.contains_points(gridpoints).reshape(glat.shape))
    
            dtcub = dtcub.copy()
            dtcub.data = dtcub.data - 100*mask[TE]
       #########################################
       # remove sections covered by previous outflow
        
        
    
        dtcub.data = filters.median_filter(dtcub.data, size = filtersize_t)
        #apply median filter to smooth field
    
        pvcub.data = filters.median_filter(pvcub.data, size = filtersize_p)
        #apply median filter to smooth field
        
        
        for TE in xrange(masks):
            dtcub.data = dtcub.data - 100*mask[TE]
            # remove mask before and after to ensure no overlap
        

        plt.subplot(thlevs, 2, 2*i+1)
    
        iplt.contourf(dtcub, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
        
        plt.text(-90, 70, timedatenow.strftime("%m%d_%H") + '_' + str(int(potemp.points[thetastart[k]+i])) + 'K', rotation = 'vertical')
        
        plt.gca().coastlines(color = [.8, .9, .8])        
        
        iplt.contour(pvcub, [2], colours = ['k'])
            
        criteria = np.logical_and(pvcub.data < 2, dtcub.data > 0)
        # from Leo's code
        
        pvcub.data = criteria.astype(int)        
        
        cs = iplt.contour(pvcub, [0.5], colors = ['g'])
        # additional requirement of pv < 2 on contour
    
        contours = trop.get_contour_verts(cs)
        # returns array of arrays of vertices for zero contours
    
        ncon = np.size(contours[0])
        #number of contours
    
        lencon = np.zeros(ncon)
        # empty array of lengths of contorus (in lat/long space)
    
        for j in xrange(ncon):
        
           lencon[j] = len_con(contours[0][j])
        
    
        imax = lencon.argmax()
        len_max = lencon[imax]
        # index of longest contour
        lcontour = contours[0][imax]
        # longest contour
        
        
        confail = False
    
        while not(trop.is_closed_contour(lcontour)):
            #Still need to check is_closed_contour works right
        
            lencon = np.delete(lencon, imax)
            del contours[0][imax]
            
            # this contour is not closed, remove & find the next one
        
            if np.size(lencon) == 0:
            
                confail = True
                # indicator that finding contour has failed
                
                break
                
            else:
        
                imax = lencon.argmax()
                len_max = lencon[imax]
                # new index of longest contour
                lcontour = contours[0][imax]
                # new longest contour

         
        
        plt.subplot(thlevs, 2, 2*i+2)
    
        iplt.contourf(dtcube[hoursafterinit[k]/3, thetastart[k]+i], V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
        
        t = hoursafterinit[k]/6        
        theta = (thetastart[k] + i)*5 + 280
        
        for TE in xrange(masks):
            
            TrEnindomain = TrEn[TE].select('air_potential_temperature', '==', theta, time = [datetime.timedelta(hours=6*t)])
            
            plt.plot(TrEnindomain.data[:, t, 0]-360, TrEnindomain.data[:, t, 1], color = [.8, .8, .8], marker = '.', ms = .8, mec = [.4, .4, .4], mfc = [.4, .4, .4])
        
        
        if confail == False:    
    
            points = increase_nodes(lcontour, resolution)
            # increase number of points on the contour such that they have minimum spacing resolution  
    
            plt.plot(points[:, 0]-360, points[:, 1], marker = '.', ms = 1, mec = 'k', mfc = 'k')
            
            points3d = np.zeros([points.shape[0], 3])
            
            points3d[:, 0:2] = points
            
            points3d[:, 2] = potemp.points[thetastart[k]+i]
            
            if i == 0:
            
                outflow_volume = points3d
            
            else:
    
                outflow_volume = np.concatenate((outflow_volume, points3d))
            
        else:
            
            plt.title('Unable to find contour around outflow at this level')
    
    #IOP = [3, 5, 6]
    
    np.save('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + timedatenow.strftime("%m%d_%H") + strapp + '.npy', outflow_volume) 
    
    plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + timedatenow.strftime("%m%d_%H") + strapp + '_image.jpg')
    
    plt.show()
    
    return
    
    
def plot_ridge_IOP3(indiv = False):
    
    basetime = datetime.datetime(2016, 9, 22, 12)
    basetime_str = basetime.strftime('%Y%m%d_%H')
    
    TrEnList = []

    for k in xrange(2):
        
        TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP3/ABsplit/{}_TrajectoryEnsemble_backward'.format(basetime_str) + str(k))
        TrF = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP3/ABsplit/{}_TrajectoryEnsemble_forward'.format(basetime_str) + str(k))
    
        TrEn = TrB.__add__(TrF)
        
        TrEnList.append(TrEn)
    
    tsteps = len(TrEnList[0].times)
    
    theta_0 = TrEnList[0].data[0, 0, 2]
    
    V = np.linspace(-10, 40, 11)
    
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
        plt.figure(figsize = (25, 40))
    
    for tht in xrange(5):
        
        theta = theta_0 + 5*tht
        
        for t in xrange(tsteps):
            
            if indiv == False:
                plt.subplot(tsteps, 5, t*5 + tht + 1)
            else:
                plt.figure(figsize = (10, 10))
                # plot frames individually
            
            dtheta = dtcube.extract(iris.Constraint(time = TrEnList[0].times[t], air_potential_temperature = theta))
            
            pvort = pvcube.extract(iris.Constraint(time = TrEnList[0].times[t], air_potential_temperature = theta))
            
            iplt.contourf(dtheta, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
            # plot contours of delta theta
            plt.gca().coastlines(color = [.6, .4, 0])
            # plot coastlines dark yellow?
            iplt.contour(pvort, [2], colors = ['m']) #, linewidths = 1.7)
            # plot 2 PVU contour, pink
            
            colour = [[0, .8, 0], 'k']
            
            hoursuntilstart = [0, 30]
            
            for k in xrange(2):
                
                if k == 1 and tht == 0:
                    
                    pass
                    # second contour not sensible on lowest theta level
                
                else:
                
                    tidx = t - hoursuntilstart[k]/6
                    #adjusted time from start of respective Trajectory Ensemble
                
                    if tidx >= 0:
            
                        TrEnindomain = TrEnList[k].select('air_potential_temperature', '==', theta, time = [datetime.timedelta(hours=t*6)]) # changed tidx to t as I actually started both ensembles from t = 0, but still don't want the mess of outflow B in the first 24 hours
                
                        if TrEnindomain.__len__() != 0:
            
                            plt.plot(TrEnindomain.data[:, t, 0]-360, TrEnindomain.data[:, t, 1], color = colour[k], marker = '.', ms = .8, mec = colour[k], mfc = colour[k])
            
            if indiv == True:
                plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP3/ABsplit/indv_plts/' + basetime_str + '_' + str(int(theta)) + 'K_t=' + str('{0:02d}'.format(t)) + '.png')
                # save individual frames            
            
            if t == 0:
                    
                plt.title(str(theta) + ' K')
                
            if tht == 0:
                
                timedatenow = TrEnList[0].times[t]
                
                plt.text(-90, 70, timedatenow.strftime("%d/%m %HUTC"), rotation = 'vertical')
                
                #label rows with dates and times
    if indiv == False:                        
        plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/IOP3/ABsplit/' + basetime_str + 'fulltraj_split.png')