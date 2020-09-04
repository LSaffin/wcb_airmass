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
import operator

from lagranto import caltra

from lagranto.trajectory import Trajectory, TrajectoryEnsemble, load
import mymodule.plot as modplt

import iris.quickplot as qplt

PVU = 2 # location of the dynamical tropopause

new_filepath = '/storage/silver/NCAS-Weather/ben/'

def compare_pv_fields(k = 0, folder = 'IOP3/T42', strapp = ''):
    
    datestr = ['0922_12', '0926_12', '0930_12', '1003_12']
    tf = ['0924_06', '0928_00', '1002_06', '1004_12']
    IOP = ['3', '5', '6', '7']
    mint = [320, 325, 310, 310]
    tstart = [42, 36, 42, 24]
    
    fileno = tstart[k] - tstart[k]%12
    
    clpvo = iris.load_cube('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_c*' + str(fileno) + '_thsfcs_5K.nc', 'ertel_potential_vorticity')
    pvo = clpvo[(tstart[k]%12)/6]
    
    clpvi = iris.load_cube('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_c000_thsfcs_5K.nc', 'ertel_potential_vorticity')
    pvi = clpvi[0]
    
    trainp_th = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/' + tf[k] + strapp + '.npy')
    #shadow_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/'+ folder + '/shadow_contour.npy')
    bound_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow_bound_3D_less_than_' + str(mint[k]) + 'K.npy')
    
    
    plt.figure(figsize=(8, 17))
    
    V = np.linspace(-1, 5, 13)
    
    for n, theta in enumerate(np.arange(mint[k]-30, mint[k], 5)):
        
        plt.subplot(6, 2, 2*n + 1)
        
        qplt.contourf(pvi.extract(iris.Constraint(air_potential_temperature = theta)), V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-5, 5))
        plt.plot(bound_cont[:, 0]-360, bound_cont[:, 1], 'k')
        plt.title(str(theta) + 'K  ts')
        plt.gca().coastlines(color = [.5, .5, .5])
        plt.subplot(6, 2, 2*n + 2)
        
        qplt.contourf(pvo.extract(iris.Constraint(air_potential_temperature = theta+20)), V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-5, 5))
        if theta+20 in trainp_th[:, 2]:
            plt.plot(trainp_th[trainp_th[:, 2] == theta+20, 0]-360, trainp_th[trainp_th[:, 2] == theta+20, 1], 'k')
#        elif theta+20 == trainp_th[0, 2] - 5:
#            plt.plot(trainp_th[trainp_th[:, 2] == trainp_th[0, 2], 0]-360, trainp_th[trainp_th[:, 2] == trainp_th[0, 2], 1], 'k')
#        elif theta+20 == trainp_th[-1, 2] + 5:
#            plt.plot(trainp_th[trainp_th[:, 2] == trainp_th[-1, 2], 0]-360, trainp_th[trainp_th[:, 2] == trainp_th[-1, 2], 1], 'k')
        plt.title(str(theta+20) + 'K  tf')
        plt.gca().coastlines(color = [.5, .5, .5])
    plt.savefig('PV_field_compare_IOP' + IOP[k] + '.jpg')
    plt.show()


def compare_outflow_size_with_time(k=0, thetastart = [4], hoursafterinit = [12], llevs = 9, times = 7):
    # script to compare smoothed contours with different filter sizes on different levels
    # repurposed to do similar but at different times

    datestr = ['0922_12', '0926_12', '0922_12', '1002_00', '1003_12']

    #thetastart = [4, 8, 8, 6, 6]#8, 9, 6

    #hoursafterinit = [12, 72, 54, 36, 30]#24, 42 #60
    
#    llevs = 9    

    V = np.linspace(-10, 40, 11)

#    times = 7
    
    filtersize = 30

    plt.figure(figsize=(times*3 + 1, llevs*2))#15

    cldt = iris.load('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_b*_thsfcs_5K.nc', 'total_minus_adv_only_theta') #  '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
    dtcube = cldt.concatenate_cube()
    # create cube of delta theta
    clpv = iris.load('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
    pvcube = clpv.concatenate_cube()
    # create cube of delta theta

    dtcub = dtcube[hoursafterinit[k]/3:hoursafterinit[k]/3+4*times, thetastart[k]:thetastart[k]+llevs]
    pvcub = pvcube[hoursafterinit[k]/3:hoursafterinit[k]/3+4*times, thetastart[k]:thetastart[k]+llevs]
    
    for i in xrange(llevs):
        
        for j in xrange(times):#xrange(fsz):#7
            
            dt = dtcub.copy()[4*j, i]
            pv = pvcub.copy()[4*j, i]
            
            pv.data = filters.median_filter(pv.data, size = filtersize)
            dt.data = filters.median_filter(dt.data, size = filtersize)
            
            #plt.subplot(fsz, 3, 3*j+i+1)#7
            plt.subplot(llevs, times, times*i + j + 1)
            
            cx = iplt.contourf(dt, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
            iplt.contour(pv, [PVU], colors = ['m'])
            iplt.contour(dt, [0], colors = ['g'])
            
            criteria = np.logical_and(pv.data < 2, dt.data > 0)
            # from Leo's code
        
            pv.data = criteria.astype(int)        
        
            cs = iplt.contour(pv, [0.5], colors = ['k'])
            # additional requirement of pv < 2 on contour
            
            if i == 0:
                
                reftime = datetime.datetime(1970, 1, 1)
                
                hourssince = dt.coord('time').points[0]
                
                timedatenow = reftime + datetime.timedelta(hours = hourssince)
                
                plt.title(timedatenow.strftime("%d%m_%H"))
                #label columns with time
                
            if j == 0:
                
                potemp = dtcube.coord('air_potential_temperature')
                
                plt.text(-90, 70, str(potemp.points[thetastart[k]+i]) + ' K', rotation = 'vertical')
                #label rows with theta level
                
    reftime = datetime.datetime(1970, 1, 1)
                
    hourssince = dtcub[0].coord('time').points[0]
                
    timedatenow = reftime + datetime.timedelta(hours = hourssince)
            
    plt.savefig('2016' + timedatenow.strftime("%d%m_%H") + '_delta_theta_' + str(PVU) + 'PVU_outflow_dtheta_time_2019.jpg')
    
    
def filter_pv():
    
    datestr = '20160922_12'
    t = 14
    thst = 8
    
    cldt = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/' + datestr +
    '/prodm_op_gl-mn_' + datestr + '_b*_thsfcs_5K.nc', 'total_minus_adv_only_theta')
    cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
    dtcube = cldt.concatenate_cube()
    # create cube of delta theta
    clpv = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/' + datestr +
    '/prodm_op_gl-mn_' + datestr + '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
    pvcube = clpv.concatenate_cube()
    # and for PV
    
    fig = plt.figure(figsize=(15, 10))
    
    dtcub = dtcube[t, thst:thst+3]
    pvcub = pvcube[t, thst:thst+3]
    
    V = np.linspace(-10, 40, 11)
    
    for i in xrange(3):
        
        dt = dtcub.copy()[i]
        pv = pvcub.copy()[i]
        
        dt.data = filters.median_filter(dt.data, size = 20)
        pv.data = filters.median_filter(pv.data, size = 15)
        
        plt.subplot(4, 3, i+1)
            
        iplt.contourf(dtcub[i], V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
        iplt.contour(pvcub[i], [PVU], colors = ['k'])
        
        potemp = dtcube.coord('air_potential_temperature')
                
        plt.title(str(potemp.points[thst+i]) + ' K')
        
        #no filter
        
        plt.subplot(4, 3, i+4)
            
        iplt.contourf(dt, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
        iplt.contour(pvcub[i], [PVU], colors = ['k'])
        iplt.contour(dt, [0], colors = ['g'])
        
        # dt filter
        
        plt.subplot(4, 3, i+7)
            
        iplt.contourf(dtcub[i], V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
        iplt.contour(pv, [PVU], colors = ['k'])
        iplt.contour(dtcub[i], [0], colors = ['g'])
        
        # pv filter
        
        plt.subplot(4, 3, i+10)
            
        iplt.contourf(dt, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
        iplt.contour(pv, [PVU], colors = ['k'])
        
        criteria = np.logical_and(pv.data < PVU, dt.data > 0)
        # from Leo's code
        
        pv.data = criteria.astype(int)        
        
        iplt.contour(pv, [0.5], colors = ['g'])
        
    plt.savefig(datestr + '_filter_20dt_15pv.jpg')
            
def increase_nodes(points, resolution):
    "Taken from l.saffin old outflow_region.py, to avoid stupid __init__ files"
    "This has since been updated to account for different weightings at different latitudes"
    "I should also update this"
    
    n = -1
    while n == -1:
        x1, y1 = points[-1]
        x2, y2 = points[0]
        dr = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        if dr > resolution:
            theta = np.arctan2(x2 - x1, y2 - y1)
            xp = x1 + resolution * np.sin(theta)
            yp = y1 + resolution * np.cos(theta)

            points = np.insert(points, 0, [xp, yp], axis=0)
        else:
            n = 0
    # to increase resolution between first and last nodes
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
            

def outflow(k=0, i=1, filtersize = 30, resolution = 0.25, pvcd = False):
    "Returns numpy array of indices on contour around WCB outflow"
    #k is an index for which case, of 0, 1, 2, 3
    #i is and index for isentropic level, of 0, 1, 2
    "Note that this function only finds outflow contours on a single isentropic"
    "Level and is not actually very useful, for caltra or otherwise"

    datestr = ['0922_12', '0926_12', '0930_12']#, '1001_00']

    thetastart = [4]#, 8, 5, 4]
    # isentropic level = 280 + thetastart*5
    hoursafterinit = [42]#, 54, 36, 54]
    # start time in hours after first time in dataset
    
    cldt = iris.load('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_b*_thsfcs_5K.nc', 'total_minus_adv_only_theta')
    cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
    dtcube = cldt.concatenate_cube()
    # create cube of delta theta
    
    # for labelling    
    potemp = dtcube.coord('air_potential_temperature')
    
    reftime = datetime.datetime(1970, 1, 1)
               
    hourssince = dtcube.coord('time').points[hoursafterinit[k]/3]
                
    timedatenow = reftime + datetime.timedelta(hours = hourssince)
    # for labelling
    
    dtcub = dtcube[hoursafterinit[k]/3, thetastart[k]+i, :, 150:].copy()
    #restrict to outflow time and selected theta level
    #the 150 artificially truncates the domain s.t. only the region we want is found
    
    dtcub.data = filters.median_filter(dtcub.data, size = filtersize)
    #apply median filter to smooth field
    
    if pvcd == True:
        
        clpv = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
        '/prodm_op_gl-mn_2016' + datestr[k] + '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
        clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
        pvcube = clpv.concatenate_cube()
        # and for PV
        
        pvcub = pvcube[hoursafterinit[k]/3, thetastart[k]+i, :, 150:].copy()
        #restrict to outflow time and selected theta level
        #the 150 artificially truncates the domain s.t. only the region we want is found
    
        pvcub.data = filters.median_filter(pvcub.data, size = filtersize)
        #apply median filter to smooth field
    
    plt.figure(1)
    
    plt.subplot(2, 1, 1)
    
    V = np.linspace(-10, 40, 11)
    
    iplt.contourf(dtcub, V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
    
    if pvcd == False:
    
        cs = iplt.contour(dtcub, [0], colors = ['g'])
        #plot zero contours
        
    elif pvcd == True:
        
        iplt.contour(pvcub, [2], colors = ['k'])
        
        criteria = np.logical_and(pvcub.data < 2, dtcub.data > 0)
        # from Leo's code
        
        pvcub.data = criteria.astype(int)        
        
        cs = iplt.contour(pvcub, [0.5], colors = ['g'])
        # additional requirement of pv < 2
    
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
    
    while not(trop.is_closed_contour(lcontour)):
        
        lencon = np.delete(lencon, imax)
        # this contour is not closed, remove & find the next one
        
        if np.size(lencon) == 0:
            
            plt.subplot(2, 1, 2)
    
            iplt.contourf(dtcube[hoursafterinit[k]/3, thetastart[k]+i], V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
            
            plt.text(-90, 70, timedatenow.strftime("%m%d_%H") + '_' + str(int(potemp.points[thetastart[k]+i])) + 'K', rotation = 'vertical')
            
            plt.savefig('outflow/' + timedatenow.strftime("%m%d_%H") + '_' + str(int(potemp.points[thetastart[k]+i])) + 'K_image.jpg')
            
            plt.show()
            
            print 'Unable to find contour around outflow at this level'
            
            return
        
        imax = lencon.argmax()
        len_max = lencon[imax]
        # new index of longest contour
        lcontour = contours[0][imax]
        # new longest contour

    
    points = increase_nodes(lcontour, resolution)
    # increase number of points on the contour such that they have minimum spacing resolution    
    
    if pvcd == False:
    
        np.save('outflow/' + timedatenow.strftime("%m%d_%H") + '_' + str(int(potemp.points[thetastart[k]+i])) + 'K.npy', points)
        
    else:
        
        np.save('outflow/' + timedatenow.strftime("%m%d_%H") + '_' + str(int(potemp.points[thetastart[k]+i])) + 'K_wpv.npy', points)
        

    plt.subplot(2, 1, 2)
    
    iplt.contourf(dtcube[hoursafterinit[k]/3, thetastart[k]+i], V, cmap = 'seismic', norm = matplotlib.colors.Normalize(-40, 40))
    
    plt.text(-90, 70, timedatenow.strftime("%m%d_%H") + '_' + str(int(potemp.points[thetastart[k]+i])) + 'K', rotation = 'vertical')
    
    plt.plot(points[:, 0]-360, points[:, 1], marker = '.', mec = 'k', mfc = 'k', ms = 0.5)
    
    plt.savefig('outflow/' + timedatenow.strftime("%m%d_%H") + '_' + str(int(potemp.points[thetastart[k]+i])) + 'K_image.jpg')
    
    plt.show()
    
    return
    
    
def outflow_th(k=0, thlevs = 3, thetastart = [8, 9, 6, 6], hoursafterinit = [42, 42, 42, 24], filtersize = 30, resolution = 0.05,
               folder = 'IOP3/contours_2019', strapp = ''):
    "Returns numpy array of indices on contour around WCB outflow"
    "k is an index for case, of 0, 1, 2, 3"
    "thlevs is the number of isentropic levels considered"
    "thetastart is an array of four initial isentropic levels considered"
    "isentropic level = 280 + thetastart*5"
    "filtersize is size for median filter"
    "resolution is minimum L2 distance between two points on returned contour"
       

    datestr = ['0922_12', '0926_12', '0930_12', '1003_12', '', '0930_12']

    #hoursafterinit = [24, 54, 42, 36, 24, 24]# 42 for all 21, 18
    # start time in hours after first time in dataset
    #filtersize = [30, 30, 30, 30, 20, 35], the 20 & 35 to deal with
    #two separate outflows incredibly close together
    
    #cldt = iris.load('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    cldt = iris.load(new_filepath + '/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_b*_thsfcs_5K.nc', 'total_minus_adv_only_theta')
    cldt[-1] = iris.util.new_axis(cldt[-1], 'time')
    dtcube = cldt.concatenate_cube()
    # create cube of delta theta
    
        
    #clpv = iris.load('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/2016' + datestr[k] +
    clpv = iris.load(new_filepath + '/nawdex/mi-ar482/2016' + datestr[k] +
    '/prodm_op_gl-mn_2016' + datestr[k] + '_c*_thsfcs_5K.nc', 'ertel_potential_vorticity')
    clpv[-1] = iris.util.new_axis(clpv[-1], 'time')
    pvcube = clpv.concatenate_cube()
    # and for PV
    
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
        
        cutoff = 150
        endoff = -1
        
        filtersize_t, filtersize_p = filtersize, filtersize
        
#        if k == 3:
#            cutoff = 50
#            endoff = -50
#        
#        if k == 4: #3rd outflow into IOP6
#            cutoff = 200
#            filtersize_p = 5 #Pv division
#            filtersize_t = 20 #narrow feature anyway
#        
        if k in [3, 5]: #outflow into IOP7, bit of a weird one
            cutoff = 0
            endoff = -150
            filtersize_p = 5 #to maintain very fine PV structure dividing the ridges
            
        dtcub = dtcube[hoursafterinit[k]/3, thetastart[k]+i, :, cutoff:endoff].copy()
        #restrict to outflow time and selected theta level
        #the cutoff artificially truncates the domain s.t. only the region we want is found
        #removes West: can't see new developing storms
        
        pvcub = pvcube[hoursafterinit[k]/3, thetastart[k]+i, :, cutoff:endoff].copy()
        #restrict to outflow time and selected theta level
        #the 150 artificially truncates the domain s.t. only the region we want is found
    
        dtcub.data = filters.median_filter(dtcub.data, size = filtersize_t)
        #apply median filter to smooth field
    
        pvcub.data = filters.median_filter(pvcub.data, size = filtersize_p)
        #apply median filter to smooth field

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
        
        if confail == False:    
    
            points = increase_nodes(lcontour, resolution)
            # increase number of points on the contour such that they have minimum spacing resolution  
    
            plt.plot(points[:, 0]-360, points[:, 1], marker = '.', ms = 0.8, mec = 'k', mfc = 'k')
            
            points3d = np.zeros([points.shape[0], 3])
            
            points3d[:, 0:2] = points
            
            points3d[:, 2] = potemp.points[thetastart[k]+i]
            
            if i == 0:
            
                outflow_volume = points3d
            
            else:
    
                outflow_volume = np.concatenate((outflow_volume, points3d))
            
        else:
            
            plt.title('Unable to find contour around outflow at this level')
            
    #np.save('outflow/T42_mfs' + str(filtersize) + '/' + timedatenow.strftime("%m%d_%H") + '.npy', outflow_volume) 
    
    #plt.savefig('outflow/T42_mfs' + str(filtersize) + '/' + timedatenow.strftime("%m%d_%H") + '_image.jpg')
    
#    np.save('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + timedatenow.strftime("%m%d_%H") + strapp + '.npy', outflow_volume) 
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/' + timedatenow.strftime("%m%d_%H") + strapp + '.npy', outflow_volume)
#    plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + timedatenow.strftime("%m%d_%H") + strapp + '_image.jpg')
    plt.savefig('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/' + timedatenow.strftime("%m%d_%H") + strapp + '_image.jpg')
    plt.show()
    
    return
    
    
def outflow_from_old(k = 1, lag = 1, folderin = 'IOP5/T42', folderout = 'IOP5/T42-6', strapp = ''):
    "define outflow volume from a trajectory ensemble already run"
    "To test consistency"
    "Taken 6*lag hours before the definition of the given trajectory"
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    TrB = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folderin + '/{}_TrajectoryEnsemble_backward'.format(basetime_str))
     
    outflow_volume = TrB.data[:, lag, :3]
    
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folderout + '/' + TrB.times[lag].strftime("%m%d_%H") + strapp + '.npy', outflow_volume) 
    
    return
    
    
def outflow_on_altitude(k = 0, folder = 'IOP3/T42', strapp = ''):
    "define outflow volume from a trajectory ensemble already run"
    "To get set of points on altitude coord to then run backwards 3d"
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    Tr2 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_2DTrajectoryEnsemble'.format(basetime_str) + strapp)
    
    llalt = np.array([True, True, False, True, False, False, False, False])
    outflow_volume = Tr2.data[:, 0, llalt]
    
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/' + basetime_str + strapp + 'initial_grid_new.npy', outflow_volume) 
    
    return
    
def try_caltra(k = 0, levs = 3, hoursafterinit = [42, 36, 42, 24], filtersize = 30, folder = 'IOP3/T42', strapp = ''):
    
    # Create a mapping from datetime objects to filenames
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    #datadir = '/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    datadir = '/storage/silver/NCAS-Weather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    timestep = datetime.timedelta(hours=6)
    
    #hoursafterinit = [42, 42, 42, 72, 96]#42
    
    times = [basetime[k] + timestep * i for i in range(hoursafterinit[k]/6 + 1)]
    print times[-1]
    #creates mapping up to and including selected outflow time    
        
    mapping = {}
    for t in times:
        leadtimehr = np.int((t - basetime[k]).total_seconds()) / 3600
        fn = [datadir + 'prodm_op_gl-mn_{0}_d{1:03d}_thgrid.pp'.\
            format(basetime_str, 12 * (leadtimehr / 12)),
              datadir + 'prodm_op_gl-mn_{0}_c{1:03d}.pp'.\
            format(basetime_str, 6 * (leadtimehr / 6))] # new addition of c to mapping
        mapping[t] = fn#datadir + fn
        
    #trainp_th = np.load('outflow/T42_mfs' + str(filtersize) + '/' + times[-1].strftime("%m%d_%H") + '.npy')
    #trainp_th = np.load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + times[-1].strftime("%m%d_%H") + strapp + '.npy')
    trainp_th = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/' + basetime_str + strapp + 'initial_grid_thlevs.npy')
    #trajectory input on theta levels
    
    levels = ('air_potential_temperature', [trainp_th[0, 2] + 5*i for i in range(levs)])
    # three isentropic levels including that of the first data point
    #levels = ('air_potential_temperature', [325])
    
    print levels
    
    tracers = ['altitude', 'x_wind', 'y_wind', 'upward_air_velocity', 'air_pressure', 'derived_pv', 'dPV_tot', 'adv_only_PV', 'dPV_LW', 'dPV_mic', 'dPV_conv', 'dPV_adv', 'dPV_SW', 'dPV_ph1', 'dPV_bl', 'dPV_cld', 'dPV_mass']
    #need these for circulation integrals
    
    traout = caltra.caltra(trainp_th, mapping, fbflag=-1, nsubs = 12, tracers = tracers, levels=levels)
    # 12 steps between files = 30 mins apart 
    
    #np.save('outflow/{}_trajectories.npy'.format(basetime_str), traout)
    #saving it as a numpy file removes functionality oof trajectory ensemble class, making it effectively useless, yay!
        
    #traout.save('outflow/T42_mfs' + str(filtersize) + '/{}_TrajectoryEnsemble_backward'.format(basetime_str))
    #traout.save('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/{}_TrajectoryEnsemble_backward'.format(basetime_str) + strapp)
    #traout.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/{}_TrajectoryEnsemble_backward'.format(basetime_str) + strapp)
    traout.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_2DTrajectoryEnsemble_dPV'.format(basetime_str) + strapp)
     
    return traout
    

def try_caltra_rev(k = 0, levs = 3, hoursafterinit = [42, 36, 42, 24], hoursuntil = [90, 78, 156, 84], filtersize = 30, folder = 'IOP3/contours_2019', strapp = '', theta_max = 315):
    
    # Create a mapping from datetime objects to filenames
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    datadir = '/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    timestep = datetime.timedelta(hours=6)
    
    #hoursafterinit = [42, 42, 42, 72, 96]#42
    #hoursuntil = [78, 78, 78, 108]
    
    times = [basetime[k] + timestep * i for i in range(hoursafterinit[k]/6, hoursuntil[k]/6 + 1)]
    print times[0]
    #creates mapping from selected outflow time until selected end time 
    print 7
    mapping = {}
    for t in times:
        leadtimehr = np.int((t - basetime[k]).total_seconds()) / 3600
        fn = 'prodm_op_gl-mn_{0}_d{1:03d}_thgrid.pp'.\
            format(basetime_str, 12 * (leadtimehr / 12))
        mapping[t] = datadir + fn
        
    #trainp_th = np.load('outflow/T42_mfs' + str(filtersize) + '/' + times[0].strftime("%m%d_%H") + '.npy')
    #trainp_th = np.load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + times[0].strftime("%m%d_%H") + strapp + '.npy')
    trainp_th = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/' + times[0].strftime("%m%d_%H") + strapp + '.npy')
    #trainp_th = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + 2)
    #trajectory input on theta levels
    
    if len(trainp_th[0]) >= 3:
        # the input file has already been given theta levels
        levels = ('air_potential_temperature', [trainp_th[0, 2] + 5*i for i in range(levs)])
        # three isentropic levels including that of the first data point
    else:
        levels = ('air_potential_temperature', [theta_max - 5*i for i in range(levs)])
        trainp_th_new = np.concatenate([trainp_th, np.transpose([theta_max*np.ones_like(trainp_th[:, 0])])], axis = 1)
        for i in range(1, levs):
            trainp_th_new = np.concatenate([trainp_th_new, np.concatenate([trainp_th, np.transpose([(theta_max - 5*i)*np.ones_like(trainp_th[:, 0])])], axis = 1)], axis = 0)
        trainp_th = trainp_th_new
    
    print levels
    
    tracers = ['altitude', 'x_wind', 'y_wind', 'upward_air_velocity', 'air_pressure']
    #need these for circulation integrals
    
    traout = caltra.caltra(trainp_th, mapping, fbflag=1, nsubs = 12, tracers = tracers, levels=levels)
    # 12 steps between files = 30 mins apart 
    
    #np.save('outflow/{}_trajectories.npy'.format(basetime_str), traout)
    #saving it as a numpy file removes functionality oof trajectory ensemble class, making it effectively useless, yay!
        
    #traout.save('outflow/T42_mfs' + str(filtersize) + '/{}_TrajectoryEnsemble_forward'.format(basetime_str))
    #traout.save('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/{}_TrajectoryEnsemble_forward'.format(basetime_str) + strapp)
    traout.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/{}_TrajectoryEnsemble_forward'.format(basetime_str) + strapp)
    
    return traout
    
    
def plot_traj(k = 0, levs = 3, filtersize = 30, field = 'total_minus_adv_only_theta', scatter3D = False, lim = 40, indiv = False, folder = 'IOP3/T42', strapp = ''):
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
    
    for tht in xrange(levs):
        
        theta = theta_0 + 5*tht
        
        TrEntheta = TrEn.select('air_potential_temperature', '==', theta, time = [-1*TrB.relative_times[-1]])
        # This gives the levels at the designated start time
        
        if scatter3D != False:
            
            Tt = Tr3.select('air_potential_temperature', '>=', theta - 2.5, time = [datetime.timedelta(hours = 0)])
            Tt = Tt.select('air_potential_temperature', '<', theta + 2.5, time = [datetime.timedelta(hours = 0)])
            Tt = Tt.select(field, '!=', -1000)
            # points starting from theta level for which field is defined
                    
        for t in xrange(tsteps):
            
            if indiv == False:
                plt.subplot(tsteps, levs, t*levs + tht + 1)
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
        pass
#        plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + basetime_str + 'fulltraj_' + field + strapp + '.png')
        #plt.savefig('easy_to_find/' + folder + strapp + '_3Dcloud.png')
    