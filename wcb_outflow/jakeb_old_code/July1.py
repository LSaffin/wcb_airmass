# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import iris
from matplotlib.path import Path
from lagranto.trajectory import load
import datetime
from mymodule import grid, convert, interpolate
import matplotlib.pyplot as plt
from lagranto import caltra, trajectory
import random

from matplotlib.dates import DayLocator, HourLocator, DateFormatter

iris.FUTURE.cell_datetime_objects=True

"I can't do 3D trajectories because I can't find the 'altitude' coordinate to start them off"
"As 3D trajectories require altitude instead of potential temperature"

ben_datadir = '/storage/silver/NCAS-Weather/ben/nawdex/mi-ar482/'




def outflow_grid(k = 0, levels = 3, hoursafterinit = [42, 36, 42, 24], thlevs = [[320, 325, 330], [325, 330, 335], [310, 315, 320], [310, 315, 320]], folder = 'IOP3/T42', strapp = ''):
    "This presently defines the grid of points from the sime time that the region is defined"
    "But by loading the forward trajectories this could easily be adapted to"
    "define the grid at any time along the trajectory"

    #save_dir = '/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/'
    save_dir = '/storage/silver/scenario/bn826011/WCB_outflow/Final/'
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    TrB = load(save_dir + folder + '/{}_TrajectoryEnsemble_backward'.format(basetime_str) + strapp)
    
    start_time = [datetime.timedelta(hours = 0)]
    # trajectory at definition time

    datadir = '/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    
    #hoursafterinit = [42, 42, 42, 72, 96]#42
    
    time = basetime[k] + datetime.timedelta(hours = hoursafterinit[k])
    
    mapping = {}
    leadtimehr = np.int((time - basetime[k]).total_seconds()) / 3600
    fn = 'prodm_op_gl-mn_{0}_d{1:03d}_thgrid.pp'.\
        format(basetime_str, 12 * (leadtimehr / 12))
    mapping[time] = datadir + fn
    
#    dexs = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/20160922_12/' + 
#            'prodm_op_gl-mn_20160922_12_d*_thsfcs_5K.nc', 'dimensionless_exner_function')
#    dexs[-1] = iris.util.new_axis(dexs[-1], 'time')
#    dex = dexs.concatenate_cube()
#    
#    temps = iris.load('/export/cloud/NCASweather/ben/nawdex/mi-ar482/20160922_12/' + 
#            'prodm_op_gl-mn_20160922_12_d*_thsfcs_5K.nc', 'air_temperature')
#    temps[-1] = iris.util.new_axis(temps[-1], 'time')
#    temp = temps.concatenate_cube()
#    
#    lvls = ('air_potential_temperature', thlevs[k])
#    
#    altd = convert.calc('altitude', iris.cube.CubeList([dex, temp]), levels = lvls)
#    # 3D caltra works on altitude not theta levels
#    # however I don't know how to get altitude?!
#    # I think leo's data must've have altitude as a coordinate
    
    cubes = iris.load(mapping[time], iris.Constraint(time=time))
    
    plt.figure(figsize = (10, 14))
    
    zc = iris.load('/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/' + basetime_str +
    '/prodm_op_gl-mn_' + basetime_str + '_d*_thsfcs_5K.nc', 'altitude')
    zc = zc[:-1].concatenate_cube()
    # the last time step has different metadata?
    
    for kl in xrange(levels):
        
        theta_level = thlevs[k][kl]
    
        trajectories = TrB.select(
            'air_potential_temperature', '==', theta_level, time = start_time)
    
        x = trajectories.x[:, 0]
        y = trajectories.y[:, 0]
        #
        
#        tlev_cstrnt = iris.Constraint(air_potential_temperature = theta_level)
#        
#        altdth = altd.extract(tlev_cstrnt)
        
        lvls = ('air_potential_temperature', [theta_level])
        
        w = convert.calc('upward_air_velocity', cubes, levels = lvls)
        # I think fairly arbitrary cube
        
#        z = grid.make_cube(w, 'altitude')
#        # altitude
#        
#        print z
#        
#        lvls = ('air_potential_temperature', [theta_level])
#        
#        coord_name, values = lvls

        # I now need some way to interpolate altitude to desired theta level
#        z = interpolate.to_level(z, **{coord_name: values})
        
#        z = convert.calc('altitude', iris.cube.CubeList([z]), levels = lvls)
    
        glon, glat = grid.get_xy_grids(w)
        gridpoints = np.array([glon.flatten(), glat.flatten()]).transpose()
    
        points = np.array([x, y]).transpose()
        pth = Path(points)

        # Mask all points that are not contained in the circuit
        mask = np.logical_not(pth.contains_points(gridpoints).reshape(glat.shape))
        
        tlev_cstrnt = iris.Constraint(air_potential_temperature = theta_level)
        time_cstrnt = iris.Constraint(time = time)        
        
        #try this for altitude        
        z = zc.extract(tlev_cstrnt & time_cstrnt)        
        
        plt.subplot(levels, 2, 2*kl+1)
        plt.contourf(mask, cmap = 'gray')
    
        masked_lon = []
        masked_lat = []
        alt_list = []
        for i, col in enumerate(mask):
            for j, point in enumerate(col):
                if point == False:
                    lat = glat[i, j]
                    lon = glon[i, j]
                    alt = z.data[i, j]
                    masked_lon.append(lon)
                    masked_lat.append(lat)
                    alt_list.append(alt)
                    
                    
        plt.subplot(levels, 2, 2*kl+2)
        plt.scatter(masked_lon, masked_lat, s = 2, c = 'k', marker = '.', edgecolor = 'k')
                
        lt = len(masked_lon)

        points3d = np.zeros([lt, 3])

        points3d[:, 0] = np.array(masked_lon)

        points3d[:, 1] = np.array(masked_lat)

        points3d[:, 2] = np.array(alt_list)
        
        pointsth = np.zeros([lt, 3])
        
        pointsth[:, 0] = np.array(masked_lon)

        pointsth[:, 1] = np.array(masked_lat)

        pointsth[:, 2] = theta_level*np.ones([lt])
        
        if kl == 0:
            
            outflow_volume = points3d
            outflow_volume_th = pointsth
            
        else:
    
            outflow_volume = np.concatenate((outflow_volume, points3d))
            outflow_volume_th = np.concatenate((outflow_volume_th, pointsth))
        
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/' + basetime_str + strapp + 'initial_grid.npy', outflow_volume) 
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/' + basetime_str + strapp + 'initial_grid_thlevs.npy', outflow_volume_th) 
    
    #plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + basetime_str + strapp + '_masks.jpg')
    
    plt.show()
        
    
def inflow_caltra(k = 0, levs = 3, hoursafterinit = [42, 36, 42, 24], folder = 'IOP3/T42', strapp = ''):
    
    # Create a mapping from datetime objects to filenames
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    #datadir = '/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    datadir = '/storage/silver/NCAS-Weather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    timestep = datetime.timedelta(hours=6)
    
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
        
    #trainp_th = np.load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/' + basetime_str + strapp + 'initial_grid.npy')
    trainp_th = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/' + basetime_str + strapp + 'initial_grid.npy')
    #trajectory input on theta levels
    
    tracers = ['air_potential_temperature', 'air_pressure', 'derived_pv', 'dPV_tot', 'adv_only_PV', 'dPV_LW', 'dPV_mic', 'dPV_conv', 'dPV_adv', 'dPV_SW', 'dPV_ph1', 'dPV_bl', 'dPV_cld', 'dPV_mass']# 'specific_humidity', 'mass_fraction_of_cloud_ice_in_air', 'mass_fraction_of_cloud_liquid_water_in_air', 'derived_pv']
    #need these for circulation integrals
    
    traout = caltra.caltra(trainp_th, mapping, fbflag=-1, nsubs = 12, tracers = tracers)
    # 12 steps between files = 30 mins apart 
    
    traout.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_dPV'.format(basetime_str) + strapp)
    
    return traout
    
def inflow_caltra_fw(k = 0, levs = 3, hoursafterinit = [42, 42, 42, 24], hoursuntil = [90, 78, 156, 84], folder = 'IOP3/T42', strapp = ''):
    
    # Create a mapping from datetime objects to filenames
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    datadir = '/export/cloud/migrated-NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    timestep = datetime.timedelta(hours=6)
    
    times = [basetime[k] + timestep * i for i in range(hoursafterinit[k]/6, hoursuntil[k]/6 + 1)]
    print times[0]
    #creates mapping up to and including selected outflow time    
    
    mapping = {}
    for t in times:
        leadtimehr = np.int((t - basetime[k]).total_seconds()) / 3600
        fn = 'prodm_op_gl-mn_{0}_d{1:03d}_thgrid.pp'.\
            format(basetime_str, 12 * (leadtimehr / 12))
        mapping[t] = datadir + fn
        
    #trainp_th = np.load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/' + basetime_str + strapp + 'initial_grid.npy')
    trainp_th = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/' + basetime_str + strapp + 'initial_grid.npy')
    #trajectory input on theta levels
    
    tracers = ['air_potential_temperature', 'air_pressure']#, 'specific_humidity', 'mass_fraction_of_cloud_ice_in_air', 'mass_fraction_of_cloud_liquid_water_in_air']
    #need these for circulation integrals
    
    traout = caltra.caltra(trainp_th, mapping, fbflag=1, nsubs = 12, tracers = tracers)
    # 12 steps between files = 30 mins apart 
    
    #traout.save('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_fw'.format(basetime_str) + strapp)
    traout.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_fw'.format(basetime_str) + strapp)    
    
    return traout
    
    
def heating_hist(k = 0, levs = 3, theta_0 = [320, 325, 310, 310], t = 0, folder = 'IOP3/T42', strapp = ''):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    #Tr3 = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble'.format(basetime_str) + strapp)
    Tr3 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble'.format(basetime_str) + strapp) 
    
    bins = np.linspace(-10, 30, 21)
    
    plt.figure(figsize = (13, 9))
    
    Tall = Tr3.select('air_potential_temperature', '!=', -1000)
    #only those which remain in domain
    
    plt.subplot(int(levs/2)+1, 2, 1)
    curve = plt.hist(Tall.data[:, 0, 3] - Tall.data[:, -(t+1), 3], bins = bins)
    #plt.title(folder + '_' + strapp)
    plt.title('Started at all surfaces')
    plt.xlabel('Delta theta')
    
    for i in xrange(levs):
        
        theta = theta_0[k] + i*5
        
        Tt = Tall.select('air_potential_temperature', '>=', theta - 2.5, time = [datetime.timedelta(hours = 0)])
        Tt = Tt.select('air_potential_temperature', '<', theta + 2.5, time = [datetime.timedelta(hours = 0)])
        # those which remain in domain that start on desired theta surface
        
        plt.subplot(int(levs/2)+1, 2, 2+i)
        plt.hist(Tt.data[:, 0, 3] - Tt.data[:, -(t+1), 3], bins = bins)
        plt.title('Started at ' + str(theta) + 'K surface')
        plt.xlabel('Delta theta')
        
    plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/delta_theta_histograms_t=' + str(t) + strapp + '.png')
    #plt.savefig('IOP3_T42_2x2_hist_t=' + str(t) + '.png')
    plt.show()
    
    return curve
    
    
def random_3D_trajs(k = 0, levs = 3, theta_0 = [320, 325, 310, 310], trajs = 50, folder = 'IOP3/T42', strapp = '', var_idxs = [8, 9, 10]):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    #Tr3 = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble'.format(basetime_str) + strapp)
    #Tr3 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_dPV'.format(basetime_str) + strapp)   
    Tr2 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_dPV'.format(basetime_str) + strapp)

    plt.figure(figsize = (12, 8))
    
    for j, threshs in enumerate([[10, 60], [5, 10], [0, 5], [-60, 0]]):
        
        thresh_lo = threshs[0]
        thresh_hi = threshs[1]
    
        Tall = Tr2.select('air_potential_temperature', '!=', -1000)
        #only those which remain in domain
        Tasc = Tall.select('air_potential_temperature', '>', thresh_lo, time=[Tall.relative_times[0], Tall.relative_times[-1]])
        Tasc = Tall.select('air_potential_temperature', '<', thresh_hi, time=[Tall.relative_times[0], Tall.relative_times[-1]])

        # select all trajectories which ascend at least 10K
    
        #plt.subplot(int(levs/2)+1, 2, 1)
        plt.subplot(2, 2, 1+j)
    
#        for i in xrange(trajs):
#            i = random.randint(0, len(Tall)-1)
#            plt.plot(Tall.times, Tall.data[i, :, 5])#3])
    
        for var_idx in var_idxs:
    
            Tmean = np.mean(Tasc.data[:, :, var_idx], axis = 0)
            #Tstdev = np.std(Tasc.data[:, :, var_idx], axis=0)
            
            plt.plot(Tasc.times, Tmean, linewidth = 3, label = Tasc.names[var_idx])
            #plt.plot(Tasc.times, Tmean + Tstdev, 'r', linewidth = 1)
            #plt.plot(Tasc.times, Tmean - Tstdev, 'r', linewidth = 1)
        
        #plt.title(folder + '_' + strapp)
        plt.title('Started at all surfaces, ' + str(thresh_lo) + ' < dtheta < ' + str(thresh_hi))
        plt.ylabel('PVU')#theta')
        plt.legend(loc = 'upper left')
        plt.ylim(-.25, .25)
    
        pg = plt.gca()
        fmt = DateFormatter('\n%m/%d')                  
        fmt2 = DateFormatter('%H')
        majorLocator = DayLocator(interval=1)
        minorLocator = HourLocator(range(0, 24, 6))
        pg.xaxis.set_major_formatter(fmt)
        pg.xaxis.set_minor_formatter(fmt2)
        pg.xaxis.set_minor_locator(minorLocator)
        pg.xaxis.set_major_locator(majorLocator)
    
#    for j in xrange(levs):
#        
#        theta = theta_0[k] + j*5
#        
#        Tt = Tall.select('air_potential_temperature', '>=', theta - 2.5, time = [datetime.timedelta(hours = 0)])
#        Tt = Tt.select('air_potential_temperature', '<', theta + 2.5, time = [datetime.timedelta(hours = 0)])
#        # those which remain in domain that start on desired theta surface
##        Tasc = Tt.select('air_potential_temperature', '>', thresh_lo, time=[Tt.relative_times[0], Tt.relative_times[-1]])
##        Tasc = Tt.select('air_potential_temperature', '<', thresh_hi, time=[Tt.relative_times[0], Tt.relative_times[-1]])
#        # select all trajectories which ascend at least 10K
#        Tasc = Tt
#        
#        plt.subplot(int(levs/2)+1, 2, 2+j)
#        
##        for i in xrange(trajs):
##            i = random.randint(0, len(Tt))
##            plt.plot(Tasc.times, Tt.data[i, :, 5])
#        
#        for var_idx in var_idxs:        
#        
#            Tmean = np.mean(Tasc.data[:, :, var_idx], axis = 0)
#            #Tstdev = np.std(Tasc.data[:, :, var_idx], axis=0)
#            
#            plt.plot(Tasc.times, Tmean, linewidth = 3, label = Tasc.names[var_idx])
#            #plt.plot(Tasc.times, Tmean + Tstdev, 'r', linewidth = 1)
#            #plt.plot(Tasc.times, Tmean - Tstdev, 'r', linewidth = 1)
#        
#        plt.title('Started at ' + str(theta) + 'K surface')
#        plt.ylabel('PVU')#theta')
#        plt.legend()
#        plt.ylim(-.25, .25)
#    
#        pg = plt.gca()
#        fmt = DateFormatter('\n%m/%d')                  
#        fmt2 = DateFormatter('%H')
#        majorLocator = DayLocator(interval=1)
#        minorLocator = HourLocator(range(0, 24, 6))
#        pg.xaxis.set_major_formatter(fmt)
#        pg.xaxis.set_minor_formatter(fmt2)
#        pg.xaxis.set_minor_locator(minorLocator)
#        pg.xaxis.set_major_locator(majorLocator)
        
    #plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/' + str(trajs) + '_random_3D-trajectories' + strapp + '.png')
    #plt.savefig('IOP3_T42_random_2x2_50_back.png')  
    plt.savefig('/home/users/bn826011/NAWDEX/From N Drive/2019_figs/dPV/' + folder[:4]
    + '_3D_PV_trajectories_all.png')# + str(var_idxs) + strapp + '.png')

    plt.show()
        
def inflow_time(k = 0, levs = 3, theta_0 = [320, 325, 310, 310], hoursafterinit = [42, 42, 42, 24], folder = 'IOP3/T42', strapp = ''):
    
    xax = []
    curve = []    
    
    for t in xrange(hoursafterinit[k]/6):
        
        hists = heating_hist(k = k, levs = levs, theta_0 = theta_0, t=t, folder = folder, strapp = strapp)
        
        xax.append((hists[1][:-1]+hists[1][1:])/2)
        curve.append(hists[0])
    
    for t in xrange(hoursafterinit[k]/6):
        
        plt.plot(xax[t], curve[t], label = str(t*6), linewidth = 2)
        
    #plt.legend()
    plt.legend(title = 'Hours since start of forecast')
    plt.xlabel('Delta theta')
    plt.title('Difference in theta between time and outflow time')
    plt.ylabel('Number of trajectories')
    
    plt.savefig('IOP3_T42_inflow_time.png')
    #plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/inflow/theta_level_distribution_all_times' + strapp + '.png')
    plt.show()
        