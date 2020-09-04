# -*- coding: utf-8 -*-
"""
Isentropic density etc.

Design in same way as circ_int_test
"""

import numpy as np
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import datetime
import mymodule.calculus as mycalc
from mymodule import convert, grid
from lagranto.trajectory import load
from mymodule import forecast, interpolate

def calc_density(k = 0, thlevs = [[320, 325, 330], [325, 330, 335], [310, 315, 320], [310, 315, 320]], plotv = True, folder = 'IOP3/T42', strapp = ''):    
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]  
    
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    datadir = '/export/cloud/NCASweather/ben/nawdex/mi-ar482/{}/'.format(basetime_str)
    
    TrB = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/{}_TrajectoryEnsemble_backward'.format(basetime_str) + strapp)
    TrF = load('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/{}_TrajectoryEnsemble_forward'.format(basetime_str) + strapp)
    
    start_time = [-1*TrB.relative_times[-1]]
    # time at which trajectories were initialised
    
    trajectories = TrB.__add__(TrF)
    
    mapping = {}
    for t in trajectories.times:
        leadtimehr = np.int((t - basetime[k]).total_seconds()) / 3600
        fn = 'prodm_op_gl-mn_{0}_d{1:03d}_thgrid.pp'.\
            format(basetime_str, 12 * (leadtimehr / 12))
        mapping[t] = datadir + fn
        
    fcast = forecast.Forecast(basetime[k], mapping)
    
    idcl = iris.cube.CubeList([])
    # create empty cube list of isentropic density
    
    pvsl = iris.cube.CubeList([])
    # and empty cube for pv substance
    
    if plotv == True:
        
        plt.figure(figsize = (15, 40))
    
    for m, theta_level in enumerate(thlevs[k]):
    
        # Select an individual theta level
        TrEntheta = trajectories.select(
        'air_potential_temperature', '==', theta_level, time = start_time)
        
        levels = ('air_potential_temperature', [theta_level])
        
        icl = iris.cube.CubeList([])
        pcl = iris.cube.CubeList([])
    
        for n, time in enumerate(TrEntheta.times):
        
            cubelist = fcast.set_time(time)
        
            #now for some reason I'm having real problems extracting the right time, so
            u = iris.unit.Unit('hours since 1970-01-01 00:00:00', calendar=iris.unit.CALENDAR_STANDARD)
            timeh = u.date2num(time)
            #this puts the time in the same format as in the cubelist
            cubes = cubelist.extract(iris.Constraint(time = time))
            # needed because each file contains two times
        
            pv = convert.calc('derived_pv', cubes, levels=levels)
        
            mass = convert.calc('mass', cubes, levels = levels)
    
            theta = convert.calc('air_potential_temperature', cubes)
    
            density = convert.calc('air_density', cubes, levels = levels)
            
            #print density
        
            coord = grid.make_coord(theta)        
            theta.add_aux_coord(coord, range(theta.ndim))
            # add theta as coordinate
        
            dtdz = mycalc.multidim(theta, 'altitude', 'z')
            
            #print dtdz
            
            coord_name, values = levels
        
            dzdtcube = interpolate.to_level(dtdz, **{coord_name: values})
        
            dzdtcube.convert_units('meter-kelvin^-1')
            
            #print dzdtcube
        
            r = dzdtcube * density
        
            r.rename('isentropic_density')
            
            #print r
            
            r.remove_coord('forecast_period')
            
            #print r
        
            idcl.append(iris.util.new_axis(r, 'time'))
            # add isentropic density to cube list
        
            pvs = pv * mass
            
            pvs.rename('pv_substance')
            
            pvs.remove_coord('forecast_period')
            
            pvsl.append(iris.util.new_axis(pvs, 'time'))
            
            #print pvsl[-1]
            
            if plotv == True:
                
                tlno = len(thlevs[k])
                
                tno = len(TrEntheta.times)
                
                plt.subplot(tno, tlno, n*tlno + m + 1)
                
                #iplt.contourf(r[0], np.linspace(0, 1000, 21), cmap = 'binary')
                iplt.contourf(pvs[0], cmap = 'binary')
                
                plt.gca().coastlines(color = [.6, .4, 0])
                
                iplt.contour(pv[0], [2], colors = ['r'])
                
                TrEnindomain = trajectories.select('air_potential_temperature', '==', theta_level, time = [datetime.timedelta(hours=6*n)])
            
                Leftdomain = TrEntheta.__len__() - TrEnindomain.__len__()
            
                if np.shape(TrEntheta) == (0,):
                
                    pass
            
                elif Leftdomain != 0:
                    #If any of the trajectories have left the domain
            
                    plt.plot(TrEnindomain.data[:, n, 0]-360, TrEnindomain.data[:, n, 1], color = [.4, .4, .4], marker = '.', ms = .8, mec = 'k', mfc = 'k')
                
                    plt.title(str(Leftdomain) + ' points have left the domain')
            
                else:
            
                    plt.plot(TrEntheta.data[:, n, 0]-360, TrEntheta.data[:, n, 1], color = [.4, .4, .4], marker = '.', ms = .8, mec = 'k', mfc = 'k')
                    # plot outflow area, black
            
                if n == 0:
                
                    if Leftdomain != 0:
                
                        plt.title(str(theta_level) + ' K \n ' + str(Leftdomain) + ' points have left the domain')
                        #label columns with isentropic levels
                    
                    else:
                    
                        plt.title(str(theta_level) + ' K')
                
                
                if m == 0:
                
                    timedatenow = trajectories.times[n]
                
                    plt.text(-90, 70, timedatenow.strftime("%d/%m %HUTC"), rotation = 'vertical')
                
                    #label rows with dates and times
                
    if plotv == True:
                            
        plt.savefig('/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' + basetime_str + '_' + 'pv_substance' + strapp + '.png')
        plt.show()
    
    #print pvsl
    
    pv_substance = pvsl.concatenate_cube()
    
    #print pv_substance
    
    isentropic_density = idcl.concatenate_cube()
    
    both = iris.cube.CubeList([isentropic_density, pv_substance])
    
    iris.save(both, '/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/' + folder + '/' +  basetime_str + '_isentropic_density' + strapp + '.nc')
    
    return isentropic_density, pv_substance
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    