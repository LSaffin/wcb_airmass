# -*- coding: utf-8 -*-
"""
Awful awful scripts to sum timeseries across multiple outflows
"""

import iris
import numpy as np

save_dir = '/glusterfs/msc/users_2018/bn826011/NAWDEX/Final/'

def Combine_IOP3_ABsplit(strapp = '_split'):
    "This may seem very stupid but I can think of no better way"

    for theta_level in [310, 315, 320, 325, 330]:

        cubes0 = iris.load(save_dir + 'IOP3/ABsplit' + '/circulation/circulations_' + str(theta_level) + 'K_' + '0' + strapp + '.nc')

        cubes1 = iris.load(save_dir + 'IOP3/ABsplit' + '/circulation/circulations_' + str(theta_level) + 'K_' + '1' + strapp + '.nc')

        results = iris.cube.CubeList([])

        for i in xrange(len(cubes0)):
            if cubes0[i].name() == cubes1[i].name():
                cubesall = iris.cube.CubeList([])
                for n, time0 in enumerate(cubes0[i].coord('time')):
                    time0 = time0.points[0]
                    for m, time1 in enumerate(cubes1[i].coord('time')):
                        time1 = time1.points[0]
                        if time0 == time1:
                            cube = cubes0[i][n] + cubes1[i][m]
                            #print cube
                coordlist = []
                for coord in cube.coords():
                    coordlist.append(coord.standard_name)
                if 'forecast_period' in coordlist:
                    cube.remove_coord('forecast_period')
                if 'time' in coordlist:
                    cube.rename(cubes0[i].name())
                    cubesall.append(iris.util.new_axis(cube, 'time'))
                        
                if cubesall == iris.cube.CubeList([]):
                    pass
                else:
                    cube = cubesall.concatenate_cube()
                    results.append(cube)
            else:
                print 'they have different names and I have no idea why'
        
        iris.save(results.merge(),
               save_dir + 'IOP3/ABsplit' + '/circulation/circulations_' + str(theta_level) + 'K_' + strapp + '.nc')
               
               
               
def Combine_all_IOP3(k = 0, thlevs = [[310, 315, 320, 325, 330], [], [310, 315, 320]], folders = ['IOP3/ABsplit', '', 'IOP6'], hoursuntil = [90, 78, 156, 84], strapp = '_split'):
    "This may seem very stupid but I can think of no better way"

    for theta_level in thlevs[k]:

        cubes0 = iris.load(save_dir + folders[k] + '/circulation/circulations_' + str(theta_level) + 'K_' + '0' + strapp + '.nc')

        cubes1 = iris.load(save_dir + folders[k] + '/circulation/circulations_' + str(theta_level) + 'K_' + '1' + strapp + '.nc')

        results = iris.cube.CubeList([])

        for i in xrange(len(cubes0)):
            if cubes0[i].name() == cubes1[i].name():
                cubesall = iris.cube.CubeList([])
                cubesn = iris.cube.CubeList([cubes1[i]])
                
                timezero = cubes0[i].coord('time')[0]
                
                for time_iter in xrange(hoursuntil[0]/6 - 1):
                    
                    timenow = timezero + 6*time_iter
                    print timenow
                    
                    for n, time in enumerate(cubes0[i].coord('time')):
                        time = time.points[0]
                        if time == timenow.points[0]:
                            cubenow = cubes0[i][n]
                
                            for cube in cubesn:
                
                                for n, time in enumerate(cube.coord('time')):
                                    time = time.points[0]
                                    if time == timenow.points[0]:
                                        cubenow += cube[n]
                                
                    #print cube
                            coordlist = []
                            for coord in cubenow.coords():
                                coordlist.append(coord.standard_name)
                            if 'forecast_period' in coordlist:
                                cubenow.remove_coord('forecast_period')
                            if 'time' in coordlist:
                                cubenow.rename(cubes0[i].name())
                                print cubenow
                                cubesall.append(iris.util.new_axis(cubenow, 'time'))
                        
                if cubesall == iris.cube.CubeList([]):
                    pass
                else:
                    print cubesall
                    cubecat = cubesall.concatenate_cube()
                    results.append(cubecat)
            else:
                print 'they have different names and I have no idea why'
        
        iris.save(results.merge(),
               save_dir + folders[k] + '/circulation/circulations_' + str(theta_level) + 'K_' + strapp + '_all.nc')
               
               
               
def Combine_IOP6(strapp = '_split'):
    "This may seem very stupid but I can think of no better way"

    for theta_level in [310, 315, 320]:

        cubes0 = iris.load(save_dir + 'IOP6' + '/circulation/circulations_' + str(theta_level) + 'K_' + '0' + strapp + '.nc')

        cubes1 = iris.load(save_dir + 'IOP6' + '/circulation/circulations_' + str(theta_level) + 'K_' + '1' + strapp + '.nc')
        
        cubes2 = iris.load(save_dir + 'IOP6' + '/circulation/circulations_' + str(theta_level) + 'K_' + '2' + strapp + '.nc')

        cubes3 = iris.load(save_dir + 'IOP6' + '/circulation/circulations_' + str(theta_level) + 'K_' + '3' + strapp + '.nc')

        results = iris.cube.CubeList([])

        for i in xrange(len(cubes0)):
            if cubes0[i].name() == cubes1[i].name() == cubes2[i].name() == cubes3[i].name():
                cubesall = iris.cube.CubeList([])
                for n, time0 in enumerate(cubes0[i].coord('time')):
                    time0 = time0.points[0]
                    for m, time1 in enumerate(cubes1[i].coord('time')):
                        time1 = time1.points[0]
                        if time0 == time1:
                            for l, time2 in enumerate(cubes2[i].coord('time')):
                                time2 = time2.points[0]
                                if time0 == time2:
                                    for k, time3 in enumerate(cubes3[i].coord('time')):
                                        time3 = time3.points[0]
                                        if time0 == time3:
                                            cube = cubes0[i][n] + cubes1[i][m] + cubes2[i][l] + cubes3[i][k]
                                            
                coordlist = []
                for coord in cube.coords():
                    coordlist.append(coord.standard_name)
                if 'forecast_period' in coordlist:
                    cube.remove_coord('forecast_period')
                if 'time' in coordlist:
                    cube.rename(cubes0[i].name())
                    cubesall.append(iris.util.new_axis(cube, 'time'))
                if cubesall == iris.cube.CubeList([]):
                    pass
                else:
                    cube = cubesall.concatenate_cube()
                    results.append(cube)
            else:
                print 'they have different names and I have no idea why'
        
        iris.save(results.merge(),
               save_dir + 'IOP6' + '/circulation/circulations_' + str(theta_level) + 'K_' + strapp + '.nc')
               
               
def Combine_all_IOP6(k = 2, thlevs = [[310, 315, 320, 325, 330], [], [310, 315, 320]], folders = ['IOP3/ABsplit', '', 'IOP6'], hoursuntil = [90, 78, 156, 84], strapp = '_split'):
    "This may seem very stupid but I can think of no better way"

    for theta_level in thlevs[k]:

        cubes0 = iris.load(save_dir + folders[k] + '/circulation/circulations_' + str(theta_level) + 'K_' + '0' + strapp + '.nc')

        cubes1 = iris.load(save_dir + folders[k] + '/circulation/circulations_' + str(theta_level) + 'K_' + '1' + strapp + '.nc')
            
        cubes2 = iris.load(save_dir + folders[k] + '/circulation/circulations_' + str(theta_level) + 'K_' + '2' + strapp + '.nc')

        cubes3 = iris.load(save_dir + folders[k] + '/circulation/circulations_' + str(theta_level) + 'K_' + '3' + strapp + '.nc')

        results = iris.cube.CubeList([])

        for i in xrange(len(cubes0)):
            if cubes0[i].name() == cubes1[i].name() == cubes2[i].name() == cubes3[i].name():
                cubesall = iris.cube.CubeList([])
                cubesn = iris.cube.CubeList([cubes1[i], cubes2[i], cubes3[i]])
                
                timezero = cubes0[i].coord('time')[0]
                
                for time_iter in xrange(hoursuntil[0]/6):
                    
                    timenow = timezero + 6*time_iter
                    
                    for n, time in enumerate(cubes0[i].coord('time')):
                        time = time.points[0]
                        if time == timenow.points[0]:
                            cubenow = cubes0[i][n]
                
                            for cube in cubesn:
                
                                for n, time in enumerate(cube.coord('time')):
                                    time = time.points[0]
                                    if time == timenow.points[0]:
                                        cubenow += cube[n]                
                                
                            #print cube
                            coordlist = []
                            for coord in cubenow.coords():
                                coordlist.append(coord.standard_name)
                            if 'forecast_period' in coordlist:
                                cubenow.remove_coord('forecast_period')
                            if 'time' in coordlist:
                                cubenow.rename(cubes0[i].name())
                                cubesall.append(iris.util.new_axis(cubenow, 'time'))
                        
                if cubesall == iris.cube.CubeList([]):
                    pass
                else:
                    cube = cubesall.concatenate_cube()
                    results.append(cube)
            else:
                print 'they have different names and I have no idea why'
        
        iris.save(results.merge(),
               save_dir + 'IOP3/ABsplit' + '/circulation/circulations_' + str(theta_level) + 'K_' + strapp + '_all.nc')