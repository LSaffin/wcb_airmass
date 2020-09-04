"""
Everything script

If you run this script with access to the appropriate modules and forecast data
And directory structure
This will do everything, from the start
"""

from outflow_area import outflow_th, try_caltra, try_caltra_rev
from outflow_maskprev import outflow_AB
from outflow_area import outflow_from_old, plot_traj

from ridge_caltra import plot_ridge_build
from outflow_maskprev import plot_ridge_IOP3
# these functions are annoyingly spacific: please generalise

from circ_int_test import main_cc, plot_timeseries_levs
from circ_split import plot_timeseries_split

from sum_circ_split import Combine_IOP3_ABsplit, Combine_IOP6, Combine_all_IOP3, Combine_all_IOP6

from density_plt import calc_density

from July1 import outflow_grid, inflow_caltra, heating_hist, random_3D_trajs, inflow_time


## Define initial outflow volumes

#outflow_th(k=0, thlevs = 3, thetastart = [8, 9, 6, 6], hoursafterinit = [42, 42, 42, 24], filtersize = 30, resolution = 0.05,
#               folder = 'IOP3/T42', strapp = ''):

outflow_th()
# IOP3 standard
outflow_th(thlevs = 5, thetastart = [6], hoursafterinit = [24], folder = 'IOP3/ABsplit', strapp = '0')
# IOP3 fist outflow of two-stage, to see cyclonic branch
outflow_th(k = 1, folder = 'IOP5/T42')
# IOP5 standard
outflow_th(k = 1, hoursafterinit = [42, 36], folder = 'IOP5/T36')
# IOP5 6 hours earlier to test min. in circulation at T+42
outflow_th(k = 1, hoursafterinit = [42, 48], folder = 'IOP5/T48')
# IOP5 6 hours later to test sensitivity to definition time
outflow_th(k = 2, folder = 'IOP6')
# IOP6 first outflow of 4


## Calculate trajectories from these initial outflow volumes
#  (as they are needed for the definition of the subsequent volumes)

#try_caltra(k = 0, levs = 3, hoursafterinit = [42, 42, 42, 24], filtersize = 30, folder = 'IOP3/T42', strapp = '')

try_caltra()
try_caltra_rev()
# IOP3 standard
try_caltra(levs = 5, hoursafterinit = [24], folder = 'IOP3/ABsplit', strapp = '0')
try_caltra_rev(levs = 5, hoursafterinit = [24], folder = 'IOP3/ABsplit', strapp = '0')
# IOP3 fist outflow of two-stage
try_caltra(k = 1, folder = 'IOP5/T42')
try_caltra_rev(k = 1, folder = 'IOP5/T42')
# IOP5 standard
try_caltra(k = 1, hoursafterinit = [42, 36], folder = 'IOP5/T36')
try_caltra_rev(k = 1, hoursafterinit = [42, 36], folder = 'IOP5/T36')
# IOP5 6 hours earlier
try_caltra(k = 1, hoursafterinit = [42, 48], folder = 'IOP5/T48')
try_caltra_rev(k = 1, hoursafterinit = [42, 48], folder = 'IOP5/T48')
# IOP5 6 hours later
try_caltra(k = 2, folder = 'IOP6', strapp = '0')
try_caltra_rev(k = 2, folder = 'IOP6', strapp = '0')
# IOP6 first outflow of 4


## Define secondary outflow volumes

#outflow_from_old(k = 1, lag = 1, folderin = 'IOP5/T42', folderout = 'IOP5/T42-6', strapp = '')
#outflow_AB(k=0, thlevs = 3, thetastart = [6, 9, 6, 6, 6, 6], hoursafterinit = [60, 42, 42, 72, 96, 108], 
#               filtersize = 30, resolution = 0.05, masks = 1, folder = 'IOP3/ABsplit', strapp = '')

outflow_AB(thlevs = 5, strapp = '1')
# IOP3 second outflow of two-stage, to see later outflow around first
outflow_th(hoursafterinit = [60], folder = 'IOP3/T60')
# IOP3 defined as same time as above second outflow, but without any prior removal
outflow_from_old()
# IOP5 6 hours earlier, but started from volume defined at T42 and
# advected back 6 horus to tesd differences with T42 & T36
outflow_AB(k = 3, folder = 'IOP6', strapp = '1')
# IOP6 second outflow of 4, completely distinct from the first


## Calculate trajectories from these secondary outflow volumes

try_caltra(levs = 5, hoursafterinit = [60], folder = 'IOP3/ABsplit', strapp = '1')
try_caltra_rev(levs = 5, hoursafterinit = [60], folder = 'IOP3/ABsplit', strapp = '1')
# IOP3 second outflow of two-stage
try_caltra(hoursafterinit = [60], folder = 'IOP3/T60')
try_caltra_rev(hoursafterinit = [60], folder = 'IOP3/T60')
# IOP3 two-in-one
try_caltra(k = 1, hoursafterinit = [42, 36], folder = 'IOP5/T42-6')
try_caltra_rev(k = 1, hoursafterinit = [42, 36], folder = 'IOP5/T42-6')
# IOP5 6 hours earlier (II)
try_caltra(k = 2, hoursafterinit = [42, 42, 72], folder = 'IOP6', strapp = '1')
try_caltra_rev(k = 2, hoursafterinit = [42, 42, 72], folder = 'IOP6', strapp = '1')
# IOP6 second outflow of 4


## Define tertiary outflow volumes
#  (I think this is the last stage of defining these bloody things)

outflow_AB(k = 4, masks = 2, folder = 'IOP6', strapp = '2')
# IOP6 third outflow of 4
outflow_th(k = 5, thetastart = [6, 6, 6, 6, 6, 6], hoursafterinit = [0, 0, 0, 0, 0, 108], folder = 'IOP6', strapp = '3')
# IOP6 fourth & final outflow: frontal wave cyclone (didn't actually need masking, just fine PV division)
outflow_th(k = 3, folder = 'IOP7')
# IOP7 same as fourth outflow in "IOP6" but from forecast started 72 hours later - can be used to assess predictability?


## Calculate trajectories from tertiary outflow volumes

try_caltra(k = 2, hoursafterinit = [42, 42, 96], folder = 'IOP6', strapp = '2')
try_caltra_rev(k = 2, hoursafterinit = [42, 42, 96], folder = 'IOP6', strapp = '2')
# IOP6 third outflow of 4
try_caltra(k = 2, hoursafterinit = [42, 42, 108], folder = 'IOP6', strapp = '3')
try_caltra_rev(k = 2, hoursafterinit = [42, 42, 108], folder = 'IOP6', strapp = '3')
# IOP6 fourth outflow of 4
try_caltra(k = 3, folder = 'IOP7')
try_caltra_rev(k = 3, folder = 'IOP7')
# IOP7


## Define another outflow
# this is getting ridiculous

outflow_AB(k = 6, thlevs = 1, thetastart = [0, 0, 0, 0, 0, 0, 8], hoursafterinit = [0, 0, 0, 0, 0, 0, 144], filtersize = 40, masks = 1, folder = 'IOP6', strapp = '4')
# IOP7, secondary outflow on 320K surface only
try_caltra(k = 2, levs = 1, hoursafterinit = [0, 0, 144], filtersize = 40, folder = 'IOP6', strapp = '4')
try_caltra_rev(k = 2, levs = 1, hoursafterinit = [0, 0, 144], filtersize = 40, folder = 'IOP6', strapp = '4')



## Plots to illustrate trajectories

#plot_traj(k = 0, levs = 3, filtersize = 30, field = 'total_minus_adv_only_theta', lim = 40, indiv = False, folder = 'IOP3/T42', strapp = '')

for idv in [True, False]:
    plot_traj(indiv = idv)
    # IOP3 standard
    plot_ridge_IOP3(indiv = idv)
    # IOP3 two-stage
    plot_traj(indiv = idv, folder = 'IOP3/T60')
    # IOP3 two-in-one
    plot_traj(k = 1, indiv = idv, folder = 'IOP5/T42')
    # IOP5 standard
    plot_traj(k = 1, indiv = idv, folder = 'IOP5/T36')
    # IOP5 6 hours earlier
    plot_traj(k = 1, indiv = idv, folder = 'IOP5/T48')
    # IOP5 6 hours later
    plot_traj(k = 1, indiv = idv, folder = 'IOP5/T42-6')
    # IOP5 6 hours earlier (II)
    plot_ridge_build(indiv = idv)
    # IOP6
    plot_traj(k = 3, indiv = idv, folder = 'IOP7')
    # IOP7


## Circulation calculations (for 10 sets of trajectories)

#main_cc(k = 0, filtersize = 30, theta_level = 325, dtheta = 1, split = False, folder = 'IOP3/T42', strapp = '')

for splt in [True, False]:
    for theta_lev in [320, 325, 330]:
        
        main_cc(theta_level = theta_lev, split = splt)
    # IOP3 standard
        
    for theta_lev in [310, 315, 320, 325, 330]:    
        
        for strpp in ['0', '1']:
            
            main_cc(theta_level = theta_lev, split = splt, folder = 'IOP3/ABsplit', strapp = strpp)
    # IOP3 two-stage
            
    for theta_lev in [320, 325, 330]:
        
        main_cc(theta_level = theta_lev, split = splt, folder = 'IOP3/T60')
    # IOP3 two-in-one

    for theta_lev in [325, 330, 335]:

        for foldr in ['IOP5/T42', 'IOP5/T36', 'IOP5/T48', 'IOP5/T42-6']:

            main_cc(k = 1, theta_level = theta_lev, split = splt, folder = foldr)
    # IOP5 all

    for theta_lev in [310, 315, 320]:

        for strpp in ['0', '1', '2', '3']:
            
            main_cc(k = 2, theta_level = theta_lev, split = splt, folder = 'IOP6', strapp = strpp)
            
        main_cc(k = 3, theta_level = theta_lev, split = splt, folder = 'IOP7')
    # IOP6 & IOP7
            
            
## Combine circulation for multi-outflow cases

for stp in ['', '_split']:
    Combine_IOP3_ABsplit(strapp = stp)
#    Combine_all_IOP3(strapp = stp) doesn't work
    # IOP3 two-stage
    Combine_IOP6(strapp = stp)
#    Combine_all_IOP6(strapp = stp) doesn't work
    # IOP6      
            
            
## Circulation plots

#plot_timeseries_levs(case = 0, levs = 3, tls = [320, 325, 310, 310], folder = 'IOP3/T42', strapp = ''):

plot_timeseries_levs()
# IOP3 standard
for strpp in ['0', '1']:
    plot_timeseries_levs(levs = 5, tls = [310], folder = 'IOP3/ABsplit', strapp = strpp)
# IOP3 two-stage
plot_timeseries_levs(folder = 'IOP3/T60')
# IOP3 two-in-one
plot_timeseries_levs(case = 1, folder = 'IOP5/T42')
# IOP5 standard
plot_timeseries_levs(case = 1, folder = 'IOP5/T36')
# IOP5 6 hours earlier
plot_timeseries_levs(case = 1, folder = 'IOP5/T48')
# IOP5 6 hours later
plot_timeseries_levs(case = 1, folder = 'IOP5/T42-6')
# IOP5 6 hours earlier (II)
for strpp in ['0', '1', '2', '3']:
    plot_timeseries_levs(case = 2, folder = 'IOP6', strapp = strpp)
# IOP6
plot_timeseries_levs(case = 3, folder = 'IOP7')
# IOP7


## Split components and circulation divided by area (= mean vorticity?)

for diva in [True, False]:
    plot_timeseries_split(divarea = diva)
    # IOP3 standard
    for strpp in ['0', '1']:
        plot_timeseries_split(levs = 5, tls = [310], divarea = diva, folder = 'IOP3/ABsplit', strapp = strpp)
    # IOP3 two-stage
    plot_timeseries_split(divarea = diva, folder = 'IOP3/T60')
    # IOP3 two-in-one
    plot_timeseries_split(case = 1, divarea = diva, folder = 'IOP5/T42')
    # IOP5 standard
    plot_timeseries_split(case = 1, divarea = diva, folder = 'IOP5/T36')
    # IOP5 6 hours earlier
    plot_timeseries_split(case = 1, divarea = diva, folder = 'IOP5/T48')
    # IOP5 6 hours later
    plot_timeseries_split(case = 1, divarea = diva, folder = 'IOP5/T42-6')
    # IOP5 6 hours earlier (II)
    for strpp in ['0', '1', '2', '3']:
        plot_timeseries_split(case = 2, divarea = diva, folder = 'IOP6', strapp = strpp)
    # IOP6
    plot_timeseries_split(case = 3, divarea = diva, folder = 'IOP7')
    # IOP7





## Calculate cubes of PV substance and isentropic density

#calc_density(k = 0, thlevs = [[320, 325, 330], [325, 330, 335], [310, 315, 320], [310, 315, 320]], plotv = True, folder = 'IOP3/T42', strapp = ''):

calc_density()
# IOP3 standard
calc_density(k = 1, folder = 'IOP5/T42')
# IOP5 standard
calc_density(k = 2, folder = 'IOP6', strapp = '0')
# IOP6 first
calc_density(k = 3, folder = 'IOP7')
# IOP7



## Define grid of points filling defined outflow region

#outflow_grid(k = 0, levels = 3, hoursafterinit = [42, 42, 42, 24], thlevs = [[320, 325, 330], [325, 330, 335], [310, 315, 320], [310, 315, 320]], folder = 'IOP3/T42', strapp = '')

outflow_grid()
# IOP3 standard
outflow_grid(levels = 5, hoursafterinit = [24], thlevs = [[310, 315, 320, 325, 330]], folder = 'IOP3/ABsplit', strapp = '0')
outflow_grid(levels = 5, hoursafterinit = [60], thlevs = [[310, 315, 320, 325, 330]], folder = 'IOP3/ABsplit', strapp = '1')
# IOP3 two-stage
outflow_grid(hoursafterinit = [60], folder = 'IOP3/T60')
# IOP3 two-in-one
outflow_grid(k = 1, folder = 'IOP5/T42')
# IOP5 standard
outflow_grid(k = 1, hoursafterinit = [42, 36], folder = 'IOP5/T36')
# IOP5 6 hours earlier
outflow_grid(k = 1, hoursafterinit = [42, 48], folder = 'IOP5/T48')
# IOP5 6 hours later
for hs in [[42, '0'], [72, '1'], [96, '2'], [108, '3']]:
    outflow_grid(k = 2, hoursafterinit = [0, 0, hs[0]], folder = 'IOP6', strapp = hs[1])
# IOP6
outflow_grid(k = 3, folder = 'IOP7')


## Calculate 3D back-trajectories

# inflow_caltra(k = 0, levs = 3, hoursafterinit = [42, 42, 42, 24], folder = 'IOP3/T42', strapp = '')
   
inflow_caltra()
# IOP3 standard
inflow_caltra(levs = 5, hoursafterinit = [24], folder = 'IOP3/ABsplit', strapp = '0')
inflow_caltra(levs = 5, hoursafterinit = [60], folder = 'IOP3/ABsplit', strapp = '1')
# IOP3 two-stage
inflow_caltra(hoursafterinit = [60], folder = 'IOP3/T60')
# IOP3 two-in-one
inflow_caltra(k = 1, folder = 'IOP5/T42')
# IOP5 standard
inflow_caltra(k = 1, hoursafterinit = [42, 36], folder = 'IOP5/T36')
# IOP5 6 hours earlier
inflow_caltra(k = 1, hoursafterinit = [42, 48], folder = 'IOP5/T48')
# IOP5 6 hours later
for hs in [[42, '0'], [72, '1'], [96, '2'], [108, '3']]:
    inflow_caltra(k = 2, hoursafterinit = [0, 0, hs[0]], folder = 'IOP6', strapp = hs[1])
# IOP6
inflow_caltra(k = 3, folder = 'IOP7')


## Plot histogram of delta theta along trajectories

# heating_hist(k = 0, levs = 3, theta_0 = [320, 325, 310, 310], folder = 'IOP3/T42', strapp = '')

heating_hist()
# IOP3 standard
heating_hist(levs = 5, theta_0 = [310], folder = 'IOP3/ABsplit', strapp = '0')
heating_hist(levs = 5, theta_0 = [310], folder = 'IOP3/ABsplit', strapp = '1')
# IOP3 two-stage
heating_hist(folder = 'IOP3/T60')
# IOP3 two-in-one
heating_hist(k = 1, folder = 'IOP5/T42')
# IOP5 standard
heating_hist(k = 1, folder = 'IOP5/T36')
# IOP5 6 hours earlier
heating_hist(k = 1, folder = 'IOP5/T48')
# IOP5 6 hours later
for s in ['0', '1', '2', '3']:
    heating_hist(k = 2, folder = 'IOP6', strapp = s)
# IOP6
heating_hist(k = 3, folder = 'IOP7')


## Plot selection of random trajectories theta against time

# random_3D_trajs(k = 0, levs = 3, theta_0 = [320, 325, 310, 310], trajs = 30, folder = 'IOP3/T42', strapp = '')

random_3D_trajs()
# IOP3 standard
random_3D_trajs(levs = 5, theta_0 = [310], folder = 'IOP3/ABsplit', strapp = '0')
random_3D_trajs(levs = 5, theta_0 = [310], folder = 'IOP3/ABsplit', strapp = '1')
# IOP3 two-stage
random_3D_trajs(folder = 'IOP3/T60')
# IOP3 two-in-one
random_3D_trajs(k = 1, folder = 'IOP5/T42')
# IOP5 standard
random_3D_trajs(k = 1, folder = 'IOP5/T36')
# IOP5 6 hours earlier
random_3D_trajs(k = 1, folder = 'IOP5/T48')
# IOP5 6 hours later
for s in ['0', '1', '2', '3']:
    random_3D_trajs(k = 2, folder = 'IOP6', strapp = s)
# IOP6
random_3D_trajs(k = 3, folder = 'IOP7')


## Plot isentropic trajectory & 3D trajectory on same figure

#plot_traj(k = 0, levs = 3, filtersize = 30, field = 'total_minus_adv_only_theta', scatter3D = False, lim = 40, indiv = False, folder = 'IOP3/T42', strapp = '')

plot_traj(field = 'air_potential_temperature', scatter3D = True)
# IOP3 standard
plot_traj(levs = 5, field = 'air_potential_temperature', scatter3D = True, folder = 'IOP3/ABsplit', strapp = '0')
plot_traj(levs = 5, field = 'air_potential_temperature', scatter3D = True, folder = 'IOP3/ABsplit', strapp = '1')
# IOP3 two-stage
plot_traj(field = 'air_potential_temperature', scatter3D = True, folder = 'IOP3/T60')
# IOP3 two-in-one
plot_traj(k = 1, field = 'air_potential_temperature', scatter3D = True, folder = 'IOP5/T42')
# IOP5 standard
plot_traj(k = 1, field = 'air_potential_temperature', scatter3D = True, folder = 'IOP5/T36')
# IOP5 6 hours earlier
plot_traj(k = 1, field = 'air_potential_temperature', scatter3D = True, folder = 'IOP5/T48')
# IOP5 6 hours later
for s in ['0', '1', '2', '3']:
    plot_traj(k = 2, field = 'air_potential_temperature', scatter3D = True, folder = 'IOP6', strapp = s)
# IOP6
plot_traj(k = 3, field = 'air_potential_temperature', scatter3D = True, folder = 'IOP7')




inflow_time()
#inflow_time(hoursafterinit = [60], folder = 'IOP3/T60')
inflow_time(k = 1, hoursafterinit = [0, 36], folder = 'IOP5/T36')
for hs in [[42, '0'], [72, '1'], [96, '2'], [108, '3']]:
    inflow_time(k = 2, hoursafterinit = [0, 0, hs[0]], folder = 'IOP6', strapp = hs[1])
inflow_time(k = 3, folder = 'IOP7')

# inflow bounds definitions are in circ_int_test
# other new things are in outflow_2019 or plots_2019