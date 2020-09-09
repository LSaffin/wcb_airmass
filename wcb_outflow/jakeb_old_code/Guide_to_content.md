General notes:
- IOP:     3 : k = 0,   5 : k = 1,    6 : k = 2,    7 : k = 3
- theta levels are referenced by indices i within files such that:
  theta = 280 + 5*i
- times are in hours from initial time (at least mostly)
- "strapp" is empty for all cases other than IOP6, for which there are 4 different outflow volumes defined indexed 0-3
- directories are mostly hard coded within functions... sorry
- files contain a lot of functions, organisation is poor, below I will list the useful functions and the files that they are in

Where to find relevant bit of code:

[Grid of points from contour](#outflow-grid)

[Fig. 4](#figure-4)

[Fig. 5](#figure-5)

[Fig. 6](#figure-6)

[Fig. 7](#figure-7)

[Fig. 8](#figure-8)

[Tabular 2](#table-2)
  
Filename: Everything_script.py

  Does what it says on the tin, this is a rundown of everything I did to produce the data and figures for the dissertation

Filename: outflow_area_2019.py

  Function name: outflow_th()
  - this is the function used to define the bounding contour for the outflow volumes on isentropic levels

    Specifically, returns numpy array of indices on contour around outflow, into file:
	  IOP#/mmdd_hh.npy,  e.g. IOP3/T42/0924_06.npy

	  Take ‘total_minus_advection_only_theta’, i.e., diabatic theta, and PV from files on theta levels, every 5K
	  Reduce size of data to predetermined outflow time, theta level and 2d domain
	  Apply median filter
	  Draw contours around regions which satisfy diabatic theta > 0 and PV < 2, i.e. tropospheric regions which have experienced cross-isentropic ascent of air into them (the points comprising the contours themselves having experienced ≈ no net heating)
	  Identify the largest contour
	  Save set of contours to file mmdd_hh.npy

  Function name: try_caltra()
  - used to calculate isentropc backwards trajectories from input grid of points
  
  Function name: try_caltra_rev()
  - as above but for trajectories forwards in time

Filename: circ_split.py
  
  Contains extensions of circulation.circuit_integrals and circulation.mass_integrals to return separate components as for the figure which we may not include
  
Filename: circ_int_test.py

  Further extensions of functions in circulation.py to use on trajectory ensembles, produce plots, and print percentage changes in a format easily pasteable into a LaTeX table
  (stop reading at "#slightly out of place tack-ons", they were very out of place and have been copied into new file WCB_inflow_functions.py)
  
  Function name: main_cc()
  - wrapper to stitch together trajectories and feed into circulation integrals
  
  Function name: calc_circulation()
  - copy of calculate.calc_circulation to facilitate the splitting of circulation into components as an option
  
###### Figure 6
  
  Function name: plot_timeseries_levs()
  - I believe this is the function to produce figure 6
  
###### Table 2
  
  Function name: print_percent_change()
  - prints in LaTeX table format change in circulation, area, volume and mass of the outflow volume V2 between the inflow time 0 and outflow time tau, 
    as a percentage of the value at outflow time tau
    
  Function name: print_mass_change()
  - Prints in LaTeX table format the change in the mass of the ridge and the change in the mass of the outflow volume between inflow and outflow times
  
Filename: produce_figs.py
  this file contains many plotting functions
  and for some reason the function I used to define ridge area?
  
  Function name: out_ridge_mass()
  
Filename: WCB_inflow_functions.py
  
  Contains functions used to define inflow volume as
  - sum of shadows of outflow volumes (shadow)
  - area bounded by 3D trajectories above the bottom of the outflow volume (bound_outflow_3D)
  - area bounded by 3D trajectories which are below the outflow volume (bound_inflow_3D)
  at inflow time
  
###### Figure 5
  
  is produced using compare_inflow_bounds()
  notes: I believe gmod refers to the use of a gaussian filter? your guess may be as good as mine
  
Filename: July1.py

  Function name: outflow_grid()
  
###### Outflow grid
  - Take contour (of outflow volume at outflow time) and define grid of points within
    This is probably the one which resulted in the incorrect start points for 3d trajectories
  
  Function name: inflow_caltra()
    Wrapper to generate 3D trajectory ensembles from grid of points, backwards in time
    
  Function name: heating_hist()
  - I believe this was used to define inflow times as times when between them & 6 hours before there was little heating, 
    and the majority of the heating occured subsequent to these, 
    to prevent excessive filamentation of the contours and make integrals easier
  
Filename: plots_2019.py

###### Figure 4
  
  Function name: trajectory_2D_3D_displacement_difference()
  - greyscale density plots of displacement at inflow time between 2D and 3D back-trajectories from outflow time
  - problem is in L76
  
###### Figure 8
  
  Function name: inflow_quantities()
  - probably
  
###### Figure 7
  
  Function name: plot_timeseries_new()
  - in here the ratios of PV, vorticity, mass, area, density etc. are defined and plotted
  - might be worth double checking the quantities we actually want to calculate and plot, but here's a skeleton
  
  Function name: this_would_do_integrals()
  - basically a script to bulk calculate circulation integrals on inflow volumes
  
  
