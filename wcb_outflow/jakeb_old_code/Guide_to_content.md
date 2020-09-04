General notes:
- IOP:     3 : k = 0,   5 : k = 1,    6 : k = 2,    7 : k = 3
- theta levels are referenced by indices i within files such that:
  theta = 280 + 5*i
- times are in hours from initial time (at least mostly)
- "strapp" is empty for all cases other than IOP6, for which there are 4 different outflow volumes defined indexed 0-3
- directories are mostly hard coded within functions... sorry
- files contain a lot of functions, organisation is poor, below I will list the useful functions and the files that they are in

Figures and tables in paper, where to find relevant bit of code:

	[Figure 6](#figure-6)
	[Table 2](#table-2)
  
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
  
  Function name: plot_timeseries_levs()
  ###### Figure 6
  - I believe this is the function to produce
  
  Function name: print_percent_change()
  ###### Table 2
  - prints in LaTeX table format change in circulation, area, volume and mass of the outflow volume V2 between the inflow time 0 and outflow time tau, 
    as a percentage of the value at outflow time tau
    
  Function name: print_mass_change()
  - Prints in LaTeX table format the change in the mass of the ridge and the change in the mass of the outflow volume between inflow and outflow times
  
  
  
  
