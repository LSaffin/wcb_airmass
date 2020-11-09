"""Code to be introduced to actual code once I'm at my desk
"""

def inflow_contours_for_isentropic_trajectories(case, resolution = 5):
  # load defined inflow contours
  boundary = np.load(str(case.data_path / "inflow_boundaries.npy"))
  
  # increase number of points around contour consistent with outflow
  increase_circuit_resolution(boundary, resolution)
  
  # concatenate over all levels below outflow layer
  for theta_level in np.arange(300, min(case.outflow_theta), 5):
    boundary3d = np.zeros([conti.shape[0], 3)
    boundary3d[:, 0:2] = boundary
    boundary3d[:, 2] = theta_level
    
    if theta == 300:
      inflow_boundaries = boundary3d
    else:
      inflow_boundaries = np.concatenate((inflow_boundaries, boundary3d))
      
  np.save(str(case.data_path / "inflow_boundaries_3d.npy"), inflow_boundaries)
  
  
def isentropic_inflow_trajectories(case):
  # get trajectories for inflow contours
  # mostly just for the first point for figure 8
  trainp = np.load(str(case.data_path / "inflow_boundaries_3d.npy")
  trainp = np.vstack(
    trainp[np.where(trainp[:, 2] == theta)] for theta in np.arange(300, min(case.outflow_theta), 5))
  
  levels = ("air_potential_temperature", np.arange(300, min(case.outflow_theta), 5))
  
  fbflag = 1
  fname = "forward"
  mapping = case.time_to_filename_winds_mapping()
  
  traout = caltra.caltra(
    trainp,
    mapping,
    nsubs = 12,
    fbflag = fbflag,
    levels = levels,
    tracers=["x_wind", "y_wind", "upward_air_velocity", "altitude"]
    )
    
  traout.save(str(case.data_path / "inflow_isentropic_trajectories.pkl")
  
  
def 
  
  
  
