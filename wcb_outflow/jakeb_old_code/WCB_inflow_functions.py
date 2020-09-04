      
 # slightly out of place tack-ons       
        
def shadow(trajectorys, fcast, theta_levels, dtheta, start_time, chosen_time_index, folder, strapp):
    """
    Shitty hack to return sum of shadows
    Chosen_time_index 
    trajectorys deliberately misspelt
    definitely move somewhere else when done
    """
    masks = []
    for theta_level in theta_levels:
    
        # Select an individual theta level
        trajectories = trajectorys.select(
            'air_potential_temperature', '==', theta_level, time = start_time)
        # JB this allows trajectories which leave the domain to be selected also
        # JB this is dealt with later (start_time addition)
        print len(trajectories)
        levels = ('air_potential_temperature', theta_levels)

        time = trajectories.times[chosen_time_index]
        n = chosen_time_index
        cubelist = fcast.set_time(time)
        #now for some reason I'm having real problems extracting the right time, so
#        u = iris.unit.Unit('hours since 1970-01-01 00:00:00', calendar=iris.unit.CALENDAR_STANDARD)
#        timeh = u.date2num(time)
        #this puts the time in the same format as in the cubelist
        cubes = cubelist.extract(iris.Constraint(time = time))
        # NOTE: THIS USED TO ONLY WORK WITH timeh, NOW ONLY WORKS WITH time ???
        # JB needed because each file contains two times
        #print cubelist
        
        if 1:
        
            # Load grid parameters
            example_cube = convert.calc('upward_air_velocity', cubes,
                                        levels=levels)

            # Create a 1d array of points for determining which gridpoints are
            # contained in the trajectory circuit when performing volume
            # integrals
            glon, glat = grid.get_xy_grids(example_cube)
            gridpoints = np.array([glon.flatten(), glat.flatten()]).transpose()
            cds = example_cube.coord_system()
            
            # JB to do integrals need dlambda, dphi, not always 0.11
            
            dlon = np.diff(glon)
            if np.diff(dlon).any():
                print 'longitudinal spacing not uniform'
                # there exists a non-zero difference between longitude spacings
            else:
                dlambda = dlon[0][0]
                #as they are all the same
                
            dlat = np.diff(glat.transpose())
            if np.diff(dlat).any():
                print 'latitudinal spacing not uniform'
                # there exists a non-zero difference between latitude spacings
            else:
                dphi = dlat[0][0]
                #as they are all the same
                

        # Load trajectory positions -(n+2) because the trajectories are
        # backwards in time. +2 to skip the analysis which is not in the
        # forecast object (i.e. n=0 corresponds to idx=-2 in the trajectories)

        #JB I'm making them forwards in time so no worries, load n
        x = trajectories.x[:, n]
        y = trajectories.y[:, n]
        z = trajectories['altitude'][:, n]
        u = trajectories['x_wind'][:, n]
        v = trajectories['y_wind'][:, n]
        w = trajectories['upward_air_velocity'][:, n]
    
        # Include points within circuit boundary
        points = np.array([x, y]).transpose()
        pth = Path(points)

        # Mask all points that are not contained in the circuit
        mask = np.logical_not(pth.contains_points(gridpoints).reshape(glat.shape)) 
    
        masks.append(mask)
        
    sumask = np.array(masks[0], dtype = int)
    for mask in masks[1:]:
        sumask = sumask + np.array(mask, dtype = int)
    value = len(masks) - 0.5  
    
    #sumooth = filters.median_filter(sumask, size = 20)
    
    cs = plt.contour(glon, glat, sumask, [value])
    
    contours = trop.get_contour_verts(cs)
    # returns array of arrays of vertices for zero contours
    
    #ncon = np.size(contours[0])
    ncon = len(contours[0])
    #number of contours
    
    lencon = np.zeros(ncon)
    # empty array of lengths of contorus (in lat/long space)
    
    print ncon
    for j in xrange(ncon):
        print j
        
        lencon[j] = len_con(contours[0][j])
    
    imax = lencon.argmax()
    # index of longest contour
    lcontour = contours[0][imax]
    # longest contour
    # don't worry about closed-ness
    
    points = increase_nodes(lcontour, resolution = .25)
    # increase number of points on the contour such that they have minimum spacing resolution    
    
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/outflow_shadow_contour.npy', points)
    
    
def bound_inflow_3D(k = 0, folder = 'IOP3/T42', in_time = 5, min_lev = 320, strapp = '', gmod = False):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    Tr3 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_new'.format(basetime_str) + strapp)
    
    Tr3_inflow = Tr3.select('air_potential_temperature', '<', min_lev - 2.5, time = [datetime.timedelta(hours = -6*in_time)])
    
    inf_area = plt.hist2d(Tr3_inflow.data[:, in_time, 0], Tr3_inflow.data[:, in_time, 1], bins = [360, 90], range = [[0, 360], [0, 90]])
    plt.show()
    if gmod:
        ia = filters.gaussian_filter(inf_area[0], gmod)
    else:
        ia = inf_area[0]
    cs = plt.contour(np.transpose(ia), [.5])
    
    contours = trop.get_contour_verts(cs)
    ncon = len(contours[0])
    #number of contours

    lencon = np.zeros(ncon)
    # empty array of lengths of contorus (in lat/long space)

    for j in xrange(ncon):

       lencon[j] = len_con(contours[0][j])

    imax = lencon.argmax()
    # index of longest contour
    lcontour = contours[0][imax]
    # longest contour
    # don't worry about closed-ness

    points = increase_nodes(lcontour, resolution = .25)
    
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow_bound_3D_less_than_' + str(min_lev) + 'K_gmod.npy', points)
    
    
def bound_outflow_3D(k = 0, folder = 'IOP3/T42', in_time = 5, min_lev = 320, strapp = '', gmod = False):
    
    basetime = [datetime.datetime(2016, 9, 22, 12), datetime.datetime(2016, 9, 26, 12), datetime.datetime(2016, 9, 30, 12), datetime.datetime(2016, 10, 03, 12)]
    basetime_str = basetime[k].strftime('%Y%m%d_%H')
    
    Tr3 = load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow/{}_3DTrajectoryEnsemble_new'.format(basetime_str) + strapp)
    
    Tr3_inflow = Tr3.select('air_potential_temperature', '>', min_lev - 2.5, time = [datetime.timedelta(hours = -6*in_time)])
    
    inf_area = plt.hist2d(Tr3_inflow.data[:, in_time, 0], Tr3_inflow.data[:, in_time, 1], bins = [360, 90], range = [[0, 360], [0, 90]])
    plt.show()
    if gmod:
        ia = filters.gaussian_filter(inf_area[0], gmod)
    else:
        ia = inf_area[0]
    cs = plt.contour(np.transpose(ia), [.5])
    
    contours = trop.get_contour_verts(cs)
    ncon = len(contours[0])
    #number of contours

    lencon = np.zeros(ncon)
    # empty array of lengths of contorus (in lat/long space)

    for j in xrange(ncon):
        print ncon
        print j

        lencon[j] = len_con(contours[0][j])

    imax = lencon.argmax()
    # index of longest contour
    lcontour = contours[0][imax]
    # longest contour
    # don't worry about closed-ness

    points = increase_nodes(lcontour, resolution = .25)
    
    np.save('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/outflow_bound_3D_less_than_' + str(min_lev) + 'K_gmod.npy', points)
   
   
   
def compare_inflow_bounds(folder = 'IOP3/T42', min_theta = 320):
    
    shadow_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/'+ folder + '/shadow_contour.npy')
    
    smooth_shadow_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/'+ folder + '/smooth_shadow_contour.npy')

    bound_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow_bound_3D_less_than_' + str(min_theta) + 'K.npy')
    
    big_cont = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/inflow_bound_3D_less_than_' + str(min_theta) + 'K_gmod.npy')
    
    bound_outf = np.load('/storage/silver/scenario/bn826011/WCB_outflow/Final/' + folder + '/outflow_bound_3D_less_than_' + str(min_theta) + 'K.npy')

    plt.plot(shadow_cont[:, 0], shadow_cont[:, 1], label = 'sum of shadows of isentropic V2 at ts')
    plt.plot(smooth_shadow_cont[:, 0], smooth_shadow_cont[:, 1], linestyle = ':', linewidth = 4, label = 'smoothed sum of shadows of isentropic V2 at ts')
    plt.plot(bound_outf[:, 0], bound_outf[:, 1],  linestyle = '--', label = 'bound around 3D trajectories at ts within V2 (should == shadow)')
    plt.plot(bound_cont[:, 0], bound_cont[:, 1], linewidth = 1.5, label = 'bound around 3D trajectories at ts below V2')
    plt.plot(big_cont[:, 0], big_cont[:, 1], linewidth = 2, label = 'bound around 3D trajectories at ts below V2, slightly larger')
    
    plt.title(folder + '  comparison of inflow lateral bounds')
    plt.legend(bbox_to_anchor=(1.2, -.1))
    
    plt.savefig('/home/users/bn826011/NAWDEX/From N Drive/2019_figs/' + folder[:4] + '_bounds_comparison2.jpg')
    plt.show()