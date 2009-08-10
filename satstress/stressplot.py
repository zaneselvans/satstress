#!/usr/bin/python
"""
Create map-projected vector and scalar plots of tidal stresses from netCDF
files output by L{gridcalc}.

Inputs:
=======
A NetCDF file containing one ore more set of scalar variables named Ttt_*,
Tpt_*, Tpp_*, representing the north-south (theta-theta), east-west (phi-phi)
and shear (phi-theta) elements of the surface membrane stress tensor.

A subset of the following stress field names to include in the plot:

  1. B{NSR} (corresponding to stresses produced by L{satstress.NSR}).

  2. B{Flat} (corresponding to stresses produced by L{satstress.Flat}).
  
  3. B{Diurnal} (corresponding to stresses produced by L{satstress.Diurnal}).

A specification of the orbital position, either as time elapsed since periapse,
or as degrees beyond periapse (must correspond to a time-slice which is present
within the NetCDF file).

A specification of the NSR period (must correspond to a value o NSR period
present within the NetCDF file).

At most one magnitude field to display, chosen from:

  * Magnitude of the most tensile principal component of the membrane stress

  * Magnitude of the most compressive principal component of the membrane stress

A flag indicating whether or not the principal components of the sum of the
two stress tensors should be displayed.

The resolution at which the stress vectors should be displayed (in degrees
separating one vector from the next)

Outputs:
========
A plot showing:

  * the chosen magnitude field as the background, in color, with an appropriate
    scale bar alongside the plot.

  * the principal components of the stress vector field that results from
    adding each of the tensor components together, with positive (tensile)
    stresses displayed as red arrows, and negative (compressive) stresses
    displayed as blue arrows, along with an appropriate reference vector for
    scale.
    
"""
import numpy as np
from matplotlib import pyplot as plt
import satstress as ss
from mpl_toolkits import basemap
import netCDF3
import sys
from optparse import OptionParser

# This old_main() function is from the version of stressplot that I wrote when
# calculating stresses wa slow enough that we really needed to be storing the
# output.  Now that it's fast, we can just do the calculation on the fly, but
# the stress plotting is going to be very similar.  This function is left here
# temporarily for reference only.

def lonlat2vector(lons=None, lats=None, times_t=[0,0,], nlons=361, nlats=181): #{{{
    """
    A short helper function that makes sure we have two equal length vectors of
    longitudes and latitudes at which to do our calculations, and only a single
    time value.

    """

    if lons is None or lats is None:
        if lons is None and nlons > 0:
            lons = np.linspace(0,np.pi,nlons)
        if lats is None and nlats > 0:
            lats = np.linspace(-np.pi/2.0,np.pi/2.0,nlats)
        # If we're doing a gridded calculation, then we need to generate the full
        # list of lats and lons, paired appropriately, to pass in to the vectorized
        # StressCalc
        calc_phis, calc_thetas = np.meshgrid(lons, np.pi/2.0-lats)
        calc_phis = np.ravel(calc_phis)
        calc_thetas = np.ravel(calc_thetas)
    else:
        assert(len(lons)==len(lats))
        calc_phis = lons
        calc_thetas = np.pi/2.0-lats
    # If we got a single value, use it:
    times_t = np.array(times_t)
    if np.shape(times_t) is ():
        time_t = np.float(times_t)
    # If we got a list or an array, just use the first one:
    else:
        time_t = times_t[0]

    return (lons, calc_phis, lats, calc_thetas, time_t)
#}}} end latlon2vector

def scalar_grid(stresscalc=None,\
                nlons=13, nlats=13,\
                lon_min=0.0, lon_max=np.pi,\
                lat_min=-np.pi/2, lat_max=np.pi/2,\
                time_t=0.0, field='tens', basemap_ax=None,\
                cmap=plt.cm.jet): #{{{
    """
    Display a rasterized scalar stress field defined by stresscalc, at the
    resolution specified by nlons, nlats, within the box bounded by lon_min,
    lon_max and lat_min, lat_max, at a time defined by time_t.  Which scalar
    field is displayed is controlled by field, which may be set to any of the
    following strings:

        'tens' -- magnitude of the more tensile principal component
        'comp' -- magnitude of the more compressive principal component
         'Ttt' -- meridional component of the stress tensor (north-south)
         'Tpt' -- off-diagonal component of the stress tensor
         'Tpp' -- east-west component of the stress tensor
    'w_stress' -- significance of the stresses, viz. magnitude, isotropy

    The time of the calculation, in seconds after periapse, is defined by
    time_t.

    The axes the stresses are plotted on is defined by basemap_ax.

    cmap defines the colormap used to color the scalar field.

    """

    calc_phis, calc_thetas = meshgrid(np.linspace(lon_min, lon_max, nlons), (np.pi/2.0)-np.linspace(lat_min, lat_max, nlats))
    calc_phis   = np.ravel(calc_phis)
    calc_thetas = np.ravel(calc_thetas)

    # some of the possible fields are easier to compute with the principal
    # components:
    if field=='tens' or field=='comp' or field=='w_stress':
        tens_mag, tens_az, comp_mag, comp_az = stresscalc.principal_components(calc_thetas, calc_phis, time_t)
        if nlons > 0 and nlats > 0:
            tens_mag = tens_mag.reshape(nlats, nlons)
            comp_mag = comp_mag.reshape(nlats, nlons)
            comp_az = comp_az.reshape(nlats, nlons)
            tens_az = tens_az.reshape(nlats, nlons)
        if field=='w_stress':
            w_stress = (tens_mag - comp_mag)/stresscalc.mean_global_stressdiff()

    # Or if people want to see the raw tensor components we can do that too:
    if field=='Ttt' or field=='Tpt' or field=='Tpp':
        Ttt, Tpt, Tpp = stresscalc.tensor(calc_thetas, calc_phis, time_t)
        if nlons > 0 and nlats > 0:
            Ttt = Ttt.reshape(nlats, nlons)
            Tpt = Tpt.reshape(nlats, nlons)
            Tpp = Tpp.reshape(nlats, nlons)
    
    # Now we need to display the results of our calculations, For a gridded
    # calculation, we can just show a raster with imshow()
    if field=='tens':
        basemap_ax.imshow(tens_mag, cmap=cmap)
    if field=='comp':
        basemap_ax.imshow(comp_mag, cmap=cmap)
    if field=='w_stress':
        basemap_ax.imshow(w_stress, cmap=cmap)
    if field=='Ttt':
        basemap_ax.imshow(Ttt, cmap=cmap)
    if field=='Tpt':
        basemap_ax.imshow(Tpt, cmap=cmap)
    if field=='Tpp':
        basemap_ax.imshow(Tpp, cmap=cmap)

    return(basemap_ax)
#}}} end scalar

def scalar_points(stresscalc=None, lons=None, lats=None, time_t=0.0,\
                  field='tens', basemap_ax=None, cmap=plt.cm.jet): #{{{
    """
    Display a set of points within the scalar stress field defined by
    stresscalc, at the locations specified by lons, lats, and a time defined by
    time_t.  Which scalar field is displayed is controlled by field, which may
    be set to any of the following strings:

        'tens' -- magnitude of the more tensile principal component
        'comp' -- magnitude of the more compressive principal component
         'Ttt' -- meridional component of the stress tensor (north-south)
         'Tpt' -- off-diagonal component of the stress tensor
         'Tpp' -- east-west component of the stress tensor
    'w_stress' -- significance of the stresses, viz. magnitude, isotropy

    The time of the calculation, in seconds after periapse, is defined by
    time_t.

    The axes the stresses are plotted on is defined by basemap_ax.

    cmap defines the colormap used to color the scalar field.

    """
    calc_phis   = lons
    calc_thetas = (np.pi/2.0) - lats

    # some of the possible fields are easier to compute with the principal
    # components:
    if field=='tens' or field=='comp' or field=='w_stress':
        tens_mag, tens_az, comp_mag, comp_az = stresscalc.principal_components(calc_thetas, calc_phis, time_t)
        if nlons > 0 and nlats > 0:
            tens_mag = tens_mag.reshape(nlats, nlons)
            comp_mag = comp_mag.reshape(nlats, nlons)
            comp_az = comp_az.reshape(nlats, nlons)
            tens_az = tens_az.reshape(nlats, nlons)
        if field=='w_stress':
            w_stress = (tens_mag - comp_mag)/stresscalc.mean_global_stressdiff()

    # Or if people want to see the raw tensor components we can do that too:
    if field=='Ttt' or field=='Tpt' or field=='Tpp':
        Ttt, Tpt, Tpp = stresscalc.tensor(calc_thetas, calc_phis, time_t)
        if nlons > 0 and nlats > 0:
            Ttt = Ttt.reshape(nlats, nlons)
            Tpt = Tpt.reshape(nlats, nlons)
            Tpp = Tpp.reshape(nlats, nlons)

    if field=='tens':
        field_norm = plt.Normalize(vmin=np.min(tens_mag),vmax=np.max(tens_mag))
        colorarray = tens_mag
    if field=='comp':
        comp_norm = plt.Normalize(vmin=np.min(comp_mag),vmax=np.max(comp_mag))
        colorarray = comp_mag
    if field=='w_stress':
        w_stress_norm = plt.Normalize(vmin=np.min(w_stress),vmax=np.max(w_stress))
        colorarray = w_stress
    if field=='Ttt':
        Ttt_norm = plt.Normalize(vmin=np.min(Ttt),vmax=np.max(Ttt))
        colorarray = Ttt
    if field=='Tpt':
        Tpt_norm = plt.Normalize(vmin=np.min(Tpt),vmax=np.max(Tpt))
        colorarray = Tpt
    if field=='Tpp':
        Tpp_norm = plt.Normalize(vmin=np.min(Tpp),vmax=np.max(Tpp))
        colorarray = Tpp

    basemap_ax.scatter(np.degrees(lons), np.degrees(lats), lw=0, s=50, c=colorarray, cmap=cmap, norm=field_norm)

#}}}

def vector_points(stresscalc=None, lons=None, lats=None, time_t=0.0, plot_all=True,\
           plot_tens=False, plot_comp=False, plot_greater=False, plot_lesser=False,\
           basemap_ax=None, lonshift=0, scale=1e8): #{{{
    """
    Display the principal components of the tidal stresses defined by the input
    stresscalc object at the points defined by lons and lats, which are equal
    length one dimensional arrays of values, and a time defined by time_t, in
    seconds after periapse.

    The stress vectors are plotted on the map axes defined by basemap_ax.

    Which stresses are plotted depends on the following flags:

       plot_tens  --  if True, plot all tensile stresses.
       plot_comp  --  if True, plot all compressive stresses.
    plot_greater  --  if True, plot the greater (more tensile) principal component
     plot_lesser  --  if True, plot the lesser (more compressive) principal component

    lonshift is a longitudinal displacement added to lons when the stresses are
    calculated, useful in simulating a backrotated shell.

    TODO: Add option of providing an array of scale factors that can be used to
    indicate how important each of the vectors are (e.g. by length of lineament
    segment, or by w_stress in thata region...

    """

    if plot_all is True:
        plot_tens = True
        plot_comp = True
        plot_lesser = True
        plot_greater = True

    calc_phis   = lons
    calc_thetas = (np.pi/2.0)-lats

    Ttt, Tpt, Tpp = stresscalc.tensor(calc_thetas, calc_phis+lonshift, time_t)

    Tau = np.array([[Ttt,Tpt],[Tpt,Tpp]])
    eigensystems = [ np.linalg.eig(Tau[:,:,N]) for N in range(len(Tau[0,0,:])) ]
    evals = np.array([ e[0] for e in eigensystems ])
    evecs = np.array([ e[1] for e in eigensystems ])

    eigval_A = evals[:,0]
    ex_A     = evecs[:,0,1]
    ey_A     = evecs[:,0,0]

    eigval_B = evals[:,1]
    ex_B     = evecs[:,1,1]
    ey_B     = evecs[:,1,0]

    if plot_greater is True or plot_lesser is True:
        mag1 = np.where(eigval_A >  eigval_B, eigval_A, eigval_B)
        ex1  = np.where(eigval_A >  eigval_B, ex_A, ex_B)
        ey1  = np.where(eigval_A >  eigval_B, ey_A, ey_B)

        mag2 = np.where(eigval_A <= eigval_B, eigval_A, eigval_B)
        ex2  = np.where(eigval_A <= eigval_B, ex_A, ex_B)
        ey2  = np.where(eigval_A <= eigval_B, ey_A, ey_B)

    mag1_comp = np.where(mag1 < 0, mag1, 0)
    mag1_tens = np.where(mag1 > 0, mag1, 0)

    mag2_comp = np.where(mag2 < 0, mag2, 0)
    mag2_tens = np.where(mag2 > 0, mag2, 0)

    # First plot all the compressional stresses
    if plot_comp is True:
        if plot_greater is True:
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas),  mag1_comp*ex1,  mag1_comp*ey1, lw=0, width=0.001, scale=scale, color='blue', pivot='tip')
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas), -mag1_comp*ex1, -mag1_comp*ey1, lw=0, width=0.001, scale=scale, color='blue', pivot='tip')

        if plot_lesser is True:
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas),  mag2_comp*ex2,  mag2_comp*ey2, lw=0, width=0.001, scale=scale, color='blue', pivot='tip')
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas), -mag2_comp*ex2, -mag2_comp*ey2, lw=0, width=0.001, scale=scale, color='blue', pivot='tip')

    # Now all the tensional stresses:
    if plot_tens is True:
        if plot_greater is True:
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas),  mag1_tens*ex1,  mag1_tens*ey1, lw=0, width=0.001, scale=scale, color='red', pivot='tail')
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas), -mag1_tens*ex1, -mag1_tens*ey1, lw=0, width=0.001, scale=scale, color='red', pivot='tail')

        if plot_lesser is True:
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas),  mag2_tens*ex2,  mag2_tens*ey2, lw=0, width=0.001, scale=scale, color='red', pivot='tail')
            basemap_ax.quiver(np.degrees(calc_phis), np.degrees(np.pi/2.0-calc_thetas), -mag2_tens*ex2, -mag2_tens*ey2, lw=0, width=0.001, scale=scale, color='red', pivot='tail')

#}}} end vector

def test_vector_grid(plot_all=True, plot_tens=False, plot_comp=False, plot_greater=False, plot_lesser=False, min_lon=0.0, max_lon=np.pi, nlons=13, min_lat=-np.pi/2.0, max_lat=np.pi/2.0, nlats=13, t=0.0): #{{{
    """
    An attempt to exercise the above routines...

    """

    satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = ss.Satellite(open(satfile,'r'))
    NSR = ss.StressCalc([ss.NSR(europa),])

    grid_min_lon = min_lon
    grid_max_lon = max_lon
    grid_min_lat = min_lat
    grid_max_lat = max_lat
    vector_grid_nlons = nlons
    vector_grid_nlats = nlats

    # Lon/Lat locations to plot the principal components:
    vector_grid_lons  = np.linspace(grid_min_lon, grid_max_lon, vector_grid_nlons)
    vector_grid_lats  = np.linspace(grid_min_lat, grid_max_lat, vector_grid_nlats)
    vector_mesh_lons, vector_mesh_lats = np.meshgrid(vector_grid_lons, vector_grid_lats)
                                                     
    vector_mesh_lons = np.ravel(vector_mesh_lons)
    vector_mesh_lats = np.ravel(vector_mesh_lats)

    vector_grid_fig = plt.figure(figsize=(10,10))
    vector_grid_ax  = vector_grid_fig.add_subplot(1,1,1)
    vector_grid_ax.set_title("Regularly gridded vectors")
    vector_grid_basemap = basemap.Basemap(llcrnrlon = np.degrees(grid_min_lon),\
                                          llcrnrlat = np.degrees(grid_min_lat),\
                                          urcrnrlon = np.degrees(grid_max_lon),\
                                          urcrnrlat = np.degrees(grid_max_lat),\
                                          ax = vector_grid_ax)

    vector_points(stresscalc=NSR, lons=vector_mesh_lons, lats=vector_mesh_lats, time_t=t,\
                  plot_greater=plot_greater, plot_lesser=plot_lesser,\
                  plot_comp=plot_comp, plot_tens=plot_tens, plot_all=plot_all,\
                  basemap_ax=vector_grid_basemap)

    vector_grid_basemap.drawmeridians(np.degrees(vector_grid_lons), labels=[1,0,0,1], linewidth=0.5, color='gray')
    vector_grid_basemap.drawparallels(np.degrees(vector_grid_lats), labels=[1,0,0,1], linewidth=0.5, color='gray')
    vector_grid_basemap.drawmapboundary()
#}}}

def test_vector_lin(lons=None, lats=None, t=0.0, plot_all=True, plot_tens=False, plot_comp=False, plot_greater=False, plot_lesser=False):
    
    import nsrhist
    import lineament

    satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = ss.Satellite(open(satfile,'r'))
    NSR = ss.StressCalc([ss.NSR(europa),])

    linzero = nsrhist.load_lins('output/lins/map_nsrfit')[0]
    best_b = linzero.best_fit()[1]
    linzero_mp_lons, linzero_mp_lats = linzero.seg_midpoints()
    min_lat = -np.pi/2.0
    max_lat = np.pi/2.0
    min_lon = 0.0
    max_lon = 2.0*np.pi

    vector_lin_fig = plt.figure()
    vector_lin_ax  = vector_lin_fig.add_subplot(1,1,1)
    vector_lin_ax.set_title("Vectors at arbitrary locations (e.g. along a feature)")

    vector_lin_basemap = basemap.Basemap(llcrnrlon = np.degrees(min_lon),\
                                         llcrnrlat = np.degrees(min_lat),\
                                         urcrnrlon = np.degrees(max_lon),\
                                         urcrnrlat = np.degrees(max_lat),\
                                         ax = vector_lin_ax)

    junk = lineament.plotlinmap([linzero,], map=vector_lin_basemap)
    vector_lin_basemap.drawmeridians(np.degrees(np.linspace(0,2*np.pi,13)), labels=[1,0,0,1], linewidth=0.5, color='gray')
    vector_lin_basemap.drawparallels(np.degrees(np.linspace(-np.pi/2., np.pi/2., 7)), labels=[1,0,0,1], linewidth=0.5, color='gray')
    vector_lin_basemap.drawmapboundary()

    vector_points(stresscalc=NSR, lons=np.mod(linzero_mp_lons,2.0*np.pi), lats=linzero_mp_lats, time_t=0.0, lonshift=best_b,\
                  plot_all=False, plot_tens=True, plot_greater=True, basemap_ax=vector_lin_basemap)
    return linzero_mp_lons, linzero_mp_lats

#}}} end testplot
    
if __name__ == "__main__":
    main(sys.argv)

def old_main(argv=sys.argv): #{{{
    # Parse the command line:
    usage = "usage: %prog [options] infile.nc"
    description = __doc__
    op = OptionParser(usage)
    op.set_defaults( diurnal     = False, \
                     flat        = False, \
                     nsr         = False, \
                     scalar_plot = False, \
                     vector_plot = False)
    
    op.add_option("--outfile", help="""write plot to OUTFILE""", metavar="OUTFILE")
    op.add_option("--diurnal", action="store_true", help="""include viscoelastic diurnal stresses""")
    op.add_option("--nsr", action="store_true", help="""include viscoelastic NSR stresses""")
    op.add_option("--flat", action="store_true", help="""include NSR flattening stresses""")
    op.add_option("--time", help="""time slice (seconds after periapse) to plot.  Cannot be used with --orbit_position.""")
    op.add_option("--time_slice", type="int", help="""time slice (seconds after periapse) to plot.  Cannot be used with --orbit_position.""")
    op.add_option("--orbital_pos", help="""orbital position (degrees after periapse) to plot.  Cannot be used with --time.""")
    op.add_option("--orbital_pos_slice", type="int", help="""orbital position slice to plot.  Cannot be used with --time.""")
    op.add_option("--nsr_period", help="""the NSR period to use in the plot (must be present in the input file)""")
    op.add_option("--nsr_period_slice", type="int", help="""the NSR period slice to use in the plot (must be present in the input file)""")
    op.add_option("--nsr_deg", help="""degrees of accumulated NSR (for flattening calculations) to use in the plot.  Cannot be used with --nsr""")
    op.add_option("--nsr_deg_slice", type="int", help="""the accumulated NSR slice (for flattening calculations).  Cannot be used with --nsr""")
    op.add_option("--scalar_plot", help="""scalar field to use as background (if any), may be either \"tens\" or \"comp\".""", choices=("tens","comp"))
    op.add_option("--vector_plot", action="store_true", help="""display the principal components of the sum of the stress fields being plotted as vectors.""")
    op.add_option("--vector_res", help="""the number of degrees separating individual vectors in the plot.""")
    
    opts, args = op.parse_args(args=argv[1:])
    
    if opts.diurnal and opts.flat:
        raise IncompatibleStressError()
    
    # Open the netCDF file containing the stress calculation:
    nc_in = netCDF3.Dataset(args[0])
    
    # Figure out which time slice we're going to include for Diurnal (if any)
    # Figure out which nsr_period slice we're going to include for NSR (if any)
    # Figure out which nsr_deg slice we're going to include for Flat (if any)
    
    lats  = nc_in.variables['latitude'][:]
    nlats = len(lats)

    lons  = nc_in.variables['longitude'][:]
    nlons = len(lons)

    times = nc_in.variables['time'][:]
    ntimes = len(times)
        
    nsr_periods = nc_in.variables['nsr_period'][:]
    nnsr_periods = len(nsr_periods)

    #nsr_degs = nc_in.variables['nsr_deg'][:]
    #nnsr_degs = len(nsr_degs)

    # Initialize to empty... but the right sizes:
    Tau_Diurnal = numpy.zeros((2,2,ntimes,nlats,nlons))
    Tau_NSR = numpy.zeros((2,2,nnsr_periods,nlats,nlons))
    #Tau_Flat= numpy.zeros((2,2,nnsr_degs,nlats,nlons))

    Tau = numpy.zeros((nlons,nlats,2,2))

    evals   = numpy.zeros((2,nlons,nlats))
    evecs   = numpy.zeros((2,2,nlons,nlats))
    pc_vecs = numpy.zeros((2,2,nlons,nlats))

    tens = numpy.zeros((nlats,nlons))
    comp = numpy.zeros((nlats,nlons))

    if opts.diurnal:
        Ttt_Diurnal = numpy.array(nc_in.variables['Ttt_Diurnal'][:,:,:])
        Tpt_Diurnal = numpy.array(nc_in.variables['Tpt_Diurnal'][:,:,:])
        Tpp_Diurnal = numpy.array(nc_in.variables['Tpp_Diurnal'][:,:,:])
        Tau_Diurnal = numpy.array([ [Ttt_Diurnal, Tpt_Diurnal],
                                    [Tpt_Diurnal, Tpp_Diurnal] ])

    if opts.nsr:
        Ttt_NSR = numpy.array(nc_in.variables['Ttt_NSR'][:,:,:])
        Tpt_NSR = numpy.array(nc_in.variables['Tpt_NSR'][:,:,:])
        Tpp_NSR = numpy.array(nc_in.variables['Tpp_NSR'][:,:,:])
        Tau_NSR = numpy.array([ [Ttt_NSR, Tpt_NSR],
                                [Tpt_NSR, Tpp_NSR] ])

    if opts.flat:
        Ttt_Flat = numpy.array(nc_in.variables['Ttt_Flat'][:,:,:])
        Tpt_Flat = numpy.array(nc_in.variables['Tpt_Flat'][:,:,:])
        Tpp_Flat = numpy.array(nc_in.variables['Tpp_Flat'][:,:,:])
        Tau_Flat = numpy.array([ [Ttt_Flat, Tpt_Flat],
                                 [Tpt_Flat, Tpp_Flat] ])

    # Now that we've got the datacubes in array form, we need to isolate the
    # slices we're actually plotting...

    if opts.nsr and opts.flat:
        # If we've got NSR and Flattening, calculate the delta between them:
        Tau = Tau_NSR[:,:,opts.nsr_period_slice,:,:] - Tau_Flat[:,:,opts.nsr_deg_slice,:,:]
    else:
        # otherwise add the stress fields together:
        Tau = Tau_NSR[:,:,opts.nsr_period_slice,:,:]# + Tau_Diurnal[:,:,opts.time_slice,:,:] + Tau_Flat[:,:,opts.nsr_deg_slice,:,:]
    
    # Now Tau contains the tensor form of whichever snapshot stress field we
    # want to plot, and needs to be turned into the principal components.  This
    # can't be done in a vectorized way because of the way the linalg.eig
    # function works... so we need to iterate over each lat-lon point.  At the
    # same time we have the opportunity to re-arrange the indicies for easier
    # x-y plotting:

    for nlon in numpy.arange(nlons):
        for nlat in numpy.arange(nlats):
            evals[:,nlon,nlat], evecs[:,:,nlon,nlat] = linalg.eig(Tau[:,:,nlat,nlon])
            #pc_vecs[:,:,nlat,nlon] = evals[:,nlon,nlat]*evecs[:,:,nlon,nlat]
            tens[nlat,nlon] = max(evals[:,nlon,nlat])
            comp[nlat,nlon] = min(evals[:,nlon,nlat])

    nx, ny = 18, 12
    i, j = 1, 1
    print np.shape(Tau)
    print np.shape(Tau[i,j])
    print Tau[i,j,nx,ny]
#    print evals[:,nx,ny]
#    print evecs[:,:,nx,ny]
#    print pc_vecs[:,:,nx,ny]

    plt.clf()

    map = basemap.Basemap(llcrnrlon = lons[0],\
                          llcrnrlat = lats[0],\
                          urcrnrlon = lons[-1],\
                          urcrnrlat = lats[-1])

    map.drawmeridians(numpy.linspace(lons[0], lons[-1], 13), labels=[1,1,1,1])
    map.drawparallels(numpy.linspace(lats[0], lats[-1], 7 ), labels=[1,1,1,1])
    map.drawmapboundary()

    if opts.scalar_plot=='tens':
        map.imshow(tens)
    elif opts.scalar_plot=='comp':
        map.imshow(comp)

    if opts.vector_plot:
        # I do the plotting twice, once with the actual values, and ones with
        # their opposites, to get two sets of arrows for each vector (two going
        # away from each other for tension, or two coming towards each other
        # for compression).

        N = 3
        print np.shape(evecs[0,0,::N,::N])
        print np.shape(evecs[0,1,::N,::N])
        print np.shape(evecs[1,0,::N,::N])
        print np.shape(evecs[1,1,::N,::N])
        print np.shape(lats[::N])
        print np.shape(lons[::N])
        map.quiver(lons[::N], lats[::N],  evecs[0,0,::N,::N], -evecs[0,1,::N,::N])#, scale=200000000, width=0.001)
        map.quiver(lons[::N], lats[::N],  evecs[1,0,::N,::N], -evecs[1,1,::N,::N])#, scale=200000000, width=0.001)
        map.quiver(lons[::N], lats[::N], -evecs[0,0,::N,::N],  evecs[0,1,::N,::N])#, scale=200000000, width=0.001)
        map.quiver(lons[::N], lats[::N], -evecs[1,0,::N,::N],  evecs[1,1,::N,::N])#, scale=200000000, width=0.001)
#}}}
