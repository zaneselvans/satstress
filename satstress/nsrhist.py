#!/usr/bin/python
"""
Routines for comparing a set of lineaments to an NSR stress field.

"""

from . import lineament
from . import satstress
from pylab import *
from matplotlib import colors, colorbar
from mpl_toolkits.basemap import Basemap
import random
import os

# These are global variables.  I know, bad.  But this file is eventually going
# to get split into two pieces.  One for general purpose stuff, and another
# which is just a script generating figures for the NSR/TPW paper and my
# thesis.

# Where to put the figures for slurping up?
figdir = os.path.join('output','figs')

# Where to store the giant tangle of pickled lineaments?
lindir = os.path.join('output','lins')

# Format to save the figures in.  If None, don't save:
save_fmt = None

###############################################################################
#  Helper functions for generating, sampling, saving, and loading lineaments  #
###############################################################################

def random_nsrlins(nsr_stresscalc=None, nlins=1000, minlen=0.0, maxlen=1.25): #{{{
    """

    Create nlins lineament objects, resulting from tensile fracture under the
    NSR stress field, having lengths between minlen and maxlen (radians of
    arc), with their locations and (maximum) lengths randomly distributed.

    Because some lineaments will self-terminate upon reaching compressive
    zones, this will naturally result in a bias towards shorter lineaments.
    This bias can be compensated for later by subsampling the synthetic
    lineament population according to whatever criteria the user desires.

    """

    if nsr_stresscalc is None:
        satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
        europa = satstress.Satellite(open(satfile,'r'))
        nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    linlist = []

    N=0
    max_lengths = minlen+(rand(2*nlins)*(maxlen-minlen))
    init_lons, init_lats = lineament.random_lonlatpoints(2*nlins)
    while len(linlist) < nlins:
        seg_len=min(0.01, max_lengths[N]/10.0)
        newlin = lineament.lingen_nsr(nsr_stresscalc, init_lon=init_lons[N], init_lat=init_lats[N], max_length=max_lengths[N], prop_dir='both', seg_len=seg_len)
        if newlin.length > minlen:
            linlist.append(newlin)
        N = N+1

    return linlist
#}}}

def random_gclins(nlins=1000, minlen=0.0, maxlen=1.25): #{{{
    """
    Create nlins lineament objects, whose paths approximate great circle
    segments on the surface of the satellite, with lengths ranging from
    minlen to maxlen (radians of arc), randomly distributed in length,
    orientation, and location on the surface.

    """

    # - pick a random endpoint v1 on the surface of the sphere
    init_lons, init_lats = lineament.random_lonlatpoints(nlins)
    # - pick a random azimuth, az
    azimuths = 2*pi*rand(nlins)
    # - pick a random length L between minlen and maxlen
    linlengths = (minlen+(rand(nlins)*(maxlen-minlen)))
    # - calculate the location of endpoint v2, L radians away
    #   from v1, along the initial heading az, from v1.
    fin_lons, fin_lats = lineament.spherical_reckon(init_lons, init_lats, azimuths, linlengths)
    # - use lingen_greatcircle() to calculate intermediary vertices.
    return([ lineament.lingen_greatcircle(init_lons[N], init_lats[N], fin_lons[N], fin_lats[N]) for N in range(nlins) ])

#}}}

def linsample(proto=None, pool=None, nbins=20, fraction=1.0): #{{{
    """
    Select a set of lineaments from a pool with a similar length distribution
    to that of a prototype set.

    nbins controls the number of discrete, evenly spaced length bins which are
    to be approximated in replicating the length distribution.  The lower edge
    of the shortest bin is the length of the shortest feature in the prototype
    set, and similarly the high edge of the longest bin is the length of the
    longest feature in the prototype set.

    """

    # define the length bins to replicate:
    proto_lengths = array([ lin.length for lin in proto ])
    proto = array(proto)
    pool_lengths  = array([ lin.length for lin in pool  ])
    pool = array(pool)

    bins = linspace(min(proto_lengths)-0.00001, max(proto_lengths)+0.00001, num=nbins+1)
    sample_lins = []

    for i in range(nbins):
        # Organize the prototypes and pool into bins by length
        proto_binned_lins    = proto[where( logical_and(proto_lengths >= bins[i], proto_lengths < bins[i+1]) )]
        proto_binned_lengths = array([lin.length for lin in proto_binned_lins])

        pool_binned_lins    = pool[where(  logical_and(pool_lengths >= bins[i], pool_lengths < bins[i+1]) )]
        pool_binned_lengths = array([lin.length for lin in pool_binned_lins])

        if len(pool_binned_lins) > 0:
            sample_binned_lins = pool_binned_lins[ np.random.randint(0, high=len(pool_binned_lins), size=int(fraction*len(proto_binned_lins))) ]
            sample_binned_lengths = array([lin.length for lin in sample_binned_lins])
            sample_lins.append(sample_binned_lins)

        #print("%3.d: %7.3g %7.3g " % (i+1, np.mean(proto_binned_lengths), np.mean(sample_binned_lengths)) )

    return(concatenate(sample_lins))
#}}}

def save_lins(lins, name="generic_lins", outdir=lindir): #{{{
    """
    A shortcut for saving a new set of analyzed lineaments.
    """
    from time import strftime
    import pickle
    import os

    outfile="%s_%s.pkl" % (name, strftime('%Y%m%d%H%M%S'))
    outpath = os.path.join(outdir,outfile)
    linkpath = os.path.join(outdir,name)
    pickle.dump(lins, open(outpath,'w'))

    try:
        os.unlink(linkpath)
    except OSError:
        pass

    os.symlink(outfile, linkpath)

#}}} end save_lins()

def load_lins(linpath): #{{{
    """
    Load a list of lineaments from either a pickle file, or a directory filled
    with pickle files (as is produced by calc_fits)

    """
    import pickle

    linlist = []
    if os.path.isdir(linpath):
        lindir = os.listdir(linpath)
        for linfilename in lindir:
            linfilepartname = os.path.join(linpath, linfilename)
            print "Loading", linfilepartname
            newpart = pickle.load(open(linfilepartname))
            print "  found %d lineaments" % (len(newpart,))
            linlist += newpart
    elif os.path.isfile(linpath):
        linlist = pickle.load(open(linpath))
    else:
        raise os.error("Path: \'%s\' is not a file or directory" % (linpath,) )

    return linlist
#}}}

def reload_nsrfits(update=False): #{{{
    """
    Loads the lineaments which I most often end up fooling around with and
    returns them as a list, ordered as follows:

    maplins    Europa's lineaments, as mapped
    synthlins  A perfect synthetic dataset
    crazylins  The mapped lineaments, with random locations/orientations
    gclins     A collection of great circle segments
    tpwlins    The mapped lineaments, transformed to PNP 80E 10N

    maplins, synthlins, crazylins, gclins, tpwlins = nsrhist.reload_nsrfits()

    if update is True, then these datasets are read in, updated to reflect the
    most recent version of the Lineament object, saved back to disk, and then
    returned.

    """

    maplins   = load_lins(os.path.join(lindir,'map_nsrfit'))
    synthlins = load_lins(os.path.join(lindir,'synth_nsrfit'))
    crazylins = load_lins(os.path.join(lindir,'crazy_nsrfit'))
    gclins    = load_lins(os.path.join(lindir,'gc_nsrfit'))
    tpwlins   = load_lins(os.path.join(lindir,'tpw_nsrfit'))

    if update is True:
        maplins   = lineament.update_lins(maplins)
        synthlins = lineament.update_lins(synthlins)
        crazylins = lineament.update_lins(crazylins)
        gclins    = lineament.update_lins(gclins)
        tpwlins   = lineament.update_lins(tpwlins)

        save_lins(synthlins, name='synth_nsrfit')
        save_lins(maplins,   name='map_nsrfit')
        save_lins(crazylins, name='crazy_nsrfit')
        save_lins(gclins,    name='gc_nsrfit')
        save_lins(tpwlins,   name='tpw_nsrfit')

    return(maplins, synthlins, crazylins, gclins, tpwlins)
#}}}

def gc(lin): #{{{
    gc_lin = lineament.Lineament(lons=lin.lons, lats=lin.lats)

    pole1_lon, pole1_lat = gc_lin.bfgc_pole()
    pole2_lon, pole2_lat = lineament.spherical_reckon(pole1_lon, pole1_lat, 0.0, pi)
    # Generate a BFGC lineament:
    gc_dop = gc_lin.doppelgen_gcseg()

    the_fig = figure()
    ax = the_fig.add_subplot(1,1,1)

    # Plot the lineament, and its BFGC doppelganger:
    ax.plot(degrees(mod(gc_lin.lons,2*pi)), degrees(gc_lin.lats), color='black', linewidth=2)
    ax.plot(degrees(mod(gc_dop.lons,2*pi)), degrees(gc_dop.lats), color='black', linewidth=2, alpha=0.5)

    # Plot the implied pole and the other points on the great circle...
    ax.scatter(degrees(mod([pole1_lon,pole2_lon],2*pi)), degrees([pole1_lat,pole2_lat]), color='black', s=50, marker='x', linewidth=2)
    far_lons, far_lats = lineament.spherical_reckon(mod(pole1_lon,2*pi), pole1_lat, linspace(0,2*pi,180), pi/2)
    ax.scatter(degrees(mod(far_lons,2*pi)), degrees(far_lats), color='black', s=3, alpha=0.3)

    print("D_bar = %f" % (gc_lin.mhd(gc_dop),))

    ax.set_xlim(0,360)
    ax.set_xticks(linspace(0,360,13))
    ax.set_ylim(-90,90)
    ax.set_yticks(linspace(-90,90,7))
    ax.set_aspect('equal')
#}}}

###############################################################################
#    Functions for calculating fits, and displaying the results in general    #
###############################################################################

def calcfits(satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
             linfile="input/GlobalLineaments",
             nb=180, pnp_lon=radians(80), pnp_lat=radians(10), nlins=0): #{{{
    """
    Starting from scratch, read in the mapped features, calculate their fits,
    and use them to generate TPW features with the given pole, and to subsample
    the two sets of synthetic features (great circle segments and perfect NSR
    lineaments)

    """
    print("Initializing satellite")
    europa = satstress.Satellite(open(satfile,'r'))
    NSR = satstress.StressCalc([satstress.NSR(europa),])

    print("Loading mapped features")
    maplins = lineament.shp2lins(linfile, stresscalc=NSR)

    if nlins > 0:
        maplins=maplins[:nlins]
    print("Transforming to pre-TPW coordinates")
    tpwlins = [ lin.poleshift(pnp_lon=pnp_lon, pnp_lat=pnp_lat) for lin in maplins ]

    print("Generating synthetic lineament pools")
    protolinlens = array([lin.length for lin in maplins])
    print("  great circle segments")
    gclins = linsample(proto=maplins, pool=random_gclins(nlins=10*len(maplins), minlen=protolinlens.min(), maxlen=protolinlens.max()))
    print("  model NSR lineaments")
    synthlins = linsample(proto=maplins, pool=random_nsrlins(nsr_stresscalc=NSR, nlins=10*len(maplins), minlen=protolinlens.min(), maxlen=protolinlens.max()))

    lins_list = (maplins, tpwlins, gclins, synthlins)
    labels = ("map_nsrfit", "tpw_nsrfit", "gc_nsrfit", "synth_nsrfit")

    # now calculate the fits...
    for lins,label in zip(lins_list, labels):
        print("Calculating fits for %s..." % (label,) )
        for lin,N in zip(lins,range(len(lins))):
            if (mod(N+1,60)==0):
                print("    N=%d" % (N+1,) )
            lin.calc_nsrfits(nb=nb, stresscalc=NSR)
        # only save if we're doing the whole dataset
        if nlins == 0:
            print("Saving %s" % (label,) )
            save_lins(lins, name=label)

    return(maplins, tpwlins, gclins, synthlins)
#}}}

def tpw_polesearch(satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                   linfile="input/GlobalLineaments",\
                   outdir="output/tpw_polesearch",\
                   np=100, nb=36, nlins=0): #{{{
    """
    Transform the lineaments read in from linfile to np different possible
    paleopoles, randomly chosen and evenly distributed over the surface of the
    sphere, fit them to the NSR stress field at nb different values of
    backrotation.  Save the results for later viewing.

    For testing purposes, set nlins to the number of features to fit.  Fits and
    transformed features will not be saved.

    """
    import os.path

    print("Initializing satellite")
    europa = satstress.Satellite(open(satfile,'r'))
    NSR = satstress.StressCalc([satstress.NSR(europa),])

    print("Loading mapped features")
    maplins = lineament.shp2lins(linfile, stresscalc=NSR)

    if nlins > 0:
        maplins = maplins[:nlins]

    # Generate the list of paleopole positions
    #pnp_lons, pnp_lats = lineament.random_lonlatpoints(np)

    # limit the search to one hemisphere because of periodicity:
    #pnp_lons = mod(pnp_lons, pi)

    # Focus on the high latitude region:
    pnp_lons, pnp_lats = fibonacci_sphere(2000)
    pnp_hilats = pnp_lats[where(degrees(pnp_lats) > 60)]
    pnp_hilons = pnp_lons[where(degrees(pnp_lats) > 60)]
    pnp_lats = pnp_hilats
    pnp_lons = pnp_hilons

    tpwlins_list = []
    for pnp_lon, pnp_lat, N in zip(pnp_lons, pnp_lats, arange(len(pnp_lons))):
        print("Fitting paleopole %d / %d (lon=%f, lat=%f)" % (N+1,len(pnp_lons), degrees(pnp_lon), degrees(pnp_lat)) )
        tpwlins = [ lin.poleshift(pnp_lon=pnp_lon, pnp_lat=pnp_lat) for lin in maplins ]

        label = "tpw_lon%.6f_lat%.6f_nsrfit" % (degrees(pnp_lon), degrees(pnp_lat))
        devnull = [ lin.calc_nsrfits(nb=nb, stresscalc=NSR) for lin in tpwlins ]
        tpwlins_list.append(tpwlins)

        # only save the features if we're doing the whole dataset:
        if nlins == 0:
            print("    Saving %s" % (os.path.join(outdir,label),) )
            save_lins(tpwlins, name=label, outdir=outdir)
    
    return(pnp_lons, pnp_lats, tpwlins_list)
#}}}

def load_tpw_results(outdir="output/tpw_polesearch", fast=True): #{{{
    """
    Look in the specified output directory for sets of lineaments which have
    been transformed to a paleopole location, and fit to the NSR field.  Parse
    their filenames to glean the paleopole information, and return the list of
    paleopoles and the amplitudes of the activity histories of the lineament
    datasets corresponding to them.

    """
    from re import match
    from os import listdir
    import os.path
    import pickle

    tpwlins_list = []
    pnp_lons = []
    pnp_lats = []

    if fast is True:
        # we've already read these things in once.... don't need to do it again.
        tpw_polesearch_results = pickle.load(open(outdir+'.pkl'))
        pnp_lons = tpw_polesearch_results['pnp_lon']
        pnp_lats = tpw_polesearch_results['pnp_lat']
        acthist_amps = tpw_polesearch_results['acthist_amp']

    else:
        # for each filename in the output directory
        for fn in os.listdir(outdir):
            # see if it matches our output filename format
            filematch = match(r'tpw_lon(?P<pnp_lon>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)_lat(?P<pnp_lat>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)_nsrfit$', fn)
            # if it does, extract information from the filename, load the lineaments:
            if filematch is not None:
                pnp_lon, pnp_lat = float(filematch.group('pnp_lon')), float(filematch.group('pnp_lat'))
                print("Found TPW dataset with pnp_lon=%.6f, pnp_lat=%.6f" % (pnp_lon, pnp_lat))
                tpwlins = load_lins(os.path.join(outdir,fn))
    
                # load the fits for the dataset:
                bs = array([ lin.bs for lin in tpwlins])
                nsrdbars = array([lin.nsrdbars for lin in tpwlins])
                nsrstresswts = array([lin.nsrstresswts for lin in tpwlins])
    
                # If it looks like we have a full fit, append it to the output:
                if shape(bs)==shape(nsrdbars)==shape(nsrstresswts) and len(bs[0]) > 0:
                    pnp_lons.append(radians(pnp_lon))
                    pnp_lats.append(radians(pnp_lat))
                    tpwlins_list.append(tpwlins)

        acthist_amps = array([ acthist_amplitude(lins) for lins in tpwlins_list ])

        # Now that we've gone to all the trouble to load this thing in, let's
        # save it for next time:
        dtype = [('pnp_lon',float),('pnp_lat',float),('acthist_amp',float)]
        pickle.dump(array([ (pnp_lons[n], pnp_lats[n], acthist_amps[n]) for n in range(len(pnp_lons))], dtype=dtype), open(outdir+'.pkl','w'))

    return(array(pnp_lons), array(pnp_lats), acthist_amps)
#}}}

def tpw_poledetails(N=2000): #{{{
    the_fig = figure(figsize=(10,10))
    np_ax = the_fig.add_subplot(1,1,1)
    np_map = Basemap(projection='npstere', boundinglat=40, lon_0=90, ax=np_ax)
    np_map.drawparallels(arange(-60.,61.,30.), labels=[1,1,0,0])
    np_map.drawmeridians(arange(-180.,181.,30.), labels=[1,1,0,1])

    pnp_lons, pnp_lats = fibonacci_sphere(N)
    np_x1_data, np_y1_data = np_map(degrees(pnp_lons), degrees(pnp_lats))
    np_map.scatter(np_x1_data, np_y1_data, s=5, alpha=0.5, color='black')

    pnp_hilats = pnp_lats[where(degrees(pnp_lats) > 60)]
    pnp_hilons = pnp_lons[where(degrees(pnp_lats) > 60)]
    print(" found %d points north of 60" % (len(pnp_hilats),) )
    np_x2_data, np_y2_data = np_map(degrees(pnp_hilons), degrees(pnp_hilats))
    np_map.scatter(np_x2_data, np_y2_data, s=5, color='red')
#}}}

def tpw_polefits(pnp_lons=None, pnp_lats=None, tpwlins_list=None, maplins=None, gclins=None, ncols=50, contour_cmap=cm.jet, fast=True): #{{{
    """
    Given a set of paleopole longitudes and latitudes defined by pnp_lons and
    pnp_lats, and a list of lineament datasets tpwlins_list corresponding to
    mapped features transformed to those paleopole locations, and fit to the
    NSR stress field there, create a 2-D plot showing how the amplitude of the
    activity history varies with paleopole location.

    """

    if pnp_lons is None and pnp_lats is None and tpwlins_list is None:
        if fast is True:
            pnp_lons, pnp_lats, acthist_amps = load_tpw_results(fast=True)
        else:
            pnp_lons, pnp_lats, tpwlins_list = load_tpw_results(fast=False)
            acthist_amps = array([ acthist_amplitude(lins) for lins in tpwlins_list ])

    if maplins is not None:
        map_acthist_amp = acthist_amplitude(maplins)
    else:
        map_acthist_amp = 0.0

    if gclins is not None:
        gc_acthist_amp = acthist_amplitude(gclins)
    else:
        gc_acthist_amp = 0.0

    levels = linspace(acthist_amps.min(), acthist_amps.max(), ncols)

    # Add another copy of the data, flipped to represent the odd-symmetry:
    pnp_lons = concatenate([pnp_lons, pnp_lons-pi])
    pnp_lats = concatenate([pnp_lats, -pnp_lats])
    acthist_amps = tile(acthist_amps,2)

    # Create a three-panel plot, one for the equatorial to mid-latitudes in a
    # simple cylindrical projection, and one each for the poles, in polar
    # stereographic underneath the other one.
    the_fig = figure(figsize=(10,10))
    lowlat_ax = the_fig.add_subplot(2,1,1)
    np_ax = the_fig.add_subplot(2,2,3)
    sp_ax = the_fig.add_subplot(2,2,4)
    lowlat_map = Basemap(llcrnrlon=-180,llcrnrlat=-60,urcrnrlon=180,urcrnrlat=60, ax=lowlat_ax)
    np_map = Basemap(projection='npstere',boundinglat=40,lon_0=90, ax=np_ax)
    sp_map = Basemap(projection='spstere',boundinglat=-40,lon_0=270, ax=sp_ax)

    # Ideally, all this interpolation would be done in spherical coordinates,
    # but that's not built-in.  Getting it right in the poles is important, and
    # can be done well by doing the interpolation in the map-projected space.
    # For the low latitude regions, we can get away with just setting up a
    # pseudo-periodic environment by copying the data a few times on either
    # side of the region we actually care about.

    # North Pole:
    np_x_data, np_y_data = np_map(degrees(pnp_lons), degrees(pnp_lats))
    np_x_grid = linspace(np_map.xmin,np_map.xmax,100)
    np_y_grid = linspace(np_map.ymin,np_map.ymax,100)
    np_pnp_amps = griddata(np_x_data, np_y_data, acthist_amps, np_x_grid, np_y_grid)
    np_x_gridmesh, np_y_gridmesh = meshgrid(np_x_grid,np_y_grid)
    np_map.contourf(np_x_gridmesh, np_y_gridmesh, np_pnp_amps, levels, cmap=cm.get_cmap('jet', ncols-1), linewidth=0)
    #np_map.contour(np_x_gridmesh, np_y_gridmesh, np_pnp_amps, [map_acthist_amp,], colors='black', linewidths=[2,])
    #np_map.contour(np_x_gridmesh, np_y_gridmesh, np_pnp_amps, [gc_acthist_amp,], colors='black', linewidths=[2,], linestyle='--')
    np_map.scatter(np_x_data, np_y_data, s=5, alpha=0.5, color='black', linewidth=0)
    np_map.drawparallels(arange(-60.,61.,30.), labels=[1,1,0,0])
    np_map.drawmeridians(arange(-180.,181.,30.), labels=[1,1,0,1])
    np_ax.set_title("North Polar Stereographic")

    # South Pole:
    sp_x_data, sp_y_data = sp_map(degrees(pnp_lons), degrees(pnp_lats))
    sp_x_grid = linspace(sp_map.xmin,sp_map.xmax,100)
    sp_y_grid = linspace(sp_map.ymin,sp_map.ymax,100)
    sp_pnp_amps = griddata(sp_x_data, sp_y_data, acthist_amps, sp_x_grid, sp_y_grid)
    sp_x_gridmesh, sp_y_gridmesh = meshgrid(sp_x_grid,sp_y_grid)
    sp_map.contourf(sp_x_gridmesh, sp_y_gridmesh, sp_pnp_amps, levels, cmap=cm.get_cmap('jet', ncols-1))
    sp_map.contour(sp_x_gridmesh, sp_y_gridmesh, sp_pnp_amps, [map_acthist_amp,], colors='black', linewidths=[2,])
    sp_map.contour(sp_x_gridmesh, sp_y_gridmesh, sp_pnp_amps, [gc_acthist_amp,], colors='black', linewidths=[2,], linestyle='--')
    sp_map.scatter(sp_x_data, sp_y_data, s=5, alpha=0.5, color='black', linewidth=0)
    sp_map.drawparallels(arange(-60.,61.,30.), labels=[1,1,0,0])
    sp_map.drawmeridians(arange(-180.,181.,30.), labels=[1,1,0,1])
    sp_ax.set_title("South Polar Stereographic")

    # Low Latitudes:
    lowlat_x_data = pnp_lons
    lowlat_y_data = pnp_lats
    lowlat_acthist_amps = acthist_amps
    lowlat_x_grid = linspace(-pi,pi,360)
    lowlat_y_grid = linspace(-pi/2,pi/2,180)
    lowlat_pnp_amps = griddata(lowlat_x_data, lowlat_y_data, lowlat_acthist_amps, lowlat_x_grid, lowlat_y_grid)
    lowlat_x_gridmesh, lowlat_y_gridmesh = meshgrid(lowlat_x_grid,lowlat_y_grid)
    lowlat_map.contourf(degrees(lowlat_x_gridmesh), degrees(lowlat_y_gridmesh), lowlat_pnp_amps, levels, cmap=cm.get_cmap('jet', ncols-1))
    lowlat_map.contour(degrees(lowlat_x_gridmesh),  degrees(lowlat_y_gridmesh), lowlat_pnp_amps, [map_acthist_amp,], colors='black', linewidths=[2,])
    lowlat_map.contour(degrees(lowlat_x_gridmesh),  degrees(lowlat_y_gridmesh), lowlat_pnp_amps, [gc_acthist_amp,], colors='black', linewidths=[2,], linestyle='--')
    lowlat_map.scatter(degrees(lowlat_x_data), degrees(lowlat_y_data), s=5, alpha=0.5, color='black', linewidth=0)
    lowlat_map.drawparallels(arange(-60.,61.,30.), labels=[1,1,0,1])
    lowlat_map.drawmeridians(arange(-180.,181.,30.), labels=[1,1,0,1])
    lowlat_ax.set_title("Equatorial and Mid-Latitudes")

    cb_ax,kw = colorbar.make_axes(lowlat_ax, orientation='horizontal')
    colorbar.ColorbarBase(cb_ax, cmap=cm.get_cmap('jet', ncols-1), norm=colors.Normalize(vmin=levels[0],vmax=levels[-1]), orientation='horizontal', format='%.4f')
    cb_ax.set_xlabel("activity amplitude")

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'PaleopoleActAmps.'+save_fmt))

#}}}

def aggregate_bestfit(lins, delta_max=radians(45), dbar_max=0.125): #{{{
    """
    An overall metric of how well a given set of features can be forced to fit
    the NSR stresses, if each feature is allowed to translate freely in
    longitude.

    """

    linlengths  = array([ lin.length for lin in lins ])
    linbestfits = array([ max(lin.nsrfits(dbar_max=dbar_max)) for lin in lins ])
    return(sum(linlengths*linbestfits)/sum(linlengths))

#}}} end aggregate_bestfit

def show_extreme_doppels(lin, dbar_max=0.125): #{{{
    mp_lon, mp_lat = lin.bfgcseg_midpoint()
    seg_len = min(0.01, lin.length/10.0)
    best_fit = min(lin.nsrdbars)
    best_b = lin.bs[where(fabs(lin.nsrdbars - best_fit) < 1e-9)[0][0]]
    best_doppel = lineament.lingen_nsr(stresscalc=lin.stresscalc, init_lon=mp_lon+best_b, init_lat=mp_lat, prop_dir="both", max_length=lin.length, seg_len=seg_len)
    best_doppel.lons -= best_b

    worst_fit = max(lin.nsrdbars)
    worst_b = lin.bs[where(fabs(lin.nsrdbars - worst_fit) < 1e-9)[0][0]]
    worst_doppel = lineament.lingen_nsr(stresscalc=lin.stresscalc, init_lon=mp_lon+worst_b, init_lat=mp_lat, prop_dir="both", max_length=lin.length, seg_len=seg_len)
    worst_doppel.lons -= worst_b

    lines, map = lineament.plotlinmap([lin,], color='black', linewidth=2.0)
    map.scatter(degrees(lin.lons), degrees(lin.lats), color='black')

    lineament.plotlinmap([best_doppel,], color='green', map=map, linewidth=2)
    #map.scatter(degrees(best_doppel.lons), degrees(best_doppel.lats), color='green')
    #map.scatter(degrees(best_doppel.lons+2*pi), degrees(best_doppel.lats), color='green')
    #map.scatter(degrees(best_doppel.lons-2*pi), degrees(best_doppel.lats), color='green')
    gca().annotate('D=%.3g' % (best_fit,), xy=(degrees(best_doppel.lons[0]),      degrees(best_doppel.lats[0])) )
    gca().annotate('D=%.3g' % (best_fit,), xy=(degrees(best_doppel.lons[0]+2*pi), degrees(best_doppel.lats[0])) )
    gca().annotate('D=%.3g' % (best_fit,), xy=(degrees(best_doppel.lons[0]-2*pi), degrees(best_doppel.lats[0])) )

    lineament.plotlinmap([worst_doppel,], color='red', map=map, linewidth=2)
    #map.scatter(degrees(worst_doppel.lons), degrees(worst_doppel.lats), color='red')
    #map.scatter(degrees(worst_doppel.lons+2*pi), degrees(worst_doppel.lats), color='red')
    #map.scatter(degrees(worst_doppel.lons-2*pi), degrees(worst_doppel.lats), color='red')
    gca().annotate('D=%.3g' % (worst_fit,), xy=(degrees(worst_doppel.lons[0]),      degrees(worst_doppel.lats[0])) )
    gca().annotate('D=%.3g' % (worst_fit,), xy=(degrees(worst_doppel.lons[0]+2*pi), degrees(worst_doppel.lats[0])) )
    gca().annotate('D=%.3g' % (worst_fit,), xy=(degrees(worst_doppel.lons[0]-2*pi), degrees(worst_doppel.lats[0])) )

    return(best_doppel, worst_doppel)
#}}}

def fibonacci_sphere(N, halfsphere=False): #{{{
    """
    Generate N points on the surface of the sphere, fairly evenly spaced from
    one another.  If halfsphere is True, return only points with longitudes
    between 0 and 180 degrees.

    """

    # if we're only doing half the sphere, need to double N first to end up
    # with the right number of points
    if halfsphere is True:
        N = 2*N

    inc = pi * (3 - sqrt(5))
    off = 2. / N
    k = arange(0,N)
    y = k*off - 1. + 0.5*off
    r = sqrt(1 - y*y)
    phi = k * inc
    x = cos(phi)*r
    z = sin(phi)*r
    theta = arctan2(sqrt(x**2+y**2),z)
    lons = arctan2(y,x)
    lats = pi/2-theta

    # Only retain those points with longtiudes between 0 and 180 degrees.
    if halfsphere is True:
        lats = lats[where(mod(lons,2*pi) < pi)]
        lons = lons[where(mod(lons,2*pi) < pi)]

    return lons, lats
#}}}

def fitcurve(lin, color='black', ax=None, dbar_max=0.125, show_dbar=True, show_w=True, show_fnsr=True, use_stress=True): #{{{
    """
    Plot delta and dbar as a function of b for lin, and show the value of the
    weighting function that results.

    """

    if ax is None:
        the_fig = figure(figsize=(9,6))
        dbar_ax = the_fig.add_subplot(1,1,1)
    else:
        dbar_ax = ax

    fit_ax  = twinx(ax=dbar_ax)

    # Plot the mapped lineament's fit curve:
    dbar_ax.set_xlabel("westward translation (b)")
    if show_dbar is True:
        dbar_ax.plot(degrees(lin.bs), lin.nsrdbars, ls='-', linewidth=2, color=color)
        dbar_ax.grid()
        dbar_ax.axis([0,179,0,0.3])
        dbar_ax.set_yticks(linspace(0,0.3,5,endpoint=False))

    if show_fnsr is True:
        fit_ax.fill_between(degrees(lin.bs), lin.nsrfits(dbar_max=dbar_max, use_stress=use_stress), 0, color=color, alpha=0.3)

    if show_w is True:
        fit_ax.plot(degrees(lin.bs), lin.nsrstresswts, ls='--', linewidth=2, color=color)

    if show_fnsr is True or show_w is True:
        fit_ax.set_yticks(linspace(0,1.5,5,endpoint=False))
        fit_ax.grid()
        fit_ax.axis([0,179,0,1.5])

    dbar_ax.invert_xaxis()
#}}}

def make_crazy(lins, tpw=True, spin=True): #{{{
    """
    Take a list of lineaments and randomize their locations and orientations on
    the sphere by applying random TPW and NSR shifts to them.

    If tpw is True, re-orient and move the features on the surface arbitrarily.
    This results in the loss of any previously calculated fit information.

    if spin is True, simply shift them in longitude.  If only this
    transformation is performed, the fit information is preserved.

    """

    newlins = []
    if spin is True:
        rand_bs = 2*pi*np.random.random(len(lins))
        for lin, b in zip(lins, rand_bs):
            newlins.append(lin.lonshift(b))
        lins = newlins

    if tpw is True:
        newlins = []
        rand_lons, rand_lats = lineament.random_lonlatpoints(len(lins))
        for lin, pnp_lon, pnp_lat in zip(lins, rand_lons, rand_lats):
            newlins.append(lin.poleshift(pnp_lon=pnp_lon, pnp_lat=pnp_lat))

    return(newlins)
#}}}

def activity_history(lins_list, the_fig=None, labels=[], colors=[], alphas=[], lin_cm=cm.jet, norm_by_all=False, dbar_max=0.125, verbose=True, titlestr="Apparent Lineament Activity History", outfile=None): #{{{
    """
    Plots apparent activity histories for one or several sets of lineaments,
    allowing visual comparison.

    lins_list is a list of lists of Lineament objects, assumed to have their
    fits already calculated.

    the_fig is the matplotlib figure in which the activity history should be
    drawn.  If None, a new plot is created.

    labels and colors are lists of strings and matplotlib color specifications
    that will be used to in legend creation and drawing (colors are
    automatically generated from the lin_cm colormap if colors is empty)

    if norm_by_sum is True, then the lineaments are treated as being subsets of
    one larger dataset, and their apparent contributions to the activity
    history are normalized by the cumulative length of the dataset, otherwise,
    they are assumed to be separate datasets, and each is normalized by its own
    cumulative length.

    """

    # Set hist_ax to the current axes, if none was supplied.
    if the_fig is None:
        the_fig = figure(figsize=(12,8))

    if len(the_fig.axes) == 0:
        hist_ax = the_fig.add_subplot(1,1,1)
    else:
        hist_ax = the_fig.axes[0]

    # Make sure that if we got labels and colors, there's the right number:
    if (len(labels) == 0):
        labels = [ 'Series '+str(N) for N in arange(len(lins_list))+1 ]
    try:
        assert(len(lins_list) == len(labels))
    except AssertionError:
        print "len(lins_list) = %g != len(labels) = %g" % (len(lins_list), len(labels))

    if (len(colors) == 0):
        colors = lin_cm(linspace(0,1,len(lins_list)))
    try:
        assert(len(lins_list) == len(colors))
    except AssertionError:
        print "len(lins_list) = %g != len(colors) = %g" % (len(lins_list), len(colors))

    if (len(alphas) == 0):
        alphas = np.ones(len(lins_list))
    try:
        assert(len(lins_list) == len(alphas))
    except AssertionError:
        print "len(lins_list) = %g != len(alphas) = %g" % (len(lins_list), len(alphas))

    # Depending on whether we're segmenting a single dataset (say, by lineament
    # length) or trying to compare different datasets, we may want to norm each
    # set of lineaments by its own length, or the total length of all the
    # lineaments.
    norm_lengths = [ sum([lin.length for lin in lins]) for lins in lins_list ]
    if norm_by_all is True:
        norm_lengths = [sum(norm_lengths),] * len(lins_list)

    for lins,label,color,alpha,norm_length in zip(lins_list,labels,colors,alphas,norm_lengths):
        if verbose is True:
            print("Calculating activity histories for %s" % (label,))
        acthist = calc_acthist(lins, dbar_max=dbar_max, norm_length=norm_length)
        hist_ax.plot(degrees(lins[0].bs), acthist, label=label, lw=3, c=color, alpha=alpha)

    hist_ax.set_ylabel("H(b)")
    hist_ax.set_xlabel("westward translation b")
    hist_ax.grid(True)
    hist_ax.set_xticks(linspace(0,170,18))
    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    hist_ax.xaxis.set_major_formatter(degree_formatter)
    hist_ax.set_xlim(179,0)
    hist_ax.legend(loc="upper left")
    if norm_by_all is True:
        titlestr += ' (Normalized)'
    hist_ax.set_title(titlestr)

    if save_fmt is not None:
        the_fig.savefig(outfile+'.'+save_fmt)
#}}}

def calc_acthist(lins, dbar_max=0.125, norm_length=None): #{{{
    if norm_length is None:
        norm_length = sum([lin.length for lin in lins])
    return(array([ lin.length*lin.nsrfits(dbar_max=dbar_max) for lin in lins ]).sum(axis=0)/norm_length)

def acthist_amplitude(lins, dbar_max=0.125, norm_length=None):
    acthist = calc_acthist(lins, dbar_max=dbar_max, norm_length=norm_length)
    return(max(acthist) - min(acthist))
#}}}

def jackstraws(N=10, nb=36, spin=True, tpw=False): #{{{
    """
    Take the mapped set of lineaments, and scatter them over the surface of
    Europa.  Then calculate their fit curve and return its amplitude.

    """

    print("Reading in mapped features")
    maplins = load_lins(os.path.join(lindir,'map_nsrfit'))
    NSR = maplins[0].stresscalc

    amps_out = []
    for n in range(N):
        if mod(n+1,100) == 0:
            print("Tossing jackstraws %d / %d" % (n+1, N,))
        crazylins = make_crazy(maplins, tpw=tpw, spin=spin)

        if tpw is True:
            print("Fitting jackstraws %d / %d" % (n+1, N,))
            [ lin.calc_nsrfits(nb=nb, stresscalc=NSR) for lin in crazylins ]

        amps_out.append(acthist_amplitude(crazylins))

    return(amps_out)
#}}}

def nsrfitquality(lins, use_stress=True): #{{{
    """
    Takes a list of lineaments which have had their NSR fits calculated and
    returns a the sum of the lineament lengths multiplied by their greatest
    value of f_{nsr}(b), divided by the overall length of the dataset.

    """

    return(sum([lin.length*lin.best_fit(use_stress=use_stress)[0] for lin in lins])/sum([lin.length for lin in lins]))

#}}}

###############################################################################
# Testing functions...
###############################################################################

def test_KDTree(linlib=None, libsize=30, N=10, d_max=1.0): #{{{
    """
    A routine for testing how well this whole idea of finding an approximate
    doppelganger in the library works.

    """
    from scipy.spatial import KDTree

    radius = 1.0
    # generate a library if we didn't get one passed in:
    if linlib is None:
        linlib = lineament.lingen_nsr_library(nlats=libsize)

    # make a list of all the longitude and latitude points in the library
    lib_lons = concatenate([lin.lons for lin in linlib])
    lib_lats = concatenate([lin.lats for lin in linlib])
    lib_x, lib_y, lib_z = lineament.sphere2xyz(radius, pi/2-lib_lats, lib_lons)
    lib_kdt = KDTree(array([lib_x, lib_y, lib_z]).T)

    test_lons, test_lats = array(lineament.random_lonlatpoints(N))
    test_x, test_y, test_z = lineament.sphere2xyz(radius, pi/2.0-test_lats, test_lons)
    dists, near_idx = lib_kdt.query(array([test_x, test_y, test_z]).T, distance_upper_bound=d_max)

    near_idx = near_idx[where(dists<=d_max)]
    dists = dists[where(dists<=d_max)]
    return(dists, near_idx)

    #near_lons = mod(lib_lons[:,near_idx],2*pi)
    #near_lats = lib_lats[:,near_idx]
    #return(near_lons, near_lats)

    #the_fig=figure(figsize=(10,5))
    #map = the_fig.add_subplot(1,1,1)
    #map.scatter(degrees(mod(lib_lons,2*pi)),  degrees(lib_lats),  c='black', linewidth=0, s=3)
    #colors = cm.jet(linspace(0,1,len(test_lons)))
    #map.scatter(degrees(mod(test_lons,2*pi)), degrees(test_lats), c=colors,  linewidth=0)
    #map.scatter(degrees(mod(near_lons,2*pi)), degrees(near_lats), c=colors,  linewidth=0, alpha=0.7)
    #map.set_xlim([0,360])
    #map.set_ylim([-90,90])
    #return(linlib, near_lons, near_lats)

    #return(nearest_pts)
#}}}

def test_fastfit(libsize=90, linlib=None, d_max=0.01): #{{{

    print("Loading and updating mapped features")
    maplins = load_lins('output/map_nsrfit')
    maplins = lineament.update_lins(maplins)
    lz = maplins[0]

    if linlib is None:
        print("Generating lineament library with N=%d" % (libsize,) )
        linlib = lineament.lingen_nsr_library(nlats=libsize)

    print("Fitting to NSR directly")
    lz.calc_nsrfits()
    plot(degrees(lz.bs), lz.nsrdbars, linewidth=2, color='black')
    print("Fitting to NSR using linlib")
    lz.calc_nsrfits(doppel_library=linlib)
    plot(degrees(lz.bs), lz.nsrdbars, linewidth=2, color='red')

#}}}

def test_fitres(lin, dbar_max=0.125, d_max=None, inter_res=20.0, nb=180): #{{{
    """
    A function for testing how the resolution of synthetic features
    (doppelgangers) affects the accuracy of our fit metric.

    """
    europa = satstress.Satellite(open("input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",'r'))
    NSR = satstress.StressCalc([satstress.NSR(europa),])

    the_fig = figure(figsize=(10,10))
    map_ax = the_fig.add_subplot(2,1,1)
    ax1 = the_fig.add_subplot(2,1,2)

    map_ax.set_aspect('equal')
    map_ax.plot(np.degrees(lin.lons), np.degrees(lin.lats), lw=2, color='black')
    map_ax.scatter(np.degrees(lin.lons), np.degrees(lin.lats), color='black', s=5)
    init_lon, init_lat = lin.bfgcseg_midpoint()
    map_ax.grid()

    if d_max is None:
        d_max = np.median(lin.seg_lengths())

    print("Performing low-resolution fits (d_max=%g)" % (d_max,))
    lowres = lin.lonshift(radians(60))
    lowres.calc_nsrfits(stresscalc=NSR, d_max=d_max, nb=nb)
    ax1.plot(np.degrees(lowres.bs), lowres.nsrfits(use_stress=False), color='red', lw=2, label="$d_{max}=%g$" % (d_max,))
    ax1.fill_between(np.degrees(lowres.bs), lowres.nsrfits(use_stress=False), 0, color='red', alpha=0.3)
    ax1.set_ylim(0,1)
    ax1.set_yticks(linspace(0,1.0,11))
    ax1.grid()

    lowres_init_lon, lowres_init_lat = lowres.bfgcseg_midpoint()
    lowres_best_dbar, lowres_best_b = lowres.best_fit(use_stress=False)
    lowres_best_doppel = lineament.lingen_nsr(stresscalc=NSR, init_lon=lowres_init_lon+lowres_best_b, init_lat=lowres_init_lat, prop_dir="both", max_length=lowres.length, seg_len=d_max)
    lrbd = lowres_best_doppel.lonshift(-lowres_best_b-radians(60))
    map_ax.plot(np.degrees(lrbd.lons), np.degrees(lrbd.lats), lw=2, color='red')
    map_ax.scatter(np.degrees(lrbd.lons), np.degrees(lrbd.lats), color='red', s=5)

    the_fig.show()

    print("Performing high-resolution fits (d_max=%g)" % (d_max/inter_res,))
    highres = lin.lonshift(radians(120))
    highres.calc_nsrfits(stresscalc=NSR, d_max=d_max/inter_res, nb=nb)
    ax1.plot(np.degrees(highres.bs), highres.nsrfits(use_stress=False), color='green', lw=2, label="$d_{max}=%g$" % (d_max/inter_res,))
    ax1.fill_between(np.degrees(highres.bs), highres.nsrfits(use_stress=False), 0, color='green', alpha=0.3)
    ax1.set_ylim(0,1)
    ax2 = ax1.twinx()
    ax2.set_yticks(linspace(0,1.0,11))
    ax1.set_title("L=%g np=%d <l>=%g" % (lin.length, len(lin.lons), np.median(lin.seg_lengths()) ) )
    ax1.legend()

    highres_init_lon, highres_init_lat = highres.bfgcseg_midpoint()
    highres_best_dbar, highres_best_b = highres.best_fit(use_stress=False)
    highres_best_doppel = lineament.lingen_nsr(stresscalc=NSR, init_lon=highres_init_lon+highres_best_b, init_lat=highres_init_lat, prop_dir="both", max_length=highres.length, seg_len=d_max/inter_res)
    hrbd = highres_best_doppel.lonshift(-highres_best_b-radians(120))
    map_ax.plot(np.degrees(hrbd.lons), np.degrees(hrbd.lats), lw=2, color='green')
    map_ax.scatter(np.degrees(hrbd.lons), np.degrees(hrbd.lats), color='green', s=5)
    map_ax.scatter([np.mod(np.degrees(init_lon),360),], [np.degrees(init_lat),], color='black', marker='x', s=50)

    the_fig.show()
#}}}

def test_initdist(lins): #{{{
    """
    Look at the correlation between the distance between the fracture
    initiation point and how good a lineament's best fit is.

    """

    init_lons, init_lats = array([ lin.bfgcseg_midpoint() for lin in lins ]).transpose()

    mp_lons_list = []
    mp_lats_list = []
    for lin in lins:
        mp_lons, mp_lats = lin.seg_midpoints()
        mp_lons_list.append(mp_lons)
        mp_lats_list.append(mp_lats)

    lin_lengths = array([ lin.length for lin in lins ])
    d_mins = array([ lineament.spherical_distance(init_lon, init_lat, mp_lons, mp_lats).min() for init_lon,init_lat,mp_lons,mp_lats in zip(init_lons,init_lats,mp_lons_list,mp_lats_list) ])
    d_mins = d_mins/lin_lengths

    best_fits, best_bs = array([ lin.best_fit(use_stress=False) for lin in lins ]).transpose()

    the_fig = figure(figsize=(10,10))
    ax1 = the_fig.add_subplot(1,1,1)
    colors = cm.jet(linspace(0,1,11))
    ax1.scatter(d_mins,best_fits, s=[100*lin.length for lin in lins],alpha=0.5,lw=0,color='black')
    ax1.set_xlabel("Normalized dist. from crack init. pt. to prototype feature")
    ax1.set_ylabel("Quality of best fit ($0 \leq \max(f_{nsr}(b)) \leq 1$)")
    ax1.grid()
    ax1.set_title("symbol size proportional to prototype length ($R^2=%g$)" % (corrcoef(d_mins, best_fits)[0][1],))
    the_fig.show()
#}}}

###############################################################################
#         Routines that generate plots in support of the TPW paper:           #
###############################################################################

def makefigs(dbar_max=0.125, maps=False, hists=False, examples=False, stats=False, stress=False, tpw=False, lindensity=False, all=False): #{{{
    """
    A wrapper that takes input parameters for the figures, and calls all of the
    individual figure plotting routines, below.

    """
    from osgeo import gdal
    from os import path

    # TODO: verify that saving PDF versions work once the stats plots are working
    # If we want to re-make all of the figures:
    if all is True:
        lindensity = True # works
        examples = True   # works
        hists = True      # works
        stats = True      # LinLengthDist colors are funny; LinLatLonStats doesn't work
        maps = True       # works
        tpw = True        # works

    maplins, synthlins, crazylins, gclins, tpwlins = reload_nsrfits()

    if maps is True: #{{{2
        print("Plotting Mapped Lineaments")
        FitMap(maplins, nbins=9, titlestr="Global Lineaments as Mapped", dbar_max=dbar_max,\
               outfile=path.join(figdir,'FitMap_Mapped'))

        print("Plotting Randomized Mapped Lineaments")
        FitMap(crazylins, nbins=9, titlestr="Randomized Global Lineaments", dbar_max=dbar_max,\
               outfile=path.join(figdir,'FitMap_Crazy'))

        print("Plotting Pre-TPW Lineaments, fit to NSR stresses")
        FitMap(tpwlins, nbins=9, titlestr="Mapped Lineaments before TPW w/ PNP=80E10N", dbar_max=dbar_max,\
               outfile=path.join(figdir,'FitMap_PreTPW'))

        print("Plotting random great circle segments")
        FitMap(gclins, nbins=9, titlestr="Random Great Circle Segments", dbar_max=dbar_max,\
               outfile=path.join(figdir,'FitMap_GreatCircle'))

        print("Plotting synthetic NSR lineaments")
        FitMap(synthlins, titlestr="Perfect Synthetic NSR Lineaments for b=0", dbar_max=dbar_max,\
               outfile=path.join(figdir,'FitMap_SyntheticNSR'))

        print("Plotting Mapped Lineaments, relative to stress field")
        FitMap(maplins, nbins=9, titlestr="Mapped Lineaments Backrotated for Best Fit", dbar_max=dbar_max, derotate=True, showbad=False,\
               outfile=path.join(figdir,'FitMap_Derotated'))

        print("Plotting Derotated Great Circle Segments")
        FitMap(gclins, nbins=9, titlestr="Random Great Circle Segments Backrotated for Best Fit", dbar_max=dbar_max, derotate=True, showbad=False,\
               outfile=path.join(figdir,'FitMap_DerotatedGC'))
    #}}}2

    if hists is True: #{{{2
        # These three figures demonstrate robustness of the metric
        ActHist_ByLength(maplins, dbar_max=dbar_max, norm_by_all=True)
        ActHist_BySin(maplins, dbar_max=dbar_max, norm_by_all=True)
        ActHist_ByDbar(maplins, norm_by_all=False)

        # Statistics of the dataset, and an attempt at synthesizing our map
        numhists = 100
        print("Calculating activity histories for %d resampled maps" % (numhists,) )
        ActHist_MapStats(maplins, N=numhists)

        print("Calculating activity histories for %d spin cycle maps" % (numhists,) )
        ActHist_SpinCycleStats(maplins, N=numhists)

        print("Calculating activity histories for %d partially derotated maps" % (numhists,) )
        ActHist_PeakStats(maplins, N=numhists, scale=np.radians(15))

        print("Calculating activity histories for %d synthetic map recreations" % (numhists,) )
        ActHist_Synthesized(maplins, N=numhists, peak_frac=0.4, scale=np.radians(15))

        # TODO:
        #print("Calculating activity histories for %d jackstraws maps" % (numhists,) )
        #ActHist_JackStrawsStats(maplins, N=numhists)

        #print("Calculating activity histories for %d resampled great circle maps" % (numhists,) )
        #ActHist_GreatCircleStats(maplins, N=numhists)

        #print("Calculating activity histories for %d derotated maps" % (numhists,) )
        #ActHist_Derotated(maplins, synthlins, crazylins, N=numhists)

        print("Plotting distribution of activity amplitudes")
        ActHist_AmpProb(maplins, tpwlins, synthlins)
    #}}}2
 
    if examples is True: #{{{2
        # List of canonical features from the map: {{{3
        good_nsr1    = maplins[1]   # 1500 km long, symmetric arcuate lineament in the N. hemisphere
        good_nsr2    = maplins[53]  #  700 km long, asymmetric arcuate lineament in the N. hemisphere
        good_nsr3    = maplins[108] # 1800 km long, nearly symmetric arcuate lineament in the S. hemisphere
        bad_nsr1     = maplins[20]  #  553 km long, north-south lineament in the S. hemisphere
        bad_nsr2     = maplins[6]   # 1222 km long, diagonal, grazing the equator
        bad_nsr3     = maplins[36]  # 1175 km long, diagonal, equator crossing
        bad_nsr4     = maplins[54]  # 1036 km long, north-south, equator crossing
        bad_nsr5     = maplins[70]  # 1061 km long, one of the SCDs, crossing the equator
        bad_nsr6     = maplins[120] #  640 km, N-S equator crossing
        bad_nsr7     = maplins[122] # 1300 km, N-S equator crossing
        cycloid1     = maplins[112] # 1132 km long cycloid, 4-5 arcs, near 30S.  Low dbar, high delta
        cycloid2     = maplins[137] #  458 km long cycloid, 5 arcs, near 60S.  Low dbar, high delta
        cycloid3     = maplins[148] #  776 km long cycloid, 5 arcs, 45-70N.
        cycloid4     = maplins[155] # 1711 km long semi-cycloidal, actually passed the test
        cycloid5     = maplins[159] # 1527 km long semi-cycloidal, actually passed the test
        sinuous1     = maplins[23]  # 1334 km long, sinuous from 30N to 75N
        sinuous2     = maplins[136] # 1189 km long, sinuous from 30S to 70S
        short_hilat1 = maplins[11]  #  450 km, east-west, just above 30N... fits perfectly
        short_hilat2 = maplins[78]  #  200 km, diagonal, 50-60N, fits perfectly
        short_hilat3 = maplins[11]  #  183 km, east-west, just above 30N... fits perfectly
        short_hilat4 = maplins[160] #  200 km, diagonal, 75S... fits very well
        short_lolat1 = maplins[26]  #  500 km, diagonal, between 0N and 30N, fits perfectly
        short_lolat2 = maplins[41]  #  177 km, diagonal, between 0N and 30N, does not fit because of dbar
        short_lolat3 = maplins[43]  #  197 km, diagonal, between 0N and 30N, fits very well
        dbar_fail1   = maplins[7]   # 1073 km, arcuate, dbar_min doesn't quite line up with delta_min, and so it gets lost.
        dbar_fail2   = maplins[62]  #  500 km, east-west almost at the equator.  Large number of "good fits" possible.
        dbar_fail3   = maplins[63]  #  500 km, east-west almost at the equator.  Large number of "good fits" possible.
        dbar_fail4   = maplins[115] # 1262 km, looks good, but dbar just doesn't quite get low enough... Grr.
        #}}}3

        eg_lins = [ good_nsr1, bad_nsr2, cycloid1, sinuous1, sinuous2, short_hilat1, short_lolat2, dbar_fail2 ]

        eg_labels = [ "Good NSR", "Bad NSR", "Cycloial", "Sinuous", "Sinuous",\
                      "Short straight mid-latitude", "Short straight low-latitude",\
                      "Poorly constrained best fit" ]

        FitCurveExamples(eg_lins, labels=eg_labels, dbar_max=dbar_max)

    #}}}2

    if stats is True: #{{{2
        LinLengthDist(maplins, label="Lineaments as Mapped",             outfile='LinLenDist_Mapped')
        LinLengthDist(gclins,  label="Random Great Circle Segments",     outfile='LinLenDist_GC')
        LinLengthDist(tpwlins, label="pre-TPW (pole=80E10N) Lineaments", outfile='LinLenDist_TPW')

        DbarLengthCorr(gclins)
        DbarSinuosityCorr(maplins)

        LinLatLonStats(maplins,   label="Lineaments as Mapped",             outfile='LinStats_Mapped')
        LinLatLonStats(gclins,    label="Great Circle Segments",            outfile='LinStats_GreatCircle')
        LinLatLonStats(tpwlins,   label="pre-TPW Lineaments (pole=80E10N)", outfile='LinStats_PreTPW')
        LinLatLonStats(synthlins, label="Synthetic NSR Lineaments",         outfile='LinStats_SyntheticNSR')
        LinLatLonStats(maplins,   label="Derotated Lineaments",             outfile='LinStats_Derotated')
        LinLatLonStats(crazylins, label="Randomized Lineaments",            outfile='LinStats_Jackstraws')
    #}}}2

    if lindensity is True: #{{{2
        # This shows the density of the mapped features, and a bunch of information
        # about how that relates to the spacecraft observations.
        LinDensityMap(maplins, maxdist=250, label="Density of Mapped Features") #}}}2

    if tpw is True: #{{{2
        # Plots having to do with the search for a best TPW pole...
        tpw_polefits()

    #}}}2

    if stress is True: #{{{
        # need to generate maps of the NSR stresses which are behind all of
        # this stuff... GridCalc here we come!
        pass

    #}}}2

# end makefigs }}}

def LinLatLonStats(goodlins=[], badlins=[], label="", outfile=None): #{{{
    """
    Shows the distribution of latitudes and longitudes of the lineament
    segments over the surface of the satellite, relative to an expected random
    distribution.

    """

    good_lat_pile = []
    good_lon_pile = []
    bad_lat_pile  = []
    bad_lon_pile  = []
    sat_radius_km = 1561

    dbar_max = max([ max(lin.nsrdbars) for lin in goodlins if max(lin.nsrdbars) > 0 ])

    good_length = sum([ lin.length for lin in goodlins ])*sat_radius_km
    bad_length = sum([ lin.length for lin in badlins ])*sat_radius_km

    total_length = good_length + bad_length

    for goodlin in goodlins:
        goodsegments  = goodlin.segments()
        goodmidpoints = goodlin.midpoints()
        good_seg_lens = [ int(round(sat_radius_km*seg.calc_length())) for seg in goodsegments ]

        for mp,seg_len in zip(goodmidpoints,good_seg_lens):
            good_lon_pile += [mp[0],]*seg_len
            good_lat_pile += [mp[1],]*seg_len

    for badlin in badlins:
        badsegments   = badlin.segments()
        badmidpoints  = badlin.midpoints()
        bad_seg_lens = [ int(round(sat_radius_km*seg.calc_length())) for seg in badsegments ]

        for mp,seg_len in zip(badmidpoints,bad_seg_lens):
            bad_lon_pile += [mp[0],]*seg_len
            bad_lat_pile += [mp[1],]*seg_len

    the_fig = figure(figsize=(12,8))
    lon_ax = the_fig.add_subplot(2,1,1)
    lon_ax.set_title(r'%s max($\delta_{rms}$)=%.0f$^\circ$ max($\bar{D}$)=%.0g' % (label, degrees(delta_max), dbar_max) )
    if len(bad_lon_pile) > 0:
        lins_n, lin_bins, lins_patches = lon_ax.hist([mod(degrees(good_lon_pile),180),mod(degrees(bad_lon_pile),180)],\
                                                     bins=90, range=(0,180), label="passed", histtype='barstacked', rwidth=1.0, facecolor='#666666')
        lins_patches[1][0].set_label("failed")
        for patch in lins_patches[1]:
            patch.set_facecolor('#CCCCCC')
    else:
        lon_ax.hist(mod(degrees(good_lon_pile),180), bins=90, range=(0,180), label="passed", histtype='barstacked', rwidth=1.0, facecolor='#666666')

    lon_ax.set_xlim(0,180)
    lon_ax.axhline(y=total_length/90.0, color='black', linewidth=3, linestyle='--', label="random")
    lon_ax.set_xticks(linspace(0,180,7))
    lon_ax.set_xlabel("longitude (modulo $\pi$)")
    lon_ax.set_ylabel("length [km]")

    lat_ax = the_fig.add_subplot(2,1,2)
    if len(bad_lat_pile) > 0:
        lins_n, lin_bins, lins_patches = lat_ax.hist([abs(degrees(good_lat_pile)),abs(degrees(bad_lat_pile))],\
                                                     bins=90, range=(0,90), label="passed", histtype='barstacked', rwidth=1.0, facecolor='#666666')
        lins_patches[1][0].set_label("failed")
        for patch in lins_patches[1]:
            patch.set_facecolor('#CCCCCC')
    else:
        lat_ax.hist(abs(degrees(good_lat_pile)), bins=90, range=(0,90), label="passed", histtype='barstacked', rwidth=1.0, facecolor='#666666')
    lat_ax.set_xlim(0,90)
    lats = linspace(0,pi/2,100)
    lat_ax.plot(degrees(lats), (total_length/90)*(cos(lats))*pi/2, color='black', linestyle='--', linewidth=3, label="random")
    lat_ax.set_xticks(linspace(0,90,7))
    lat_ax.set_xlabel("abs(latitude)")
    lat_ax.set_ylabel("length [km]")
    lat_ax.legend()

    if outfile is None:
        the_fig.show()
    else:
        the_fig.savefig(outfile)

# end LinLatLonStats }}}

def LinLengthDist(lins, label="", outfile=None): #{{{
    """
    Plot the overall distribution of feature lengths, color coding retained and
    rejected features differently.

    """
    # Need to generate two arrays of values pertaining to the combined set of
    # good and bad lineaments.
    #  - one of floats: lineament lengths
    #  - one of booleans: True/False depending on whether it was kept/rejected
    # Then I need to sort them both in order of lineament length.

    # The plot I want to draw is a curve, the Y-value of which is lineament
    # length, and the x-value is just N, the number of the lineament, ordered by length.
    # For each point in the curve, if the lineament was kept, it should be black beneath
    # and if it was rejected, it should be white.

    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000

    # generate a list containing the lineament lengths in km:
    linlengths = [ l.length*sat_radius_km for l in lins ]

    # and a list of boolean values saying whether or not it passed:
    best_fit = [ l.best_fit()[0] for l in lins ]

    lengthskept = array([ (linlengths[n],best_fit[n]) for n in range(len(lins)) ], dtype=[('length',float),('best_fit',float)])
    lengthskept.sort(order='length')

    the_fig = figure(figsize=(9,3))
    plot_ax = the_fig.add_subplot(1,1,1)
    plot_ax.plot(lengthskept['length'], color='black', linewidth=2)
    plot_ax.fill_between(range(len(lengthskept)), zeros(len(lengthskept)), lengthskept['length'], where=lengthskept['best_fit'], facecolor=lengthskept['best_fit'])
    plot_ax.set_title(label + ' Length Distribution')
    plot_ax.set_xlabel('N')
    plot_ax.set_ylabel('lineament length [km]')

    if outfile is None:
        the_fig.show()
    else:
        the_fig.savefig(outfile)

#}}}

def ActHist_ByLength(lins, dbar_max=0.125, norm_by_all=True): #{{{
    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000
    labels_by_length = []
    lins_by_length = []
    lengths = [0,150,300,500,1000,2500]
    numlens = len(lengths)-1
    for N in range(numlens):
        lins_by_length.append([])
        labels_by_length.append("%d km < L < %d km" % (lengths[N], lengths[N+1]) )
        for lin in lins:
            lin_length_km = lin.length*sat_radius_km
            if (lin_length_km > lengths[N]) and (lin_length_km < lengths[N+1]):
                lins_by_length[N].append(lin)
                continue

    activity_history(lins_by_length, dbar_max=dbar_max, labels=labels_by_length, norm_by_all=norm_by_all, outfile='ActHist_ByLength')
#}}}

def ActHist_BySin(lins, dbar_max=0.125, norm_by_all=True): #{{{
    labels_by_sin = []
    lins_by_sin = []
    sins = [1.0, 1.02, 1.04, 1.06, 1.08, 1.10, 1.20]
    numsins = len(sins)-1
    for N in range(numsins):
        lins_by_sin.append([])
        labels_by_sin.append("%.3g < S < %.3g" % (sins[N], sins[N+1]) )
        for lin in lins:
            S = lin.sinuosity()
            if (S > sins[N]) and (S < sins[N+1]):
                lins_by_sin[N].append(lin)
                continue

    activity_history(lins_by_sin, dbar_max=dbar_max, labels=labels_by_sin, norm_by_all=norm_by_all, outfile='ActHist_BySin')
#}}}

def ActHist_ByDbar(maplins, norm_by_all=False): #{{{
    """
    Show that the time variability is not sensitive to the particular value
    chosen for dbar_max.

    """
    the_fig = figure(figsize=(12,8))

    dbars = (0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25)
    colors = cm.jet(linspace(0,1,len(dbars)))

    for dbar_max,N in zip(dbars,range(len(dbars))):
        activity_history([maplins,], labels=[r'$\bar{D}_{max}=%.3g$'%(dbar_max,),], norm_by_all=norm_by_all, the_fig=the_fig, dbar_max=dbar_max, colors=[colors[N],], outfile="ActHist_ByDbar")
#}}}

def ActHist_MapStats(maplins, nbins=20, N=100): #{{{
    """
    Do a Monte Carlo subsampling of our lineament map and compare the activity
    history of the true map with the synthetic ones, to see how self-consistent
    it is.

    """

    # pseudo-maps in transparent gray
    maps = [ linsample(proto=maplins, pool=maplins, nbins=nbins) for i in range(N) ]
    alphas = [10./N,]*N
    colors = ['gray',]*N
    labels = [None,]*(N-1)
    labels.append('Synthetic Subsamples (N=%d)' % (N,))

    # map minus E15 mapping swath in solid blue, as an example of geographic exclusion
    labels.append('Map minus E15 Swath')
    noE15lins = [ lin for lin in maplins if (degrees(mod(array(lin.lons),2*pi)) > 290).all() or (degrees(mod(array(lin.lons),2*pi)) < 240).all() ]
    maps.append(noE15lins)
    colors.append('blue')
    alphas.append(1.0)

    # Real map in solid black
    labels.append('Mapped Lineaments')
    maps.append(maplins)
    colors.append('black')
    alphas.append(1.0)

    the_fig = figure(figsize=(12,8))
    activity_history(maps, labels=labels, colors=colors, alphas=alphas, the_fig=the_fig, verbose=False, titlestr="Natural Variability in the Mapped Dataset", outfile="ActHist_MapStats")
#}}}

def ActHist_SpinCycleStats(maplins, nbins=20, N=100): #{{{
    """
    Do a Monte Carlo subsampling of our lineament map, randomizing the
    longitudes of the subsampled data, and then compare the activity history of
    the true map with the synthetic ones, to see how far outside of the normal
    range our data is.

    """

    maps = [ make_crazy(linsample(proto=maplins, pool=maplins, nbins=nbins),tpw=False) for i in range(N) ]
    maps.append(maplins)

    # pseudo-maps in gray
    colors = ['gray',]*N
    # real map in black
    colors.append('black')

    # pseudo-maps nearly transparent
    alphas = [15./N,]*N
    # real map opaque
    alphas.append(1.0)

    labels = [None,]*(N-1)
    labels.append('Spun Subsamples (N=%d)' % (N,))
    labels.append('Mapped Lineaments')

    the_fig = figure(figsize=(12,8))
    activity_history(maps, labels=labels, colors=colors, alphas=alphas, the_fig=the_fig, verbose=False, outfile="ActHist_SpinCycleStats")
#}}}

def ActHist_PeakStats(maplins, nbins=20, N=100, scale=np.radians(15)): #{{{
    """
    Create a large number of synthetic datasets composed partly of a
    longitudinally randomized sampling of maplins, and partly from a de-rotated
    set of maplins, and see how their activity histories compare to that of the
    mapped features in their real locations.

    """

    map_acthist = calc_acthist(maplins)
    mu = int(where(np.abs(map_acthist-np.max(map_acthist)) < 1e-6)[0])*(np.pi/len(map_acthist))
    maps = []
    for i in range(N):
        peak_dist = linsample(proto=maplins, pool=maplins, nbins=nbins)
        peak_dist = [ lin.lonshift(-b+lin.best_fit()[1]) for lin,b in zip(peak_dist,np.random.normal(loc=mu, scale=scale, size=len(peak_dist))) ]
        maps.append(peak_dist)


    maps.append(maplins)

    # pseudo-maps in gray
    colors = ['gray',]*N
    # real map in black
    colors.append('black')

    # pseudo-maps nearly transparent
    alphas = [10./N,]*N
    # real map opaque
    alphas.append(1.0)

    labels = [None,]*(N-1)
    labels.append('Synthetic Peak (N=%d)' % (N,))
    labels.append('Mapped Lineaments')

    the_fig = figure(figsize=(12,8))
    titlestr='Peak subsampled from derotated map $\mu=%g^\circ$, $\sigma=%g^\circ$' % (degrees(mu),degrees(scale))
    activity_history(maps, labels=labels, colors=colors, alphas=alphas, the_fig=the_fig, verbose=False, titlestr=titlestr, outfile='ActHist_PeakStats')
#}}}

def ActHist_Synthesized(maplins, nbins=20, N=100, peak_frac=0.4, scale=np.radians(15)): #{{{
    """
    Create a large number of synthetic datasets composed partly of a
    longitudinally randomized sampling of maplins, and partly from a de-rotated
    set of maplins, and see how their activity histories compare to that of the
    mapped features in their real locations.

    """

    map_acthist = calc_acthist(maplins)
    mu = int(where(np.abs(map_acthist-np.max(map_acthist)) < 1e-6)[0])*(np.pi/len(map_acthist))
    maps = []
    for i in range(N):
        random_part = make_crazy(linsample(proto=maplins, pool=maplins, nbins=nbins, fraction=(1.0-peak_frac)),tpw=False)
        derotated_part = linsample(proto=maplins, pool=maplins, fraction=peak_frac, nbins=nbins)
        derotated_part = [ lin.lonshift(-b+lin.best_fit()[1]) for lin,b in zip(derotated_part,np.random.normal(loc=mu, scale=scale, size=len(derotated_part))) ]
        maps.append(concatenate([random_part,derotated_part]))


    maps.append(maplins)

    # pseudo-maps in gray
    colors = ['gray',]*N
    # real map in black
    colors.append('black')

    # pseudo-maps nearly transparent
    alphas = [10./N,]*N
    # real map opaque
    alphas.append(1.0)

    labels = [None,]*(N-1)
    labels.append('%d%% Rand + %d%% Peak (N=%d)' % (int(100*(1-peak_frac)),int(100*peak_frac),N))
    labels.append('Mapped Lineaments')

    the_fig = figure(figsize=(12,8))
    titlestr='Mix of random and peaked features w/ $\mu=%g^\circ$, $\sigma=%g^\circ$' % (degrees(mu),degrees(scale))
    activity_history(maps, labels=labels, colors=colors, alphas=alphas, the_fig=the_fig, verbose=False, outfile='ActHist_Synthesized' )
#}}}

def ActHist_AmpProb(maplins, tpwlins, synthlins): #{{{
    """
    Create a figure showing the likelihood of getting various activity history
    amplitudes:
      - map completely randomized (TPW + spin: XXX < |H(b)| < YYY)
      - map randomized in longitude (spin only: 0.007 < |H(b)| < 0.057)
      - features as mapped (|H(b)| = 0.0936)
      - greatest amplitude possible with mapped features (|H(b)| = 0.427)
      - greatest amplitude possible with synthetic features (|H(b)| = 0.7)

    """
    import pickle

    spin_amps = pickle.load(open('output/spin_amps.pkl'))
    tpw_amps = pickle.load(open('output/tpw_spin_amps.pkl'))
    pnp_lons, pnp_lats, pnp_amps = load_tpw_results(fast=True)

    map_amp = acthist_amplitude(maplins)
    synth_amp = acthist_amplitude(synthlins)
    tpw80E10N_amp = acthist_amplitude(tpwlins)
    bestfits = [ lin.best_fit()[1] for lin in maplins ]
    derotated = [ lin.lonshift(b) for lin,b in zip(maplins, bestfits) ]
    derotated_amp = acthist_amplitude(derotated)

    the_fig = figure(figsize=(12,8))
    hist_ax = the_fig.add_subplot(1,1,1)
    spin_counts, spin_bins, spin_patches = hist_ax.hist(spin_amps, bins=200, range=(0,0.125), color='red',   lw=2, normed=True, histtype='step', label='spun')
    tpw_counts,  tpw_bins,  tpw_patches  = hist_ax.hist(tpw_amps,  bins=50,  range=(0,0.125), color='green', lw=2, normed=True, histtype='step', label='jackstraws')
    pnp_counts,  pnp_bins,  pnp_patches  = hist_ax.hist(pnp_amps,  bins=25,  range=(0,0.125), color='blue',  lw=2, normed=True, histtype='step', label='map + TPW')
    hist_ax.set_xlim(0,0.11)
    hist_ax.set_xticks(linspace(0,0.11,12))

    # add a vertical line to show where the actual mapped features fall:
    hist_ax.axvline(x=median(spin_amps), linewidth=3, linestyle=':', color='red')
    hist_ax.annotate('median: %.2g' % (median(spin_amps),), xy=(median(spin_amps)-0.001, max(spin_counts)/5.0), rotation='vertical', ha='right', color='red')

    hist_ax.axvline(x=median(tpw_amps), linewidth=3, linestyle=':', color='green')
    hist_ax.annotate('median: %.2g' % (median(tpw_amps),), xy=(median(tpw_amps)+0.001, max(tpw_counts)/3.0), rotation=270, ha='left', color='green')

    hist_ax.axvline(x=tpw80E10N_amp, linewidth=3, linestyle='--', color='gray')
    hist_ax.annotate('PNP 80E10N: %.2g' % (tpw80E10N_amp,), xy=(tpw80E10N_amp-0.001, 0.7*max(spin_counts)), rotation='vertical', ha='right', color='gray')

    hist_ax.axvline(x=map_amp, linewidth=3, linestyle='--', color='black')
    hist_ax.annotate('mapped: %.2g' % (map_amp,), xy=(map_amp-0.001, max(spin_counts)/5.0), rotation='vertical', ha='right', color='black')

    hist_ax.set_title(r'$|H(b)|$ distributions for several lineament populations')
    hist_ax.set_xlabel(r'$|H(b)|$')
    hist_ax.set_ylabel('N')
    hist_ax.legend(loc='upper right')
    hist_ax.grid()

    if save_fmt is not None:
        the_fig.savefig('output/ActHist_AmpProb.'+save_fmt)

#}}}

def DbarLengthCorr(lins): #{{{
    """
    Several scatter plots showing the correlation (or lack thereof) between
    quality of best fit and lineament length, in equatorial and non-equatorial
    regions.

    """

    the_fig = figure(figsize=(12,8))
    ax1 = the_fig.add_subplot(1,1,1)

    hilatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) > 30 ]
    hilat_best_fits = array([ lin.nsrdbars.min() for lin in hilatlins ])
    hilat_lengths = array([ lin.length for lin in hilatlins ])*1561
    hilat_r2 = corrcoef(hilat_lengths, hilat_best_fits)[0,1]**2
    hilat_symbs = scatter(hilat_lengths, hilat_best_fits, c='blue', marker='s', label=r'$|\theta_{min}|>30^\circ R^2=%.3g$' % ( hilat_r2, ) )

    lolatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) <= 30 ]
    lolat_best_fits = array([ lin.nsrdbars.min() for lin in lolatlins ])
    lolat_lengths = array([ lin.length for lin in lolatlins ])*1561
    lolat_r2 = corrcoef(lolat_lengths,lolat_best_fits)[0,1]**2
    lolat_symbs = scatter(lolat_lengths, lolat_best_fits, c='red', marker='^', label=r'$\|\theta_{min}|\leq30^\circ R^2=%.3g$' % ( lolat_r2, ) )

    ax1.legend(loc='lower right')
    ax1.set_xlabel('lineament length [km]')
    ax1.set_ylabel(r'$\min(\bar{D}(b))$')
    ax1.set_ylim(0,0.3)
    ax1.set_xlim( min((hilat_lengths.min(),lolat_lengths.min())),\
                  max((hilat_lengths.max(),lolat_lengths.max())) )
    ax1.set_title('Effects of length and latitude on fit to NSR for great circle segments, N=%d'\
                   % (len(hilatlins)+len(lolatlins)) )

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'DbarLengthCorr.'+save_fmt))

# end DbarLengthCorr}}}

def DbarSinuosityCorr(lins): #{{{
    """
    Scatter plot showing the correlation (or lack thereof) between quality of
    best fit and lineament sinuosity in equatorial and non-equatorial regions.

    """

    the_fig = figure(figsize=(12,8))
    ax1 = the_fig.add_subplot(1,1,1)

    hilatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) > 30 ]
    hilat_best_fits = array([ lin.nsrdbars.min() for lin in hilatlins ])
    hilat_sins = array([ lin.sinuosity() for lin in hilatlins ])
    hilat_r2 = corrcoef(hilat_sins,hilat_best_fits)[0,1]**2
    hilat_symbs = scatter(hilat_sins, hilat_best_fits, c='blue', marker='s', label=r'$|\theta_{min}|>30^\circ R^2=%.3g$' % ( hilat_r2, ) )

    lolatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) <= 30 ]
    lolat_best_fits= array([ lin.nsrdbars.min() for lin in lolatlins ])
    lolat_sins = array([ lin.sinuosity() for lin in lolatlins ])
    lolat_r2 = corrcoef(lolat_sins,lolat_best_fits)[0,1]**2
    lolat_symbs = scatter(lolat_sins, lolat_best_fits, c='red', marker='^', label=r'$\|\theta_{min}|\leq30^\circ R^2=%.3g$' % ( lolat_r2, ) )

    ax1.legend(loc='lower right')
    ax1.set_xlabel('sinuosity')
    ax1.set_ylabel(r'$\min(\bar{D}(b))$')
    ax1.set_xlim(1,1.1)
    ax1.set_ylim(0,0.3)
    ax1.set_title('Effects of sinuosity and latitude on fit to NSR for all mapped features, N=%d'\
                   % (len(hilatlins)+len(lolatlins)) )

    if save_fmt is not None:
        the_fig.savefig(os.path.join(figdir,'DbarLengthCorr.'+save_fmt))

# end DbarSinuosityCorr }}}

def LinDensityMap(lins, maxdist=500, N=0, label="", cmap=cm.jet): #{{{
    """
    Calculate an interpolated grid showing the density of lineaments on the
    surface, as reflected by the sum of the lengths of the lineament segments
    within the given distance maxdist, of N sample points on the surface of the
    satellite.

    Compare that lineament density to the resolution of coverage we have from
    the USGS mosaic that was used to make the lineament map, determining to
    what degree the resolution of coverage determines the density of mapped
    lineaments.

    """
    import matplotlib.colors as colors
    from osgeo import gdal
    import re

    # The number of points we need to sample the distribution at scales
    # proportional to the surface area that each point is sampling.  This
    # number works well for Europa anyway.
    if N==0:
        N = 5.0e8/(maxdist**2)

    #randlons, randlats = lineament.random_lonlatpoints(N)
    randlons, randlats = fibonacci_sphere(N)
    randlons = mod(randlons,2*pi)
    reglats = linspace(-90,90,180)
    reglons = linspace(0,360,360)

    # adding these corner points makes sure that the griddata goes to the edges of the map
    edge_N = 10
    toplats = array(pi/1.99).repeat(edge_N)
    toplons = linspace(-0.01,2.01*pi,edge_N)
    bottomlats = -toplats
    bottomlons = toplons
    westlats = linspace(-pi/1.98,pi/1.98,edge_N)
    westlons = array(-0.02).repeat(edge_N)
    eastlats = westlats
    eastlons = array(2.01*pi).repeat(edge_N)
    randlons = hstack([randlons, toplons, bottomlons, eastlons, westlons])
    randlats = hstack([randlats, toplats, bottomlats, eastlats, westlats])

    seglons = array([])
    seglats = array([])
    seglens = array([])
    for lin in lins:
        newlons, newlats = lin.seg_midpoints()
        newlens = lin.seg_lengths()
        seglons = concatenate([seglons, newlons])
        seglats = concatenate([seglats, newlats])
        seglens = concatenate([seglens, newlens])

    # For each point (lon,lat) defined by randlons, randlats, calculate the sum
    # of the lengths of the segments closer than dist radians of arc away:
    nsegs = len(seglens)
    lensums = array([])

    print("Calculating lineament density map with d=%d km and N=%d" % (maxdist, N) )
    for lon,lat in zip(randlons, randlats):
        lon_arr = array(lon)
        lat_arr = array(lat)
        newsum = sum(where(lineament.spherical_distance(lon_arr.repeat(nsegs),lat_arr.repeat(nsegs),seglons,seglats) < maxdist/1561.0, seglens, 0.0))
        lensums = hstack([lensums,newsum])

    # convert these values of radians per footprint, into m/km^2
    lindensity = lensums*1561*1000/(pi*maxdist**2)

    lindensity_grid = griddata(degrees(randlons), degrees(randlats), lindensity, reglons, reglats)

    # get rid of the out-of-bounds points we added so that the interpolated grid would stretch to the edges of the map
    randlons, randlats, lindensity = randlons[:N], randlats[:N], lindensity[:N]

    lindens_fig = figure(figsize=(12,8))
    lindens_ax = lindens_fig.add_subplot(1,1,1)
    lindens_ax.contourf(reglons, reglats, lindensity_grid, 64, cmap=cmap, linewidth=0)
    lindens_ax.scatter(degrees(randlons), degrees(randlats), marker='o', color='white', s=2, edgecolor='white', alpha=0.5, linewidth=0)
    lindens_ax.scatter(mod(degrees(seglons),360), degrees(seglats), marker='o', color='black', s=200*seglens, edgecolor='black', alpha=0.375, linewidth=0)

    lindens_ax.set_xlim(0,360)
    lindens_ax.set_xticks(linspace(0,360,13))
    lindens_ax.set_ylim(-90,90)
    lindens_ax.set_yticks(linspace(-90,90,7))
    lindens_ax.set_title(label+" N=%d, d=%g km" % (N,maxdist) )

    lindensnorm = colors.Normalize(vmin=0,vmax=lindensity.max())
    cb_ax,kw = colorbar.make_axes(lindens_ax, pad=0.05, orientation='horizontal', shrink=0.5)
    colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=lindensnorm, orientation='horizontal')
    cb_ax.set_xlabel(r"mapped lineament density [ m/km$^2$]")

    if save_fmt is not None:
       lindens_fig.savefig(os.path.join(figdir,'LinDensity_Grid.'+save_fmt))

    print("Reading in USGS mosaic raster data")
    # Read in the resolution/coverage map from USGS:
    # The data cube has 4 bands:
    # 1. resolution [km/px]
    # 2. emission angle [degrees]
    # 3. incidence angle [degrees]
    # 4. phase angle [degrees]
    EuropaGrid = gdal.Open("input/europa_coverage_simpcyl.cub")
    raster_numlons = EuropaGrid.RasterXSize
    raster_numlats = EuropaGrid.RasterYSize
    raster_lonidx = raster_numlons*(randlons/(2*pi))-1.0
    raster_lonidx = raster_lonidx.round().astype(int)
    raster_latidx = (raster_numlats*(pi/2-(randlats))/pi)-1.0
    raster_latidx = raster_latidx.round().astype(int)

    # changing resolution to [px/km]:
    resolution_raster = 1.0 / EuropaGrid.GetRasterBand(1).ReadAsArray()
    # mask the raster to remove extreme polar distortions, nodata values
    resolution_raster = ma.masked_outside(resolution_raster, 0.0, 6.0)

    emission_angle_raster = EuropaGrid.GetRasterBand(2).ReadAsArray()
    emission_angle_raster = ma.masked_where(emission_angle_raster < 0, emission_angle_raster)
    incidence_angle_raster = EuropaGrid.GetRasterBand(3).ReadAsArray()
    incidence_angle_raster = ma.masked_where(incidence_angle_raster < 0, incidence_angle_raster)
    phase_angle_raster    = EuropaGrid.GetRasterBand(4).ReadAsArray()
    phase_angle_raster    = ma.masked_where(phase_angle_raster < 0, phase_angle_raster)

    rasters   = [resolution_raster, emission_angle_raster, incidence_angle_raster, phase_angle_raster]
    rastnames = ['Resolution', 'Emission Angle', 'Incidence Angle', 'Phase Angle']
    rastunits = ['[km/px]', '[degrees]', '[degrees]', '[degrees]']
    rastfigs = [ figure(figsize=(6,12)), figure(figsize=(6,12)), figure(figsize=(6,12)), figure(figsize=(6,12)) ]
    rastoutfiles = [ 'LinDensity_Resolution', 'LinDensity_Emission', 'LinDensity_Incidence', 'LinDensity_Phase' ]

    for raster,rastname,rastunit,rastfig,outfile in zip(rasters, rastnames, rastunits, rastfigs, rastoutfiles):
        rast_ax = rastfig.add_subplot(2,1,1)
        rast_ax.imshow(raster, extent=(0,360,-90,90))
        rast_ax.scatter(mod(degrees(seglons),360), degrees(seglats), marker='o', color='black', s=200*seglens, edgecolor='black', alpha=0.375, linewidth=0)

        rast_ax.set_xlim(0,360)
        rast_ax.set_xticks(linspace(0,360,13))
        rast_ax.set_ylim(-90,90)
        rast_ax.set_yticks(linspace(-90,90,7))
        rast_ax.set_title("USGS Europa Mosaic: %s" % (rastname,) )
        rast_norm = colors.Normalize(vmin=0,vmax=raster.max())
        cb_ax,kw = colorbar.make_axes(rast_ax, pad=0.1, orientation='horizontal', shrink=0.5)
        colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=rast_norm, orientation='horizontal')
        cb_ax.set_xlabel(rastunit)

        # The array lindensity contains values of lineament density for the N
        # (lon,lat) points defined by (randlons, randlats).  I need a similar
        # set of values for the same (lon,lat) points, but corresponding to the
        # values stored within various rasters we're comparing to:
        randsamples = ma.masked_less(raster[raster_latidx, raster_lonidx],0)
        # Calculate the correlation between lindensity and the raster in question:
        lindensity_raster_corrcoef = ma.corrcoef(randsamples, lindensity)[0,1] 
        print(r"lindensity v. %s: $R^2$=%g" % (rastname, lindensity_raster_corrcoef**2) )

        # Make a scatter plot showing the correlation (or lack thereof):
        lin_rast_corr_ax = rastfig.add_subplot(2,1,2)
        lin_rast_corr_ax.scatter(randsamples, lindensity, s=10, linewidth=0, marker='o', color='black', alpha=0.375)
        lin_rast_corr_ax.set_xlabel("USGS Mosaic %s, %s" % (rastname, rastunit) )
        lin_rast_corr_ax.set_ylabel(r'Lineament density [m/km$^2$]')
        lin_rast_corr_ax.set_title(r'd=%g km, N=%d, $R^2$=%.4g' % (maxdist, N, lindensity_raster_corrcoef**2) )
        lin_rast_corr_ax.set_ylim(0,lindensity.max())
        lin_rast_corr_ax.set_xlim(0,randsamples.max())

        if save_fmt is not None:
            rastfig.savefig(os.path.join(figdir,outfile+'.'+save_fmt))
#}}} end LinDensityMap

def FitCurveExamples(lins, labels=[], outfile=None, dbar_max=0.125): #{{{
    # Create a full page plot, with the top half consisting of a map of the example lineaments,
    # with each one having a different color.  In the bottom half, show on individual subplots
    # the fitcurves for each of the features, color coded similarly.
    the_fig= figure(figsize=(9,12))

    colors = cm.jet(linspace(0,1,len(lins)))

    eg_ax1 = the_fig.add_axes((0.1,0.4625,0.4,0.1375))
    eg_ax2 = the_fig.add_axes((0.5,0.4625,0.4,0.1375))
    eg_ax3 = the_fig.add_axes((0.1,0.3250,0.4,0.1375))
    eg_ax4 = the_fig.add_axes((0.5,0.3250,0.4,0.1375))
    eg_ax5 = the_fig.add_axes((0.1,0.1875,0.4,0.1375))
    eg_ax6 = the_fig.add_axes((0.5,0.1875,0.4,0.1375))
    eg_ax7 = the_fig.add_axes((0.1,0.0500,0.4,0.1375))
    eg_ax8 = the_fig.add_axes((0.5,0.0500,0.4,0.1375))

    eg_axes = [eg_ax1,eg_ax2,eg_ax3,eg_ax4,eg_ax5,eg_ax6,eg_ax7,eg_ax8]

    # this is the map
    eg_map_ax = the_fig.add_axes((0.1,0.6,0.8,0.4))
    llcrnrlon=0
    llcrnrlat=0
    urcrnrlon=180
    urcrnrlat=90
    eg_map = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat, ax=eg_map_ax)
    eg_map.drawmeridians(range(llcrnrlon,urcrnrlon+1,15), labels=[1,0,0,1])
    eg_map.drawparallels(range(llcrnrlat,urcrnrlat+1,15), labels=[1,0,0,1])

    for lin,eg_ax,N in zip(lins,eg_axes,range(len(lins))):
        # Plot the lineament, and color it
        lineament.plotlinmap([lin,], map=eg_map, lon_cyc=radians(180), lat_mirror=True, color=colors[N], linewidth=2)
        fitcurve(lin, ax=eg_ax, color=colors[N], dbar_max=dbar_max)

    # clean up the massively multiple axes:
    ys_to_hide = [ the_fig.axes[N] for N in (1,3,5,7,9,11,13,15) ]
    [ ax.set_ylabel('') for ax in ys_to_hide ]
    [ setp(ax.get_yticklabels(),visible=False) for ax in ys_to_hide ]
    [ setp(ax.get_xticklabels(),visible=False) for ax in the_fig.axes[0:6] ]
    [ setp(ax.get_xticklabels(),visible=False) for ax in the_fig.axes[8:17] ]

    eg_map_ax.set_title(r'Example Lineaments and Fit Probabilities')

    if save_fmt is not None:
        the_fig.savefig(outfile)
#}}}
