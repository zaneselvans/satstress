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

###############################################################################
#  Helper functions for generating, sampling, saving, and loading lineaments  #
###############################################################################

def random_nsrlins(nsr_stresscalc=None, nlins=1000, minlen=0.0, maxlen=1.0): #{{{
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
        NSR = satstress.StressCalc([satstress.NSR(europa),])

    linlist = []

    while len(linlist) < nlins:
        max_length = minlen+(rand()*(maxlen-minlen))
        init_lon, init_lat = lineament.random_lonlatpoints(nlins)
        newlin = lineament.lingen_nsr(NSR, init_lon=init_lon, init_lat=init_lat, max_length=max_length, propagation='both')
        if newlin.length >= minlen:
            linlist.append(newlin)

    return linlist
#}}}

def random_gclins(nlins=1000, minlen=0.0, maxlen=1.0): #{{{
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
    return([ lineament.lingen_greatcircle(init_lons, init_lats, fin_lons, fin_lats, seg_len=0.01) for N in range(nlins) ])

#}}}

def regular_nsr_lingen(nlats=36): #{{{
    """
    Create a regularaly spaced "grid" of synthetic NSR lineaments, against
    which to compare the orientations of de-rotated mapped features - allowing
    an intuitive visualization of a lineament's delta_rms value.

    """
    satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    NSR = satstress.StressCalc([satstress.NSR(europa),])

    linlist = []
    for lat in linspace(-pi/2,pi/2,nlats+2):
        if lat != pi/2 and lat != -pi/2 and lat != 0:
            linlist.append(lineament.lingen_nsr(NSR, init_lon=0.0, init_lat=lat, max_length=2*pi, propagation='both'))
            linlist.append(lineament.lingen_nsr(NSR, init_lon=-pi, init_lat=lat, max_length=2*pi, propagation='east'))
            linlist.append(lineament.lingen_nsr(NSR, init_lon=pi, init_lat=lat, max_length=2*pi, propagation='west'))
    return linlist
#}}}

def linsample(proto=None, pool=None, nbins=20): #{{{
    """
    Select lineaments from pool that approximate the length distribution of
    those in proto.  Approximate the distribution using nbins bins, ranging
    from the minimum to the maximum length of th lineaments in proto.

    """

    # define the length bins to replicate:
    Ns, bins = histogram([ lin.length for lin in proto ], bins=nbins, new=True)

    # randomize lins2fit to avoid habitually selecting from the same pool:
    shuffle(pool)

    keepers = []
    numlins = len(proto)
    for N,i in zip(Ns, range(nbins)):
        lins_by_bin = []
        for lin in pool:
            linlen = lin.length
            if linlen > bins[i] and linlen < bins[i+1]:
                lins_by_bin.append(lin)

        keepers += lins_by_bin[:N]

    # avoid ordering the lineaments by length:
    shuffle(keepers)
    return(keepers)
#}}}

def save_lins(lins, name="global_lins"): #{{{
    """
    A shortcut for saving a new set of analyzed lineaments.
    """
    from time import strftime
    import pickle
    import os

    outdir="output"
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

    print("Loading synthetic lineament pools")
    gc_pool = load_lins("output/gc_pool")
    synth_pool = load_lins("output/synth_pool")

    print("Sampling synthetic lineament pools")
    gclins = linsample(proto=maplins, pool=gc_pool)
    synthlins = linsample(proto=maplins, pool=synth_pool)

    lins_list = (maplins, tpwlins, gclins, synthlins)
    labels = ("map_nsrfit", "tpw_nsrfit", "gc_nsrfit", "synth_nsrfit")

    # now calculate the fits...
    for lins,label in zip(lins_list, labels):
        print("Calculating fits for %s..." % (label,) )
        for lin,N in zip(lins,range(len(lins))):
            if (mod(N+1,10)==0):
                print("    N=%d" % (N+1,) )
            lin.calc_nsrfits(nb=nb, stresscalc=NSR)
        # only save if we're doing the whole dataset
        if nlins == 0:
            print("Saving %s" % (label,) )
            save_lins(lins, name=label)

    return(maplins, tpwlins, gclins, synthlins)
#}}}

def aggregate_fitmetric(lins, delta_max=radians(45), dbar_max=0.125): #{{{
    """
    An overall metric of how well a given set of features can be forced to fit
    the NSR stresses, if each feature is allowed to translate freely in
    longitude.

    """

    linlengths  = array([ lin.length for lin in lins ])
    linbestfits = array([ max(lin.nsrfits(delta_max=delta_max, dbar_max=dbar_max)) for lin in lins ])
    return(sum(linlengths*linbestfits)/sum(linlengths))

#}}}

def fitcurve(lin, delta_max=radians(45), dbar_max=0.125, color='black', ax=None): #{{{
    """
    Plot delta_rms as a function of b for lin, and show the value of the
    weighting function that results.

    """

    if ax is None:
        the_fig = figure(figsize=(9,6))
        delta_ax = the_fig.add_subplot(1,1,1)
    else:
        delta_ax = ax

    sat_radius_km = lin.stresscalc.stresses[0].satellite.radius()/1000

    # Plot the mapped lineament's fit curve:
    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    delta_ax.axis([0,179,0,100])
    delta_ax.invert_xaxis()
    delta_ax.set_yticks(unique(sort(hstack([arange(0,91,15),round(degrees(delta_max))]))))
    delta_ax.set_xticks(arange(0,179,15))
    delta_ax.set_xlabel("westward translation (b)")
    delta_ax.xaxis.set_major_formatter(degree_formatter)
    delta_ax.set_ylabel(r'$\delta_{rms}$')
    wt_ax = twinx(ax=delta_ax)
    wt_ax.set_ylabel('length weighting')
    wt_ax.set_yticks(linspace(0,0.9,7))
    delta_ax.yaxis.set_major_formatter(degree_formatter)
    delta_ax.plot(degrees(lin.bs), degrees(lin.nsrdeltas), ls='-', linewidth=2, color=color)

    wt_ax.plot(degrees(lin.bs), lin.nsrstresswts/2.0, ls=':', linewidth=2, color=color)
    wt_ax.plot(degrees(lin.bs), lin.nsrdbars*3.6, ls='--', linewidth=2, color=color)

    wt_ax.fill_between(degrees(lin.bs), lin.nsrfits(delta_max=delta_max, dbar_max=dbar_max), 0, color=color, alpha=0.4)
    delta_ax.set_ylim(0,100)
    wt_ax.set_ylim(0,1.0)

    if ax is None:
        the_fig.show()
#}}}

def activity_history(lins_list, the_fig=None, labels=[], colors=[], lin_cm=cm.jet, norm_by_all=False, outfile=None, delta_max=radians(45), dbar_max=0.125): #{{{
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

    print("delta_max = %.2g degrees" % (degrees(delta_max),))

    # Depending on whether we're segmenting a single dataset (say, by lineament
    # length) or trying to compare different datasets, we may want to norm each
    # set of lineaments by its own length, or the total length of all the
    # lineaments.
    norm_lengths = [ sum([lin.length for lin in lins]) for lins in lins_list ]
    if norm_by_all is True:
        norm_lengths = [sum(norm_lengths),] * len(lins_list)

    for lins,label,color,norm_length in zip(lins_list,labels,colors,norm_lengths):
        print("Calculating activity histories for %s" % (label,))
        acthist = array([ lin.length*lin.nsrfits(delta_max=delta_max, dbar_max=dbar_max) for lin in lins ]).sum(axis=0)/norm_length
        hist_ax.plot(degrees(lins[0].bs), acthist, label=label, lw=3, c=color)

    hist_ax.set_ylabel("retained length fraction")
    hist_ax.grid(True)
    hist_ax.set_xticks(linspace(0,170,18))
    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    hist_ax.xaxis.set_major_formatter(degree_formatter)
    hist_ax.set_xlim(179,0)
    hist_ax.legend(loc="upper left")
    hist_ax.set_title('Apparent Lineament Activity History')

    if outfile is not None:
        the_fig.savefig(outfile)
    else:
        the_fig.show()

#}}}

###############################################################################
#         Routines that generate plots in support of the TPW paper:           #
###############################################################################

def makefigs(delta_max=radians(45), dbar_max=0.125, maps=False, hists=False, examples=False, stats=False, stress=False, save_format=None): #{{{
    """
    A wrapper that takes input parameters for the figures, and calls all of the
    individual figure plotting routines, below.

    """
    from osgeo import gdal
    from os import path

    outdir = 'output'
    figure_names = ['FitMap_Mapped',\
                    'FitMap_PreTPW',\
                    'FitMap_Derotated',\
                    'FitMap_GreatCircle',\
                    'FitMap_SyntheticNSR',\
                    'ActHist_ByLength',\
                    'ActHist_BySin',\
                    'ActHist_Compare',\
                    'LinStats_DeltaLenCorr',\
                    'LinStats_DeltaSinCorr',\
                    'LinStats_Mapped',\
                    'LinStats_Derotated',\
                    'LinStats_GreatCircle',\
                    'LinStats_PreTPW',\
                    'LinStats_SyntheticNSR',\
                    'LinLengthDist_Mapped',\
                    'LinLengthDist_GreatCircle',\
                    'LinLengthDist_PreTPW',\
                    'LinDensity_Grid',\
                    'LinDensity_Resolution',\
                    'LinDensity_Emission',\
                    'LinDensity_Incidence',\
                    'LinDensity_Phase',\
                    'FitCurveExamples']

    figure_outfiles = {}

    for figname in figure_names:
        if save_format is None:
            figure_outfiles[figname] = None
        else:
            figure_outfiles[figname] = path.join(outdir,figname) + '.' + save_format

    # These are the lineaments in their modern locations, compared to an
    # elastic NSR stress field, with Delta ~ 0.001.
    nsrlins = load_lins(path.join(outdir,'map_nsrfit'))
    # All the mapped lineaments, transformed to a pre-TPW coordinate system,
    # with the paleo-north pole at 80E, 10N in current lat/lon.
    tpwlins = load_lins(path.join(outdir,'tpw_nsrfit'))
    # A set of 661 synthetic lineaments, each of which is a portion of a great
    # circle, oriented randomly on the surface, having a random location, in a
    # variety of lengths, allowing a very basic null hypothesis test.  They've
    # had their fits calculated relative to an elastic NSR stress field.
    gclins = load_lins(path.join(outdir,'gc_nsrfit'))
    # A set of 661 perfect synthetic NSR lineaments, as a control dataset:
    synlins = load_lins(path.join(outdir,'synth_nsrfit'))

    if maps is True: #{{{2
        print("Plotting Mapped Lineaments, fit to NSR stresses")
        FitMap(nsrlins, nbins=9, titlestr="global lins, fit to NSR", delta_max=delta_max, dbar_max=dbar_max, outfile=figure_outfiles['FitMap_Mapped'])
        print("Plotting Pre-TPW Lineaments, fit to NSR stresses")
        FitMap(tpwlins, nbins=9, titlestr="pre-TPW lins, fit to NSR", delta_max=delta_max, dbar_max=dbar_max, outfile=figure_outfiles['FitMap_PreTPW'])
        print("Plotting random great circle segments")
        FitMap(gclins, nbins=9, titlestr="Great Circle Segments fit to NSR", delta_max=delta_max, dbar_max=dbar_max, outfile=figure_outfiles['FitMap_GreatCircle'])
        print("Plotting synthetic NSR lineaments")
        FitMap(synlins, titlestr="Perfect Synthetic NSR Lineaments", delta_max=delta_max, dbar_max=dbar_max, outfile=figure_outfiles['FitMap_SyntheticNSR'])
        print("Plotting Mapped Lineaments, relative to stress field")
        FitMap(nsrlins, nbins=9, titlestr="global lins, in NSR formation locations", delta_max=delta_max, dbar_max=dbar_max, derotate=True, stresscentric=True, outfile=figure_outfiles['FitMap_Derotated'], showbad=False)
    #}}}2

    if hists is True: #{{{2
        ActHist_ByLength(nsrlins, delta_max=delta_max, dbar_max=dbar_max, outfile=figure_outfiles['ActHist_ByLength'], norm_by_all=True)
        ActHist_BySin(nsrlins, delta_max=delta_max, dbar_max=dbar_max, outfile=figure_outfiles['ActHist_BySin'], norm_by_all=True)
        ActHist_Compare(nsrlins, gclins, tpwlins, delta_max=delta_max, dbar_max=dbar_max, outfile=figure_outfiles['ActHist_Compare'], norm_by_all=False)
    #}}}2
 
    if examples is True: #{{{2
        # List of canonical features from the map: {{{3
        good_nsr1    = nsrlins[1]   # 1500 km long, symmetric arcuate lineament in the N. hemisphere
        good_nsr2    = nsrlins[53]  #  700 km long, asymmetric arcuate lineament in the N. hemisphere
        good_nsr3    = nsrlins[108] # 1800 km long, nearly symmetric arcuate lineament in the S. hemisphere
        bad_nsr1     = nsrlins[20]  #  553 km long, north-south lineament in the S. hemisphere
        bad_nsr2     = nsrlins[6]   # 1222 km long, diagonal, grazing the equator
        bad_nsr3     = nsrlins[36]  # 1175 km long, diagonal, equator crossing
        bad_nsr4     = nsrlins[54]  # 1036 km long, north-south, equator crossing
        bad_nsr5     = nsrlins[70]  # 1061 km long, one of the SCDs, crossing the equator
        bad_nsr6     = nsrlins[120] #  640 km, N-S equator crossing
        bad_nsr7     = nsrlins[122] # 1300 km, N-S equator crossing
        cycloid1     = nsrlins[112] # 1132 km long cycloid, 4-5 arcs, near 30S.  Low dbar, high delta
        cycloid2     = nsrlins[137] #  458 km long cycloid, 5 arcs, near 60S.  Low dbar, high delta
        cycloid3     = nsrlins[148] #  776 km long cycloid, 5 arcs, 45-70N.
        cycloid4     = nsrlins[155] # 1711 km long semi-cycloidal, actually passed the test
        cycloid5     = nsrlins[159] # 1527 km long semi-cycloidal, actually passed the test
        sinuous1     = nsrlins[23]  # 1334 km long, sinuous from 30N to 75N
        sinuous2     = nsrlins[136] # 1189 km long, sinuous from 30S to 70S
        short_hilat1 = nsrlins[11]  #  450 km, east-west, just above 30N... fits perfectly
        short_hilat2 = nsrlins[78]  #  200 km, diagonal, 50-60N, fits perfectly
        short_hilat3 = nsrlins[11]  #  183 km, east-west, just above 30N... fits perfectly
        short_hilat4 = nsrlins[160] #  200 km, diagonal, 75S... fits very well
        short_lolat1 = nsrlins[26]  #  500 km, diagonal, between 0N and 30N, fits perfectly
        short_lolat2 = nsrlins[41]  #  177 km, diagonal, between 0N and 30N, does not fit because of dbar
        short_lolat3 = nsrlins[43]  #  197 km, diagonal, between 0N and 30N, fits very well
        dbar_fail1   = nsrlins[7]   # 1073 km, arcuate, dbar_min doesn't quite line up with delta_min, and so it gets lost.
        dbar_fail2   = nsrlins[62]  #  500 km, east-west almost at the equator.  Large number of "good fits" possible.
        dbar_fail3   = nsrlins[63]  #  500 km, east-west almost at the equator.  Large number of "good fits" possible.
        dbar_fail4   = nsrlins[115] # 1262 km, looks good, but dbar just doesn't quite get low enough... Grr.
        #}}}3

        eg_lins = [ good_nsr1, bad_nsr2, cycloid1, sinuous1, sinuous2, short_hilat1, short_lolat2, dbar_fail2 ]

        eg_labels = [ "Good NSR", "Bad NSR", "Cycloial", "Sinuous", "Sinuous",\
                      "Short straight mid-latitude", "Short straight low-latitude",\
                      "Poorly constrained best fit" ]

        FitCurveExamples(eg_lins, labels=eg_labels, delta_max=delta_max, dbar_max=dbar_max, outfile=figure_outfiles['FitCurveExamples'])

    #}}}2

    if stats is True: #{{{2
        # Uncomment if we need to demonstrate the variation of orientations with latitude.
        #NSRFailDirByLat(nsrlins[0].stresscalc, lats=linspace(20,36,9))

        LinLengthDist(nsrlins, label="Mapped Lineament", outfile=figure_outfiles['LinLengthDist_Mapped'])
        LinLengthDist(gclins,  label="Great Circle Segment", outfile=figure_outfiles['LinLengthDist_GreatCircle'])
        LinLengthDist(tpwlins, label="pre-TPW (pole=80E10N) Lineament", outfile=figure_outfiles['LinLengthDist_PreTPW'])

        DeltaSinuosityCorr(nsrlins, outfile=figure_outfiles['LinStats_DeltaSinCorr'])
        DeltaLengthCorr(gclins, outfile=figure_outfiles['LinStats_DeltaLenCorr'])

        LinLatLonStats(nsrlins, label="Mapped Lineaments", outfile=figure_outfiles['LinStats_Mapped'])
        LinLatLonStats(gclins,  label="Great Circle Segments", outfile=figure_outfiles['LinStats_GreatCircle'])
        LinLatLonStats(tpwlins, label="pre-TPW Lineaments (pole=80E10N)", outfile=figure_outfiles['LinStats_PreTPW'])
        LinLatLonStats(synlins, label="Synthetic NSR Lineaments", outfile=figure_outfiles['LinStats_SyntheticNSR'])
        LinLatLonStats(nsrlins, label="Derotated Lineaments", outfile=figure_outfiles['LinStats_Derotated'])

        # This shows the density of the mapped features, and a bunch of information
        # about how that relates to the spacecraft observations.
        LinDensityMap(nsrlins, maxdist=250, label="Density of Mapped Features",\
                      grid_outfile       = figure_outfiles['LinDensity_Grid'],\
                      resolution_outfile = figure_outfiles['LinDensity_Resolution'],\
                      emission_outfile   = figure_outfiles['LinDensity_Emission'],\
                      incidence_outfile  = figure_outfiles['LinDensity_Incidence'],\
                      phase_outfile      = figure_outfiles['LinDensity_Phase'])

    #}}}2

    if stress is True: #{{{2
        # need to generate maps of the NSR stresses which are behind all of
        # this stuff... GridCalc here we come!
        pass

    #}}}2

# end makefigs }}}

def ActHist_ByLength(lins, delta_max=radians(45), dbar_max=0.125, outfile=None, norm_by_all=True): #{{{
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

    activity_history(lins_by_length, delta_max=delta_max, dbar_max=dbar_max, labels=labels_by_length, outfile=outfile, norm_by_all=norm_by_all)
    #}}}

def ActHist_BySin(lins, delta_max=radians(45), dbar_max=0.125, outfile=None, norm_by_all=True): #{{{
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

    activity_history(lins_by_sin, delta_max=delta_max, dbar_max=dbar_max, labels=labels_by_sin, outfile=outfile, norm_by_all=norm_by_all)
    #}}}

def ActHist_DeltaRMS(nsrlins, dbar_max=0.125, outfile=None, norm_by_all=False): #{{{
    """
    Show that the time variability is not sensitive to the particular value
    chosen for delta_max.

    """
    the_fig = figure(figsize=(12,8))

    deltas = (10,15,20,30,45)
    colors = cm.jet(linspace(0,1,len(deltas)))

    for delta_max,N in zip(deltas,range(len(deltas))):
        activity_history([nsrlins,], labels=[r'$\max(\delta_{rms})=%d^\circ$'%(delta_max,),], outfile=outfile, norm_by_all=norm_by_all, the_fig=the_fig, delta_max=radians(delta_max), dbar_max=dbar_max, colors=[colors[N],])
    #}}}

def ActHist_Compare(nsrlins, gclins, tpwlins, delta_max=radians(45), dbar_max=0.125, outfile=None, norm_by_all=False): #{{{
    # the global lineaments, with anything in the E15REGMAP01 observation
    # excluded: i.e. any feature with a portion of its length between
    # 250E and 290E longitude.
    noE15lins = [ lin for lin in nsrlins if (degrees(mod(array(lin.lons),2*pi)) > 290).all() or \
                                            (degrees(mod(array(lin.lons),2*pi)) < 240).all() ]

    activity_history([nsrlins, noE15lins, gclins, tpwlins], delta_max=delta_max, dbar_max=dbar_max,\
                     labels=['Mapped','Mapped Minus E15','Great Circles','Pre-TPW'], outfile=outfile, norm_by_all=norm_by_all)
    #}}}

def FitMap(lins, titlestr="Features colored by fit", lin_cm=cm.jet, nbins=9, stresscentric=False, outfile=None, delta_max=radians(45), dbar_max=0.125, showbad=True, derotate=False): #{{{
    """
    Creates a global map of the lineaments, color coding them by what amount of
    backrotation (b) minimizes delta_rms(b) when compared to the NSR stress
    field.  Lineament width indicates the value of min(delta_rms(b)), with
    wider lineaments agreeing better with NSR stresses.

    Those lineaments which did not meet the criteria for inclusion in the
    analysis are plotted thin and black.

    """

    # We can't both derotate and show the bad features, since bad features
    # don't necessarily have a best fit backrotation
    if derotate is True:
        showbad = False

    lin_cm = colors.ListedColormap(lin_cm(linspace(0,1,nbins)))
    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000

    fig = figure(figsize=(12,8))
    if stresscentric:
        llcrnrlon=0
        llcrnrlat=0
        urcrnrlon=180
        urcrnrlat=90
        lat_mirror=True
        gridspace=15
    else:
        llcrnrlon=0
        llcrnrlat=-90
        urcrnrlon=360
        urcrnrlat=90
        lat_mirror=False
        gridspace=30

    lon_cyc=abs(radians(llcrnrlon-urcrnrlon))
    linfitmap = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat)
    linfitmap.drawmeridians(range(llcrnrlon,urcrnrlon+1,gridspace), labels=[1,0,0,1])
    linfitmap.drawparallels(range(llcrnrlat,urcrnrlat+1,gridspace), labels=[1,0,0,1])
    linfitmap.drawmapboundary(fill_color="white")
    map_ax = fig.axes[0]

    for lin in lins:
        if len(lin.bs) > 0:
            backrot = 0
            # load some information about the lineament's best fit:
            nsrfits = lin.nsrfits(delta_max=delta_max, dbar_max=dbar_max, use_stress=False)
            best_fit = max(nsrfits)

            if best_fit > 0.0:
                best_b = float(lin.bs[where(fabs(nsrfits-best_fit) < 1e-9)[0]])
                if derotate:
                    backrot = best_b
                # Map the color of the lineament to its best_b
                lin_color = lin_cm(int((lin_cm.N)*(best_b/pi)))
                # use line width to indicate goodness of best_fit
                lin_width = 1.0 + 5.0*best_fit

            elif showbad:
                lin_width = 1.0
                lin_color = 'black'

            else:
                continue

            if backrot == 0:
                lin2plot = lin
            else:
                lin2plot = lin.lonshift(backrot)

            newline, linfitmap = lineament.plotlinmap([lin2plot,], map=linfitmap, linewidth=lin_width, color=lin_color, lon_cyc=lon_cyc, lat_mirror=lat_mirror)

    map_ax.set_title(titlestr + " " + 'max($\delta_{rms}$)=%.2g$^\circ$' % (degrees(delta_max),))
    cb_ax,kw = colorbar.make_axes(map_ax, orientation="horizontal", pad=0.05, shrink=0.5)
    colorbar.ColorbarBase(cb_ax, cmap=lin_cm, orientation="horizontal", norm=colors.BoundaryNorm(linspace(0,180,nbins+1),nbins), format=r'%.0f$^\circ$')
    # Fix up the colorbar a bit:
    cb_ax.invert_xaxis()
    cb_ax.set_xlabel("backrotation b")
    if outfile is None:
        fig.show()
    else:
        fig.savefig(outfile)
# end FitMap }}}

def DeltaLengthCorr(lins, outfile=None): #{{{
    """
    Several scatter plots showing the correlation (or lack thereof) between
    min(delta_rms) and lineament length, in equatorial and non-equatorial
    regions.

    """

    the_fig = figure(figsize=(12,8))
    ax1 = the_fig.add_subplot(1,1,1)

    hilatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) > 30 ]
    hilat_best_deltas = array([ lin.nsrdeltas.min() for lin in hilatlins ])
    hilat_lengths = array([ lin.length for lin in hilatlins ])*1561
    hilat_r2 = corrcoef(hilat_lengths, hilat_best_deltas)[0,1]**2
    hilat_symbs = scatter(hilat_lengths, degrees(hilat_best_deltas), c='blue', marker='s', label=r'$|\theta_{min}|>30^\circ R^2=%.3g$' % ( hilat_r2, ) )

    lolatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) <= 30 ]
    lolat_best_deltas = array([ lin.nsrdeltas.min() for lin in lolatlins ])
    lolat_lengths = array([ lin.length for lin in lolatlins ])*1561
    lolat_r2 = corrcoef(lolat_lengths,lolat_best_deltas)[0,1]**2
    lolat_symbs = scatter(lolat_lengths, degrees(lolat_best_deltas), c='red', marker='^', label=r'$\|\theta_{min}|\leq30^\circ R^2=%.3g$' % ( lolat_r2, ) )

    ax1.legend(loc='lower right')
    ax1.set_xlabel('lineament length [km]')
    ax1.set_ylabel(r'$\min(\delta_{rms}(b))$')
    ax1.set_ylim(0,30)
    ax1.set_xlim( min((hilat_lengths.min(),lolat_lengths.min())),\
                  max((hilat_lengths.max(),lolat_lengths.max())) )
    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    ax1.yaxis.set_major_formatter(degree_formatter)
    ax1.set_title('Effects of length and latitude on fit to NSR for great circle segments, N=%d'\
                   % (len(hilatlins)+len(lolatlins)) )

    if outfile is None:
        the_fig.show()
    else:
        the_fig.savefig(outfile)

# end DeltaLengthCorr}}}

def DeltaSinuosityCorr(lins, outfile=None): #{{{
    """
    Several scatter plots showing the correlation (or lack thereof) between
    min(delta_rms) and other lineament characteristics, including:

    A: lineament.sinuosity()
    B: min(abs(lineament.latitudes()))
    """

    the_fig = figure(figsize=(12,8))
    ax1 = the_fig.add_subplot(1,1,1)

    hilatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) > 30 ]
    hilat_best_deltas = array([ lin.nsrdeltas.min() for lin in hilatlins ])
    hilat_sins = array([ lin.sinuosity() for lin in hilatlins ])
    hilat_r2 = corrcoef(hilat_sins,hilat_best_deltas)[0,1]**2
    hilat_symbs = scatter(hilat_sins, degrees(hilat_best_deltas), c='blue', marker='s', label=r'$|\theta_{min}|>30^\circ R^2=%.3g$' % ( hilat_r2, ) )

    lolatlins = [ lin for lin in lins if degrees(min(fabs(lin.lats))) <= 30 ]
    lolat_best_deltas = array([ lin.nsrdeltas.min() for lin in lolatlins ])
    lolat_sins = array([ lin.sinuosity() for lin in lolatlins ])
    lolat_r2 = corrcoef(lolat_sins,lolat_best_deltas)[0,1]**2
    lolat_symbs = scatter(lolat_sins, degrees(lolat_best_deltas), c='red', marker='^', label=r'$\|\theta_{min}|\leq30^\circ R^2=%.3g$' % ( lolat_r2, ) )

    ax1.legend(loc='lower right')
    ax1.set_xlabel('sinuosity')
    ax1.set_ylabel(r'$\min(\delta_{rms}(b))$')
    ax1.set_ylim(0,30)
    ax1.set_xlim(1,1.10)
    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    ax1.yaxis.set_major_formatter(degree_formatter)
    ax1.set_title('Effects of sinuosity and latitude on fit to NSR for all mapped features, N=%d'\
                   % (len(hilatlins)+len(lolatlins)) )

    if outfile is None:
        the_fig.show()
    else:
        the_fig.savefig(outfile)

# end DeltaSinuosityCorr }}}

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

    delta_max = max([ max(lin.good_fits()[:,1]) for lin in goodlins if lin.good_fits().any() ])
    dbar_max = max([ max(lin.good_fits()[:,2]) for lin in goodlins if lin.good_fits().any() ])

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
    linkept = [ l.good_fits().any() for l in lins ]

    lengthskept = array([ (linlengths[n],linkept[n]) for n in range(len(lins)) ], dtype=[('length',float),('kept',bool)])
    lengthskept.sort(order='length')

    the_fig = figure(figsize=(12,4))
    plot_ax = the_fig.add_subplot(1,1,1)
    plot_ax.plot(lengthskept['length'], color='black', linewidth=2)
    plot_ax.fill_between(range(len(lengthskept)), zeros(len(lengthskept)), lengthskept['length'], where=lengthskept['kept'], facecolor='black')
    plot_ax.set_title(label + ' Length Distribution for ' + 'max($\delta_{rms}$)=%.2g$^\circ$, max($\\bar{D}$)=%.2g' % (degrees(delta_max),dbar_max))
    plot_ax.set_xlabel('N')
    plot_ax.set_ylabel('lineament length [km]')

    if outfile is None:
        the_fig.show()
    else:
        the_fig.savefig(outfile)

#}}}

def LinDensityMap(lins, maxdist=500, N=0, label="", cmap=cm.jet,\
                  grid_outfile=None, resolution_outfile=None,\
                  emission_outfile=None, incidence_outfile=None,\
                  phase_outfile=None): #{{{
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

    randlons, randlats = lineament.random_lonlatpoints(N)
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
    lindens_ax.contourf(reglons, reglats, lindensity_grid, 64, cmap=cmap)
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
    if grid_outfile is None:
        lindens_fig.show()
    else:
       lindens_fig.savefig(grid_outfile)

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
    rastoutfiles = [ resolution_outfile, emission_outfile, incidence_outfile, phase_outfile ]

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
        print("lindensity v. %s: R=%g" % (rastname, lindensity_raster_corrcoef) )

        # Make a scatter plot showing the correlation (or lack thereof):
        lin_rast_corr_ax = rastfig.add_subplot(2,1,2)
        lin_rast_corr_ax.scatter(randsamples, lindensity, s=10, linewidth=0, marker='o', color='black', alpha=0.375)
        lin_rast_corr_ax.set_xlabel("USGS Mosaic %s, %s" % (rastname, rastunit) )
        lin_rast_corr_ax.set_ylabel(r'Lineament density [m/km$^2$]')
        lin_rast_corr_ax.set_title(r'd=%g km, N=%d, r=%.4g' % (maxdist, N, lindensity_raster_corrcoef) )
        lin_rast_corr_ax.set_ylim(0,lindensity.max())
        lin_rast_corr_ax.set_xlim(0,randsamples.max())

        if outfile is None:
            rastfig.show()
        else:
            rastfig.savefig(outfile)
#}}}

def FitCurveExamples(lins, labels=[], delta_max=radians(45), dbar_max=0.125, outfile=None): #{{{
    # Create a full page plot, with the top half consisting of a map of the example lineaments,
    # with each one having a different color.  In the bottom half, show on individual subplots
    # the fitcurves for each of the features, color coded similarly.
    the_fig= figure(figsize=(12,15))

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
    eg_map.drawmapboundary(fill_color="white")

    for lin,eg_ax,N in zip(lins,eg_axes,range(len(lins))):
        # Plot the lineament, and color it
        lineament.plotlinmap([lin,], map=eg_map, lon_cyc=radians(180), lat_mirror=True, color=colors[N], linewidth=2)
        fitcurve(lin, ax=eg_ax, delta_max=delta_max, dbar_max=dbar_max, color=colors[N])

    # clean up the massively multiple axes:
    ys_to_hide = [ the_fig.axes[N] for N in (1,3,5,7,9,11,13,15) ]
    [ ax.set_ylabel('') for ax in ys_to_hide ]
    [ setp(ax.get_yticklabels(),visible=False) for ax in ys_to_hide ]
    [ setp(ax.get_xticklabels(),visible=False) for ax in the_fig.axes[0:6] ]
    [ setp(ax.get_xticklabels(),visible=False) for ax in the_fig.axes[8:17] ]

    eg_map_ax.set_title(r'Example Lineaments and Fit Curves for $\delta_{rms}<%d^\circ$'%(round(degrees(delta_max)),))

    if outfile is None:
        the_fig.show()
    else:
        the_fig.savefig(outfile)
#}}}
