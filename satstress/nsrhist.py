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

def nsrfits(nb=180, nlins=0, save=100, name="linfits", nbins=20,\
            satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
            mapped_lins=None, synthetic_lins=None,\
            pnp_lon=None, pnp_lat=None): #{{{
    """
    Calculate the fit curves for a set of lineaments.

    nb controls the resolution of the fit curves.  It is the number of b values
    for which fits are calculated.  They are spread evenly throughout 180
    degrees of longitudinal translation.  The default, 180, is one sample per
    degree.

    nlins controls what proportion of the features to calculate fits for.  If
    nlins=0, all lineaments are used, otherwise, only the first nlins of them
    are used.  Mostly useful for testing before starting a long run.

    save defines how many fits should be calculated before output is done.  If
    save is 0, then no output is done.  Otherwise, for each save lineaments, a
    new pickle file is dumped into the directory created using the name and
    datestamp provided.

    name defines a useful identifying prefix for the saved file.

    nbins is the number of bins to divide the lengths of the prototype
    lineaments into for replication.

    satfile is a SatStress runfile, defining the satellite on which the features
    reside, and what stresses they are subject to.

    mapped_lins and synthetic lins can both be either lists of Lineament
    objects or a file containing a pickled list of lineaments.

    pnp_lon, pnp_lat define the location in current coordinates of a paleo-
    north pole, which should be used to transform the proffered lineaments to
    pre-polar wander locations before calculating the fits.

    At least one of mapped_lins and synthetic_lins must be defined.  If
    synthetic_lins is None, calculate the first for the mapped_lins.  If both
    mapped_lins and synthetic_lins are defined, use the mapped_lins as a
    prototype dataset, and subsample synthetic_lins so as to statistically
    resemble mapped_lins, before calculating fits for synthetic_lins.  If only
    synthetic_lins is not None, calculate fits for all of them.
    
    """
    import pickle
    from time import strftime

    europa = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    # If both sets of lineaments are defined, we use the mapped as prototypes,
    # and fit the synthetics:
    if mapped_lins is not None and synthetic_lins is not None:
        proto_lins = mapped_lins
        lins2fit   = synthetic_lins

    # if we only have mapped_lins, then we fit them:
    elif mapped_lins is not None:
        proto_lins = None
        lins2fit = mapped_lins

    # if we only have synthetic_lins, fit them all:
    elif synthetic_lins is not None:
        proto_lins = None
        lins2fit = synthetic_lins

    # otherwise, both are None, and we can't do anything:
    else:
        return(False)

    # Now we need to turn proto_lins and lins2fit into actual lists of
    # Lineament objects, reading them in from their files if need be.
    if type(lins2fit) is str:
        lins2fit = pickle.load(open(lins2fit))
    if type(proto_lins) is str:
        proto_lins = pickle.load(open(proto_lins))

    # if we've been told to truncate the input for testing...
    if nlins > 0:
        if proto_lins is None:
            # here we shorten the list of lineament we're fitting...
            lins2fit = lins2fit[:nlins]
        else:
            # here we shorten the prototypes, which thereby shortens
            # the list of lineaments we ultimately end up fitting.
            proto_lins = proto_lins[:nlins]

    # if we have proto_lins we need to subsample lins2fit
    if proto_lins is not None:
        lins2fit = linsample(proto=proto_lins, pool=lins2fit, nbins=nbins)

    # If we need to transform for a paleopole, do so here:
    if(pnp_lon is not None and pnp_lat is not None):
        lins2fit = [ lin.poleshift(radians(pnp_lon), radians(pnp_lat)) for lin in lins2fit ]

    # If save is non-zero, we need to create a new and uniquely named directory
    # into which we will save the lineaments and their fits, in several bunches as
    # the calculations are completed.
    if save > 0:
        # create the uniquely named directory
        savedir = "output/%s_%s" % (name, strftime('%Y%m%d%H%M%S')) 
        print "Creating output directory", savedir, "\n"
        os.mkdir(savedir)
        partcount=1

    numlins = len(lins2fit)
    # Now we calculate the fits:
    for n,lin in zip(range(numlins),lins2fit):
        print("    Fitting lin %d / %d from %s" % (n+1,numlins,name))
        lin.stresscalc = nsr_stresscalc
        fits = lin.fit(bs=linspace(0,pi,nb, endpoint=False))
        # If we're going to be saving our results, do it periodically so we don't 
        # lose a bunch of work when something goes wrong elsewhere, and so we can
        # play with the results as they come out.
        if save > 0 and ((mod(n+1,save) == 0) or (n == numlins-1)):
            savepartfilename = savedir+"/%s_part_%d_of_%d.pkl" % (strftime('%Y%m%d%H%M%S'), partcount, numlins/save + 1)
            print "Saving lins %d-%d of %d in %s\n" % ((partcount-1)*save+1,n+1,numlins,savepartfilename)
            if mod(n+1,save) == 0:
                pickle.dump(lins2fit[(partcount-1)*save:partcount*save],open(savepartfilename,'w'))
            else:
                pickle.dump(lins2fit[(partcount-1)*save:],open(savepartfilename,'w'))
            partcount = partcount + 1

    return lins2fit
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

def browse_fits(lins, delta_max=radians(20), dbar_max=0.05): #{{{
    lin_num=0

    the_fig = figure()

    stresscalc = lins[0].stresscalc
    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000

    print("""Caching good fits for modern mapped features:
    delta_max = %.2g degrees
    dbar_max  = %.2g""" % (degrees(delta_max), dbar_max))
    for lin in lins:
        lin.cache_good_fits(delta_max=delta_max, dbar_max=dbar_max)

    for lin in lins:
        # Draw mapped lineament in a Basemap map:
        map_ax = the_fig.add_subplot(2,1,1)
        lin_length_km = sat_radius_km*lin.length
        linmap = Basemap(rsphere=sat_radius_km*1000, ax=map_ax, projection='cyl')
        map_ax.set_title("lin=%d, length=%g km" % (lin_num, lin_length_km) )
        linmap.drawmeridians(range(-180,180,30), labels=[1,0,0,1])
        linmap.drawparallels(range(-90,90,30), labels=[1,0,0,1])

        lin.plot(map=linmap, linewidth=2.0)

        good_fits = lin._Lineament__good_fits
        b_vals     = lin.bs()
        delta_vals = degrees(lin.delta_rms(bs=b_vals))
        dbar_vals  = lin.dbar(bs=b_vals)
        b_vals = degrees(b_vals)

        # Plot the mapped lineament's fit curve:
        delta_ax = the_fig.add_subplot(2,1,2)
        dbar_ax = twinx()

        degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')

        delta_ax.set_title(r'fit curves (solid=$\delta_{rms}$ and dashed=$\bar{D}$)')
        delta_ax.axis([0,180,0,90])
        delta_ax.invert_xaxis()
        delta_ax.set_yticks(unique(sort(hstack([arange(0,91,15),round(degrees(delta_max))]))))
        delta_ax.set_xticks(arange(0,180,15))
        delta_ax.set_xlabel("westward translation (b)")
        delta_ax.xaxis.set_major_formatter(degree_formatter)
        delta_ax.set_ylabel(r'$\delta_{rms}$')
        delta_ax.yaxis.set_major_formatter(degree_formatter)
        delta_ax.plot(b_vals, delta_vals, c='black', ls='-', lw=2)
        delta_ax.plot(b_vals, [degrees(delta_max),]*len(b_vals), color='green', lw=2)
        delta_ax.fill_between(b_vals, 90, 0, where=delta_vals<=degrees(delta_max), color='green', alpha=0.5)
        delta_ax.set_ylim(0,90)

        dbar_ax.axis([0,180,0,0.6])
        dbar_ax.invert_xaxis()
        dbar_ax.plot(b_vals, dbar_vals, c='black', ls='--', lw=2)
        dbar_ax.plot(b_vals, [dbar_max,]*len(b_vals), color='green', lw=2, ls='--')
        dbar_ax.fill_between(b_vals, 0.6, 0, where=dbar_vals<=dbar_max, color='green', alpha=0.5)
        dbar_ax.set_ylabel(r'$\bar{D}$')
        dbar_ax.set_yticks(unique(sort(hstack([arange(0,0.61,0.1),round(dbar_max,2)]))))
        dbar_ax.set_ylim(0,0.6001)
        dbar_ax.grid(True)

        the_fig.show()

        print("""\nLineament #%d:
=========================
length    = %g km
delta_min = %.2g deg (< %.2g required)
dbar_min  = %.3g (< %.2g required)
sinuosity = %g
lon_range = (%g, %g)
lat_range = (%g, %g)
fit_corr  = %g
good_fits: %d / %d (%.2g%%)""" % (lin_num,\
                                 lin_length_km,\
                                 (min(delta_vals)), delta_max,\
                                 min(dbar_vals), dbar_max,\
                                 lin.sinuosity(),\
                                 degrees(min(lin.longitudes())), degrees(max(lin.longitudes())),\
                                 degrees(min(lin.latitudes())),  degrees(max(lin.latitudes())),\
                                 corrcoef(delta_vals, dbar_vals)[0][1],\
                                 len(good_fits), len(b_vals), 100*float(len(good_fits))/float(len(b_vals)) ))
        for fit in good_fits:
            print("    [ %g, %g, %g ]" % (degrees(fit[0]), degrees(fit[1]), fit[2]))
        if len(good_fits) == 0:
            print("    None")
        print("\nclick map for next fit")
        x = ginput(n=1, timeout=-1, show_clicks=True)
        clf()
        lin_num += 1
#}}}

def calc_fit_hist(lins): #{{{
    """
    For a list of Lineaments, lins, generate a histogram describing the number
    of kilometers of lineament fitting at each value of b (0-180).

    Assumes that the lineaments have already had their fits calculated, and
    their good_fits cached, according to some predefined values for
    delta_max and dbar_max.

    """

    linhist = []
    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000

    tot_lin_length_km = 0
    fit_rate_by_N = 0
    fit_rate_by_km = 0

    for lin in lins:
        # Select the fits which meet our criteria for this lineament:
        good_fits = lin.good_fits()
        # if there weren't any good fits, skip to the next lineament:
        if len(good_fits) == 0:
            continue

        lin_length_km = lin.length*sat_radius_km
        fit_rate = float(len(good_fits)) / float(len(lin.bs()))
        tot_lin_length_km  += lin_length_km
        fit_rate_by_N  += fit_rate
        fit_rate_by_km += fit_rate * lin_length_km

        # We generate the histogram dataset, adding one "count" for each km of
        # lineament that fit best at that value of b(ackrotation).  Instead of
        # just looking at the local minima within the fit curves though, we 
        # want to include *all* of the b values for which the two fit metrics
        # meet our criteria, dividing up the 'length' of the lineament among all
        # the b values for which the fit is good.  This will smooth the histogram,
        # and more accurately reflect the uncertainty in the historical activity.
        for b, delta, dbar in good_fits:
            # we need to turn this into an integer in order to use it to generate
            # histogram counts by multiplying the [b,] array.  It shouldn't create
            # much truncation error though, because the length is in meters, and
            # the lineaments are often millions of meters long...  Actually, this
            # might be overkill.  We'll see
            split_length_km = int(lin_length_km/len(good_fits))
            linhist += split_length_km * [b,]

    fit_rate_by_N  /= len(lins)
    fit_rate_by_km /= tot_lin_length_km

    print("    fit rate by N:  %.2g%%\n    fit rate by km: %.2g%%" % ( 100.0*fit_rate_by_N, 100.0*fit_rate_by_km) )

    return(linhist)
# }}}

def save_lins(lins, name="global_lins"): #{{{
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

def activity_histogram(lins_list, hist_ax=None, nbins=18, labels=[], colors=[], lin_cm=cm.jet): #{{{
    """
    Plots activity histograms for several sets of lineaments, allowing visual
    comparison.

    lins_list is a list of lists of Lineament objects, assumed to have their
    fits already calculated.

    hist_ax is the matplotlib axes in which the histograms should be drawn.  If
    None, the current axes is used.

    nbins is the number of bins to divide the histograms into.  Defaults to 18,
    each 10 degrees wide.

    labels and colors are lists of strings and matplotlib color specifications
    that will be used to in legend creation and drawing.

    """

    # Set hist_ax to the current axes, if none was supplied.
    if hist_ax is None:
        hist_ax = gca()

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

    # Check to see what the required values of dbar_max and delta_max are
    delta_max = max( [ max([ max(lin._Lineament__good_fits[:,1]) for lin in lins if len(lin._Lineament__good_fits) > 0 ]) for lins in lins_list ])
    dbar_max = max( [ max([ max(lin._Lineament__good_fits[:,2]) for lin in lins if len(lin._Lineament__good_fits) > 0 ]) for lins in lins_list ])

    print("delta_max = %.2g degrees" % (degrees(delta_max),))
    print("dbar_max  = %.2g\n" % (dbar_max,))

    # Calculate the histogram data:
    lins_histdata = []
    for lins,label in zip(lins_list,labels):
        print("Generating activity histograms for %s" % (label,))
        lins_histdata.append(calc_fit_hist(lins))

    # Alas, need to special case this, because multi-histograms behave a
    # little differently than normal ones.
    if len(lins_histdata) == 1:
        lins_n, lin_bins, lins_patches = hist_ax.hist( degrees(lins_histdata).transpose(), bins=nbins, range=(0,180), rwidth=0.8 )
        lins_n = [ lins_n, ]
        lins_patches = [ lins_patches, ]
    else:
        assert(len(lins_histdata) > 1)
        lins_n, lin_bins, lins_patches = hist_ax.hist( [ degrees(histdata) for histdata in lins_histdata ], bins=nbins, range=(0,180) )

    # Assign the colors and labels to the histogram bars:
    for label,color,patches in zip(labels,colors,lins_patches):
        patches[0].set_label(label)
        for patch in patches:
            patch.set_facecolor(color)

    hist_ax.set_ylabel("cumulative feature length [km]")
    hist_ax.grid(True)
    hist_ax.set_xticks(linspace(0,180,nbins+1))
    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    hist_ax.xaxis.set_major_formatter(degree_formatter)
    
    hist_ax.invert_xaxis()
    hist_ax.legend(loc="upper left")
    for textobj in hist_ax.get_legend().get_texts():
        textobj.set_size('x-small')

    hist_ax.set_title('Lineament Formation Activity, max($\delta_{rms}$)=%.2g$^\circ$, max($\\bar{D}$)=%.2g' % (degrees(delta_max),dbar_max))

#}}}

def fitcurves(lins_list, fit_ax=None, labels=[], colors=[], lin_cm=cm.jet): #{{{
    """
    Plots delta_rms(b) and dbar(b) for several sets of lineaments.  If only one
    set containing only a single lineament is supplied, this reduces to simply
    plotting that feature's fit curves.

    lins_list is a list of lists of Lineament objects, assumed to have their
    fits calculated already.

    fit_ax is the matplotlib axes in which the curves are plotted.

    labels and colors are lists of strings and matplotlib color specifications
    that will be used to in legend creation and drawing.

    """

    if fit_ax is None:
        fit_ax = gca()

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

    # Calculate the aggregate fit curves:
    lins_agg_delta_rms = []
    lins_agg_delta_rms_derotated = []
    lins_agg_length = []

    for lins,label in zip(lins_list,labels):
        print("Generating fit curves for %s" % (label,))
        agg_delta_rms = zeros(len(lins[0].bs()))
        agg_delta_rms_derotated = zeros(len(lins[0].bs()))
        agg_length = 0
        for lin,N in zip(lins,range(len(lins))):
            linlength = lin.length
            agg_length += linlength
            # de-rotation can only work if we know how much to de-rotate a lineament...
            # which means having a best_fit, which means having good_fits.  If we don't
            # have any fits, add zeros instead.
            if len(lin._Lineament__good_fits > 0):
                # Find the index corresponding to the minimum delta_rms, while
                # still satisfying the dbar_max requirement:
                best_idx = mod(int( find(lin._Lineament__fits[:,1]==lin.best_fit()[1])), len(lin._Lineament__fits[:,0]))

                # assemble an array containing the same values as the delta_rms
                # portion of the fits array, but shifted by best_idx positions:
                derotated_fits = hstack([lin._Lineament__fits[best_idx:,1],lin._Lineament__fits[:best_idx,1]])
                agg_delta_rms_derotated += derotated_fits*linlength

            fits_to_add = lin._Lineament__fits[:,1]
            agg_delta_rms += fits_to_add*linlength

        agg_delta_rms /= agg_length
        lins_agg_delta_rms.append(agg_delta_rms)

        agg_delta_rms_derotated /= agg_length
        lins_agg_delta_rms_derotated.append(agg_delta_rms_derotated)

        lins_agg_length.append(agg_length)

    for lins,delta_rms,label,color in zip(lins_list,lins_agg_delta_rms,labels,colors):
        nb = len(lins[0].bs())
        fit_ax.plot(linspace(0,180,nb), degrees(delta_rms), color=color, linewidth=2, linestyle='-', label=label)

    for lins,delta_rms,label,color in zip(lins_list,lins_agg_delta_rms_derotated,labels,colors):
        nb = len(lins[0].bs())
        fit_ax.plot(linspace(0,180,nb), degrees(delta_rms), color=color, linewidth=2, linestyle='--', label="")

    fit_ax.set_ylabel(r'mean $\delta_{rms}$')
    fit_ax.set_xlabel("westward translation (b)")

    leg_x = arange(10)
    agg_line   = Line2D(leg_x, leg_x, linewidth=2, color='black', linestyle='-')
    derot_line = Line2D(leg_x, leg_x, linewidth=2, color='black', linestyle='--')
    
    fit_ax.legend([agg_line, derot_line], ('geographic', 'tide-centric'), loc='lower center')

    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    fit_ax.xaxis.set_major_formatter(degree_formatter)
    fit_ax.yaxis.set_major_formatter(degree_formatter)

    fit_ax.set_xticks(arange(0,181,10))
    fit_ax.invert_xaxis()
    fit_ax.grid(True)
    fit_ax.set_ylim(0,90)

    delta_max = max( [ max([ max(lin._Lineament__good_fits[:,1]) for lin in lins ]) for lins in lins_list ])
    dbar_max  = max( [ max([ max(lin._Lineament__good_fits[:,2]) for lin in lins ]) for lins in lins_list ])
    fit_ax.set_title('aggregate fit curves, max($\delta_{rms}$)=%.2g$^\circ$, max($\\bar{D}$)=%.2g' % (degrees(delta_max),dbar_max))
#}}}

def regular_nsr_lingen(nlats=36): #{{{
    """
    Create a regularaly spaced "grid" of synthetic NSR lineaments, against
    which to compare the orientations of de-rotated mapped features - allowing
    an intuitive visualization of a lineament's delta_rms value.

    """
    satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    linlist = []
    for lat in linspace(-pi/2,pi/2,nlats+2):
        if lat != pi/2 and lat != -pi/2 and lat != 0:
            linlist.append(lineament.lingen(nsr_stresscalc, init_lon=0.0, init_lat=lat, max_length=2*pi, propagation='both'))
            linlist.append(lineament.lingen(nsr_stresscalc, init_lon=-pi, init_lat=lat, max_length=2*pi, propagation='east'))
            linlist.append(lineament.lingen(nsr_stresscalc, init_lon=pi, init_lat=lat, max_length=2*pi, propagation='west'))
    return linlist
    #}}}

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
        nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    linlist = []

    while len(linlist) < nlins:
        # Choose a target length, and divide by two... since the "both" direction
        # propagation applies the target length to each half.
        max_length = (minlen+(rand()*(maxlen-minlen)))/2.0
        newlin = lineament.lingen(nsr_stresscalc, max_length=max_length, propagation='both')
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
    initial_points = array(zip(2*pi*rand(nlins), arccos(2*rand(nlins)-1)-pi/2))
    # - pick a random azimuth, az
    azimuths = 2*pi*rand(nlins)
    # - pick a random length L between minlen and maxlen
    linlengths = (minlen+(rand(nlins)*(maxlen-minlen)))
    # - calculate the location of endpoint v2, L radians away
    #   from v1, along the initial heading az, from v1.
    final_points = [ lineament.spherical_reckon(initial_points[N][0], initial_points[N][1], azimuths[N], linlengths[N]) for N in range(nlins) ]
    # - use lingen_greatcircle() to calculate intermediary vertices.

    return([ lineament.lingen_greatcircle(initial_points[N][0], initial_points[N][1],\
                                          final_points[N][0],   final_points[N][1], seg_len=0.01) for N in range(nlins) ])

#}}}

def w_stress_map(nlats=90, nlons=180): #{{{
    """
    Plot a grid showing the value of w_stress for an NSR stress field:

    """

    satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    mgsd = nsr_stresscalc.mean_global_stressdiff()

    stress_weights = zeros((nlats,nlons))

    for N,latitude in zip(range(nlats),linspace(-pi/2,pi/2,nlats)):
        for M,longitude in zip(range(nlons),linspace(-pi,pi,nlons)):
            (tens_mag, tens_az, comp_mag, comp_az) = nsr_stresscalc.principal_components(theta = (pi/2.0)-latitude,\
                                                                                           phi = longitude,\
                                                                                             t = 0 )
            stress_weights[N,M] = (tens_mag - comp_mag)/mgsd

    the_fig = figure()
    plot_ax = the_fig.add_subplot(1,1,1)
    plot_ax.imshow(stress_weights, extent=[-180,180,-90,90], origin='lower')
    plot_ax.set_xticks(linspace(-180,180,9))
    plot_ax.set_yticks(linspace(-90,90,5))
    cb_ax, kw = colorbar.make_axes(plot_ax, orientation="horizontal", pad=0.05)
    colorbar.ColorbarBase(cb_ax, norm=colors.Normalize(vmin=0.0, vmax=stress_weights.max()), orientation="horizontal")
    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    plot_ax.yaxis.set_major_formatter(degree_formatter)
    plot_ax.xaxis.set_major_formatter(degree_formatter)
    cb_ax.set_xlabel(r'$w_{stress}$')
    plot_ax.set_title(r'global map of $w_{stress}$ values for the NSR stress field')
    return(stress_weights)
#}}}

def cache_linstats(lins, name="", delta_max=radians(20), dbar_max=0.10): #{{{
    """
    Segregate a list of lineaments by whether or not they've got fits of the
    specified goodness, and also calculate some statistics on the population,
    to be used later in plotting.

    """

    print("""Caching good fits for %s features:
    delta_max = %.2g degrees
    dbar_max  = %.2g""" % (name, degrees(delta_max), dbar_max))

    stresscalc = lins[0].stresscalc
    sat_radius_km = stresscalc.stresses[0].satellite.radius()/1000.0

    for lin in lins:
        lin.cache_good_fits(delta_max=delta_max, dbar_max=dbar_max)

    # How many of the lineaments pass the fit test, by number:
    goodlins = [ lin for lin in lins if len(lin.good_fits()) > 0 ]
    badlins  = [ lin for lin in lins if lin not in goodlins ]

    # What proportion of their overall length passes the fit tests:
    lins_lengths     = array([ lin.length for lin in lins ])
    goodlins_lengths = sum( [lin.length for lin in goodlins ]) * sat_radius_km
    badlins_lengths  = sum( [lin.length for lin in badlins  ]) * sat_radius_km

    # How unique are the fits for these values of delta_max and dbar_max?
    # The fit degeneracy 
    fit_degeneracy   = array([ float(len(lin.good_fits()))/float(len(lin._Lineament__fits)) for lin in lins ])
    mean_fdn = mean( [ x for x in fit_degeneracy if x > 0 ])

    print("    N kept: %d / %d (%.0f%%)" % (len(goodlins), len(lins), 100.0*len(goodlins)/float(len(lins))) )
    print("    km kept: %d / %d (%.0f%%)" % (goodlins_lengths, badlins_lengths+goodlins_lengths, 100.0*goodlins_lengths/float(goodlins_lengths+badlins_lengths)) )
    print("    mean fit degeneracy:  %.2f%% " % (100.0*mean_fdn) )

    return(goodlins, badlins, goodlins_lengths, badlins_lengths, fit_degeneracy)

#}}}

###############################################################################
#         Routines that generate plots in support of the TPW paper:           #
###############################################################################

def makefigs(delta_max=radians(20), dbar_max=0.10, faildir=False, maps=False, hists=False, examples=False, stats=False): #{{{
    """
    A wrapper that takes input parameters for the figures, and calls all of the
    individual figure plotting routines, below.

    """
    # These are the lineaments in their modern locations, compared to an
    # elastic NSR stress field, with Delta ~ 0.001.
    nsrlins = load_lins("output/map_nsrfit")
    # All the mapped lineaments, transformed to a pre-TPW coordinate system,
    # with the paleo-north pole at 80E, 10N in current lat/lon.
    tpwlins = load_lins("output/tpw_nsrfit")
    # A set of 661 synthetic lineaments, each of which is a portion of a great
    # circle, oriented randomly on the surface, having a random location, in a
    # variety of lengths, allowing a very basic null hypothesis test.  They've
    # had their fits calculated relative to an elastic NSR stress field.
    gclins = load_lins("output/gc_nsrfit")
    # A set of 661 perfect synthetic NSR lineaments, as a control dataset:
    synlins = load_lins("output/synth_nsrfit")

    # All the mapped lineaments, compared to various NSR stress fields, for
    # different amounts of viscous relaxation (Delta=[0.1,1,10,100,1000])
    #nsrD01lins   = load_lins("output/Delta01_nsrfit")
    #nsrD1lins    = load_lins("output/Delta1_nsrfit")
    #nsrD10lins   = load_lins("output/Delta10_nsrfit")
    #nsrD100lins  = load_lins("output/Delta100_nsrfit")
    #nsrD1000lins = load_lins("output/Delta1000_nsrfit")

    # The implicit assumption here is that we've got the same StressCalc for all of the
    # lists of lineaments... should really check to make sure that it's true:
    nsr_stresscalc = nsrlins[0].stresscalc
    sat_radius_km = nsr_stresscalc.stresses[0].satellite.radius()/1000.0

    nsr_goodlins, nsr_badlins, nsr_goodlins_lengths, nsr_badlins_lengths, nsr_fitprop = cache_linstats(nsrlins, name="Modern Mapped",  delta_max=delta_max, dbar_max=dbar_max)
    tpw_goodlins, tpw_badlins, tpw_goodlins_lengths, tpw_badlins_lengths, tpw_fitprop = cache_linstats(tpwlins, name="Pre-TPW 80E10N", delta_max=delta_max, dbar_max=dbar_max)
    gc_goodlins,   gc_badlins,  gc_goodlins_lengths,  gc_badlins_lengths,  gc_fitprop = cache_linstats(gclins,  name="Great Circle",   delta_max=delta_max, dbar_max=dbar_max)
    syn_goodlins, syn_badlins, syn_goodlins_lengths, syn_badlins_lengths, syn_fitprop = cache_linstats(synlins, name="Synthetic NSR",  delta_max=delta_max, dbar_max=dbar_max)

    if faildir is True:
        # the stresscalc object has information about the satellite buried in it.
        NSRFailDirByLat(nsr_stresscalc, lats=linspace(20,36,9))

    if maps is True:
        FitMap(goodlins=nsr_goodlins, badlins=nsr_badlins, nbins=9, titlestr="global lins, fit to NSR")
        FitMap(goodlins=tpw_goodlins, badlins=tpw_badlins, nbins=9, titlestr="pre-TPW lins, fit to NSR")
        FitMap(goodlins=nsr_goodlins, badlins=regular_nsr_lingen(nlats=60), nbins=9, titlestr="global lins, in NSR formation locations", derotate=True, stresscentric=True)
        # Map of the fake (great circle) NSR lineaments:
        FitMap(goodlins=gc_goodlins, badlins=gc_badlins, nbins=9, titlestr="Great Circle Segments fit to NSR")
        # Map of the perfect synthetic NSR lineaments:
        FitMap(goodlins=synlins, titlestr="Perfect Synthetic NSR Lineaments", nbins=9)

        # This shows the resolution of coverage we've got to work with:
        EuropaCoverageMap()

    if hists is True:
        ActHist_NSRvTPW(nsr_goodlins, nsr_badlins, tpw_goodlins, tpw_badlins, nbins=18)
        ActHist_ByLength(nsr_goodlins, nsr_badlins, lengths=[0,150,300,500,1000,2000], nbins=18)
        ActHist_MappedGC(nsr_goodlins, nsr_badlins, gc_goodlins, gc_badlins, nbins=18)
        #ActHist_Viscoelastic([nsrlins, nsrD01lins, nsrD1lins, nsrD10lins, nsrD100lins, nsrD1000lins])
 
    if examples is True: #{{{2
        good_nsr1    = nsrlins[2]   # 1500 km long, symmetric arcuate lineament in the N. hemisphere
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

        example_lins = [ good_nsr, bad_nsr, cycloid, sinuous, short_hilat, short_lolat ]

        FitCurveExamples(example_lins)
        DoppelgangerExamples(example_lins)
#}}}2

    if stats is True:
        MinDeltaRMSCorr(nsrlins)
        #DbarDeltaStats(nsrlins)
        LinLatLonStats(goodlins=nsr_goodlins, badlins=nsr_badlins, label="Mapped Lineaments")
        LinLatLonStats(goodlins=gc_goodlins,  badlins=gc_badlins,  label="Great Circle Segments")
        LinLatLonStats(goodlins=tpw_goodlins, badlins=tpw_badlins, label="pre-TPW Lineaments (pole=80E10N)")
        LinLatLonStats(goodlins=syn_goodlins, badlins=syn_badlins, label="Synthetic NSR Lineaments")
        # This is the pool of 12,000 synthetic NSR lineaments from which the
        # above control dataset was drawn.  We can use it to show what dense
        # "perfect" coverage of the surface would look like, 
        synlins_pool = load_lins("output/synth_pool")
        LinLatLonStats(goodlins=synlins_pool, badlins=[],          label="All Synthetic NSR Lineaments")

# end makefigs }}}

def ActHist_NSRvTPW(nsr_goodlins, nsr_badlins, tpw_goodlins, tpw_badlins, nbins=18): #{{{
    """
    A plot showing the activity histograms and the aggregate fit curves for the
    set of lineaments included in each of the analyses.  One for the lineaments
    in their mapped locations/orientations, and one for the lineaments
    transformed to their pre-TPW locations/orientations.

    """

    the_fig = figure()
    # Two subplots, one above the other:
    act_hist_axes = the_fig.add_subplot(2,1,1)
    fit_curv_axes = the_fig.add_subplot(2,1,2)

    nsr_goodlen = sum([lin.length for lin in nsr_goodlins])
    nsr_badlen  = sum([lin.length for lin in nsr_badlins])
    nsr_pct_kept = 100.0*nsr_goodlen/(nsr_goodlen + nsr_badlen)

    tpw_goodlen = sum([lin.length for lin in tpw_goodlins])
    tpw_badlen  = sum([lin.length for lin in tpw_badlins])
    tpw_pct_kept = 100.0*tpw_goodlen/(tpw_goodlen + tpw_badlen)

    lin_labels = ['Modern Mapped (%.0f%% retained)'  % (nsr_pct_kept),\
                  'Pre-TPW 80E10N (%.0f%% retained)' % (tpw_pct_kept) ]

    activity_histogram([nsr_goodlins, tpw_goodlins], hist_ax=act_hist_axes,nbins=nbins, labels=lin_labels)
    fitcurves([nsr_goodlins, tpw_goodlins], fit_ax=fit_curv_axes, labels=lin_labels)

    the_fig.show()

# end ActHist_NSRvTPW }}}

def ActHist_MappedGC(nsr_goodlins, nsr_badlins, gc_goodlins, gc_badlins, nbins=18): #{{{
    """
    A plot comparing the mapped features to random great circle segments having
    the same overall length distribution.

    """

    the_fig = figure()
    # Two subplots, one above the other:
    act_hist_axes = the_fig.add_subplot(2,1,1)
    fit_curv_axes = the_fig.add_subplot(2,1,2)

    nsr_goodlen = sum([lin.length for lin in nsr_goodlins])
    nsr_badlen  = sum([lin.length for lin in nsr_badlins])
    nsr_pct_kept = 100.0*nsr_goodlen/(nsr_goodlen + nsr_badlen)

    gc_goodlen = sum([lin.length for lin in gc_goodlins])
    gc_badlen  = sum([lin.length for lin in gc_badlins])
    gc_pct_kept = 100.0*gc_goodlen/(gc_goodlen + gc_badlen)

    lin_labels = ['Modern Mapped (%.0f%% retained)'  % (nsr_pct_kept),\
                  'Great Circles (%.0f%% retained)' % (gc_pct_kept) ]

    activity_histogram([nsr_goodlins, gc_goodlins], hist_ax=act_hist_axes,nbins=nbins, labels=lin_labels)
    fitcurves([nsr_goodlins, gc_goodlins], fit_ax=fit_curv_axes, labels=lin_labels)

    the_fig.show()
#}}}

def ActHist_ByLength(goodlins, badlins, lengths=[0,150,300,500,1000,2000], nbins=18): #{{{
    """
    Activity histogram and aggregate fit curves for the lineaments compared to
    the NSR stress field, but with the lineaments broken into several different
    length categories, showing the degree to which the result is supported by
    all the different length categories.

    """

    the_fig = figure()
    # Two subplots, one above the other:
    act_hist_axes = the_fig.add_subplot(211)
    fit_curv_axes = the_fig.add_subplot(212)

    numlens = len(lengths)-1
    goodlins_by_length = []
    goodlengths_by_length = [0,]*numlens
    badlins_by_length = []
    badlengths_by_length = [0,]*numlens

    for N in range(numlens):
        goodlins_by_length.append([])
        badlins_by_length.append([])

    sat_radius_km = goodlins[0].stresscalc.stresses[0].satellite.radius()/1000
    for N in range(numlens):
        for lin in goodlins:
            lin_length_km = lin.length*sat_radius_km
            if (lin_length_km > lengths[N]) and (lin_length_km < lengths[N+1]):
                goodlins_by_length[N].append(lin)
                goodlengths_by_length[N] += lin_length_km
                continue

        for lin in badlins:
            lin_length_km = lin.length*sat_radius_km
            if (lin_length_km > lengths[N]) and (lin_length_km < lengths[N+1]):
                badlins_by_length[N].append(lin)
                badlengths_by_length[N] += lin_length_km
                continue

    pct_kept = []
    lin_labels = []
    for N in range(numlens):
        pct_kept.append(100.0*goodlengths_by_length[N]/(goodlengths_by_length[N]+badlengths_by_length[N]))
        lin_labels.append("%d - %d km (%.0f%% retained)" % (lengths[N], lengths[N+1], pct_kept[N]))

    activity_histogram(goodlins_by_length, hist_ax=act_hist_axes, nbins=nbins, labels=lin_labels)
    fitcurves(goodlins_by_length, fit_ax=fit_curv_axes, labels=lin_labels)

    show()

# end ActHist_ByLength }}}

def ActHist_Viscoelastic(D01, D1, D10, D100, delta_max, dbar_max): #{{{
    """
    Activity histogram and aggregate fit curves for the lineaments compared to NSR
    stress fields for three values of Delta, representing four different viscous
    relaxation scenarios, ranging from elastic, to viscous: Delta={0.1,1,10,100}

    """
# end ActHist_Viscoelastic }}}

def FitMap(goodlins=[], badlins=[], titlestr="Features colored by fit", lin_cm=cm.jet, nbins=9, derotate=False, stresscentric=False): #{{{
    """
    Creates a global map of the lineaments, color coding them by what amount of
    backrotation (b) minimizes delta_rms(b) when compared to the NSR stress
    field.  Lineament width indicates the value of min(delta_rms(b)), with
    wider lineaments agreeing better with NSR stresses.

    Those lineaments which did not meet the criteria for inclusion in the
    analysis are plotted thin and black.

    The background might be a grayscale map of the resolution of the images
    that went into the creation of the 500 meter resolution USGS mosaic of
    Europa, from which the lineaments were mapped.

    """

    lin_cm = colors.ListedColormap(lin_cm(linspace(0,1,nbins)))

    sat_radius_km = goodlins[0].stresscalc.stresses[0].satellite.radius()/1000
    delta_max = max([ max(lin._Lineament__good_fits[:,1]) for lin in goodlins ])
    dbar_max = max([ max(lin._Lineament__good_fits[:,2]) for lin in goodlins ])
    badlen  = sum( [ lin.length for lin in badlins  ] ) * sat_radius_km
    goodlen = sum( [ lin.length for lin in goodlins ] ) * sat_radius_km

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

    if len(badlins) > 0:
        badlines,linfitmap = lineament.plotlinmap(badlins, map=linfitmap, color='black', alpha=1.0, linewidth=1.0, lon_cyc=lon_cyc, lat_mirror=lat_mirror)

    goodlines = []
    for lin in goodlins:
        # load some information about the lineament's best fits...
        best_b, best_delta = lin.best_fit()[0:2]

        # Map the color of the lineament to its best_b
        lin_color = lin_cm(int((lin_cm.N)*(best_b/pi)))

        # use line width to indicate goodness of best_fit
        lin_width = 5*(delta_max - best_delta)/delta_max
        lin_alpha = 1.0

        backrot = 0
        if derotate:
            backrot = best_b

        newline, linfitmap = lineament.plotlinmap([lin.lonshift(backrot)], map=linfitmap, linewidth=lin_width, color=lin_color, alpha=lin_alpha, lon_cyc=lon_cyc, lat_mirror=lat_mirror)
        goodlines.append(newline)

    map_ax.set_title(titlestr + " " + 'max($\delta_{rms}$)=%.2g$^\circ$, max($\\bar{D}$)=%.2g' % (degrees(delta_max),dbar_max))
    cb_ax,kw = colorbar.make_axes(map_ax, orientation="horizontal", pad=0.05, shrink=0.5)
    colorbar.ColorbarBase(cb_ax, cmap=lin_cm, orientation="horizontal", norm=colors.BoundaryNorm(linspace(0,180,nbins+1),nbins), format=r'%.0f$^\circ$')
    # Fix up the colorbar a bit:
    cb_ax.invert_xaxis()
    cb_ax.set_xlabel("backrotation b")

    show()

# end FitMap }}}

def MinDeltaRMSCorr(lins): #{{{
    """
    Several scatter plots showing the correlation (or lack thereof) between
    min(delta_rms) and other lineament characteristics, including:

    A: lineament.sinuosity()
    B: lineament.length
    C: min(lineament.latitudes()), color coded by length?
    """

    goodlins = [ lin for lin in lins if len(lin.good_fits()) > 0 ]
    delta_max = max([lin.best_fit()[1] for lin in goodlins])
    dbar_max  = max([lin.best_fit()[2] for lin in goodlins])

    linlens  = [ lin.length for lin in goodlins ]
    meanlats = [ mean(fabs(lin.latitudes())) for lin in goodlins ]
    best_deltas = [ lin.best_fit()[1] for lin in goodlins ]
    sins = [ lin.sinuosity() for lin in goodlins ]

    the_fig = figure()
    
    ax1 = the_fig.add_subplot(1,1,1)

    lat_cols = [ cm.jet(lat/max(meanlats)) for lat in meanlats ]

    symbols = scatter( sins, degrees(best_deltas), c=lat_cols, s=300*array(linlens) )
    symbols.set_alpha(0.75)
    
    ax1.set_xlabel('sinuosity (S)')
    ax1.set_ylabel(r'best fit ($\delta_{rms}$)')
    ax1.set_ylim(0,degrees(delta_max))
    ax1.set_xlim(1,1.1)

    degree_formatter = matplotlib.ticker.FormatStrFormatter(r'%.0f$^\circ$')
    ax1.yaxis.set_major_formatter(degree_formatter)

    cb_ax,kw = colorbar.make_axes(ax1, pad=0.05)
    colorbar.ColorbarBase(cb_ax, cmap=cm.jet, norm=colors.Normalize(vmin=0,vmax=degrees(max(meanlats))), format=r'%.0f$^\circ$')
    cb_ax.set_ylabel(r'mean latitude ($\bar{\theta}$)')

    ax1.set_title(r'N=%d, $\bar{D}\leq$%.3g, $R^2(S-\delta_{rms})$=%.3g, $R^2(\bar{\theta}-\delta_{rms})$=%.3g, $R^2(L-\delta_{rms})$=%.3g'\
                   % (len(goodlins),dbar_max,corrcoef(sins,best_deltas)[0,1]**2,corrcoef(meanlats,best_deltas)[0,1]**2,corrcoef(linlens,best_deltas)[0,1]**2) )

# end MinDeltaRMSCorr }}}

def NSRFailDirByLat(nsr_stresscalc, lats=linspace(0,80,9), lons=linspace(0,180,1000)): #{{{
    """
    Graphic display of how the range of possible failure orientations changes
    with latitude in the NSR stress field, demonstrating that for any short
    lineament above 30 degrees latitude, there will always be a good fit.

    """
    lats = radians(lats)
    lons = radians(lons)

    N = len(lats)
    fail_orient = [ [] for n in range(N) ]
    for lat,n in zip(lats,range(N)):
        for lon in lons:
            nsr_pcs = nsr_stresscalc.principal_components(pi/2-lat,lon,0)
            if nsr_pcs[0] > 0:
                fail_orient[n].append(nsr_pcs[3])

    roses = figure()
    rose_axes = [ roses.add_subplot(3,3,n, polar=True) for n in arange(9)+1 ]

    for ax,n in zip(rose_axes,arange(len(rose_axes))):
        ax.hist(fail_orient[n], range=(0,pi), bins=180, facecolor='gray', edgecolor='gray')
        ax.set_title("lat = %g          " % (degrees(lats[n]),) )

# end NSRFailDirByLat }}}

def EuropaCoverageMap(worst_res=2000): #{{{
    from osgeo import gdal
    import matplotlib.colors as colors
    from mpl_toolkits import basemap

    # Incredibly, I can just read in the whole raster dataset like this!
    EuropaGrid = gdal.Open("input/europa_coverage_simpcyl.cub")

    # The data cube has 4 bands:
    # 1. resolution [km/px]
    # 2. emission angle [degrees]
    # 3. incident angle [degrees]
    # 4. phase angle [degrees]

    # changing it to m/px here
    res_arr = 1000*EuropaGrid.GetRasterBand(1).ReadAsArray()

    the_fig = figure(figsize=(12,8))

    palette = cm.jet
    # "bad" values are where there is no coverage.
    palette.set_bad('black')
    # "over" values are where the resolution is worse than worst_res
    palette.set_over('gray')
    palette.set_under(palette(0))

    # mask the array so that the "bad" values disappear
    res_masked = ma.masked_where(res_arr < 0, res_arr)
    res_norm = colors.Normalize(vmin=500,vmax=worst_res)

    #res_ax = the_fig.add_subplot(1,1,1)
    res_map = Basemap(llcrnrlon=0,llcrnrlat=-90,urcrnrlon=360,urcrnrlat=90)
    res_ax = the_fig.add_subplot(1,1,1)
    res_ax.set_title("Europa mosaic resolution and coverage")
    res_map.imshow(res_masked, cmap=palette, norm=res_norm, extent=[0,360,-90,90], origin='upper')
    res_map.drawmeridians(linspace(0,360,13), labels=[1,0,0,1])
    res_map.drawparallels(linspace(-90,90,7), labels=[1,0,0,1])

    cb_ax,kw = colorbar.make_axes(res_ax, pad=0.05, orientation='horizontal', shrink=0.5)
    colorbar.ColorbarBase(cb_ax, cmap=palette, norm=res_norm, orientation='horizontal')
    cb_ax.set_xlabel('resolution [m/px]')
#}}}

def DbarDeltaStats(lins_list, labels, deltas=[], dbars=[]): #{{{
    """
    A routine for plotting the behavior of several quantities 
    which depend on both delta_max and dbar_max

      - proportion of lineaments retained (those with at least one good fit)
      - degeneracy of fits (how many good fits?)
    
    Do these things behave differently depending on what kind of dataset I'm
    fitting?

      - mapped lineaments
      - mapped TPW 80E10N lineaments
      - great circle segments
      - synthetic NSR lineaments

    """

    numsets = len(lins_list)
    scale_factor=3
    the_fig = figure(figsize=(numsets*scale_factor,3*scale_factor))

    # how many stats are we going to collect?
    numstats=3

    # The (number) fraction of lineaments kept
    num_frac = zeros((numsets,len(deltas),len(dbars)))
    # The (length) fraction of lineaments kept
    len_frac = zeros((numsets,len(deltas),len(dbars)))
    # The fraction of B values that are "good enough"
    fit_prop = zeros((numsets,len(deltas),len(dbars)))

    # iterate over the lineament sets, and the space of possible delta_max, and
    # dbar_max values
    for lins,i in zip(lins_list,range(numsets)):
        for delta,j in zip(deltas,range(len(deltas))):
            for dbar,k in zip(dbars,range(len(dbars))):

                goodlins, badlins, goodlengths, badlengths, fitprop = cache_linstats(lins, name=labels[i], delta_max=delta, dbar_max=dbar)

                num_frac[i,j,k] = float(len(goodlins))/float((len(badlins)+len(goodlins)))
                len_frac[i,j,k] = sum(goodlengths)/(sum(goodlengths+badlengths))
                fit_prop[i,j,k] = mean(fitprop)

        num_frac_ax = the_fig.add_subplot(numstats, numsets, (i+1)) # i=(0,1) => (1,2) 
        num_frac_ax.imshow(num_frac[i,:,:], interpolation='nearest', vmin=0.0, vmax=1.0, origin='lower', cmap=cm.gray)
        num_frac_ax.set_title("N frac kept")
        num_frac_ax.set_xlabel(r'max($\bar{D}$)')
        num_frac_ax.set_ylabel(r'max($\delta_{rms}$)')
        # need to effectively label x and y axes with good ticklabels

        len_frac_ax = the_fig.add_subplot(numstats, numsets, numsets+(i+1)) # i=(0,1) => (3,4)
        len_frac_ax.imshow(len_frac[i,:,:], interpolation='nearest', vmin=0.0, vmax=1.0, origin='lower', cmap=cm.gray)
        len_frac_ax.set_title("km frac kept")
        len_frac_ax.set_xlabel(r'max($\bar{D}$)')
        len_frac_ax.set_ylabel(r'max($\delta_{rms}$)')
        # need to effectively label x and y axes with good ticklabels

        fit_prop_ax = the_fig.add_subplot(numstats, numsets, 2*numsets+(i+1)) # i(0,1) => (5,6)
        fit_prop_ax.imshow(fit_prop[i,:,:], interpolation='nearest', vmin=0.0, vmax=1.0, origin='lower', cmap=cm.gray)
        fit_prop_ax.set_title("fit frac")
        fit_prop_ax.set_xlabel(r'max($\bar{D}$)')
        fit_prop_ax.set_ylabel(r'max($\delta_{rms}$)')
        # need to effectively label x and y axes with good ticklabels

        the_fig.show()

    return(num_frac,len_frac,fit_prop)

# end DbarDeltaCorr }}}

def LinLatLonStats(goodlins=[], badlins=[], label=""): #{{{
    """
    Shows the distribution of latitudes and longitudes of the lineament
    segments over the surface of the satellite, relative to an expected random
    distribution.

    """

    # TODO:
    #   - compare derotated "good" lins from each dataset vs. the synthetics?

    good_lat_pile = []
    good_lon_pile = []
    bad_lat_pile  = []
    bad_lon_pile  = []
    sat_radius_km = 1561

    try:
        delta_max = max([ max(lin.good_fits()[:,1]) for lin in goodlins if lin.good_fits() ])
        dbar_max = max([ max(lin.good_fits()[:,2]) for lin in goodlins if lin.good_fits() ])
    except ValueError:
        delta_max = pi/2
        dbar_max = 0.5

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

    the_fig = figure()
    lon_ax = the_fig.add_subplot(2,1,1)
    lon_ax.set_title(r'%s max($\delta_{rms}$)=%.0f$^\circ$ max($\bar{D}$)=%.0g' % (label, degrees(delta_max), dbar_max) )
    if len(bad_lon_pile) > 0:
        lins_n, lin_bins, lins_patches = lon_ax.hist([mod(degrees(good_lon_pile),180),mod(degrees(bad_lon_pile),180)],\
                                                     bins=180, range=(0,180), label="passed", histtype='barstacked', rwidth=1.0)
        lins_patches[1][0].set_label("failed")
    else:
        lon_ax.hist(mod(degrees(good_lon_pile),180), bins=180, range=(0,180), label="passed", histtype='barstacked', rwidth=1.0)

    lon_ax.set_xlim(0,180)
    lon_ax.axhline(y=total_length/180.0, color='red', linewidth=2, label="random")
    lon_ax.set_xticks(linspace(0,180,7))
    lon_ax.set_xlabel("longitude (modulo $\pi$)")
    lon_ax.set_ylabel("length [km]")
    lon_ax.legend()

    lat_ax = the_fig.add_subplot(2,1,2)
    if len(bad_lat_pile) > 0:
        lins_n, lin_bins, lins_patches = lat_ax.hist([abs(degrees(good_lat_pile)),abs(degrees(bad_lat_pile))],\
                                                     bins=90, range=(0,90), label="passed", histtype='barstacked', rwidth=1.0)
        lins_patches[1][0].set_label("failed")
    else:
        lat_ax.hist(abs(degrees(good_lat_pile)), bins=90, range=(0,90), label="passed", histtype='barstacked', rwidth=1.0)
    lat_ax.set_xlim(0,90)
    lats = linspace(0,pi/2,100)
    lat_ax.plot(degrees(lats), (total_length/90)*(cos(lats))*pi/2, color='red', linewidth=2, label="random")
    lat_ax.set_xticks(linspace(0,90,7))
    lat_ax.set_xlabel("abs(latitude)")
    lat_ax.set_ylabel("length [km]")
    lat_ax.legend()

    the_fig.show()

# end LinStats }}}

########################################
#         Not Yet Implemented          #
########################################

def LinLatLonStats2D(lins, label=""): #{{{
    """
    Create a map of lineament km per unit area over the surface of the
    satellite, and the deviation that represents from a random distribution
    having the same cumulative length.
    
    A correlation between this deviation and the resolution of coverage map can
    be calculated, showing how much of the spatial distribution is the result
    of our uneven coverage.

    """
    # Need the lats and lons of every lineament segment midpoint, and the
    # length of each segment, but we don't care what order it's all in,
    # or whether lineament

    pass
#}}}

def LinStats(lins, label=""): #{{{
    """
    Creates several histograms, showing information about the statistical
    distribution and morphologies of the lineaments used in the analysis.
    Several of these can be shown in each individual bar chart, showing
    side-by-side comparisons.

    Possibilities include:

    A1: Number weighted length distribution of all included lineaments
    A2: Number weighted length distribution of all excluded lineaments
        (barstacked histogram, so that total bar heights show the distribution
        of the entire mapped dataset)

    B1: Length weighted length distribution of all included lineaments
    B2: Length weighted length distribution of all excluded lineaments
        (barstacked histogram, so that total bar heights show the distribution
        of the entire mapped dataset)

    C1: Length weighted sinuosity distribution of included lineaments
    C2: Length weighted sinuosity distribution of excluded lineaments
        (barstacked histogram, so that total bar heights show the distribution
        of the entire mapped dataset)
    C3: Length weighted sinuosity distribution of a synthetic population with
        the same length distribution as the included features.

    D1: Latitude and longitude biases present in the resolution of images that
        went into the 500 meter global mosaic, as captured by the number of image
        pixels present at each latitude and longitude, compared to the number that
        there would have been if the entire surface had been imaged at 500 m/px.
    D2: Correlation between lat/lon biases in mapped lineaments, and resolution
        of images.

    """
#}}}

def FitCurveExamples(example_lins, delta_max, dbar_max): #{{{
    """
    Several archetypal fit curves (delta_rms(b) and dbar(b)), showing a
    selection of different behaviors, for different lineament morphologies.

    Possibly to include:
    
    A: An "obviously" NSR related lineament (long, arcuate, good best fit).
    B: An "obviously" non-NSR lineament (long E-W near equator?).
    C: Cycloidal lineament.
    D: High sinuosity, non-cycloidal lineament.
    E: Short lineament at high latitude.
    F: Short lineament at low latitude.

    """
# end FitCurveExamples }}}

def DoppelgangerExamples(example_lins): #{{{
    """
    A close up of the example lineaments shown in the fit curve examples, with
    a whole range of their doppelganger lineaments plotted, color coded by the
    value of b that was used to generate them, labeled with the value of dbar
    that they result in.

    """
# end DoppelgangerExamples }}}


