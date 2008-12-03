#!/usr/bin/python
"""
Routines for comparing a set of lineaments to an NSR stress field.

"""

from . import lineament
from . import satstress
from pylab import *
from matplotlib import colors, colorbar
from mpl_toolkits.basemap import Basemap
from scipy.stats.mmorestats import idealfourths
import random

def visco_lincomp(nb=180, delta_max=radians(20), dbar_max=0.05,\
                  satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                  shapefile = "input/GlobalLineaments",\
                  nlins=0, save=True):
#{{{
    europa = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(the_sat),])

    europa_D01 = satstress.Satellite(open(satfile,'r'))
    europa_D01.nsr_period *= 10

    europa_D1 = satstress.Satellite(open(satfile,'r'))
    europa_D10 = satstress.Satellite(open(satfile,'r'))
    europa_D100 = satstress.Satellite(open(satfile,'r'))

    lins = lineament.shp2lins(shapefile, stresscalc=nsr_stresscalc)

    if nlins:
        lins = lins[:nlins]

    # For TPW:
    if(pnp_lon is not None and pnp_lat is not None):
        lins = [ lin.poleshift(radians(pnp_lon), radians(pnp_lat)) for lin in lins ]
        filetitle = title+"_TPW_%dE_%dN" % (int(pnp_lon), int(pnp_lat))
        plotlabel = "PNP: %dE, %dN" % (pnp_lon,pnp_lat)
    else:
        filetitle = title
        plotlabel = "global lins"

    for n,lin in zip(range(len(lins)),lins):
        print("Calculating fits for lineament %d / %d" % (n+1,len(lins)))
        fits = lin.fit(bs=linspace(0,pi,nb, endpoint=False))
        print("    Caching the best fits")
        lin.cache_good_fits(delta_max=delta_max, dbar_max=dbar_max)

    if save:
        save_lins(lins, name=filetitle)

    # Here's where I'd do output... if I had the script written!

    return(lins)
#}}}

def mapped_lincomp(lins=None, nb=180, delta_max=radians(20), dbar_max=0.05,\
                  satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                  shapefile = "input/GlobalLineaments",\
                  pnp_lon = None, pnp_lat = None, nlins=0, title="global", save=True):
#{{{
    the_sat = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(the_sat),])

    if lins is None:
        lins = lineament.shp2lins(shapefile, stresscalc=nsr_stresscalc)

    if nlins:
        lins = lins[:nlins]

    # For TPW:
    if(pnp_lon is not None and pnp_lat is not None):
        lins = [ lin.poleshift(radians(pnp_lon), radians(pnp_lat)) for lin in lins ]
        filetitle = title+"_TPW_%dE_%dN" % (int(pnp_lon), int(pnp_lat))
        plotlabel = "PNP: %dE, %dN" % (pnp_lon,pnp_lat)
    else:
        filetitle = title
        plotlabel = "global lins"

    for n,lin in zip(range(len(lins)),lins):
        print("Calculating fits for lineament %d / %d" % (n+1,len(lins)))
        fits = lin.fit(bs=linspace(0,pi,nb, endpoint=False))
        print("    Caching the best fits")
        lin.cache_good_fits(delta_max=delta_max, dbar_max=dbar_max)

    if save:
        save_lins(lins, name=filetitle)

    # Here's where I'd do output... if I had the script written!

    return(lins)
#}}}

def synth_nsr_lincomp(lins=None, nb=180,):
#{{{
    # - read in the mapped NSR lineaments
    # - read in the synthetic NSR lineaments
    # - subsample the synthetic lineaments so as to mimic the length
    #   distribution of the mapped lineaments
    # - introduce longitude shifts in the synthetic features, attempting
    #   to recreate the observed peak in activity, with a maximum around
    #   b=30
    # - calculate the fit metrics for the synthetic lineaments
    # - save the fitted synthetic features for later analysis.
    pass
#}}}

def synth_gc_lincomp(lins=None, nb=180):
#{{{
    # - read in the mapped NSR lineaments
    # - read in the fake great circle lineaments
    # - subsample the great circle lineaments so as to mimic the length
    #   distribution of the mapped features
    # - calculate fit metrics for the synthetic features.
    # - save fitted synthetic features for later analysis.
    pass
#}}}

def browse_fits(lins, delta_max=radians(20), dbar_max=0.05, fixlon=True): #{{{
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
        lin_length_km = sat_radius_km*lin.length()
        linmap = Basemap(rsphere=sat_radius_km*1000, ax=map_ax, projection='cyl')
        map_ax.set_title("lin=%d, length=%g km" % (lin_num, lin_length_km) )
        linmap.drawmeridians(range(-180,180,30), labels=[1,0,0,1])
        linmap.drawparallels(range(-90,90,30), labels=[1,0,0,1])

        lin.plot(map=linmap, fixlon=fixlon, lw=2.0)

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

        lin_length_km = lin.length()*sat_radius_km
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
    from pickle import dump

    dump(lins, open("output/%s_%s.pkl" % (name, strftime('%Y%m%d%H%M%S')),'w'))
#}}} end save_lins()

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

# mean/stddev aren't *really* meaningful here, because there's no reason to
# expect this to have a normal distribution.
#    bin_means   = [ mean(N) for N in lins_n ]
#    bin_stddevs = [ std(N) for N in lins_n ]
#    for bin_mean,bin_stddev,label,color in zip(bin_means, bin_stddevs,labels,colors):
#        print("%s bin mean count = %.0f +/- %.0f" % (label, bin_mean, bin_stddev) )
#        b_vals = linspace(0,180,len(lins_n))
#
#        hist_ax.fill_between(b_vals, bin_mean+bin_stddev, max([0,bin_mean-bin_stddev]), color=color, alpha=0.2)
#        hist_ax.fill_between(b_vals, bin_mean+2*bin_stddev, max([0,bin_mean-2*bin_stddev]), color=color, alpha=0.2)

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
        for lin in lins:
            linlength = lin.length()
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
                derotated = hstack([lin._Lineament__fits[best_idx:,1],lin._Lineament__fits[:best_idx,1]])
                agg_delta_rms_derotated += derotated*linlength

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

def cleanfits(lin): #{{{
    """
    Remove all fits with b values corresponding to non-integer numbers of degrees

    Remove all doppelganger lineaments from the lineament as well.

    """

    lin._Lineament__doppelgangers=[]
    new_fits = array([])
    for fit in lin.fit():
        if fit[0] >= pi:
            continue
        if mod(degrees(fit[0]),1) < 1e-6 or mod(degrees(fit[0]),1) > 1-1e-4:
            if len(new_fits) == 0:
                new_fits = array([fit,])
            else:
                new_fits = vstack([new_fits,fit])

    lin._Lineament__fits = new_fits
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

def random_nsrlins(nsr_stresscalc=None, nlins=1000, minlen=0.1, maxlen=1.0): #{{{
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
        if newlin.length() >= minlen:
            linlist.append(newlin)

    return linlist
    #}}}

def random_gclins(nlins=1000, minlen=0.1, maxlen=1.0): #{{{
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

# Routines that generate the plots I have listed in the paper outline:

def makefigs(delta_max=radians(20), dbar_max=0.10,\
             make_faildir=False, make_maps=True, make_hists=True,\
             make_examples=False, make_stats=False): #{{{
    """
    A wrapper that takes input parameters for the figures, and calls all of the
    individual figure plotting routines, below.

    """
    import pickle

    # These are the lineaments in their modern locations, compared to an
    # elastic NSR stress field, with Delta ~ 0.001.
    nsrlinfile = "output/nsr_delta20_dbar05_20081124174547.pkl"

    # All the mapped lineaments, transformed to a pre-TPW coordinate system,
    # with the paleo-north pole at 80E, 10N in current lat/lon.
    tpwlinfile = "output/tpw_80E10N_delta20_dbar05_20081124174727.pkl"

    # All the mapped lineaments, compared to various NSR stress fields, for
    # different amounts of viscous relaxation:
    # Delta = 0.1
    nsrD01linfile = "output/nsr_D01_200811.pkl"
    # Delta = 1.0
    nsrD1linfile = "output/nsr_D1_200811.pkl"
    # Delta = 10
    nsrD10linfile = "output/nsr_D10_200811.pkl"
    # Delta = 100
    nsrD100linfile = "output/nsr_D100_200811.pkl"
    # Delta = 1000
    nsrD1000linfile = "output/nsr_D1000_200811.pkl"

    # A set of synthetic lineaments, each of which is a portion of a great
    # circle, oriented randomly on the surface, having a random location, in a
    # variety of lengths, allowing a very basic null hypothesis test.  They've
    # had their fits calculated relative to an elastic NSR stress field.
    nsr_SynthGCL_linfile = "output/nsr_synth_gcl.pkl"

    # A set of 1000 perfect synthetic NSR lineaments, as a control dataset:
    nsr_SynthPerfect_linfile = "output/nsr_synth_perfect.pkl"

    # Read in the NSR lineaments, and cache their good fits, it necessary:
    nsrlins = pickle.load(open(nsrlinfile))
    nsr_stresscalc = nsrlins[0].stresscalc

    # going to assume we're all using the same satellite here...
    sat_radius_km = nsr_stresscalc.stresses[0].satellite.radius()/1000.0

    print("""Caching good fits for modern mapped features:
    delta_max = %.2g degrees
    dbar_max  = %.2g""" % (degrees(delta_max), dbar_max))
    for lin in nsrlins:
        lin.cache_good_fits(delta_max=delta_max, dbar_max=dbar_max)
    nsr_goodlins = [ lin for lin in nsrlins if len(lin.good_fits()) > 0 ]
    nsr_badlins  = [ lin for lin in nsrlins if lin not in nsr_goodlins ]
    nsr_goodlength = sum( [lin.length() for lin in nsr_goodlins ]) * sat_radius_km
    nsr_badlength  = sum( [lin.length() for lin in nsr_badlins  ]) * sat_radius_km
    print("""    lineaments kept: %d / %d (%.0f%%)
    kilometers kept: %.0f / %.0f (%.0f%%)""" % (len(nsr_goodlins), len(nsrlins), 100.0*len(nsr_goodlins)/float(len(nsrlins)),\
           nsr_goodlength, nsr_goodlength+nsr_badlength, 100.0*nsr_goodlength/float(nsr_goodlength+nsr_badlength)) )

    # Read in the pre-TPW lineaments, and cache their good fits:
    tpwlins = pickle.load(open(tpwlinfile))
    print("""Caching good fits for pre-TPW features:
    delta_max = %.2g degrees
    dbar_max  = %.2g""" % (degrees(delta_max), dbar_max))
    for lin in tpwlins:
        lin.cache_good_fits(delta_max=delta_max, dbar_max=dbar_max)
    tpw_goodlins = [ lin for lin in tpwlins if len(lin.good_fits()) > 0 ]
    tpw_badlins  = [ lin for lin in tpwlins if lin not in tpw_goodlins ]
    tpw_goodlength = sum( [lin.length() for lin in tpw_goodlins ]) * sat_radius_km
    tpw_badlength  = sum( [lin.length() for lin in tpw_badlins  ]) * sat_radius_km
    print("""    lineaments kept: %d / %d (%.0f%%)
    kilometers kept: %.0f / %.0f (%.0f%%)""" % (len(tpw_goodlins), len(tpwlins), 100.0*len(tpw_goodlins)/float(len(tpwlins)),\
           tpw_goodlength, tpw_goodlength+tpw_badlength, 100.0*tpw_goodlength/float(tpw_goodlength+tpw_badlength)) )

    if make_faildir is True:
        # the stresscalc object has information about the satellite buried in it.
        NSRFailDirByLat(nsr_stresscalc, lats=linspace(20,36,9))

    if make_maps is True:
        FitMap(goodlins=nsr_goodlins, badlins=nsr_badlins, nbins=9, titlestr="global lins, fit to NSR")
        #FitMap(goodlins=tpw_goodlins, badlins=tpw_badlins, nbins=9, titlestr="pre-TPW lins, fit to NSR", fixlon=True)
        FitMap(goodlins=nsr_goodlins, badlins=regular_nsr_lingen(nlats=60), nbins=9, titlestr="global lins, in NSR formation locations", derotate=True)

        # Map of the fake synthetic NSR lineaments:
        #FitMap()
        # Map of the perfect synthetic NSR lineaments:
        #Fitmap()

    if make_hists is True:
        ActHist_NSRvTPW(nsr_goodlins, tpw_goodlins, nbins=18)
        ActHist_ByLength(nsr_goodlins, lengths=[0,150,300,500,1000,2000], nbins=18)

        #ActHist_SynthMap()

        # Lineaments for doing the viscous relaxation comparison:
        #D01  = pickle.load(open(nsrD01linfile))
        #D1   = pickle.load(open(nsrD1linfile))
        #D10  = pickle.load(open(nsrD010linfile))
        #D100 = pickle.load(open(nsrD0100linfile))
        #ViscoActHist([nsrlins, D01, D1, D10, D100, D1000])
 
    if make_examples is True: #{{{2
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

    if make_stats is True:
        LinStats(nsrlins, delta_max, dbar_max)
        MinDeltaRMSCorr(nsrlins)
        DbarDeltaCorr(nsrlins)

# end makefigs }}}

def ActHist_NSRvTPW(nsrlins, tpwlins, nbins=18): #{{{
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

    sat_radius_km = nsrlins[0].stresscalc.stresses[0].satellite.radius()/1000
    lin_labels = ['modern (%d km)'  % (sum([lin.length() for lin in nsrlins])*sat_radius_km),\
                  'pre-TPW (%d km)' % (sum([lin.length() for lin in tpwlins])*sat_radius_km) ]

    activity_histogram([nsrlins, tpwlins], hist_ax=act_hist_axes,nbins=nbins, labels=lin_labels)
    fitcurves([nsrlins, tpwlins], fit_ax=fit_curv_axes, labels=lin_labels)

    the_fig.show()

# end ActHist_NSRvTPW }}}

def ActHist_ByLength(nsrlins, lengths=[0,150,300,500,1000,2000], nbins=18): #{{{
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
    lins_by_length = []

    for N in range(numlens):
        lins_by_length.append([])

    sat_radius_km = nsrlins[0].stresscalc.stresses[0].satellite.radius()/1000
    for lin in nsrlins:
        lin_length_km = lin.length()*sat_radius_km
        for N in range(numlens):
            if (lin_length_km > lengths[N]) and (lin_length_km < lengths[N+1]):
                lins_by_length[N].append(lin)
                continue

    lin_labels = ["all lengths (%d km)" % (sum([lin.length()*sat_radius_km for lin in nsrlins]))]

    for N in range(numlens):
        lin_labels.append("%d - %d km (%d km)" % (lengths[N], lengths[N+1], sum([ lin.length()*sat_radius_km for lin in lins_by_length[N] ])))

    activity_histogram([nsrlins,]+lins_by_length, hist_ax=act_hist_axes, nbins=nbins, labels=lin_labels)
    fitcurves([nsrlins,]+lins_by_length, fit_ax=fit_curv_axes, labels=lin_labels)

    show()

# end ActHist_ByLength }}}

def ViscoActHist(D01, D1, D10, D100, delta_max, dbar_max): #{{{
    """
    Activity histogram and aggregate fit curves for the lineaments compared to NSR
    stress fields for three values of Delta, representing four different viscous
    relaxation scenarios, ranging from elastic, to viscous: Delta={0.1,1,10,100}

    """
# end ViscoActHist }}}

def FitMap(goodlins=[], badlins=[], titlestr="Features colored by fit", lin_cm=cm.jet, nbins=18, derotate=False, fixlon=False): #{{{
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
    badlen  = sum( [ lin.length() for lin in badlins  ] ) * sat_radius_km
    goodlen = sum( [ lin.length() for lin in goodlins ] ) * sat_radius_km

    fig = figure()
    linfitmap = Basemap()
    linfitmap.drawmapboundary(fill_color="white")
    linfitmap.drawmeridians(range(-180,181,30), labels=[1,0,0,1])
    linfitmap.drawparallels(range(-90,91,30), labels=[1,0,0,1])
    map_ax = fig.axes[0]

    if len(badlins) > 0:
        badlines,linfitmap = lineament.plotlinmap(badlins, map=linfitmap, color='black', alpha=1.0, lw=1.0, fixlon=fixlon)

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
            fixlon = True

        newline, linfitmap = lineament.plotlinmap([lin.lonshift(backrot)], map=linfitmap, lw=lin_width, color=lin_color, alpha=lin_alpha, fixlon=fixlon)
        goodlines.append(newline)

    map_ax.set_title(titlestr + " " + 'max($\delta_{rms}$)=%.2g$^\circ$, max($\\bar{D}$)=%.2g' % (degrees(delta_max),dbar_max))
    cb_ax,kw = colorbar.make_axes(map_ax, orientation="horizontal", pad=0.05, shrink=0.5)
    colorbar.ColorbarBase(cb_ax, cmap=lin_cm, orientation="horizontal", norm=colors.BoundaryNorm(linspace(0,180,nbins+1),nbins), format=r'%.0f$^\circ$')
    # Fix up the colorbar a bit:
    cb_ax.invert_xaxis()
    cb_ax.set_xlabel("backrotation b")

    show()

# end FitMap }}}

def LinStats(nsrlins, delta_max, dbar_max): #{{{
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

    D1: Latitude of mapped lineament segment midpoints, scaled by the length of
        the segments that they are part of, normalized by the factor
        corresponding to the decrease in surface area at higher latitudes.
        (Showing the deviation from a normal distribution across the surface)
    D2: Longitude of lineament segment midpoints, scaled by the length of the
        segments that they are part of.
    D3: Latitude and longitude biases present in the resolution of images that
        went into the 500 meter global mosaic, as captured by the number of image
        pixels present at each latitude and longitude, compared to the number that
        there would have been if the entire surface had been imaged at 500 m/px.
    D4: Correlation between lat/lon biases in mapped lineaments, and resolution
        of images.

    """
# end LinStats }}}

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

def MinDeltaRMSCorr(lins): #{{{
    """
    Several scatter plots showing the correlation (or lack thereof) between
    min(delta_rms) and other lineament characteristics, including:

    A: lineament.sinuosity()
    B: lineament.length()
    C: min(lineament.latitudes()), color coded by length?
    """

    goodlins = [ lin for lin in lins if len(lin.good_fits()) > 0 ]
    delta_max = max([lin.best_fit()[1] for lin in goodlins])
    dbar_max  = max([lin.best_fit()[2] for lin in goodlins])

    linlens  = [ lin.length() for lin in goodlins ]
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

def DbarDeltaCorr(nsrlins): #{{{
    """
    Scatter plots exploring the correlation between dbar and delta, as a
    function of several lineament variables, including:

    A: lineament.sinuosity()
    B: lineament.length()
    C: min(lin.delta_rms())

    """
# end DbarDeltaCorr }}}
 
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
