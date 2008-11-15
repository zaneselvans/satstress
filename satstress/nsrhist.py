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

def overnight_run(lins=None, nb=61,\
                  satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                  shapefile = "input/GlobalLineaments",\
                  pnp_lon = None, pnp_lat = None, nlins=0, title="global"):
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
        filetitle = "tpw_%dE_%dN" % (int(pnp_lon), int(pnp_lat))
        plotlabel = "PNP: %dE, %dN" % (pnp_lon,pnp_lat)
    else:
        filetitle = title
        plotlabel = "global lins"

    for n,lin in zip(range(len(lins)),lins):
        print("Calculating fits for lineament %d / %d" % (n+1,len(lins)))
        fits = lin.fit(bs=linspace(0,pi,nb))

    save_lins(lins, name=filetitle)

    plot_histagg(lins, label=plotlabel, bins=18)

    return(lins)
#}}}

def plot_fits(lins, fixlon=False, delta_max=radians(20), dbar_max=0.05): #{{{
    lin_num=0

    clf()

    stresscalc = lins[0].stresscalc
    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000

    for lin in lins:
        # Draw mapped lineament in a Basemap map:
        subplot(2,1,1)
        title("lin=%d, length=%g km" % (lin_num, sat_radius_km*lin.length()))
        linmap = Basemap(rsphere=sat_radius_km*1000, projection='cyl')
        linmap.drawmeridians(range(-180,180,30), labels=[1,0,0,1])
        linmap.drawparallels(range(-90,90,30), labels=[1,0,0,1])

        lin.plot(map=linmap, fixlon=fixlon, lw=2.0)

        # Plot the mapped lineament's fit curve:
        delta_ax = subplot(2,1,2)
        delta_ax.set_title(r'$\delta_{rms}$')
        delta_ax.axis([0,180,0,90])
        delta_ax.set_yticks(arange(0,91,15))
        delta_ax.set_xticks(arange(0,181,15))
        delta_ax.grid(True)
        delta_ax.set_xlabel("westward translation, b [degrees]")
        delta_ax.set_ylabel(r'$\delta_{rms}$ [degrees]')
        delta_ax.plot(degrees(lin.bs()), degrees(lin.delta_rms()), c='black', ls='-', lw=2)

        # We're going to assume at this point that we've already calculated the fits:
        good_bs = lin.bs(delta_ismin=True, delta_max=delta_max, dbar_max=dbar_max, winwidth=18)
        if len(good_bs) > 0:
            scatter(degrees(good_bs), degrees(lin.delta_rms(bs=good_bs)),s=150,c='green',marker='o')

        # Mark on the plot where we had good fits that generated bad doppelgangers
        okay_bs = lin.bs(delta_ismin=True, winwidth=18)
        bad_bs = [ b for b in okay_bs if b not in good_bs ]
        if len(bad_bs) > 0:
            scatter(degrees(bad_bs), degrees(lin.delta_rms(bs=bad_bs)), s=150, c='red',marker='o')

        dbar_ax = twinx()
        dbar_ax.plot(degrees(lin.bs()), lin.dbar(), c='blue', ls='-', lw=2)
        dbar_ax.set_ylabel(r'$\bar{D}$')
        dbar_ax.set_yticks(arange(0,0.61,0.1))
        dbar_ax.set_xticks(arange(0,181,15))

        if len(good_bs) > 0:
            subplot(2,1,1)
            lin.plotdoppels(bs=good_bs, map=linmap, backrot=True, fixlon=fixlon, alpha=1.0, color='green')
        if len(bad_bs) > 0:
            subplot(2,1,1)
            lin.plotdoppels(bs=bad_bs, map=linmap, backrot=True, fixlon=fixlon, alpha=1.0, color='red')

        # wait for the user to hit return before we continue drawing...
        #x = raw_input("press return for next fit: ")
        print("click for next fit")
        x = ginput(n=1, timeout=-1, show_clicks=True)
        clf()
        lin_num += 1
#}}}

def calc_fit_hist(lins, delta_max=radians(15), dbar_max=0.05, winwidth=18): #{{{
    """
    For a list of Lineaments, lins, generate a histogram describing the number
    of kilometers of lineament at each value of b (0-180) having delta_rms <
    delta_max and dbar < dbar_max.

    Assumes that the lineaments have already had their fits calculated.

    """

    linhist = []
    for lin in lins:
        # Select the fits which meet our criteria for this lineament:
        good_fits = lin.fit(bs=lin.bs(delta_ismin=True, delta_max=delta_max,\
                                      dbar_max=dbar_max, winwidth=winwidth))

        # if there weren't any good fits, skip to the next lineament:
        if len(good_fits) == 0:
            continue

        # We generate the histogram dataset, adding one "count" for each km of
        # lineament that fit best at that value of b(ackrotation)
        for b, delta, dbar in good_fits:
            linlenkm = int(lin.length()*lin.stresscalc.stresses[0].satellite.radius()/1000)
            print("[ %g, %g, %g ]" % (degrees(b),degrees(delta), dbar) )
            linhist += linlenkm * [b,]

    return(linhist)
# }}}

def plot_fit_hist(fithist, bins=36, hist_title="Lineament formation vs. backrotation", color="green", rwidth=0.8, label=None): #{{{
    n, bins, patches = hist(degrees(fithist), bins=bins, facecolor=color, rwidth=rwidth, range=(0,180), label=label)
    grid(True)
    ylabel("total lineaments formed [km]")
    title(hist_title)
    xticks(arange(0,181,10))
    axis(xmin=0, xmax=180)
    return n, bins, patches
#}}}

def plot_aggregate_fits(lins, color="green", label=None): #{{{
    """
    For a given set of lineaments (lins), calculate the quality of their fits,
    over the range of their backrotation values, and plot it, weighting the fit
    by lineament length.

    """
    nb = len(lins[0].bs())
    aggregate_fits = zeros(nb)

    #filt_lins = lineament.filter(lins, dbar_max=dbar_max, delta_max=delta_max)

    for lin in lins:
        aggregate_fits += degrees(lin.delta_rms())*(lin.length())

    total_lin_length = array([ lin.length() for lin in lins ]).sum()

    aggregate_fits /= total_lin_length

    plot(linspace(0,180,nb),aggregate_fits, linewidth="2", color=color, label=label)
    ylabel("mean $\delta_{rms}$ [degrees]")
    xlabel("backrotation [degrees]")
    axis(ymin=0,ymax=90,xmin=0,xmax=18)
    xticks(range(0,181,10))
    grid(True)
    return(aggregate_fits)
#}}}

def plot_histagg(lins, color="green", rwidth=0.8, bins=18, label=None, dbar_max=0.05, delta_max=15): #{{{
    subplot(2,1,1)
    n, bins, patches = plot_fit_hist(calc_fit_hist(lins, dbar_max=dbar_max, delta_max=delta_max), color=color, rwidth=rwidth, bins=bins, label=label)
    legend()

    subplot(2,1,2)
    plot_aggregate_fits(lins, color=color, label=label, dbar_max=dbar_max, delta_max=delta_max)
    return bins
#}}}

def save_lins(lins, name="global_lins"): #{{{
    from time import strftime
    from pickle import dump

    dump(lins, open("output/%s_%s.pkl" % (name, strftime('%Y%m%d%H%M%S')),'w'))
#}}} end save_lins()

def fitlinmap(lins, dbar_max=0.05, delta_max=radians(15), lin_cm=cm.jet): #{{{
    """
    Plots the lineaments in lins, with those whose average MHD to their
    doppelgangers is greater than max_dop_mhd (as a fraction of their overall
    length) in gray.

    """
    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000

    goodlins = lineament.filter(lins, dbar_max=0.05, delta_max=delta_max)
    goodlen = sum( [ lin.length() for lin in goodlins ] ) * sat_radius_km

    badlins = [ lin for lin in lins if lin not in goodlins ]
    badlen = sum( [ lin.length() for lin in badlins ] ) * sat_radius_km

    print("found %d good NSR lins, (%d km)" % (len(goodlins),int(goodlen)))
    print("found %d bad NSR lins, (%d km)" % (len(badlins),int(badlen)))

    lines,fitmap = lineament.plotlinmap(badlins, color='black', alpha=1.0, lw=1.0)

    for lin in goodlins:
        # load some information about the lineament's best fits...
        best_delta = min(lin.delta_rms())
        best_b     = lin.bs(delta_rms=best_delta)

        # Map the color of the lineament to its best_b
        lin_color = lin_cm(int((lin_cm.N)*(best_b/pi)))

        # use line width to indicate goodness of best_fit
        lin_width = 5*(delta_max - best_delta)/delta_max
        lin_alpha = 1.0

        newline, fitmap = lineament.plotlinmap([lin], map=fitmap, lw=lin_width, color=lin_color, alpha=lin_alpha)
        lines.append(newline)

    cbax,kw = colorbar.make_axes(gca(), orientation="horizontal")
    colorbar.ColorbarBase(cbax, cmap=lin_cm, orientation="horizontal", norm=colors.Normalize(0,180))
    return lines, fitmap
    #}}}

def failure_orientations(lats=linspace(0,80,9), lons=linspace(0,180,1000)): #{{{
    """
    Plots rose diagrams showing the proportion of failure orientations that are
    taken on at various latitudes.

    """
    satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    shapefile = "input/GlobalLineaments"

    the_sat = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(the_sat),])

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

    return roses

#}}}

# Routines that generate the plots I have listed in the paper outline:

def makefigs(delta_max=radians(20), dbar_max=0.05): #{{{
    """
    A wrapper that takes input parameters for the figures, and calls all of the
    individual figure plotting routines, below.

    """

    # Assumes that all the stress comparisons have already been done, and that
    # the lineaments simply need to be loaded from pickled files.

    # These are the lineaments in their modern locations, compared to an
    # elastic NSR stress field:
    nsrlinfile = "output/global_20081112210635.pkl"
    nsrlins = pickle.load(open(nsrlinfile))

    # the stresscalc object has information about the satellite buried in it.
    nsr_stresscalc = nsrlins[0].stresscalc

    # All the mapped lineaments, transformed to a pre-TPW coordinate system,
    # with the paleo-north pole at 80E, 10N in current lat/lon.
    tpwlinfile = "output/tpw_80E_10N_20081112211815.pkl"
    tpwlins = pickle.load(open(tpwlinfile))

    # Lineaments for doing the viscous relaxation comparison:
#    D01  = pickle.load(open(nsrD01linfile))
#    D1   = pickle.load(open(nsrD1linfile))
#    D10  = pickle.load(open(nsrD010linfile))
#    D100 = pickle.load(open(nsrD0100linfile))
 
#    good_nsr = nsrlins[]
#    bad_nsr  = nsrlins[]
#    cycloid  = nsrlins[]
#    sinuous  = nsrlins[]
#    short_hilat = nsrlins[]
#    short_lolat = nsrlins[]

#    example_lins = [ good_nsr, bad_nsr, cycloid, sinuous, short_hilat, short_lolat ]

#    FitMap(nsrlins, delta_max, dbar_max)
#    TPWFitMap(tpwlins)

#    LinStats(nsrlins, delta_max, dbar_max)

#    FitCurveExamples(example_lins)
#    DoppelgangerExamples(example_lins)

    NSRvTPWActHist(nsrlins, tpwlins, delta_max, dbar_max)
#    ViscoActHist(D01, D1, D10, D100)
#    ActHistByLength(nsrlins)
#    MinDeltaRMSCorr(nsrlins)
#    DbarDeltaCorr(nsrlins)
#    NSRFailDirByLat(nsr_stresscalc)

# end makefigs }}}

def FitMap(nsrlins, delta_max, dbar_max): #{{{
    """
    Creates a global map of the lineaments, color coding them by what amount of
    backrotation (b) minimizes delta_rms(b) when compared to the NSR stress
    field.  Lineament width indicates the value of min(delta_rms(b)), with
    wider lineaments agreeing better with NSR stresses.

    Those lineaments which did not meet the criteria for inclusion in the
    analysis are plotted thin and black.

    The background is a grayscale map of the resolution of the images that went
    into the creation of the 500 meter resolution USGS mosaic of Europa, from
    which the lineaments were mapped.

    """
# end FitMap }}}

def TPWFitMap(tpwlins, delta_max, dbar_max): #{{{
    """
    Creates a global map of the lineaments, transformed into their pre-polar
    wander locations and orientations, color coded by what backrotation (b)
    minimizes delta_rms(b) when compared to the NSR stress field.  Lineament
    width indicates the value of min(delta_rms(b)), with wider lineaments
    indicating better agreement with the NSR stress field.

    Those lineaments which did not meet the criteria for inclusion in the
    analysis are plotted thin and black.

    The background is a grayscale map of the resolution of the images that went
    into the creation of the 500 meter resolution USGS mosaic of Europa, from
    which the lineaments were mapped, transformed so as to represent the
    pre-polar wander location of the surface.

    """
# end TPWFitMap }}}

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

def NSRvTPWActHist(nsrlins, tpwlins, delta_max, dbar_max, bins=18): #{{{
    """
    A plot showing the activity histograms and the aggregate fit curves for the
    set of lineaments included in each of the analyses.  One for the lineaments
    in their mapped locations/orientations, and one for the lineaments
    transformed to their pre-TPW locations/orientations.

    """

    # Calculate the histograms for the two sets of lineaments:
    nsrhistdata = calc_fit_hist(nsrlins, delta_max=delta_max, dbar_max=dbar_max, winwidth=18)
    tpwhistdata = calc_fit_hist(tpwlins, delta_max=delta_max, dbar_max=dbar_max, winwidth=18)

    # make that histogram data into the first half of the figure:
    fig = figure(1)
    ax1 = fig.add_subplot(2,1,1)
    (nsr_n,tpw_n) (nsr_bins,tpw_bins), (nsr_patches,tpw_patches) = ax1.hist([degrees(nsrhistdata), degrees(tpwhistdata)], bins=bins, range=(0,180))
    nsr_patches[0].set_label("modern")
    tpw_patches[0].set_label("pre-TPW")
    for nsr_bar,tpw_bar in zip(nsr_patches, tpw_patches):
        nsr_bar.set_facecolor((.5,.5,.5,1)) # gray
        tpw_bar.set_facecolor((0,0,0,1))    # black
    ax1.set_xticks(arange(0,181,10))
    ax1.axis(xmin=0, xmax=180)
    ax1.grid(True)
    ax1.legend()
    ax1.set_ylabel('aggregate feature length [km]')
    ax1.set_xlabel('backrotation b [degrees]')
    ax1.set_title('Lineament Formation Activity, $\max(\delta_{rms})=$%g$^\circ$, $\max(\\bar{D})=$%g' % (degrees(delta_max),dbar_max))

    # Filter the input lineaments to keep only the good fits:
    nsrkeepers = lineament.filter(nsrlins, delta_max=delta_max, dbar_max=dbar_max)
    tpwkeepers = lineament.filter(tpwlins, delta_max=delta_max, dbar_max=dbar_max)

    # Calculate the aggregate fits:
    nsr_nb = len(nsrkeepers[0].bs())
    nsr_agg_fits = zeros(nsr_nb)
    total_nsr_linlen = 0
    for lin in nsrkeepers:
        linlen = lin.length()
        nsr_agg_fits += lin.delta_rms()*linlen
        total_nsr_linlen += linlen
    nsr_agg_fits /= total_nsr_linlen

    tpw_nb = len(tpwkeepers[0].bs())
    tpw_agg_fits = zeros(tpw_nb)
    total_tpw_linlen = 0
    for lin in tpwkeepers:
        linlen = lin.length()
        tpw_agg_fits += lin.delta_rms()*linlen
        total_tpw_linlen += linlen
    tpw_agg_fits /= total_tpw_linlen

    ax2 = fig.add_subplot(2,1,2)
    nsr_line, tpw_line = ax2.plot(linspace(0,180,nsr_nb), degrees(nsr_agg_fits), linspace(0,180,tpw_nb), degrees(tpw_agg_fits))
    nsr_line.set(color=(0,0,0,1), linewidth=2, linestyle='--', label='modern')
    tpw_line.set(color=(.5,.5,.5,1), linewidth=2, linestyle='-', label='pre-TPW')
    ax2.set_ylabel("mean $\delta_{rms}$ [degrees]")
    ax2.set_xlabel("backrotation b [degrees]")
    ax2.set_title("aggregate $\delta_{rms}$ of all features")
    ax2.set_xticks(arange(0,181,10))
    ax2.legend(loc='lower right')
    ax2.grid(True)

    show()
# end NSRvTPWActHist }}}

def ViscoActHist(D01, D1, D10, D100, delta_max, dbar_max): #{{{
    """
    Activity histogram and aggregate fit curves for the lineaments compared to NSR
    stress fields for three values of Delta, representing four different viscous
    relaxation scenarios, ranging from elastic, to viscous: Delta={0.1,1,10,100}

    """
# end ViscoActHist }}}

def ActHistByLength(nsrlins, delta_max, dbar_max): #{{{
    """
    Activity histogram and aggregate fit curves for the lineaments compared to
    the NSR stress field, but with the lineaments broken into several different
    length categories, showing the degree to which the result is supported by
    all the different length categories.

    """
# end ActHistByLength }}}

def MinDeltaRMSCorr(nsrlins): #{{{
    """
    Several scatter plots showing the correlation (or lack thereof) between
    min(delta_rms) and other lineament characteristics, including:

    A: lineament.sinuosity()
    B: lineament.length()
    C: min(lineament.latitudes()), color coded by length?
    D: 
    """
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
 
def NSRFailDirByLat(nsr_stresscalc): #{{{
    """
    Graphic display of how the range of possible failure orientations changes
    with latitude in the NSR stress field, demonstrating that for any short
    lineament above 30 degrees latitude, there will always be a good fit.

    """
# end NSRFailDirByLat }}}
