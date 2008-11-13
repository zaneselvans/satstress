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

def plot_fits(lins, fixlon=True, delta_max=radians(20), dbar_max=0.05): #{{{
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
            lin.plotdoppels(bs=good_bs, map=linmap, backrot=True, fixlon=fixlon, alpha=0.5)


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
    max_delta and dbar < max_dbar.

    Assumes that the lineaments have already had their fits calculated.

    """

    linhist = []
    for lin in lins:
        # Select the fits which meet our criteria for this lineament:
        good_fits = lin.fit(bs=lin.bs(delta_ismin=True, delta_max=delta_max,\
                                      dbar_max=dbar_max, winwidth=winwidth))

        # if there were any sufficiently good fits:
        if len(good_fits) == 0:
            continue

        # We generate the histogram dataset, adding one "count" for each km of
        # lineament that fit best at that value of b(ackrotation)
        for b, delta, dbar in good_fits:
            linlenkm = int(lin.length()*lin.stresscalc.stresses[0].satellite.radius()/1000)
            print("     b = %g" % (degrees(b),))
            print(" delta = %g" % (degrees(delta),))
            print("  dbar = %g" % (dbar))
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

def plot_aggregate_fits(lins, color="green", label=None, dbar_max=0.05, delta_max=15): #{{{
    """
    For a given set of lineaments (lins), calculate the quality of their fits,
    over the range of their backrotation values, and plot it, weighting the fit
    by lineament length.

    """
    nb = len(lins[0].bs())
    aggregate_fits = zeros(nb)

    filt_lins = lineament.filter(lins, dbar_max=dbar_max, delta_max=delta_max)

    for lin in filt_lins:
        aggregate_fits += degrees(lin.delta_rms())*(lin.length())

    total_lin_length = array([ lin.length() for lin in filt_lins ]).sum()

    aggregate_fits /= total_lin_length

    plot(linspace(0,180,nb),aggregate_fits, linewidth="2", color=color, label=label)
    ylabel("mean RMS delta [degrees]")
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
