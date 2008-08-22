#!/usr/bin/python
"""
Routines for comparing a set of lineaments to an NSR stress field.

"""

from . import lineament
from . import satstress
from pylab import *
from mpl_toolkits.basemap import Basemap
from numpy.ma.mstats import idealfourths
import random

def overnight_run(lins=None, nb=19,\
                  satfile   = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                  shapefile = "input/GlobalLineaments",\
                  doppel_fits = False):
#{{{
    the_sat = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(the_sat),])

    if lins is None:
        lins = lineament.shp2lins(shapefile)

    pnp_lon=170
    pnp_lat=60
    lins = [ lin.poleshift(radians(pnp_lon), radians(pnp_lat)) for lin in lins ]

    # Just for testing that this works...
    # lins = lins[:5]

    for n in range(len(lins)):
        print("Calculating fits for lineament %d / %d" % (n+1,len(lins)))
        lins[n].calc_fits(nsr_stresscalc,nb=nb)

    save_lins(lins)
    plot_histagg(lins, label="PNP: %dE, %dN"%(pnp_lon,pnp_lat), bins=18)

    return(lins)
#}}}

def calc_fit_hist(lins): #{{{
    linhist = []
    for lin,N in zip(lins,range(len(lins))):

        # Just in case we haven't already calculated the better fits:
        if lin.better_fits is None or lin.better_bs is None:
            print("Calculating histogram contribution for lineament %d / %d" % (N+1,len(lins)))
            lin.good_doppel_fits(fit_thresh=0.1, max_dop_mhd=0.1, window=15)

        # We generate the histogram dataset, adding one "count" for each km of
        # lineament that fit best at that value of b(ackrotation)
        for b,fit in zip(lin.better_bs,lin.better_fits):
            print("    b   = %g" % (degrees(b),))
            print("    fit = %g" % (degrees(fit),))
            linlenkm = int(lin.length()*lin.stresscalc.stresses[0].satellite.radius()/1000)
            linhist += linlenkm * [b,]
        N+=1

    return(linhist)
# }}}

def plot_fits(lins, fixlons=True): #{{{
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
        if fixlons:
            linmap.plot(degrees(lin.fixed_longitudes()),\
                        degrees(lin.latitudes()),\
                        c='black', ls='-')
        else:
            linmap.plot(degrees(lin.longitudes()),\
                        degrees(lin.latitudes()),\
                        c='black', ls='-')

        # Plot the mapped lineament's fit curve:
        subplot(2,1,2)
        title("best_fit=%g, best_b=%g" % (degrees(lin.best_fit()), degrees(lin.best_b())))
        axis([0,180,0,90])
        xticks( range(0,181,15) )
        yticks( range(0,91,15) )
        grid(True)
        xlabel("degrees of backrotation")
        ylabel("RMS weighted misfit (degrees)")
        plot(linspace(0,180,num=len(lin.fits)), degrees(lin.fits), c='black', ls='-', lw=2)
        if lin.better_fits is None or lin.better_bs is None:
            lin.good_doppel_fits()

        if len(lin.better_bs) > 0:
            scatter(degrees(lin.better_bs), degrees(lin.better_fits),s=150,c='green',marker='o')

        # Mark on the map the band of acceptably good fits:
        fit_min = degrees(min(lin.fits))
        fit_max = degrees(max(lin.fits))
        fit_amp = fit_max - fit_min
        fitband_xs = [0,180,180,0]
        fitband_ys = [fit_min, fit_min, fit_min+fit_amp*0.1, fit_min+fit_amp*0.1]
        fill(fitband_xs, fitband_ys, facecolor='green', alpha=0.20)

        # Mark on the plot where we had good fits that generated bad doppelgangers
        okay_bs, okay_fits = lin.good_fits()
        for b,fit in zip(okay_bs, okay_fits):
            if b not in lin.better_bs:
                scatter(degrees([b,]), degrees([fit,]), s=150, c='red',marker='o')

        # Now deal with the doppelgangers, if they exist:
        if len(lin.doppels) > 0:
            for dop in lin.doppels:
                if dop.backrot in lin.better_bs:
                    # Map view:
                    subplot(2,1,1)
                    if fixlons:
                        linmap.plot(degrees(dop.lonshift(-dop.backrot).fixed_longitudes()),\
                                    degrees(dop.lonshift(-dop.backrot).latitudes()),\
                                    ls='-', c='black', alpha=0.5)
                    else:
                        linmap.plot(degrees(dop.lonshift(-dop.backrot).longitudes()),\
                                    degrees(dop.lonshift(-dop.backrot).latitudes()),\
                                    ls='-', c='black', alpha=0.5)
                    # Fit curve, if they've got one:
                    if len(dop.fits) > 0:
                        subplot(2,1,2)
                        plot(linspace(0,180,num=len(dop.fits)), degrees(dop.fits), ls='--', lw=2)
        #show()

        # wait for the user to hit return before we continue drawing...
        x = raw_input("press return for next fit: ")
        clf()
        lin_num += 1
#}}}

def plot_fit_hist(fithist, bins=36, hist_title="Lineament formation vs. backrotation", color="green", rwidth=0.8, label=None): #{{{
    hist(fithist, bins=bins, facecolor=color, rwidth=rwidth, range=(0,180), label=label)
    grid(True)
    ylabel("total lineaments formed [km]")
    title(hist_title)
    xticks(arange(0,181,5))
    axis(xmin=0, xmax=180)
#}}}

def plot_aggregate_fits(lins, color="green", label=None): #{{{
    """
    For a given set of lineaments (lins), calculate the quality of their fits,
    over the range of their backrotation values, and plot it, weighting the fit
    by lineament length.

    """
    nb = len(lins[0].fits)
    aggregate_fits = zeros(nb)
    for lin in lins:
        aggregate_fits += degrees(array(lin.fits))*(lin.length())

    total_lin_length = array([ lin.length() for lin in lins ]).sum()

    aggregate_fits /= total_lin_length

    plot(linspace(0,180,nb),aggregate_fits, linewidth="2", color=color, label=label)
    ylabel("mean RMS delta [degrees]")
    xlabel("backrotation [degrees]")
    axis(ymin=0,ymax=90,xmin=0,xmax=18)
    xticks(range(0,181,10))
    grid(True)
    return(aggregate_fits)
#}}}

def plot_histagg(lins, color="green", rwidth=0.8, bins=36, label=None): #{{{
    subplot(2,1,1)
    plot_fit_hist(degrees(calc_fit_hist(lins)), color=color, rwidth=rwidth, bins=bins, label=label)
    legend()

    subplot(2,1,2)
    plot_aggregate_fits(lins, color=color, label=label)
#}}}

def save_lins(lins, name="global_lins"): #{{{
    from time import strftime
    from pickle import dump

    dump(lins, open("output/%s_%s.pkl" % (name, strftime('%Y%m%d%H%M%S')),'w'))
#}}} end save_lins()

