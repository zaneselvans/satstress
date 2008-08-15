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

    # Just for testing that this works...
    #lins = lins[:5]

    for n in range(len(lins)):
        print("Calculating fits for lineament %d..." % (n,))
        lins[n].calc_fits(nsr_stresscalc,nb=nb)

    save_lins(lins)
    #plot_histagg(lins)

    return(lins)
#}}}

def calc_fit_hist(lins): #{{{
    linhist = []
    for lin,N in zip(lins,range(len(lins))):
        print("Calculating histogram contribution for lineament %d" % (N,))
        # We generate the histogram dataset, adding one "count" for each km of
        # lineament that fit best at that value of b(ackrotation)
        good_bs, good_fits = lin.good_doppel_fits(fit_thresh=0.1, max_dop_mhd=0.1, window=15)
        for b,fit in zip(good_bs,good_fits):
            print("    b   = %g" % (degrees(b),))
            print("    fit = %g" % (degrees(fit),))
            linlenkm = int(lin.length()*lin.stresscalc.stresses[0].satellite.radius()/1000)
            linhist += linlenkm * [b,]
        N+=1

    return(linhist)
# }}}

def plot_fits(lins): #{{{
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
        good_bs, good_fits = degrees(lin.good_doppel_fits(fit_thresh=0.1, window=15, max_dop_mhd=0.1))
        if len(good_bs) > 0:
            scatter(good_bs, good_fits)

        # Now deal with the doppelgangers, if they exist:
        if len(lin.doppels) > 0:
            for dop in lin.doppels:
                if degrees(dop.backrot) in good_bs:
                    # Map view:
                    subplot(2,1,1)
                    linmap.plot(degrees(dop.lonshift(-dop.backrot).longitudes()),\
                                degrees(dop.lonshift(-dop.backrot).latitudes()),\
                                ls='--')
                    # Fit curve, if they've got one:
                    if len(dop.fits) > 0:
                        subplot(2,1,2)
                        plot(linspace(0,180,num=len(dop.fits)), degrees(dop.fits), ls='--', lw=2)
        show()

        # wait for the user to hit return before we continue drawing...
        x = raw_input("press return for next fit: ")
        clf()
        lin_num += 1
#}}}

def plot_fit_hist(fithist, bins=36, hist_title="Lineament formation vs. backrotation", color="green", rwidth=0.8): #{{{
    hist(fithist, bins=bins, facecolor=color, rwidth=rwidth, align="mid")
    grid(True)
    ylabel("total lineaments formed [km]")
    title(hist_title)
    xticks(arange(0,181,5))
    axis(xmin=0, xmax=180)
#}}}

def plot_aggregate_fits(lins, color="green"): #{{{
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

    # this doubles the waveform so we can see more clearly its width...
    #plot(linspace(0,360,nb*2-1),hstack([aggregate_fits,aggregate_fits[1:]]))

    plot(linspace(0,180,nb),aggregate_fits, linewidth="2", color=color)
    ylabel("mean RMS delta [degrees]")
    xlabel("backrotation [degrees]")
    axis(ymin=0,ymax=90,xmin=0,xmax=18)
    xticks(range(0,181,10))
    grid(True)
    return(aggregate_fits)
#}}}

def plot_histagg(lins, color="green", rwidth=0.8, bins=36): #{{{
    subplot(2,1,1)
    plot_fit_hist(calc_fit_hist(lins), color=color, rwidth=rwidth, bins=bins)
    subplot(2,1,2)
    plot_aggregate_fits(lins, color=color)
#}}}

def split_by_smhd(lins, smhd_ratio=0.1): #{{{
    """
    Within a set of lineaments lins, discriminate between those whose East and
    West doppelgangers are similar in shape (having a symmetric mean Hausdorff
    distance from each other of less than smhd_ratio*lin.length()), vs.  those
    whose doppelgangers are different from each other.  If one of the
    doppelganger lineaments is None, the lineament is put in the "bad" bin.
    Returns two lists of lineaments (good_lins, bad_lins).

    """

    good_lins = []
    bad_lins = []

    sat_rad_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000

    for lin in lins:
        if lin.doppel_E is not None and lin.doppel_W is not None and lin.doppel_E.smhd(lin.doppel_W)/lin.length() < smhd_ratio:
            good_lins.append(lin)
        else:
            bad_lins.append(lin)

    return(good_lins, bad_lins)

#}}}

def save_lins(lins, name="global_lins"): #{{{
    from time import strftime
    from pickle import dump

    dump(lins, open("output/%s_%s.pkl" % (name, strftime('%Y%m%d%H%M%S')),'w'))
#}}} end save_lins()

def test_leastsq(p0): #{{{
    from scipy.optimize import leastsq  
    x = arange(0,2*pi,pi/100)
    A,theta = 90.0, 0.0
    y_true = A*sin(x+theta)  
    y_meas = y_true + (A/10)*randn(len(x))
    
    plsq = leastsq(residuals, p0, args=(y_meas, x))  

    plot(degrees(x),peval(x,plsq[0]),'g--',degrees(x),y_meas,'bo',degrees(x),y_true,'r-')  
    xlim(0,360)
    xticks(arange(0,361,45))
    grid()
    title("Least squares fit to noisy data")

    return(plsq)
 
def residuals(p, y, x):
    A,theta = p  
    err = y-A*sin(x+theta)  
    return err  
 
def peval(x, p):
    A,theta = p
    return A*sin(x+theta)  
#}}}

def fit2hist(fits, bins=20):

    # for a given range of fit values, how many points are within it?
    #thresholds = linspace(min(fits),max(fits),bins)
    thresholds = linspace(0,pi/2,bins)
    fithist=[]

    for n_t in range(len(thresholds)-1):
        fithist.append(len([ x for x in fits if x > thresholds[n_t] and x < thresholds[n_t+1] ]))

    return(fithist)

def showfits(lins):
    for lin in lins:
        clf()
        hist(fit2hist(lin.fits), bins=20)
        raw_input("press return: ")
