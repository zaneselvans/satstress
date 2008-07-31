#!/usr/bin/python
"""
Routines for comparing a set of lineaments to an NSR stress field.

"""

from satstress import *
from pylab import *

def calc_nsr_fits(satfile = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                  shapefile = "input/GlobalLineaments",\
                  nb = 180,\
                  numlins = "all",\
                  w_length = True,\
                  w_anisotropy = False,\
                  w_stressmag = False): #{{{
    
    # Define our satellite and stress field based on the input files:
    the_sat = satstress.Satellite(open(satfile,'r'))
    nsr = satstress.StressCalc([satstress.NSR(the_sat),])

    # This is useful just to keep things easy to read
    sat_radius_km = nsr.stresses[0].satellite.radius()/1000.0

    # Turn the shapefile into a list of lineaments:
    lins = lineament.shp2lins(shapefile)
    # By default do all the lineaments... but allow us to shorten for testing
    if numlins != "all":
        lins = lins[:numlins]

    # Do the NSR-Lineament comparison and store the results in the Lineaments
    # This can be time consuming depending on how many lineaments we've got,
    # and how detailed of a calculation we're doing, so it's polite to print
    # out some progress indication... also build up the fit histogram.
    lin_num = 0
    linhist=[]
    for lin in lins:
        print "lineament %d (%d segments, %g km)" % (lin_num, len(lin.vertices), lin.length()*sat_radius_km)
        lin.calc_fits(nsr, nb, failuremode="tensilefracture", w_length=w_length, w_anisotropy=w_anisotropy, w_stressmag=w_stressmag)
        lin.sat_radius_km = sat_radius_km
        lin_num += 1

    return (lins)
# }}}

def calc_fit_hist(lins): #{{{
    linhist = []
    for lin in lins:
        # We generate the histogram dataset, adding one "count" for each km of
        # lineament that fit best at that value of b(ackrotation)
        linhist += int(lin.length()*lin.sat_radius_km)*[degrees(lin.best_b()),]

    return(linhist)

def draw_fit_hist(linhist, bins=18):
    hist(linhist, bins=bins)
    axis(xmin=0, xmax=180)
    xlabel("degrees of backrotation")
    ylabel("km of lineaments")
    title("Lineament formation vs. best fit backrotation")
    show()
# }}}

def draw_fits(lins, sat_radius_km, stresscalc, doppels_E=None, doppels_W=None): #{{{
    lin_num=0

    # Make sure we got the same number of doppelganger lineaments as
    # mapped lineaments, if we're not generating them on the fly
    if doppels_E is not None:
        assert len(doppels_E) == len(lins)
    if doppels_W is not None:
        assert len(doppels_W) == len(lins)

    for lin in lins:
        # Plot the lineament itself:
        mapped_lons = degrees([ v[0] for v in lin.vertices ])
        mapped_lats = degrees([ v[1] for v in lin.vertices ])

        subplot(211)
        if doppels_E is None:
            doppel_E = lin.doppelganger(stresscalc, propagation="east")

        if doppel_E is not None:
            doppel_E = doppel_E.lonshift(-lin.best_b())
            doppel_E_lons = degrees( [ v[0] for v in doppel_E.vertices ] )
            doppel_E_lats = degrees( [ v[1] for v in doppel_E.vertices ] )
            plot(doppel_E_lons, doppel_E_lats, 'k--', linewidth=1)

        doppel_W      = lin.doppelganger(stresscalc, propagation="west")
        if doppel_W is not None:
            doppel_W = doppel_W.lonshift(-lin.best_b())
            doppel_W_lons = degrees( [ v[0] for v in doppel_W.vertices ] )
            doppel_W_lats = degrees( [ v[1] for v in doppel_W.vertices ] )
            plot(doppel_W_lons, doppel_W_lats, 'k--', linewidth=1)

        plot(mapped_lons, mapped_lats, 'k-', linewidth=2)
        title("lin=%d, length=%g km" % (lin_num, sat_radius_km*lin.length()))
        xlabel("longitude")
        ylabel("latitude")
        grid(True)

        subplot(212)
        plot(linspace(0,180,num=len(lin.fits)), degrees(lin.fits), 'b-', linewidth=2)
        plot(linspace(0,181,10), 10*[degrees(lin.best_fit()),], 'r--', linewidth=2)
        plot(10*[degrees(lin.best_b()),], linspace(0,91,10), 'r--', linewidth=2)
        grid(True)

        axis([0,180,0,90])
        xticks( range(0,181,15) )
        yticks( range(0,91,15) )
        xlabel("degrees of backrotation")
        ylabel("average misfit (degrees)")
        title("best_fit=%g, best_b=%g, iqr=%g" % (degrees(lin.best_fit()), degrees(lin.best_b()), degrees(lin.fits_width()) ))
        show()
        x = raw_input("press return for next fit: ")
        clf()
        lin_num += 1
# }}}

def draw_nsr_fit_hist(fithist, bins=18, hist_title="Lineament formation vs. backrotation"):
    clf()
    hist(fithist, bins=bins, facecolor="green")
    grid(True)
    xlabel("degrees of backrotation")
    ylabel("km of lineaments formed")
    title(hist_title)
    axis(xmin=0, xmax=180)


# What do I want this tool to do, anyway?
# - assume for the moment that storing all of the stresses for the entire set
#   of backrotated fits isn't something I want to do... meaning that I can't
#   store all of the possible fits and their associated weights.
#
# - Given (lineaments, stress field):
#   - calculate fits to stress field over 180 degrees of backrotation
#   - for the best fit backrotation, create both and east and a west doppelganger
#   - calculate fits for the doppelgangers
#   - display all three lineaments and their fit curves
#     - use 2 panel plot (fits, lineaments)
#     - use Basemap for the lineaments, so we can have all the mappy goodness.
#   - save all the fit data for later reference
#
# TODO:
# =====
# - figure out what to do about weighting functions, play around with it
# - devise measure of statistical significance:
#   - monte carlo method
#   - analysis of fit-curve (IQR, etc)
# - Better binning of "best fits" in histogram
#   - add lineament length to all bins in which fit is better than X?
#   - or add length/N (where N is the # of bins in which it's better than X)
