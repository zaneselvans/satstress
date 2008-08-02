#!/usr/bin/python
"""
Routines for comparing a set of lineaments to an NSR stress field.

"""

from satstress import *
from pylab import *
from mpl_toolkits.basemap import Basemap

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

def draw_fits(lins, sat_radius_km, stresscalc, doppels_E=[], doppels_W=[]): #{{{
    lin_num=0

    clf()

    # Make sure we got the same number of doppelganger lineaments as
    # mapped lineaments, if we're not generating them on the fly
    if len(doppels_E) > 0:
        assert len(doppels_E) == len(lins)
    else: # Create the array of doppelgangers:
        doppels_E = [ lin.doppelganger(stresscalc, propagation="east") for lin in lins ]

    if len(doppels_W) > 0:
        assert len(doppels_W) == len(lins)
    else: # Create the array of doppelgangers:
        doppels_W = [ lin.doppelganger(stresscalc, propagation="west") for lin in lins ]

    for (lin,doppel_E,doppel_W) in zip(lins,doppels_E,doppels_W):
        subplot(2,1,1)
        linmap = Basemap(rsphere=sat_radius_km*1000, projection='cyl')
        linmap.drawmeridians(range(-180,180,30), labels=[1,0,0,1])
        linmap.drawparallels(range(-90,90,30), labels=[1,0,0,1])
        # Draw lineament in a Basemap map
        linmap.plot(degrees(lin.longitudes()),\
                    degrees(lin.latitudes()))

        # Draw Doppelgangers in the Basemap map
        if doppel_E is not None:
            linmap.plot(degrees(doppel_E.lonshift(-lin.best_b()).longitudes()),\
                        degrees(doppel_E.lonshift(-lin.best_b()).latitudes()),\
                        'r--')

        if doppel_W is not None:
            linmap.plot(degrees(doppel_W.lonshift(-lin.best_b()).longitudes()),\
                        degrees(doppel_W.lonshift(-lin.best_b()).latitudes()),\
                        'g--')

        title("lin=%d, length=%g km" % (lin_num, sat_radius_km*lin.length()))

        # Now plot the fits for both the lineaments and the doppelgangers
        subplot(2,1,2)
        plot(linspace(0,180,num=len(lin.fits)), degrees(lin.fits), 'b-', linewidth=2)
        #plot(linspace(0,181,10), 10*[degrees(lin.best_fit()),], 'r--')
        #plot(10*[degrees(lin.best_b()),], linspace(0,91,10), 'r--')

        if doppel_E is not None:
            if len(doppel_E.fits) > 0:
                plot(linspace(0,180,num=len(doppel_E.fits)), degrees(doppel_E.fits), 'r--', linewidth=2)
        if doppel_W is not None:
            if len(doppel_W.fits) > 0:
                plot(linspace(0,180,num=len(doppel_W.fits)), degrees(doppel_W.fits), 'g--', linewidth=2)

        grid(True)

        axis([0,180,0,90])
        xticks( range(0,181,15) )
        yticks( range(0,91,15) )
        xlabel("degrees of backrotation")
        ylabel("average misfit (degrees)")
        title("best_fit=%g, best_b=%g, iqr=%g" % (degrees(lin.best_fit()), degrees(lin.best_b()), degrees(lin.fits_width()) ))
        show()

        # wait for the user to hit return before we continue drawing...
        x = raw_input("press return for next fit: ")
        clf()
        lin_num += 1
# }}}

def draw_nsr_fit_hist(fithist, bins=18, hist_title="Lineament formation vs. backrotation"): #{{{
    clf()
    hist(fithist, bins=bins, facecolor="green")
    grid(True)
    xlabel("degrees of backrotation")
    ylabel("km of lineaments formed")
    title(hist_title)
    axis(xmin=0, xmax=180)
    #}}}

def doppel_fits(numlins=10, nb=18,\
         satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
         shapefile="input/GlobalLineaments"):

    """
    Given a list of lineaments:
       - calculate their fits to an elastic NSR stress field
       - create East and West doppelgangers for them
       - calculate the fits for those doppelgangers
       - return the three lists of lineaments for later plotting

    """

    europa = satstress.Satellite(open(satfile,'r'))

    elastic_nsr = satstress.StressCalc([satstress.NSR(europa),])
    global_lins = lineament.shp2lins(shapefile)
    global_lins = [lin for lin in global_lins if len(lin.vertices) > 1]
    if numlins == 0:
        test_lins = global_lins
        numlins = len(global_lins)
    else:
        test_lins = global_lins[:numlins]

    N=1
    for lin in test_lins:
        print "fitting lin %d of %d: length=%g km; segs=%d; nb=%d" % (N,numlins, lin.length()*europa.radius()/1000,len(lin.segments()), nb)
        lin.calc_fits(elastic_nsr, nb)
        print "    best_fit = %g deg" % (degrees(lin.best_fit()))
        print "    best_b   = %g deg" % (degrees(lin.best_b()))
        N+=1

    doppels_E = [ lin.doppelganger(elastic_nsr, propagation="east") for lin in test_lins ]
    doppels_W = [ lin.doppelganger(elastic_nsr, propagation="west") for lin in test_lins ]

    N=1
    for lin, dop_E, dop_W in zip(test_lins, doppels_E, doppels_W):
        print "fitting doppelgangers for lin %d of %d" % (N, numlins)
        if dop_E is not None:
            dop_E_shifted = dop_E.lonshift(-lin.best_b())
            dop_E_shifted.calc_fits(elastic_nsr, nb)
            dop_E.fits = dop_E_shifted.fits
        else:
            print "    dop_E was None"
        if dop_W is not None:
            dop_W_shifted = dop_W.lonshift(-lin.best_b())
            dop_W_shifted.calc_fits(elastic_nsr, nb)
            dop_W.fits = dop_W_shifted.fits
        else:
            print "    dop_W was None"

        lin.doppel_E = dop_E
        lin.doppel_W = dop_W
        N+=1

    #draw_fits(test_lins, europa.radius()/1000, elastic_nsr, doppels_E, doppels_W)
    return(test_lins)

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
# - improve Lineament.draw() to use Basemap, and a "doppel" style
# - figure out what to do about weighting functions, play around with it
# - devise measure of statistical significance:
#   - monte carlo method
#   - analysis of fit-curve (IQR, etc)
# - Better binning of "best fits" in histogram
#   - add lineament length to all bins in which fit is better than X?
#   - or add length/N (where N is the # of bins in which it's better than X)
