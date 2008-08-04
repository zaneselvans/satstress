#!/usr/bin/python
"""
Routines for comparing a set of lineaments to an NSR stress field.

"""

from satstress import *
from pylab import *
from mpl_toolkits.basemap import Basemap
import random

def calc_nsr_fits(lins=None, numlins=0, nb=19, metric="rms", satfile="input/ConvectingEuropa_GlobalDiurnalNSR.ssrun", shapefile="input/GlobalLineaments", w_length=True, w_stress=False, min_length=0.3, doppel_fits=False): #{{{

    """
    Given a list of lineaments:
       - calculate their fits to an elastic NSR stress field
       - create East and West doppelgangers for them
       - calculate the fits for those doppelgangers
       - return the original list of lineaments, with its fits and doppels

    """

    the_sat = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(the_sat),])

    if lins is None:
        # If we didn't get a list of lineaments passed in, read them from the shapefile
        lins = lineament.shp2lins(shapefile)

    # Make sure we don't have any lineaments that are actually just points
    lins = [ lin for lin in lins if len(lin.vertices) > 1 ]
    # make sure they're long enough to be interesting...
    lins = [ lin for lin in lins if lin.length() > min_length ]

    # if numlins is zero, that means do all the lineaments
    if numlins==0:
        numlins=len(lins)
    # if it's not zero, we select numlins random lineaments from the list of lineaments
    else:
        lins = random.sample(lins,numlins)

    for lin,N in zip(lins,range(len(lins))):
        print "fitting lin %d of %d: length=%g km; segs=%d; nb=%d" % (N+1,numlins, lin.length()*the_sat.radius()/1000,len(lin.segments()), nb)
        lin.calc_fits(nsr_stresscalc, nb, metric=metric, w_length=w_length, w_stress=w_stress)

        print "    best_fit = %g deg" % (degrees(lin.best_fit()))
        print "    best_b   = %g deg" % (degrees(lin.best_b()))
        N+=1

    doppels_E = [ lin.doppelganger(nsr_stresscalc, propagation="east") for lin in lins ]
    doppels_W = [ lin.doppelganger(nsr_stresscalc, propagation="west") for lin in lins ]

    for lin, dop_E, dop_W, N in zip(lins, doppels_E, doppels_W, range(len(lins))):
        print "making doppelgangers for lin %d of %d" % (N+1, numlins)
        if dop_E is not None:
            dop_E_shifted = dop_E.lonshift(-lin.best_b())
            if doppel_fits is True:
                dop_E_shifted.calc_fits(nsr_stresscalc, nb, metric=metric, w_length=w_length, w_stress=w_stress)
                dop_E.fits = dop_E_shifted.fits
        else:
            print "    dop_E was None"
        if dop_W is not None:
            dop_W_shifted = dop_W.lonshift(-lin.best_b())
            if doppel_fits is True:
                dop_W_shifted.calc_fits(nsr_stresscalc, nb, metric=metric, w_length=w_length, w_stress=w_stress)
                dop_W.fits = dop_W_shifted.fits
        else:
            print "    dop_W was None"

        lin.doppel_E = dop_E
        lin.doppel_W = dop_W
        N+=1

    return(lins)
# }}} end doppel_fits

def calc_fit_hist(lins): #{{{
    linhist = []
    for lin in lins:
        # We generate the histogram dataset, adding one "count" for each km of
        # lineament that fit best at that value of b(ackrotation)
        linhist += int(lin.length()*lin.stresscalc.stresses[0].satellite.radius()/1000)*[degrees(lin.best_b()),]

    return(linhist)

def draw_fit_hist(linhist, bins=18):
    hist(linhist, bins=bins)
    axis(xmin=0, xmax=180)
    xlabel("degrees of backrotation")
    ylabel("km of lineaments")
    title("Lineament formation vs. best fit backrotation")
    show()
# }}}

def draw_fits(lins): #{{{
    lin_num=0
    doppels_E = [ lin.doppel_E for lin in lins ]
    doppels_W = [ lin.doppel_W for lin in lins ]

    clf()

    stresscalc = lins[0].stresscalc
    sat_radius_km = lins[0].stresscalc.stresses[0].satellite.radius()/1000
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
        ylabel("%s weighted misfit (degrees)" % (lin.metric,) )
        title("best_fit=%g, best_b=%g" % (degrees(lin.best_fit()), degrees(lin.best_b())))
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


def overnight_run():
    import pickle
    global_lins=lineament.shp2lins("input/GlobalLineaments")

    linear_stress = calc_nsr_fits(lins=global_lins, metric="linear", w_stress=True)
    pickle.dump(linear_stress,open("linear_stress.pkl",'w'))
    subplot(2,2,1)
    title("metric=linear; w_stress=True")
    draw_fit_hist(calc_fit_hist(linear_stress))

    rms_stress = calc_nsr_fits(lins=global_lins, metric="rms", w_stress=True)
    pickle.dump(rms_stress,open("rms_stress.pkl",'w'))
    subplot(2,2,2)
    title("metric=RMS; w_stress=True")
    draw_fit_hist(calc_fit_hist(rms_stress))

    linear_nostress = calc_nsr_fits(lins=global_lins, metric="linear", w_stress=False)
    pickle.dump(linear_nostress,open("linear_nostress.pkl",'w'))
    subplot(2,2,3)
    title("metric=linear; w_stress=False")
    draw_fit_hist(calc_fit_hist(linear_nostress))

    rms_nostress = calc_nsr_fits(lins=global_lins, metric="rms", w_stress=False)
    pickle.dump(rms_nostress,open("rms_nostress.pkl",'w'))
    subplot(2,2,4)
    title("metric=RMS; w_stress=False")
    draw_fit_hist(calc_fit_hist(rms_nostress))

    return(linear_stress, rms_stress, linear_nostress, rms_nostress)

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
# - Better binning of "best fits" in histogram
#   - add lineament length to all bins in which fit is better than X?
#   - or add length/N (where N is the # of bins in which it's better than X)
#   - develop better ways to discriminate which fits are real
#     - SMHD of lin vs. doppels for local minima
#     - MC method
