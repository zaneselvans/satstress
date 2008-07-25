#!/usr/bin/python
"""
A script to compare a set of lineaments to a given NSR stress field.

"""

from satstress import *
from pylab import *

def calc_nsr_fit_hist(satfile = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                      shapefile = "input/GlobalLineaments",\
                      nb = 180,\
                      numlins = "all"):
    
    # By default do all the lineaments... but allow us to shorten for testing
    if numlins != "all":
        lins = lins[:numlins]

    # Define our satellite and stress field based on the input files:
    the_sat = satstress.Satellite(open(satfile,'r'))
    nsr = satstress.StressCalc([satstress.NSR(the_sat),])

    # This is useful just to keep things easy to read
    sat_radius_km = nsr.stresses[0].satellite.radius()/1000.0

    # Turn the shapefile into a list of lineaments:
    lins = lineament.shp2lins(shapefile)

    # Do the NSR-Lineament comparison and store the results in the Lineaments
    # This can be time consuming depending on how many lineaments we've got,
    # and how detailed of a calculation we're doing, so it's polite to print
    # out some progress indication... also build up the fit histogram.
    lin_num = 0
    linhist=[]
    for lin in lins:
        print "lineament %d (%d segments, %g km)" % (lin_num, len(lin.vertices), lin.length()*sat_radius_km)
        lin.calc_fits(nsr, nb)
        lin_num += 1

        # We generate the histogram dataset, adding one "count" for each km of
        # lineament that fit best at that value of b(ackrotation)
        linhist += int(lin.length()*sat_radius_km)*[degrees(lin.best_b()),]

    drawnsrhist(histdata=linhist, radius=sat_radius_km)

    # Output the shapefile with its new attributes
    # TODO

    return (lins, linhist, sat_radius_km)

    axis([0,180,0,90])
    plot(degrees(backrots), degrees(lindict['fits']))
    grid(True)
    show()
# }}}

def draw_nsr_fit_hist(linhist, bins=18): #{{{
     hist(linhist, bins=bins)
     axis(xmin=0, xmax=180)
     xlabel("degrees of backrotation")
     ylabel("km of lineaments")
     title("Lineament formation vs. best fit backrotation")
     show()
 # }}}
