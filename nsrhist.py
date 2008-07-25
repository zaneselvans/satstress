#!/usr/bin/python
"""
Routines for comparing a set of lineaments to an NSR stress field.

"""

from satstress import *
from pylab import *

def calc_nsr_fits(satfile = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun",\
                  shapefile = "input/GlobalLineaments",\
                  nb = 180,\
                  numlins = "all"): #{{{
    
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

def draw_fits(lins): #{{{
    lin_num=0
    for lin in lins:
        # Plot the lineament itself:
        lons = [ v[0] for v in lin.vertices ]
        lats = [ v[1] for v in lin.vertices ]
        plot(mod(degrees(lons),180), abs(degrees(lats)), 'k-', linewidth=2)
        plot(linspace(0,180,num=len(lin.fits)), degrees(lin.fits))
        plot(linspace(0,180,10), 10*[degrees(lin.best_fit()),], 'r--')
        plot(10*[degrees(lin.best_b()),], linspace(0,90,10), 'r--')
        #grid(True)
        axis([0,180,0,91])
        xticks( range(0,181,15) )
        yticks( range(0,91,15) )
        xlabel("degrees of backrotation")
        ylabel("average misfit (degrees)")
        title("lin=%d, length=%gkm, best_fit=%g, best_b=%g, iqr=%g" % (lin_num, lin.sat_radius_km*lin.length(), degrees(lin.best_fit()), degrees(lin.best_b()), degrees(lin.fits_width()) ))
        show()
        x = raw_input("press return for next fit: ")
        clf()
        lin_num += 1
# }}}
