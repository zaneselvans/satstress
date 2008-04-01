"""Check to see that SatStress gives the expected output from a series
of known calculations.

Calculates the stresses due to the L{NSR} and L{Diurnal} forcings at a
series of lat lon points on Europa, over the course of most of an orbit, and
also at a variety of different amounts of viscous relaxation.  Compares the
calculated values to those listed in the C{pickle} file passed in via the
command line.

C{test} is called from the C{SatStress Makefile}, when one does C{make test},
with the appropriate C{pickle}d input to compare against (it is provided with
the source code).

This function also acts as a short demonstration of how to use the SatStress
module.
"""

import pickle
import sys
from optparse import OptionParser
import SatStress as SS

usage = "usage: %prog satfile [testoutput.pkl]"
op = OptionParser(usage)
(options, args) = op.parse_args(sys.argv[1:])

if len(args) < 1 or len(args) > 2:
    op.error("incorrect number of arguments")

# Create a new satellite object, as defined by the input file:
the_sat = SS.Satellite(open(args[0],'r'))

# Create a StressCalc object, that calculates both the NSR and Diurnal
# stresses on the Satellite just instantiated:
the_stresses = SS.StressCalc([SS.Diurnal(the_sat), SS.NSR(the_sat)])

# do a series of calculations, varying all the inputs at each iteration
theta = 0.0
phi   = 0.0
t     = 0.0

Taulist = []

while the_stresses.stresses[1].Delta() > 0.01 and t < 2.0*scipy.pi/the_sat.mean_motion():
    print("""Tau(theta = %g, phi = %g, time = %g, Delta = %g) = """ % (scipy.degrees(theta), scipy.degrees(phi), t, the_stresses.stresses[1].Delta()))
    Taulist.append(the_stresses.tensor(theta=theta, phi=phi, t=0))
    print Taulist[-1], "\n"
    theta = theta + scipy.pi/23.0
    phi   = phi   + scipy.pi/11.0
    t     = t     + (scipy.pi/11.0)/the_sat.mean_motion()

    # Re-construct the NSR stress (and the_stresses StressCalc object)
    # changing the forcing period, so that we can have the NSR stresses
    # vary with each iteration too.
    the_sat.nsr_period = the_sat.nsr_period / 1.5
    the_stresses = SS.StressCalc([SS.Diurnal(the_sat), SS.NSR(the_sat)])

# Now we want to compare this against reference output, just to make sure
# that we're getting the right numbers...

pickledTau = pickle.load(open(args[1]))
for (tau1, tau2) in zip(Taulist, pickledTau):
    if tau1.all() != tau2.all():
        print("\nTest failed.  :(\n")
        sys.exit(1)

print("\nTest passed! :)\n")
sys.exit(0)
