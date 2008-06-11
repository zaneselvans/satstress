"""
Create map-projected vector and scalar plots of tidal stresses from netCDF
files output by L{gridcalc}.

Inputs:
=======
A NetCDF file containing one ore more set of scalar variables named Ttt_*,
Tpt_*, Tpp_*, representing the north-south (theta-theta), east-west (phi-phi)
and shear (phi-theta) elements of the surface membrane stress tensor.

A subset of the following stress field names to include in the plot:

  1. B{NSR} (corresponding to stresses produced by L{satstress.NSR}).

  2. B{Flat} (corresponding to stresses produced by L{satstress.Flat}).
  
  3. B{Diurnal} (corresponding to stresses produced by L{satstress.Diurnal}).

A specification of the orbital position, either as time elapsed since periapse,
or as degrees beyond periapse (must correspond to a time-slice which is present
within the NetCDF file).

A specification of the NSR period (must correspond to a value o NSR period
present within the NetCDF file).

At most one magnitude field to display, chosen from:

  * Magnitude of the most tensile principal component of the membrane stress

  * Magnitude of the most compressive principal component of the membrane stress

A flag indicating whether or not the principal components of the sum of the
two stress tensors should be displayed.

The resolution at which the stress vectors should be displayed (in degrees
separating one vector from the next)

Outputs:
========
A plot showing:

  * the chosen magnitude field as the background, in color, with an appropriate
    scale bar alongside the plot.

  * the principal components of the stress vector field that results from
    adding each of the tensor components together, with positive (tensile)
    stresses displayed as red arrows, and negative (compressive) stresses
    displayed as blue arrows, along with an appropriate reference vector for
    scale.
    
"""

import satstress as ss
import netCDF3
import numpy
import pylab
from mpl_toolkits import basemap
import sys
from numpy import linalg
from optparse import OptionParser

def main(argv=sys.argv):
    # Parse the command line:
    usage = "usage: %prog [options] infile.nc"
    description = __doc__
    op = OptionParser(usage)
    op.set_defaults( diurnal     = False, \
                     flat        = False, \
                     nsr         = False, \
                     scalar_plot = False, \
                     vector_plot = False)
    
    op.add_option("--outfile", help="""write plot to OUTFILE""", metavar="OUTFILE")
    op.add_option("--diurnal", action="store_true", help="""include viscoelastic diurnal stresses""")
    op.add_option("--nsr", action="store_true", help="""include viscoelastic NSR stresses""")
    op.add_option("--flat", action="store_true", help="""include NSR flattening stresses""")
    op.add_option("--time", help="""time slice (seconds after periapse) to plot.  Cannot be used with --orbit_position.""")
    op.add_option("--time_slice", type="int", help="""time slice (seconds after periapse) to plot.  Cannot be used with --orbit_position.""")
    op.add_option("--orbital_pos", help="""orbital position (degrees after periapse) to plot.  Cannot be used with --time.""")
    op.add_option("--orbital_pos_slice", type="int", help="""orbital position slice to plot.  Cannot be used with --time.""")
    op.add_option("--nsr_period", help="""the NSR period to use in the plot (must be present in the input file)""")
    op.add_option("--nsr_period_slice", type="int", help="""the NSR period slice to use in the plot (must be present in the input file)""")
    op.add_option("--nsr_deg", help="""degrees of accumulated NSR (for flattening calculations) to use in the plot.  Cannot be used with --nsr""")
    op.add_option("--nsr_deg_slice", type="int", help="""the accumulated NSR slice (for flattening calculations).  Cannot be used with --nsr""")
    op.add_option("--scalar_plot", help="""scalar field to use as background (if any), may be either \"tens\" or \"comp\".""", choices=("tens","comp"))
    op.add_option("--vector_plot", action="store_true", help="""display the principal components of the sum of the stress fields being plotted as vectors.""")
    op.add_option("--vector_res", help="""the number of degrees separating individual vectors in the plot.""")
    
    opts, args = op.parse_args(args=argv[1:])
    
    if opts.diurnal and opts.flat:
        raise IncompatibleStressError()
    
    # Open the netCDF file containing the stress calculation:
    nc_in = netCDF3.Dataset(args[0])
    
    # Figure out which time slice we're going to include for Diurnal (if any)
    # Figure out which nsr_period slice we're going to include for NSR (if any)
    # Figure out which nsr_deg slice we're going to include for Flat (if any)
    
    lats  = nc_in.variables['latitude'][:]
    nlats = len(lats)

    lons  = nc_in.variables['longitude'][:]
    nlons = len(lons)

    times = nc_in.variables['time'][:]
    ntimes = len(times)
        
    nsr_periods = nc_in.variables['nsr_period'][:]
    nnsr_periods = len(nsr_periods)

    #nsr_degs = nc_in.variables['nsr_deg'][:]
    #nnsr_degs = len(nsr_degs)

    # Initialize to empty... but the right sizes:
    Tau_Diurnal = numpy.zeros((2,2,ntimes,nlats,nlons))
    Tau_NSR = numpy.zeros((2,2,nnsr_periods,nlats,nlons))
    #Tau_Flat= numpy.zeros((2,2,nnsr_degs,nlats,nlons))

    Tau = numpy.zeros((nlons,nlats,2,2))

    evals   = numpy.zeros((2,nlons,nlats))
    evecs   = numpy.zeros((2,2,nlons,nlats))
    pc_vecs = numpy.zeros((2,2,nlons,nlats))

    tens = numpy.zeros((nlats,nlons))
    comp = numpy.zeros((nlats,nlons))

    if opts.diurnal:
        Ttt_Diurnal = numpy.array(nc_in.variables['Ttt_Diurnal'][:,:,:])
        Tpt_Diurnal = numpy.array(nc_in.variables['Tpt_Diurnal'][:,:,:])
        Tpp_Diurnal = numpy.array(nc_in.variables['Tpp_Diurnal'][:,:,:])
        Tau_Diurnal = numpy.array([ [Ttt_Diurnal, Tpt_Diurnal],
                                    [Tpt_Diurnal, Tpp_Diurnal] ])

    if opts.nsr:
        Ttt_NSR = numpy.array(nc_in.variables['Ttt_NSR'][:,:,:])
        Tpt_NSR = numpy.array(nc_in.variables['Tpt_NSR'][:,:,:])
        Tpp_NSR = numpy.array(nc_in.variables['Tpp_NSR'][:,:,:])
        Tau_NSR = numpy.array([ [Ttt_NSR, Tpt_NSR],
                                [Tpt_NSR, Tpp_NSR] ])

    if opts.flat:
        Ttt_Flat = numpy.array(nc_in.variables['Ttt_Flat'][:,:,:])
        Tpt_Flat = numpy.array(nc_in.variables['Tpt_Flat'][:,:,:])
        Tpp_Flat = numpy.array(nc_in.variables['Tpp_Flat'][:,:,:])
        Tau_Flat = numpy.array([ [Ttt_Flat, Tpt_Flat],
                                 [Tpt_Flat, Tpp_Flat] ])

    # Now that we've got the datacubes in array form, we need to isolate the
    # slices we're actually plotting...

    if opts.nsr and opts.flat:
        # If we've got NSR and Flattening, calculate the delta between them:
        Tau = Tau_NSR[:,:,opts.nsr_period_slice,:,:] - Tau_Flat[:,:,opts.nsr_deg_slice,:,:]
    else:
        # otherwise add the stress fields together:
        Tau = Tau_NSR[:,:,opts.nsr_period_slice,:,:]# + Tau_Diurnal[:,:,opts.time_slice,:,:] + Tau_Flat[:,:,opts.nsr_deg_slice,:,:]
    
    # Now Tau contains the tensor form of whichever snapshot stress field we
    # want to plot, and needs to be turned into the principal components.  This
    # can't be done in a vectorized way because of the way the linalg.eig
    # function works... so we need to iterate over each lat-lon point.  At the
    # same time we have the opportunity to re-arrange the indicies for easier
    # x-y plotting:

    for nlon in numpy.arange(nlons):
        for nlat in numpy.arange(nlats):
            evals[:,nlon,nlat], evecs[:,:,nlon,nlat] = linalg.eig(Tau[:,:,nlat,nlon])
            #pc_vecs[:,:,nlat,nlon] = evals[:,nlon,nlat]*evecs[:,:,nlon,nlat]
            tens[nlat,nlon] = max(evals[:,nlon,nlat])
            comp[nlat,nlon] = min(evals[:,nlon,nlat])

    nx, ny = 18, 12
    i, j = 1, 1
    print pylab.shape(Tau)
    print pylab.shape(Tau[i,j])
    print Tau[i,j,nx,ny]
#    print evals[:,nx,ny]
#    print evecs[:,:,nx,ny]
#    print pc_vecs[:,:,nx,ny]

    pylab.clf()

    map = basemap.Basemap(llcrnrlon = lons[0],\
                          llcrnrlat = lats[0],\
                          urcrnrlon = lons[-1],\
                          urcrnrlat = lats[-1])

    map.drawmeridians(numpy.linspace(lons[0], lons[-1], 13), labels=[1,1,1,1])
    map.drawparallels(numpy.linspace(lats[0], lats[-1], 7 ), labels=[1,1,1,1])
    map.drawmapboundary()

    if opts.scalar_plot=='tens':
        map.imshow(tens)
    elif opts.scalar_plot=='comp':
        map.imshow(comp)

    if opts.vector_plot:
        # I do the plotting twice, once with the actual values, and ones with
        # their opposites, to get two sets of arrows for each vector (two going
        # away from each other for tension, or two coming towards each other
        # for compression).

        N = 3
        print pylab.shape(evecs[0,0,::N,::N])
        print pylab.shape(evecs[0,1,::N,::N])
        print pylab.shape(evecs[1,0,::N,::N])
        print pylab.shape(evecs[1,1,::N,::N])
        print pylab.shape(lats[::N])
        print pylab.shape(lons[::N])
        map.quiver(lons[::N], lats[::N],  evecs[0,0,::N,::N], -evecs[0,1,::N,::N])#, scale=200000000, width=0.001)
        map.quiver(lons[::N], lats[::N],  evecs[1,0,::N,::N], -evecs[1,1,::N,::N])#, scale=200000000, width=0.001)
        map.quiver(lons[::N], lats[::N], -evecs[0,0,::N,::N],  evecs[0,1,::N,::N])#, scale=200000000, width=0.001)
        map.quiver(lons[::N], lats[::N], -evecs[1,0,::N,::N],  evecs[1,1,::N,::N])#, scale=200000000, width=0.001)

if __name__ == "__main__":
    main(argv)

        
