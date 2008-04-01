"""
Calculate stresses on an icy satellite over a rectangular geographic region on
a regularly spaced lat-lon grid.

The datacube containing the results of the calculation are output as a
U{Unidata netCDF <http://www.unidata.ucar.edu/software/netcdf>} (.nc) file,
which can be displayed using a wide variety of visualization software.

"""

import SatStress as SS
import time
import netCDF3
import physcon as pc
import numpy
from numpy import linalg
import scipy
from optparse import OptionParser


def main():
    """
    Calculate satellite stresses on a regular grid within a lat-lon window.

    """

    usage = "usage: %prog [options] satfile gridfile outfile" 
    description = __doc__
    op = OptionParser(usage)

    (options, args) = op.parse_args()

    # Do a little error checking on the arguments:
    if len(args) != 3:
        op.error("incorrect number of arguments")
 
    # The satfile defines the satellite and forcings it is subject to:
    satfile = open(args[0], 'r')
    the_sat = SS.Satellite(satfile)
    satfile.close()

    # The gridfile defines the scope and resolution of the calculation:
    gridfile = open(args[1], 'r')
    the_grid = Grid(gridfile, satellite=the_sat)
    gridfile.close()

    # Initialize a StressCalc object for the diurnal and NSR stresses
    diurnal_stress = SS.StressCalc([SS.Diurnal(the_sat),])
    nsr_stress = SS.StressCalc([SS.NSR(the_sat),])

    # Make a netCDF object, and set its metadata
    # Create a netCDF file object to stick the calculation results in:
    nc_out = netCDF3.Dataset(args[2], 'w')

    # Set metadata fields of nc_out appropriate to the calculation at hand.

    # - Make sure that the file complies with either COARDS or CF.
    # - Include all of the parameters that go into doing the calculation (whatever
    #   is listed within the .satellite and .grid files).
    # - Make sure that the units for those parameters are clear.
    # - A decent description of the contents of the file in English.

    # Set the global data attributes of the netCDF file describing the run
    # This allows us to re-construct the entire model run if we want to, given
    # the SatStress library and a netCDF output file.
    nc_out.description = "Testing pySatStress on a regular grid"
    nc_out.history     = "Created: %s using pySatStress" % ( time.ctime(time.time()) )
    nc_out.Conventions = "COARDS"

    # Specify the size and shape of the output datacube  
    
    n_deltas = the_grid.nsr_delta_numsteps
    nc_out.createDimension('nsr_delta', n_deltas)
    nsr_deltas = nc_out.createVariable('nsr_delta', 'f4', ('nsr_delta',))
    nsr_deltas.units = ""
    nsr_deltas.long_name = "NSR Surface Delta (mu/(eta*omega))"
    nsr_deltas[:] = numpy.logspace(the_grid.nsr_delta_min, the_grid.nsr_delta_max, num=n_deltas)

    n_times = len(numpy.arange(the_grid.time_min, the_grid.time_max, the_grid.time_step))
    nc_out.createDimension('time', n_times)
    times = nc_out.createVariable('time', 'f4', ('time',))

    # Check to see what kind of units we're using for time,
    # and name the variables and their units appropriately
    if the_grid.orbit_min is None:
        times.units = "seconds"
        times.long_name = "time after periapse"
        times[:] = numpy.arange(the_grid.time_min, the_grid.time_max, the_grid.time_step)
    else:
        times.units = "degrees"
        times.long_name = "degrees after periapse"
        times[:] = numpy.arange(the_grid.orbit_min, the_grid.orbit_max, the_grid.orbit_step)

    n_lats = len(numpy.arange(the_grid.lat_min, the_grid.lat_max+the_grid.latlon_step, the_grid.latlon_step))
    nc_out.createDimension('latitude',  n_lats)
    lats = nc_out.createVariable('latitude',  'f4', ('latitude',))
    lats.units = "degrees_north"
    lats.long_name = "latitude"
    lats[:] = numpy.arange(the_grid.lat_min, the_grid.lat_max  + (the_grid.latlon_step/2.0), the_grid.latlon_step)

    n_lons = len(numpy.arange(the_grid.lon_min, the_grid.lon_max+the_grid.latlon_step, the_grid.latlon_step))
    nc_out.createDimension('longitude', n_lons)
    lons = nc_out.createVariable('longitude',  'f4', ('longitude',))
    lons.units = "degrees_east"
    lons.long_name = "longitude"
    lons[:] = numpy.arange(the_grid.lon_min, the_grid.lon_max  + (the_grid.latlon_step/2.0), the_grid.latlon_step)
    #  end dimension and coordinate variable definitions

    # Define the data variables 
    Ttt_D = nc_out.createVariable('Ttt_D', 'f4', ('time', 'latitude', 'longitude',))
    Ttt_D.units = "Pa"
    Ttt_D.long_name = "north-south component stress of Diurnal stresses"

    Tpt_D = nc_out.createVariable('Tpt_D', 'f4', ('time', 'latitude', 'longitude',))
    Tpt_D.units = "Pa"
    Tpt_D.long_name = "shear component of Diurnal stresses"

    Tpp_D = nc_out.createVariable('Tpp_D', 'f4', ('time', 'latitude', 'longitude',))
    Tpp_D.units = "Pa"
    Tpp_D.long_name = "east-west component of Diurnal stresses"

    Ttt_N = nc_out.createVariable('Ttt_N', 'f4', ('nsr_delta', 'latitude', 'longitude',))
    Ttt_N.units = "Pa"
    Ttt_N.long_name = "north-south component of NSR stresses"

    Tpt_N = nc_out.createVariable('Tpt_N', 'f4', ('nsr_delta', 'latitude', 'longitude',))
    Tpt_N.units = "Pa"
    Tpt_N.long_name = "shear component of NSR stresses"

    Tpp_N = nc_out.createVariable('Tpp_N', 'f4', ('nsr_delta', 'latitude', 'longitude',))
    Tpp_N.units = "Pa"
    Tpp_N.long_name = "east-west component of NSR stresses"
    #  end data variable definitions

    # Loop over the time variable, doing diurnal calculations over an orbit  
    for t in range(n_times):
        # We need some kind of progress update, and we need to make sure that
        # we have a representation of the time coordinate in seconds, because
        # that's what the SatStress library expects - even if we're ultimately
        # communicating time to the user in terms of "degrees after periapse"
        if the_grid.orbit_min is None:
            time_sec = times[t]
        else:
            time_sec = the_sat.orbit_period()*(times[t]/360.0)

        print "Calculating diurnal stresses at", times[t], times.long_name
        for lon in range(n_lons):
            for lat in range(n_lats):

                Tau_D = diurnal_stress.tensor(theta = scipy.radians(90.0-lats[lat]),\
                                                phi = scipy.radians(lons[lon]),\
                                                  t = time_sec )

                nc_out.variables['Ttt_D'][t,lat,lon] = Tau_D[0,0]
                nc_out.variables['Tpt_D'][t,lat,lon] = Tau_D[1,0]
                nc_out.variables['Tpp_D'][t,lat,lon] = Tau_D[1,1] 
                # Not doing the principal components in this app... but save for later.
                #(magA, magB), eigenvecs = linalg.eigh(Tau_D)
                #mag1_D[t,lat,lon], mag2_D[t,lat,lon] = max(magA, magB), min(magA, magB)

    # Make sure everything gets written out to the file.
    nc_out.sync()
    #  end Diurnal calculations

    # Begin the NSR calculations: 
    for delta in range(n_deltas):
        # Set the nsr_period such that it corresponds to the current value of Delta:
        the_sat.nsr_period = 4*numpy.pi*nsr_deltas[delta]*the_sat.layers[-1].viscosity/the_sat.layers[-1].lame_mu
        # Re-set the StressCalc object based on the new nsr_period:
        nsr_stress = SS.StressCalc([SS.NSR(the_sat),])
        # We need some kind of progress update...
        print "Calculating NSR stresses at Delta =", nsr_deltas[delta], "(P_nsr = %g yr)" % (the_sat.nsr_period/31556926)
        for lon in range(n_lons):
            for lat in range(n_lats):

                # Note that t is always zero, because in the NSR potential, all time
                # does is translate surface features through the field, and here we're
                # actually not interested in what the surface features are doing - we
                # just want to know what the stresses look like relative to a fixed 
                # longitude (relative to the planet's location in the satellite's sky)
                # There's more discussion of this issue in Wahr et al. (2008).
                Tau_N = nsr_stress.tensor(theta = scipy.radians(90.0-lats[lat]),\
                                            phi = scipy.radians(lons[lon]),\
                                              t = 0)

                nc_out.variables['Ttt_N'][delta,lat,lon] = Tau_N[0,0]
                nc_out.variables['Tpt_N'][delta,lat,lon] = Tau_N[1,0]
                nc_out.variables['Tpp_N'][delta,lat,lon] = Tau_N[1,1] 

                # Not doing the principal components in this app... but save for later.
                #(magA, magB), eigenvecs = linalg.eigh(Tau_N)
                #mag1_N[delta,lat,lon], mag2_N[delta,lat,lon] = max(magA, magB), min(magA, magB)

        # re-set the NSR period to correspond to the next value of Delta...

    #  end NSR calculation loop.
    
    # Make sure everything gets written out to the file.
    nc_out.sync()

# end main()

class Grid(): #
    """
    A container class defining the temporal and geographic range and resolution
    of the calculation.

    The parameters defining the calculation grid are read in from a name value
    file, parsed into a Python dictionary using L{SatStress.nvf2dict}, and used
    to set the data attributes of the L{Grid} object.

    The geographic extent of the calculation is specified by minimum and maximum
    values for latitude and longitude.

    The geographic resolution of the calculation is defined by an angular
    separation between calculations.  This angular separation is the same in
    the north-south and the east-west direction.

    The temporal range and resolution of the calculation can be specified either
    in terms of actual time units (seconds) or in terms of the satellite's orbital
    position (in degrees).  In both cases, time=0 is taken to occur at periapse.

    S{Delta} is a measure of how viscous or elastic the response of the body
    is.  It's equal to (S{mu})/(S{eta}S{omega}) where S{mu} and S{eta} are the
    shear modulus and viscosity of the surface layer, respectively, and
    S{omega} is the forcing frequency to which the body is subjected (see Wahr
    et al. (2008) for a detailed discussion).  It is a logarithmic parameter, so
    its bounds are specified as powers of 10, e.g. if the minimum value is -3,
    the initial S{Delta} is 10^-3 = 0.001.

    @ivar grid_id: A string identifying the grid
    @type grid_id: str
    @ivar lat_min: Southern bound, degrees (north positive).
    @type lat_min: float
    @ivar lat_max: Northern bound, degrees (north positive).
    @type lat_max: float
    @ivar lon_min: Western bound, degrees (east positive).
    @type lon_min: float
    @ivar lon_max: Eastern bound, degrees (east positive).
    @type lon_max: float
    @ivar latlon_step: Angular separation between calculations.
    @type latlon_step: float
    @ivar time_min: Initial time at which calculation begins (0 = periapse).
    @type time_min: float
    @ivar time_max: Final time at which calculation ends.
    @type time_max: float
    @ivar time_step: Seconds between subsequent calculations.
    @type time_step: float
    @ivar orbit_min: Initial orbital position in degrees (0 = periapse)
    @type orbit_min: float
    @ivar orbit_max: Final orbital position in degrees (0 = periapse)
    @type orbit_max: float
    @ivar orbit_step: Orbital angular separation between calculations in degrees
    @type orbit_step: float
    @ivar nsr_delta_min: Initial S{Delta} = 10^(nsr_delta_min)
    @type nsr_delta_min: float
    @ivar nsr_delta_max: Final S{Delta} = 10^(nsr_delta_max)
    @type nsr_delta_max: float
    @ivar nsr_delta_numsteps: How many S{Delta} values to calculate total
    @type nsr_delta_numsteps: int

    """

    def __init__(self, gridFile, satellite=None):
        """Initialize the Grid object from a gridFile."""

        # Slurp up the proffered parameters:
        gridParams = SS.nvf2dict(gridFile, comment='#')
        
        self.grid_id = gridParams['GRID_ID']

        self.lat_min     = float(gridParams['LAT_MIN'])
        self.lat_max     = float(gridParams['LAT_MAX'])
        self.lon_min     = float(gridParams['LON_MIN'])
        self.lon_max     = float(gridParams['LON_MAX'])
        self.latlon_step = float(gridParams['LATLON_STEP'])

        self.nsr_delta_numsteps = float(gridParams['NSR_DELTA_NUMSTEPS'])
        self.nsr_delta_min  = float(gridParams['NSR_DELTA_MIN'])
        self.nsr_delta_max  = float(gridParams['NSR_DELTA_MAX'])

        # We MUST have all three of either of these and NONE of the other.
        if gridParams.has_key('TIME_MIN'):
            self.time_min  = float(gridParams['TIME_MIN'])
            self.time_max  = float(gridParams['TIME_MAX'])
            self.time_step = float(gridParams['TIME_STEP'])

            self.orbit_min  = None
            self.orbit_max  = None
            self.orbit_step = None

        elif gridParams.has_key('ORBIT_MIN'):
            self.orbit_min  = float(gridParams['ORBIT_MIN'])
            self.orbit_max  = float(gridParams['ORBIT_MAX'])
            self.orbit_step = float(gridParams['ORBIT_STEP'])

            self.time_min  = satellite.orbit_period()*(self.orbit_min/360.0)
            self.time_max  = satellite.orbit_period()*(self.orbit_max/360.0)
            self.time_step = satellite.orbit_period()*(self.orbit_step/360.0)

        else:
            # Die with an error about bad input file...
            pass
    
#  end class Grid

if __name__ == "__main__":
    main()
