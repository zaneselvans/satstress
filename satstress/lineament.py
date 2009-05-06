from numpy import *
from pylab import *
from mpl_toolkits.basemap import Basemap
import pickle

# Open Source Geospatial libraries:
from osgeo import ogr

class Lineament(object): #{{{1
    """
    A one dimensional feature on the surface of a spherical satellite.

    """

    def __init__(self, lons=None, lats=None, stresscalc=None, bs=None, nsrdbars=None, nsrstresswts=None): #{{{2
        """
        Create a lineament from a given list of (lon,lat) points.

        lons and lats are arrays of longitudes and latitudes, having equal
        length, defining the vertices making up a polyline.  East and North are
        taken to be positive.  Units are radians
        
        stresscalc is a L{satstress.StressCalc} object defining the field we
        are going to be comparing the lineament to.
        
        If bs, nsrdbars, and nsrstresswts are set and have the same length, then we have
        a complete NSR fit and we can go ahead and set those values for the
        feature.  Otherwise, set these attributes to None.

        """

        # make sure we've got good input points...
        assert (len(lons) == len(lats))
        assert (len(lons) > 0)

        self.lons = array(lons)
        self.lats = array(lats)
        self.length = self.calc_length()
        self.stresscalc = stresscalc

        if bs is None or nsrdbars is None or nsrstresswts is None:
            self.bs = None
            self.nsrdbars = None
            self.nsrstresswts = None
        else:
            assert(len(bs) == len(nsrdbars) == len(nsrstresswts))
            self.bs = array(bs)
            self.nsrdbars = array(nsrdbars)
            self.nsrstresswts = array(nsrstresswts)
    #}}}2

    def __str__(): #{{{2
        """
        """
        pass
    #}}}2

    def calc_length(self): #{{{2
        """
        Return the total length of the L{Lineament} in radians of arc, along
        the surface of the sphere.

        """
        if len(self.lons) > 1:
            return(sum(self.seg_lengths()))
        else:
            return(0.0)

    #}}}2 end calc_length
    
    def sinuosity(self): #{{{2
        """
        Return the sinuosity of the lineament, using spherical distances.

        Sinuosity is the ratio of the sum of the lengths of all the line
        segments making up the lineament, to the distance separating the
        endpoints.

        """
        # If the two endpoints are the same, we'll get a divide by zero error...
        assert( (self.lons[0],self.lats[0]) != (self.lons[-1], self.lats[-1]))
        return (self.length / spherical_distance(self.lons[0], self.lats[0], self.lons[-1], self.lats[-1]))

    #}}}2

    def midpoint(self): #{{{2
        """
        Returns the (lon,lat) point halfway along the lineament.

        """

        seg_lengths = self.seg_lengths()
        half_length = self.length/2.0

        cumulative_lengths = concatenate([ [0.0,], array([ seg_lengths[:N+1].sum() for N in range(len(seg_lengths)) ])])
        just_past_half = where(cumulative_lengths > half_length)[0][0]
        passed_by = cumulative_lengths[just_past_half] - half_length
        back_az = spherical_azimuth(self.lons[just_past_half], self.lats[just_past_half], self.lons[just_past_half-1], self.lats[just_past_half-1])
        mp_lon, mp_lat = spherical_reckon(self.lons[just_past_half], self.lats[just_past_half], back_az, passed_by)

        return(mp_lon, mp_lat)

    #}}}2

    def seg_midpoints(self): #{{{2
        """
        Return a list of the midpoints of each line segment in the Lineament.

        """
        if len(self.lons) == 1:
            mp_lons, mp_lats = self.lons, self.lats
        else:
            mp_lons,mp_lats = spherical_midpoint(self.lons[:-1], self.lats[:-1], self.lons[1:], self.lats[1:])

        return(mp_lons, mp_lats)
    #}}}2

    def seg_azimuths(self): #{{{2
        """
        Calculate the lineament azimuths (orientations) at the midpoints of the
        segments making up the lineament.  Because these segments are not
        directional, this will be an angle between 0 and pi.

        """

        mplons, mplats = self.seg_midpoints()
        return(mod(spherical_azimuth(mplons, mplats, self.lons[:-1], self.lats[:-1]), pi))
    # }}}2 end midpoint_azimuths

    def seg_lengths(self): #{{{2
        """
        Calculate the lengths (in radians of arc) of the line segments making
        up the lineament.

        """

        return(spherical_distance(self.lons[:-1], self.lats[:-1], self.lons[1:], self.lats[1:]))

    #}}}2 end seg_lengths

    def lonshift(self, b): #{{{2
        """
        Return the lineament shifted in longitude by 'b' radians on the
        surface, with its fits as they would have been if they'd been
        calculated from that location.

        """
        # Add b to each of self.lons
        nb = len(self.bs)
        shift_lons = self.lons+b
        shift_bs = mod(self.bs-b,pi)
        dtype = [('bs',float),('nsrdbars',float),('nsrstresswts',float)]
        shift_fits = zeros(nb, dtype=dtype)
        shift_fits['bs'] = shift_bs
        shift_fits['nsrdbars'] = self.nsrdbars
        shift_fits['nsrstresswts'] = self.nsrstresswts
        shift_fits.sort(order='bs')
        return(Lineament(lons=shift_lons, lats=self.lats, stresscalc=self.stresscalc,\
                         bs=shift_fits['bs'], nsrdbars=shift_fits['nsrdbars'],\
                         nsrstresswts=shift_fits['nsrstresswts']))
    #}}}2

    def poleshift(self, pnp_lon=0.0, pnp_lat=pi/2): #{{{2
        """
        Return a Lineament object representing the location and orientation of
        self, when the point defined in modern lon,lat terms by pnp_lon pnp_lat
        was the north pole.

        """
        tpw_lons, tpw_lats = paleopole_transform(lon_in=self.lons, lat_in=self.lats, pnp_lon=pnp_lon, pnp_lat=pnp_lat)
        return(Lineament(lons=tpw_lons, lats=tpw_lats, stresscalc=self.stresscalc))
    #}}}2

    def d_min(self, linB): #{{{2
        """
        Return an array containing the distances from the midpoints of each
        segment in self to the nearest midpoint of a segment in linB.

        """
        
        return(d_min(self, linB))
    #}}}2

    def mhd(self, linB): #{{{2
        """
        MHD form A to B.

        """
        
        return(mhd(self, linB))
    #}}}2

    def doppelgen_midpoint_nsr(self, stresscalc=None): #{{{2
        """
        Find a midpoint, representative of the lineament, and generate a synthetic
        NSR feature there which approximates the feature, for each value of
        backrotation (b, longitudinal shift) stored in self.bs.  Return those
        doppelgangers as a list.

        """
        if stresscalc is None:
            stresscalc=self.stresscalc

        # Figure out where to initiate the formation of the doppelganger:
        ep1_lon, ep1_lat, ep2_lon, ep2_lat = self.bfgcseg_endpoints()
        mp_lon, mp_lat = spherical_midpoint(ep1_lon, ep1_lat, ep2_lon, ep2_lat)

        # Set the max_length of the doppelganger to be the length of the best
        # fit great circle, and not the prototype, since NSR features are almost
        # perfectly straight.
        max_length = spherical_distance(ep1_lon, ep1_lat, ep2_lon, ep2_lat)

        # Create an array of all the midpoints at which we need to generate
        # doppelgangers, given the our b values.
        init_lons = self.bs + mp_lon

        seg_len = min(0.01, self.length/10.0)

        return([ lingen_nsr(stresscalc=stresscalc, init_lon=init_lon, init_lat=mp_lat,\
                            max_length=self.length, prop_dir="both", seg_len=seg_len) for init_lon in init_lons ])

    #}}}2 end doppelgen_midpoint_nsr

    def doppelgen_gcseg(self): #{{{2
        """
        Generate a great circle segment approximating the lineament.

        """
        init_lon, init_lat, fin_lon, fin_lat = self.bfgcseg_endpoints()
        return(lingen_greatcircle(init_lon, init_lat, fin_lon, fin_lat))
    #}}}2

    def bfgcseg_endpoints(self): #{{{2
        """
        Find the points on the lineament's best fit great circle which are
        closest to the endpoints of the lineament.

        """

        mplons, mplats = self.seg_midpoints()
        bfgcpole_lon, bfgcpole_lat = self.bfgc_pole()

        # find points on the great circle nearest to the lineament endpoints.
        # Go mod((pi/2)-(distance to gcpole),pi/2) radians away from the pole if closer than pi/2
        # and toward the pole if further away than pi/2.
        dist_to_bfgcpole1 = spherical_distance(self.lons[0], self.lats[0], bfgcpole_lon, bfgcpole_lat)
        az_to_bfgcpole1   = spherical_azimuth(self.lons[0], self.lats[0], bfgcpole_lon, bfgcpole_lat)

        # if further than pi/2 from the gc_pole, go mod(dist,pi/2) toward it.
        if(dist_to_bfgcpole1 > pi/2):
            bfgcseg_lon1, bfgcseg_lat1 = spherical_reckon(self.lons[0], self.lats[0], az_to_bfgcpole1, mod(dist_to_bfgcpole1, pi/2))
        # otherwise, go pi/2-dist away from it.
        else:
            bfgcseg_lon1, bfgcseg_lat1 = spherical_reckon(self.lons[0], self.lats[0], az_to_bfgcpole1+pi, pi/2-dist_to_bfgcpole1)

        dist_to_bfgcpole2 = spherical_distance(self.lons[-1], self.lats[-1], bfgcpole_lon, bfgcpole_lat)
        az_to_bfgcpole2   = spherical_azimuth(self.lons[-1], self.lats[-1], bfgcpole_lon, bfgcpole_lat)

        # if further than pi/2 from the gc_pole, go mod(dist,pi/2) toward it.
        if(dist_to_bfgcpole2 > pi/2):
            bfgcseg_lon2, bfgcseg_lat2 = spherical_reckon(self.lons[-1], self.lats[-1], az_to_bfgcpole2, mod(dist_to_bfgcpole2, pi/2))
        # otherwise, go pi/2-dist away from it.
        else:
            bfgcseg_lon2, bfgcseg_lat2 = spherical_reckon(self.lons[-1], self.lats[-1], az_to_bfgcpole2+pi, pi/2-dist_to_bfgcpole2)

        return(bfgcseg_lon1, bfgcseg_lat1, bfgcseg_lon2, bfgcseg_lat2)
    #}}}2

    def bfgcseg_midpoint(self): #{{{2
        """
        Find the lon,lat point which lies on the lineament's best fit great
        circle, and is halfway between points on the great circle which are
        closest to the endpoints of the lineament.  This point represents a
        kind of geometric center of the feature, which can be used to
        approximate its location for a variety of purposes, including
        generating doppelgangers.

        """

        bfgcseg_lon1, bfgcseg_lat1, bfgcseg_lon2, bfgcseg_lat2 = self.bfgcseg_endpoints()
        bfgc_mp_lon, bfgc_mp_lat = spherical_midpoint(bfgcseg_lon1, bfgcseg_lat1, bfgcseg_lon2, bfgcseg_lat2)
        return(bfgc_mp_lon, bfgc_mp_lat)

    #}}}2

    def calc_nsrfits(self, nb=180, stresscalc=None): #{{{2
        """
        For nb evenly spaced values of longitudinal translation, b, ranging
        from 0 to pi, calculate the fit metric (dbar) for the lineament,
        using the stresscalc object associated with the feature to generate
        the hypothesized synthetic features.  Use the stresscalc that was
        passed in if there is none associated with the feature, or 

        Store the resulting values of b and dbar in the lineament's
        attributes bs, and nsrdbars for later reference.

        This only works for a purely NSR stress field.

        Calculate a metric of the significance of the feature's location within
        the stress field, based on its magnitude and anisotropy relative to the
        global average, and store it for later use in weighting the feature's
        contribution to the overall results.  This is self.nsrstresswts

        Generate a synthetic lineament meant to replicate the feature (a
        doppelganger) by using the NSR stress field, and a starting point near
        the center of the feature (at the midpoint of a great circle segment
        approximating the feature).  Calculate the mean Hausdorff distance
        (MHD) between the prototype (self) and this doppelganger for every
        value of b, and store the results in self.nsrdbars.

        """

        # we have to have at least one stresscalc or this is pointless:
        assert(self.stresscalc is not None or stresscalc is not None)

        if stresscalc is None and self.stresscalc is not None:
            stresscalc = self.stresscalc
        elif self.stresscalc is None and stresscalc is not None:
            self.stresscalc = stresscalc

        # set the b values first, so that the fit metrics can refer to them.
        self.bs = linspace(0,pi,nb,endpoint=False)

        # Create vectors of non-stress and non-b dependent values:
        mp_lons, mp_lats = self.seg_midpoints()
        w_length = self.seg_lengths()/self.length

        # number of segments in this feature
        nsegs = len(w_length)

        # Create vectors containing all of the lons and lats at which stresses
        # will need to be calculated, based on the locations of the vertices
        # making up the lineament and the values of b at which we are doing
        # calculations
        calc_thetas = tile((pi/2)-mp_lats, nb)
        calc_phis = repeat(linspace(0,pi,nb,endpoint=False),nsegs) + tile(mp_lons,nb)

        # use SatStress to perform the stress calculations at those locations
        tens_mag, tens_az, comp_mag, comp_az = stresscalc.principal_components(calc_thetas, calc_phis, 0.0)

        # Create an (nsegs*nb) length array of w_stress values
        w_stress = (tens_mag - comp_mag)/stresscalc.mean_global_stressdiff()

        # generate NSR doppelgangers for all b in self.bs:
        doppels = self.doppelgen_midpoint_nsr(stresscalc)

        # Shift all of the doppelgangers such that they are superimposed upon
        # the prototype lineament:
        for doppel,b in zip(doppels,self.bs):
            doppel.lons -= b

        # for each point in self (the prototype) find the minimum distance to
        # any point in each doppelganger - this measures the similarity in
        # their shape.  It's normalized by self.length so as to be unitless,
        # and scaled linearly to the overall length of the feature.
        d_min = ravel(array([ self.d_min(doppel) for doppel in doppels ]))/self.length

        # note that this is the RMS minimum distance...
        self.nsrdbars = sqrt(((tile(w_length,nb)*(d_min**2))).reshape(nb,nsegs).sum(axis=1))

        self.nsrstresswts = ((tile(w_length,nb)*w_stress)).reshape(nb,nsegs).sum(axis=1)
    #}}}2 end calc_nsrfits

    def nsrfits(self, dbar_max=0.125): #{{{2
        """
        Return a heuristic metric of the probability that the feature fits the
        chosen stress field at each value of b, multiplied by the significance
        of such a fit, as measured by the lack of perturbability of the fit in
        this location within the stress field.

        """

        assert(self.bs is not None and self.nsrdbars is not None and self.nsrstresswts is not None)

        return(self.nsrstresswts*(1.0 - where(self.nsrdbars/dbar_max < 1.0, self.nsrdbars/dbar_max, 1.0))**2)

    #}}}2

    def best_fit(self, dbar_max=0.125): #{{{2
        """
        Return best fit, and the value of b at which it occurs.

        """

        nsrfits = self.nsrfits()
        best_fit = max(nsrfits)
        best_b = float(self.bs[where(fabs(nsrfits-best_fit) < 1e-9)[0][0]])

        return(best_fit, best_b)
    #}}}2

    def bfgc_pole(self): #{{{2
        """
        Return the (lon, lat) point on the sphere defining the great circle
        which minimizes the sum of the squares of the distance between it and
        the lineament.

        """

        weights = self.seg_lengths()/self.length
        mp_lons, mp_lats = self.seg_midpoints()

        # special case for features consisting of a single segment... in which
        # case all the eigenstuff is both broken and unnecessary:
        if len(mp_lons) == 1:
            best_pole_phi, bp_lat = spherical_reckon(mp_lons[0], mp_lats[0], spherical_azimuth(self.lons[0], self.lats[0], self.lons[-1], self.lats[-1])+pi/2.0, pi/2.0)
            best_pole_theta = pi/2.0 - bp_lat
        else:
            mp_x, mp_y, mp_z = sphere2xyz(weights, pi/2 - mp_lats, mp_lons)
            A = dot(array([mp_x, mp_y, mp_z]),array([mp_x, mp_y, mp_z]).transpose())
            eigenvals, eigenvecs = eig(A)
            best_pole_x, best_pole_y, best_pole_z = eigenvecs[:,where(fabs(eigenvals - min(eigenvals)) < 1e-9)[0][0] ]
            best_pole_r, best_pole_theta, best_pole_phi = xyz2sphere(best_pole_x, best_pole_y, best_pole_z)

        return(best_pole_phi, pi/2.0 - best_pole_theta)
    #}}}2 end bfgc_pole

    def plot(self, map, lon_cyc=2*pi, lat_mirror=False, color='black', alpha=1.0, linewidth=1.0): #{{{2
        """
        Plot the lineament on the provided Basemap object, map, with the
        associated color, alpha, and width.

        if lon_cyclic is non-zero, modulo the longitudes of the feature by the
        cyclic value, until any that can ever fall in the displayed area of the
        map do so, plotting the feature once for each distinct time it is the
        case.  The default value, 2*pi, means that features which cross the
        edge of the map will have that portion that runs out of the display
        plotted on the far side.  If lon_cyclic=pi, then all the features will
        be placed in the same non-repeating pi-wide portion of the stress
        field, for easy comparison therein.

        If lat_mirror is True, use the absolute values of the latitudes,
        instead of the latitudes, in order to map the lineament into the
        northern hemisphere.  This allows more compact comparison of lineaments
        in stress-equivalent locations.
        
        """

        plotted_lines = []

        # What does it mean to have a longitude discontinuity?  It means
        # that two adjacent points in the lineament have longitude values
        # which differ by more than they need to... that is, for two points
        # ptA, ptB, having longitudes lonA and lonB, there exists some integer
        # N, such that:
        # abs(lonA - (2*pi*N + lonB)) < abs(lonA - lonB) 
        # this can only be true if:
        # abs(lonA - lonB) > pi
        lons = fixlons(self.lons)

        # In order to create the illusion that the lineaments are wrapping
        # around in longitude, we plot one copy at each longitude where a
        # portion of them will show up on the map...
        nrange = unique1d(floor((lons-radians(map.llcrnrlon))/lon_cyc))
        for N in nrange:
            (x,y) = map(degrees(lons-N*lon_cyc), degrees(self.lats))
            plotted_lines.append(map.plot(x, y, color=color, alpha=alpha, linewidth=linewidth))

            # If we're mirroring across the equator, we need to plot all of the
            # lineaments with their anti-latitudes as well:
            if lat_mirror:
                (x,y) = map(degrees(lons-N*lon_cyc), -degrees(self.lats))
                plotted_lines.append(map.plot(x, y, color=color, alpha=alpha, linewidth=linewidth))

        return(plotted_lines)
    #}}}2

#}}}1 end of the Lineament class

################################################################################
# Helpers having to do with fit metrics or lineament generation.
################################################################################
def lingen_nsr(stresscalc, init_lon=None, init_lat=None, max_length=2*pi, prop_dir="both", seg_len=0.01): # {{{
    """
    Generate a synthetic NSR feature, given a starting location, maximum
    length, propagation direction, and a step size on the surface.

    Assumes tensile fracture, perpendictular to the most tensile principal
    component of the stresses.

    """


    # Make a note of the fact that we've got to do the other half too
    if prop_dir == "both":
        max_length = max_length/2.0
        prop_dir = "east"
        done = False
    else:
        done = True

    lons = array([init_lon,])
    lats = array([init_lat,])
    # Calculate the stresses at the given time and initial location
    (tens_mag, tens_az, comp_mag, comp_az) = stresscalc.principal_components(theta=(pi/2.0)-lats[0], phi=lons[0], t=0.0)

    lin_length = 0.0
    ice_strength = stresscalc.stresses[0].satellite.layers[-1].tensile_str

    # This can't be vectorized because we don't know where we'll end up until
    # we get there.
    while lin_length < max_length and tens_mag > ice_strength:
        # Choose which direction to go based on prop_dir, knowing that comp_az
        # is always an angle clockwise from north, between 0 and pi
        if prop_dir == "east":
            prop_az = comp_az
        else:
            assert(prop_dir == "west")
            prop_az = comp_az + pi

        # this introduces irregularity into the length of the segments, to
        # avoid having problems with exact integer multiples in doppelganger
        # generation, etc.
        newseglen = 2.0*seg_len*random_sample()

        newlon, newlat = spherical_reckon(lons[-1], lats[-1], prop_az, newseglen)

        # Make sure that our new longitude is within 2*pi in longitude of the
        # previous point in the feature, i.e. don't allow any big
        # discontinuities (this is a hack to deal with longitude cyclicity)
        while (abs(newlon - lons[-1]) > abs((newlon - 2*pi) - lons[-1])):
            newlon = newlon - 2*pi
        while (abs(newlon - lons[-1]) > abs((newlon + 2*pi) - lons[-1])):
            newlon = newlon + 2*pi
        
        lons = concatenate([lons,array([newlon,])])
        lats = concatenate([lats,array([newlat,])])

        lin_length += newseglen

        # Calculate the stresses at the new location
        (tens_mag, tens_az, comp_mag, comp_az) = stresscalc.principal_components(theta=(pi/2.0)-lats[-1], phi=lons[-1], t=0.0)

    # if we only got a single point, then we failed to initiate a fracture, and
    # should not even try to make the second half.
    if not done and len(lons) > 1:
        second_half = lingen_nsr(stresscalc, init_lon=init_lon, init_lat=init_lat, max_length=max_length, prop_dir="west", seg_len=seg_len)
        # the second_half lat/lon arrays are in reversed order here to preserve
        # the overall directionality of the vertices (see Numpy fancy slicing):
        lons = concatenate([second_half.lons[:0:-1], lons])
        lats = concatenate([second_half.lats[:0:-1], lats])

    return(Lineament(lons=lons, lats=lats, stresscalc=stresscalc))

#}}} end lingen_nsr

def lingen_greatcircle(init_lon, init_lat, fin_lon, fin_lat, seg_len=0.01): #{{{
    """
    Return a L{Lineament} object closely approximating the shortest great
    circle route between (lon1,lat1) and (lon2,lat2).

    """

    gc_length  = spherical_distance(init_lon, init_lat, fin_lon, fin_lat)

    # this forces each feature to have at least 10 segments
    num_segs = max(10, gc_length/seg_len)

    # randomizing the lengths of the segments in order to prevent creation of
    # artifacts in the fit curves.
    shape(gc_length)
    shape(gc_length*random(num_segs))
    seg_dists = concatenate([[0,],sort(gc_length*random(num_segs)),ravel(array([gc_length,]))])

    init_az    = spherical_azimuth(init_lon, init_lat, fin_lon, fin_lat)
    lons, lats = spherical_reckon(init_lon, init_lat, init_az, seg_dists)
    return(Lineament(lons=lons, lats=lats))

#}}} end lingen_greatcircle

def d_min(linA, linB): #{{{
    """
    For each point in linA, find the minimum distance to a point in linB.

    """
    mp_lonsA, mp_latsA = linA.seg_midpoints()
    mp_lonsB, mp_latsB = linB.seg_midpoints()

    # Create a matrix of lat/lon points such that pairwise calculations in
    # spherical_distance do the entire grid of possible distances between
    # points in linA and linB:
    distmatrix = spherical_distance(repeat(mp_lonsA,len(mp_lonsB)), repeat(mp_latsA,len(mp_latsB)),\
                                      tile(mp_lonsB,len(mp_lonsA)),   tile(mp_latsB,len(mp_latsA)))

    # reshape distmatrix such that each column corresponds to a point in linA and find the
    # minimum distances for each point in A:
    return(distmatrix.reshape((len(mp_lonsA), len(mp_lonsB))).min(axis=1))

#}}} end d_min

def mhd(linA, linB): #{{{
    """
    Calculate the MHD from linA to linB.  Duh.

    """

    return(sum(d_min(linA, linB)*(linA.seg_lengths()/linA.length))/linA.length)
#}}}

################################################################################
# Helpers having to do with spherical geometry
################################################################################
def spherical_distance(lon1, lat1, lon2, lat2): #{{{
    """
    Calculate the distance between two points on a sphere, in radians.

    """

    sph_len = 0.0
    sph_len += sin((lat1-lat2)/2)**2 
    sph_len += cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2)**2)
    sph_len =  2.0*arcsin(sqrt(sph_len))
    return(sph_len)
# }}} end spherical_distance

def spherical_azimuth(lon1, lat1, lon2, lat2): # {{{
    """
    Calculate the azimuth between one point and another.

    The returned azimuth angle is the initial heading (angle clockwise from
    north) to set out on along a great-circle route in order to travel from
    point 1 to point 2.

    """
    lon1 = 2.0*pi - lon1
    lon2 = 2.0*pi - lon2

    return(mod(arctan2(sin(lon1-lon2)*cos(lat2), cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon1-lon2)), 2.0*pi))

# }}} end spherical_azimuth

def spherical_reckon(lon1, lat1, az12, ang_dist): # {{{
    """
    Calculate current location given a starting point, heading, and distance
    traveled.

    """
    # This formula uses WEST positive longitude... so we need to adjust.
    lon1 = 2.0*pi - lon1

    lat2 = arcsin( sin(lat1)*cos(ang_dist) + cos(lat1)*sin(ang_dist)*cos(az12) )
    dlon = arctan2( sin(az12)*sin(ang_dist)*cos(lat1), cos(ang_dist) - sin(lat1)*sin(lat2) )
    lon2 = mod( lon1 - dlon + pi, 2.0*pi ) - pi

    # This formula used WEST positive longitude... so we need to re-adjust.
    lon2 = 2.0*pi - lon2
    return(array([lon2,lat2]))

# }}} end spherical_reckon

def spherical_midpoint(lon1, lat1, lon2, lat2): #{{{
    """
    Given two points on the surface of a sphere, return the point that lies
    halfway between them on a great circle route.

    """

    ang_dist = spherical_distance(lon1, lat1, lon2, lat2)
    az12     = spherical_azimuth(lon1, lat1, lon2, lat2)

    return(spherical_reckon(lon1, lat1, az12, ang_dist/2.0))

# }}}

def random_lonlatpoints(N): #{{{
    """
    Generate N evenly distributed random points on a sphere.

    """

    return(2*pi*random(N), (pi/2)-arccos(2*random(N)-1))
#}}}

def sphere2xyz(r, theta, phi): #{{{
    """
    Takes a point in spherical coordinates and returns its cartesian
    equivalent.

    r is the distance from the origin.

    theta is south positive co-latitude (zero at the north pole)

    phi is the east positive longitude.

    Assumes that the z-axis is the rotation axis.

    """

    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)

    return(x, y, z)
#}}}

def xyz2sphere(x, y, z): #{{{
    """
    Takes a Cartesian (x,y,z) point and returns the equivalent point in terms
    of spherical (r, theta, phi), where:

    r is the radial distance from the origin.

    theta is co-latitude.

    phi is east-positive longitude.

    And the z-axis is taken to be the same as the rotation axis.

    """

    r     = sqrt(x**2 + y**2 + z**2)
    theta = arctan2(sqrt(x**2 + y**2), z)
    phi   = arctan2(y, x)

    return(r, theta, phi)
#}}}

def paleopole_transform(pnp_lon, pnp_lat, lon_in, lat_in): #{{{
    """
    Transforms the location of a point on the surface of a sphere, defined by
    (lon_in,lat_in) in east-positive longitude, returning the (lon,lat)
    location that point would be located at if the point defined by
    (pnp_lon,pnp_lat) is moved directly north until it is at the north pole.

    """

    # First we convert the point to be transformed into Cartesian coords.
    # Since we only care about the lon,lat position in the end, we'll treat the
    # the body as a unit sphere.  Remember that sphere2xyz needs CO-latitude:
    colat_in = pi/2 - lat_in
    xyz_in = array(sphere2xyz(1.0, colat_in, lon_in))

    # Now, remember that what we're doing is bringing a wayward pole back to
    # the top of the sphere... which means the angles we're interested in are
    # actually -colat, -lon:
    alpha = pi/2 + pnp_lon 
    beta  = pi/2 - pnp_lat 
    gamma = 0

    # The definition of the X-Z rotation matrix:
    rot_mat = array([ [ cos(alpha)*cos(gamma)-sin(alpha)*cos(beta)*sin(gamma), -cos(alpha)*sin(gamma)-sin(alpha)*cos(beta)*cos(gamma),  sin(beta)*sin(alpha) ],\
                      [ sin(alpha)*cos(gamma)+cos(alpha)*cos(beta)*sin(gamma), -sin(alpha)*sin(gamma)+cos(alpha)*cos(beta)*cos(gamma), -sin(beta)*cos(alpha) ],\
                      [                 sin(beta)*sin(gamma),                                     sin(beta)*cos(gamma),                       cos(beta)      ] ])

    # Do the matrix multiplication using dot():
    xyz_out = dot(xyz_in.transpose(), rot_mat)

    # Transform back to spherical coordinates:
    r_out, theta_out, phi_out = xyz2sphere(xyz_out[:,0], xyz_out[:,1], xyz_out[:,2])

    # and back to lon/lat from colon/lat:
    lon_out = mod(phi_out + alpha,2*pi)
    lat_out = pi/2 - theta_out # co-lat to lat

    return(lon_out, lat_out)
#}}}

def fixlons(lons): #{{{
    """
    Takes a set of longitudes, and forces it to be within a continuous range
    (i.e. no skips at "international date line")

    """

    lons = mod(lons, 2*pi)
    while len(where(abs(lons[:-1]-lons[1:]) > pi)[0] > 0):
        lons[1:] = where(abs(lons[:-1]-lons[1:]) > pi, lons[1:]+2*pi*sign(lons[:-1]-lons[1:]), lons[1:])
    return(lons)
#}}}

################################################################################
# Helpers for input/output/saving/updating: (TODO: portable output...)
################################################################################
def update_lins(lins): # {{{
    """
    Just a quick and easy way to remember what all has to be done in order to
    update a set of lineaments to include all the newest and bestest
    functionality from the module.

    """

    newlins = []
    linlen=len(lins)
    for lin,N in zip(lins,range(linlen)):
        newlin = Lineament(lons=lin.lons, lats=lin.lats, stresscalc=lin.stresscalc, bs=lin.bs, nsrdbars=lin.nsrdbars, nsrstresswts=lin.nsrstresswts)
        newlins.append(newlin)

    return(newlins)
#}}} end update_lins

def plotlinmap(lins, map=None, color='black', alpha=1.0, linewidth=1.0, lon_cyc=2*pi, lat_mirror=False): #{{{
    """
    Plot a map of the lineaments listed in 'lins'.  Plot it to 'map' if given.

    """
    if map is None:
        thefig = figure()
        map = Basemap()
        map.drawmapboundary(fill_color="white")
        map.drawmeridians(range(-180,181,30), labels=[1,0,0,1])
        map.drawparallels(range(-90,91,30), labels=[1,0,0,1])

    wasinteractive = isinteractive()
    interactive(False)
    lines = []

    lines += [ lin.plot(map, color=color, alpha=alpha, linewidth=linewidth, lon_cyc=lon_cyc, lat_mirror=lat_mirror) for lin in lins ]

    interactive(wasinteractive)

    return([line for line in flatten(lines)], map)
#}}} end plotlinmap

def shp2lins(shapefile, stresscalc=None): #{{{
    """
    Create a list of L{Lineament} objects from an ESRI shapefile.

    The shapefile must contain one or more linear features, defined in geographic
    coordinates using decimal degrees.  The shapefile is read in using the GDAL/OGR
    geospatial library.

    """

    # This is the list of lineaments we are going to return eventually:
    linlist = []

    # First read in the shapefile as an OGR "data source"
    data_source = ogr.Open(shapefile, update=0)

    # OGR data sources can in general have many data layers, but ours will only have
    # one, containing linear features.  Get that layer:
    lineaments = data_source.GetLayer(0)

    # Within that layer, there will be potentially many features: individual
    # lineaments which should be extracted and turned into Lineament objects
    # independently.  Each one should be extracted and made into a Lineament object
    # independently:

    ogr_lin_feat = lineaments.GetNextFeature()

    while ogr_lin_feat is not None:
        ogr_lin_geom = ogr_lin_feat.GetGeometryRef()
        pointlist = array([ ogr_lin_geom.GetPoint(i)[0:2] for i in range(ogr_lin_geom.GetPointCount()) ])

        if(len(pointlist) > 0):
            newlats = radians(pointlist[:,1])
            newlons = radians(pointlist[:,0])
            linlist.append(Lineament(lons=newlons, lats=newlats, stresscalc=stresscalc))

        ogr_lin_feat = lineaments.GetNextFeature()

    return linlist
# }}} end shp2lins

def lins2kml(lins=[], kmlfile=None): #{{{ TODO: WRITE IT!
    """
    Create a KML file from a list of L{Lineament} objects, defining their
    shapes.

    """
    pass
#}}} end lins2kml

def lins2shp(lins=[], shapefile=None): #{{{ TODO: WRITE IT!
    """
    Create a shapefile from a list of L{Lineament} objects, including additions
    to the attribute tables that have been calculated so that those attributes
    can be displayed easily within ArcGIS.

    """
    pass
#}}} end lins2shp
