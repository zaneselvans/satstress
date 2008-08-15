# Numeric Python routines, like arcsin2:
from numpy import *
from numpy.ma.mstats import idealfourths
from pylab import *
from mpl_toolkits.basemap import Basemap

# Open Source Geospatial libraries:
from osgeo import ogr

class Lineament(object): #{{{1
    """
    A one dimensional feature on the surface of a spherical satellite.

    """
    
    def __init__(self, vertices, fits=[], stresscalc=None, doppels=[], is_doppel=False, backrot=0): #{{{2
        """
        Create a lineament from a given list of (lon,lat) points.

        vertices is a list of (lon,lat) tuples in radians, representing the
        vertices making up a polyline.  East and North are taken to be
        positive.
        
        fits is a list of fits to some stress field (defined by stresscalc), as
        the lineament is rotated back through pi radians of longitude.

        doppels is a list of synthetic lineaments that are meant to approximate
        the mapped lineament.

        """

        self.vertices = vertices
        self.fits = fits
        self.stresscalc = stresscalc
        self.backrot = backrot
        self.is_doppel = is_doppel
        if self.is_doppel:
            self.doppels = []
        else:
            self.doppels = doppels
    #}}}2

    def __str__(self): #{{{2
        """
        Output the lineament as list of (lon, lat) tuples.

        """
        if self.stresscalc is not None:
            stressnames = ", ".join( [ stress.__name__ for stress in self.stresscalc.stresses ])
            my_sat = self.stresscalc.stresses[0].satellite.system_id
            linlen = "%g %s" % (self.length()*self.stresscalc.stresses[0].satellite.radius()/1000, 'km')
        else:
            stressnames = None
            my_sat = None
            linlen = "%g %s" % (degrees(self.length()), 'degrees')

        if len(self.fits) > 0:
            has_fits = True
        else:
            has_fits = False

        if len(self.doppels) > 0:
            has_doppels = True
        
        return("""
linhash:   %d
system_id: %s
stresses:  %s
fits[]:    %s
doppels:   %s
length():  %s
vertices:  %d
""" % (self.__hash__(), str(my_sat), str(stressnames), str(has_fits), str(has_doppels), linlen, len(self.vertices)))

    #}}}2

    def __hash__(self): #{{{2
        """
        Return the hash of the list of the lineament's vertices, so we can
        have a unique ID that corresponds to the geometry of the feature.

        """
        return self.vertices.__hash__()

    #}}}2

    def fixlon(self): # {{{
        """
        Ensures that the longitudes of all points making up the L{Lineament} lie
        within the rage -pi < lon < pi, by shifting the lineament east or west in
        increments of pi.  This is acceptable because, for stress comparisons, a
        lineament's location only matters to modulo pi.
        """

        self.vertices = zip(self.fixed_longitudes(), self.latitudes())
        return(self.vertices)

# }}}2 end fixlon

    def length(self): #{{{2
        """
        Return the total length of the L{Lineament} in radians of arc, along
        the surface of the sphere.

        """
        sph_len = 0
        for i in range(len(self.vertices)-1):
            lon1, lat1 = self.vertices[i]
            lon2, lat2 = self.vertices[i+1]
            sph_len += spherical_distance(lon1, lat1, lon2, lat2)

        return sph_len
    #}}}2 end spherical_length

    def segments(self): #{{{2
        """
        Return a list of line segments making up the L{Lineament}

        """
        seglist = []
        for i in range(len(self.vertices)-1):
            lon1, lat1 = self.vertices[i]
            lon2, lat2 = self.vertices[i+1]
            seglist.append(Segment(lon1, lat1, lon2, lat2))

        return(seglist)
    #}}}2 end segments

    def midpoints(self): #{{{2
        return([ seg.midpoint() for seg in self.segments() ])
    #}}}2

    def longitudes(self): # {{{2
        """
        Return a list of the latitude values from the lineament's vertices

        """
        return( [ v[0] for v in self.vertices ] )
# }}}2

    def fixed_longitudes(self): # {{{2
        """
        Return a list of the latitude values from the lineament's vertices,
        but shifted so as to lie within -pi < lon < pi.

        """
        fixed_lons = self.longitudes()
        while max(fixed_lons) > pi or min(fixed_lons) < -pi:
            if max(fixed_lons) > pi:
                pi_sign = -1
            else:
                pi_sign = +1

            fixed_lons = [ lon+(pi_sign*pi) for lon in fixed_lons ]

        return(fixed_lons)
# }}}2

    def latitudes(self): # {{{2
        """
        Return a list of the latitude values from the lineament's vertices"

        """
        return( [ v[1] for v in self.vertices ] )
    #}}}2

    def lonshift(self, b): #{{{2
        """
        Return the lineament shifted in longitude by 'b'.

        """
        return(Lineament([(v[0]+b, v[1]) for v in self.vertices], stresscalc=self.stresscalc))
    #}}}2

    def eastend(self): #{{{2
        """
        Return a tuple representing the easternmost endpoint of the L{Lineament}.

        """

        (lon1, lat1) = self.vertices[0]
        (lon2, lat2) = self.vertices[-1]

        if spherical_azimuth(lon1, lat1, lon2, lat2) < pi:
            return(lon2, lat2)
        else:
            return(lon1, lat1)
    # }}}2

    def westend(self): #{{{2
        """
        Return a tuple representing the westernmost endpoint of the L{Lineament}.

        """

        (lon1, lat1) = self.vertices[0]
        (lon2, lat2) = self.vertices[-1]

        if spherical_azimuth(lon1, lat1, lon2, lat2) < pi:
            return(lon1, lat1)
        else:
            return(lon2, lat2)
    #}}}2

    def doppelganger_endpoint(self, stresscalc=None, time_sec=0.0, propagation="east", failure="tens_frac", lonshift=0.0): #{{{2
        """
        Synthesize and return a new lineament consistent with the given
        stresscalc object and failure mode, of the same length, and having one
        of the same endpoints as self (the same west endpoint if propagation is
        "east", and the same east endpoint if propagation is "west")

        """

        if stresscalc is None:
            if self.stresscalc is not None:
                stresscalc = self.stresscalc
            else:
                raise(LinMissingStressCalcError())

        if propagation=="east":
            (init_lon, init_lat) = self.westend()
        else:
            assert (propagation=="west")
            (init_lon, init_lat) = self.eastend()

        init_lon += lonshift 

        doppel = lingen(stresscalc, init_lon=init_lon, init_lat=init_lat, max_length=self.length(), time_sec=time_sec, propagation=propagation, failure=failure)
        if doppel is not None:
            doppel.backrot = lonshift
        return(doppel)

    #}}}2

    def stresscomp(self, stresscalc, time_sec=0.0, failure = "tens_frac", w_length=True, w_stress=True): # {{{2
        """
        Return a metric of how likely it is this lineament resulted from a
        given failure mode and stress field.
        
        Calculates the stresses at the midpoint of each line segment making up
        the lineament, and compares the expected orientation of failure given
        that stress and failure mode, to the observed orientation of the line
        segment.

        The difference between expected and observed orientations (an angle)
        for each line segment may be weighted according to several factors:
        
            1. the proportion of the overall length of the lineament which that
               line segment represents.

            2. the magnitude of the stress, relative to the strength of the
               material being stressed, or relative to the mean global stress
               field on the body.

            3. the anisotropy of the stress field (since the principal
               components of an isotropic field can be re-oriented
               significantly with only a small change in their relative
               magnitudes)

            4. the applicability of the stress to the failure mode being
               considered (e.g. a compressive field is not applicable if we
               are considering tensile fracture)
        """
        seglist = self.segments()
        mean_global_stressdiff = stresscalc.mean_global_stressdiff()

        for segment in seglist:
            segment.delta = pi/2.0
            segment.w_length = w_length
            segment.w_stress = w_stress

            (tens_mag, tens_az, comp_mag, comp_az) = stresscalc.principal_components(theta = (pi/2.0)-segment.midpoint()[1],\
                                                                                       phi = segment.midpoint()[0],\
                                                                                         t = time_sec )

            # Delta = angle between observed and expected fractures
            #   - this is the raw material for the fit-metric
            #   - set to pi/2 if stress is inapplicable to failure mode
            #   - will be scaled according to several weights below...
            #   - will be different, depending on the failure mode
            if (failure == "tens_frac"):
                # If both PCs are compressive (negative), tensile fracture
                # does not apply, and we assign the largest possible delta:
                if comp_mag < 0 and tens_mag < 0:
                    segment.delta = pi/2
                # Otherwise, the angle between the segment at its midpoint,
                # and the orientation of the most compressive stress is the
                # difference between the observed and expected orientations
                # for tensile fractures (both comp_az and midpoint_azimuth()
                # should always be less than pi... and the biggest angle 
                # you can get between two intersecting lines (not having
                # directionality) is pi/2... hence the modulo)
                else:
                    mid_az = segment.midpoint_azimuth()
                    assert(comp_az <= pi)
                    assert(mid_az <= pi)
                    segment.delta = mod(abs(comp_az - segment.midpoint_azimuth()), pi/2)
            # Haven't worked this out for non-tensile failure yet.
            else:
                raise(UnimplementedFeatureError("Failure mode not defined"))

            # W_{length}: 
            #   - want the impact of a particular segment on the overall fit
            #     to be linearly related to its length, as a fraction of the
            #     overall length of the lineament
            #   - thus, this is simply the fractional length of the segment.

            if segment.w_length:
                segment.w_length = segment.length()/self.length()
        
            # W_{stress}:
            #   - want to reduce the weight of segments experiencing nearly
            #     isotropic stresses, since small background perturbations can
            #     significantly re-orient the principal compoents
            #   - because the background perturbations have some finite amplitude,
            #     presumably unrelated to the magnitude of the stresses we're
            #     calculating, we can have this related linearly to the difference
            #     between the magnitudes of the two principal components.
            #   - May want to saturate at some threshold stress - possibly the
            #     strength of the ice or some fraction thereof?
            #   - should be linear with the magnitude of the stress relative to the
            #     strength of the ice.
            #   - Should possibly allow saturation at some threshold above the
            #     strength of the ice... say 1x, or 2x?  Dunno.
            #   - If we want to ignore the question of ice strength... since we
            #     actually have no idea what it is, we can just use Francis'
            #     suggestion:
            #                    tens_mag - comp_mag
            #     w_stress = ----------------------------------
            #                  < tens_mag - comp_mag > (global)

            if segment.w_stress:
                segment.w_stress = (tens_mag - comp_mag)/mean_global_stressdiff

        # Now that we have all of the error weighting functions set for each segment:
        #    - segment.delta
        #    - segment.w_length
        #    - segment.w_stress
        # we can run through the entire list of segments making up the
        # lineament, and calculate the aggregate fit-metric.  This is done
        # outside the loop that's calculating the per-segment statistics so
        # that we can play around with various ways of combining them...

        fit = 0
        for segment in seglist:

            segfit = segment.delta**2

            if segment.w_length:
                segfit *= segment.w_length
            if segment.w_stress:
                segfit *= segment.w_stress

            fit += segfit
        
        return(sqrt(fit))

    # }}}2 end stresscomp

    def calc_fits(self, stresscalc, nb=19, time_sec=0.0, failure="tens_frac", w_length=True, w_stress=True): #{{{2
        """
        Given a stresscalc object, calculate the fits between that stress and
        the lineament, at a number of backrotations (nb), evenly distributed
        throughout pi radians of longitude.

        """
        # store the stresscalc object that was used in the comparison, so we know where the fits came from.
        self.stresscalc=stresscalc

        self.fits=[ self.lonshift(b).stresscomp(stresscalc, time_sec=time_sec, failure=failure, w_length=w_length, w_stress=w_stress) for b in linspace(0,pi,num=nb) ]

    #}}}2 end calc_fits

    def best_fit(self): #{{{2
        """
        Return the minimum value found within the array of fits.

        """
        return(min(self.fits))
    #}}}2

    def good_fits(self, fit_thresh=0.1, window=15): #{{{2
        """
        Return a list of backrotation values and the corresponding fits, which
        are close to being the best fit.

        fit_thresh is the proportion of the entire fit curve's amplitude away from
        the best fit that a fit can lie, and still be considered "good".  Only
        local minima are returned, meaning that those fits which are
        immediately adjacent to the best fit (or another good fit) will be
        ignored.

        """

        min_indices, min_fits = local_minima(self.fits, window=window)
        fit_amplitude = max(self.fits) - min(self.fits)
        best_fit = self.best_fit()

        good_fits = []
        for fit in min_fits:
            if (fit-best_fit) < fit_thresh*fit_amplitude:
                good_fits.append(fit)

        good_fits.sort()
        good_indices = [ self.fits.index(fit) for fit in good_fits ]
        good_bs = [ linspace(0,pi,num=len(self.fits))[idx] for idx in good_indices ]



        return(good_bs, good_fits)
    #}}}2

    def good_doppel_fits(self, fit_thresh=0.1, max_dop_mhd=0.1, window=15): #{{{2
        """
        Returns a lsit of backrotation values and the corresponding fits, which
        are close to being the best fit, similar to L{good_fits}, but with the
        additional constraint that for a particular backrotation result in a
        "good fit", it must also result in doppelgangers that are similar to
        the mapped lineament in question.

        max_dop_mhd is the greatest average MHD that a lineament can have to
        all of its doppelgangers generated at a particular value of
        backrotation, and still have that backrotation be considered a "good
        fit".  It is measured as a proportion of the mapped lineament's length.

        If the lineament in question does not have doppelgangers, then they are
        generated in order to do the average MHD comparison.

        """
        # make sure that the lineament has fits:
        if len(self.fits) == 0:
            raise(MissingFitsError())

        # calculate the good_fits:
        good_bs, good_fits = self.good_fits(fit_thresh=fit_thresh, window=window)

        # generate doppelgangers for the lineament:
        generate_doppels([self,], fit_thresh=fit_thresh, window=window)

        better_bs = []
        better_fits = []

        # restrict those good_fits results based on the MHDs:
        #   - for each good_b, find all the doppelgangers created there
        #   - calculate the MHD from self to those doppelgangers
        #   - average them
        #   - if the average is less than max_dop_mhd, retain that (b,fit)
        for b,fit in zip(good_bs,good_fits):
            if len(self.doppels) > 0:
                dops = [ dop for dop in self.doppels if dop.backrot == b ]

                avg_dop_mhd = sum([ self.mhd(dop.lonshift(-dop.backrot))/self.length() for dop in dops ])/len(dops)

                if avg_dop_mhd < max_dop_mhd:
                    better_bs.append(b)
                    better_fits.append(fit)

        return(better_bs, better_fits)
    #}}}2

    def best_b(self): #{{{2
        """
        Return the amount of backrotation at which the best fit occurs.

        """
        return(linspace(0,pi,num=len(self.fits))[self.fits.index(self.best_fit())])

    # }}}2

    def mhd(self, linB): #{{{2
        """
        Calculate the mean Hausdorff Distance between linA and linB.

        Use the midpoints of the line segments making up the lineament, instead of
        the verticies, and weight by the lengths of the line segments.

        (See: http://en.wikipedia.org/wiki/Hausdorff_distance)

        """
        return(mhd(self, linB))
    #}}}2 end mhd()

    def smhd(self, linB): #{{{2
        """
        Calculate the symmetric mean Hausdorff Distance between linA and
        linB - this is just (mhd(linA, linB) + mhd(linB, linA))/2

        """
        return( (mhd(self, linB) + mhd(linB, self))/2.0 )
    #}}}2 end smhd()

################################################################################
# NOT YET IMPLEMENTED:
################################################################################
    def bfgc(self): #{{{2
        """
        Returns two (lon,lat) tuples defining the great circle that best fits the
        L{Lineament} object.

        """
        pass
    #}}}2 end bfgc()

    def power_spectrum(self): #{{{2
        """
        Returns the power spectrum of the L{Lineament} object's distances from the
        great circle which best fits it.

        """
        pass
    #}}}2 end power_spectrum

    def centroid(self): #{{{2
        """
        Returns the spherical centroid of the lineament, defined as the lon,lat
        point which minimizes the sum of the distances from itself to the midpoints
        of all the line segments making up the lineament, weighted by the lengths
        of those line segments.

        """
        pass
    #}}}2 end centroid()

    def fit_percentile(self, stress=None, nb=19, num_lins=100): #{{{2
        """
        Calculate and return the proportion of synthetic lineaments similar to self
        (generated using doppelganger_power_spectrum) which self does better at
        fitting to the given stress.

        I.e. a fit_percentile of 0.95 would mean that self's best fit to the given
        stress is better than the best fits of 95% of the synthetic lineaments, and
        would ghus be a very robustly significant fit.

        This also sets self.fit_percentile to this value, for later reference.

        """
        pass
    #}}}2 end fit_percentile()

    def doppelganger_power_spectrum(self): #{{{2
        """
        Generate and return a lineament whose centroid has the same latitude as self,
        and which has the same power spectrum relative to its best-fit great circle
        as self, but random overall orientation, and random phase information.

        Such synthetic lineaments are used to test the statistical significance of
        a given lineament's fit to a given stress field.

        """
        pass
    #}}}2 end doppelganger_power_spectrum()

    def doppelganger_centroid(self, stress=None, seg_len=0.05): #{{{2
        """
        Generate and return a L{Lineament} object whose midpoint is
        self.centroid(), which has the same length as self, which is consistent with
        the given stress.

        """
        pass
    #}}}2 end doppelganger_centroid()

    def doppelganger_bfgc(): #{{{2
        """
        Generate and return a L{Lineament} object which is a segment of self's best
        fit great circle, has the same length as self, and whose midpoint is the
        point closest to self.centroid() on the great circle path.

        """
        pass
    #}}}2 end doppelganger_bfgc()

#}}}1 end of the Lineament class

class Segment(object): #{{{
    """
    A great-circle path line defined by two end points on the surface of a sphere.

    """

    def __init__(self, lon1, lat1, lon2, lat2, parent_lin=None):
        self.lon1 = lon1
        self.lat1 = lat1
        self.lon2 = lon2
        self.lat2 = lat2

    def length(self):
        return(spherical_distance(self.lon1, self.lat1, self.lon2, self.lat2))

    def midpoint(self):
        """
        Find the point halfway between the segment's endpoints, along a great
        circle route.

        """
        return(spherical_midpoint(self.lon1, self.lat1, self.lon2, self.lat2))

    def midpoint_azimuth(self):
        """
        Calculate the azimuth (orientation) at the midpoint of the segment.

        Because the segment is not directional, this will be an angle between 0 and pi.
        """

        return(mod(spherical_azimuth(self.midpoint()[0], self.midpoint()[1], self.lon2, self.lat2), pi))
    
    # }}} end segment

################################################################################
# Helper functions that aren't part of any object class:
################################################################################

def plotlinmap(lins, map=None, fixlon=False): #{{{
    if map is None:
        map = Basemap()

    for lin in lins:
        if lin is not None:
            if fixlon:
                map.plot(degrees(lin.fixed_longitudes()),degrees(lin.latitudes()))
            else:
                map.plot(degrees(lin.longitudes()),degrees(lin.latitudes()))

    map.drawmeridians(range(-180,181,30), labels=[1,0,0,1])
    map.drawparallels(range(-90,91,30), labels=[1,0,0,1])

    return(map)
#}}} end plotlinmap

def shp2lins(shapefile): #{{{
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
        pointlist = [ ogr_lin_geom.GetPoint(i)[0:2] for i in range(ogr_lin_geom.GetPointCount()) ]
        linlist.append(Lineament(radians(pointlist)))
        ogr_lin_feat = lineaments.GetNextFeature()

    linlist = [ x for x in linlist if len(x.vertices) > 0 ]
    return linlist
# }}} end shp2lins

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
    return((lon2,lat2))

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

def lingen(stresscalc, init_lon="rand", init_lat="rand", max_length="rand", time_sec=0.0, propagation="rand", failure="tens_frac", segment_length=0.01, lonshift=0.0): # {{{
    """
    Generate a "perfect" lineament, given initial crack location, final length,
    and the stress field and failure mode.

    """
    lin_length = 0.0

    # If we're choosing a random point on the surface... let's make sure that
    # the points are evenly distributed all over the sphere:

    if init_lon == "rand":
        init_lon = 2*pi*rand()
    if init_lat == "rand":
        init_lat = arccos(2*rand()-1)-pi/2
    if max_length == "rand":
        max_length = pi*rand()

    vertices = [(init_lon, init_lat),]

    ice_strength = stresscalc.stresses[0].satellite.layers[-1].tensile_str

    # Calculate the stresses at the given time and initial location
    (tens_mag, tens_az, comp_mag, comp_az) = stresscalc.principal_components(theta = (pi/2.0)-vertices[-1][1],\
                                                                               phi = vertices[-1][0],\
                                                                                 t = time_sec )

    # Convert the stress tensor into principal components:

    if propagation == "rand":
        if rand() > 0.5:
            propagation="east"
        else:
            propagation="west"

    # Make a note of the fact that we've got to do the other half too
    if propagation == "both":
        done = False
        propagation = "east"
    else:
        done = True

    while lin_length < max_length and tens_mag > ice_strength:
        # Choose which direction to go based on "propagation", and which
        # half of the lineament we are synthesizing
        if propagation == "east":
            prop_az = comp_az
        else:
            assert(propagation == "west")
            prop_az = comp_az + pi

        newlon, newlat = spherical_reckon(vertices[-1][0], vertices[-1][1], prop_az, segment_length)

        while (abs(newlon - vertices[-1][0]) > abs((newlon - 2*pi) - vertices[-1][0])):
            newlon = newlon - 2*pi

        while (abs(newlon - vertices[-1][0]) > abs((newlon + 2*pi) - vertices[-1][0])):
            newlon = newlon + 2*pi
        
        vertices.append((newlon, newlat))

        lin_length += segment_length

        # Calculate the stresses at the given time and initial location
        (tens_mag, tens_az, comp_mag, comp_az) = stresscalc.principal_components(theta = (pi/2.0)-vertices[-1][1],\
                                                                                   phi = vertices[-1][0],\
                                                                                     t = time_sec )

    if len(vertices) > 1:
        if not done:
            second_half = lingen(stresscalc, init_lon=init_lon, init_lat=init_lat, max_length=max_length, time_sec=time_sec, propagation="west", failure=failure)
            second_half_reversed = second_half.vertices[:0:-1]
            vertices = vstack([second_half_reversed, array(vertices)])

        newlin = Lineament(vertices, stresscalc=stresscalc)

        if lonshift:
            if lonshift == "rand":
                lonshift = 2*pi*(rand()-0.5)

            newlin = newlin.lonshift(lonshift)

        return newlin

    else:
        return None

#}}} end lingen

def update_lins(lins): # {{{
    """
    Just a quick and easy way to remember what all has to be done in order to
    update a set of lineaments to include all the newest and bestest
    functionality from the module.

    """

    newlins = []
    for lin in lins:
        new_doppels = []
        if len(lin.doppels) > 0:
            for dop in lin.doppels:
                new_doppels.append(Lineament(dop.vertices, fits=dop.fits, stresscalc=dop.stresscalc, is_doppel=True))

        newlin = Lineament(lin.vertices, fits=lin.fits, stresscalc=lin.stresscalc, doppels=new_doppels, is_doppel=False, backrot=0.0)

        newlins.append(newlin)
        
    return(newlins)
#}}} end update_lins

def aggregate_length(lins): #{{{
    """
    Return the sum of the lengths of the L{Lineaments} in lins, in radians of
    arc on the surface of the body.

    """
    return(sum([lin.length() for lin in lins]))
#}}} end aggregate_length()

def mhd(linA, linB): #{{{
    """
    Calculate the mean Hausdorff Distance between linA and linB.

    Use the midpoints of the line segments making up the lineament, instead of
    the verticies, and weight by the lengths of the line segments, because
    we're really trying to compare linear features, not sets of points.

    (See: http://en.wikipedia.org/wiki/Hausdorff_distance)

    """
    if linA is None or linB is None:
        return None

    midpoints = linA.midpoints()
    linlen = linA.length()

    # the proportion of the overall metric that will come from each point...
    length_weights = [ seg.length()/linlen for seg in linA.segments() ]

    mhd = 0.0
    for pt,wt in zip(midpoints, length_weights):
        mhd += wt*min( [ spherical_distance(pt[0], pt[1], v[0], v[1]) for v in linB.midpoints() ] )

    return(mhd)

#}}} end mhd

def smhd(linA, linB): #{{{
    """
    Calculate the symmetric mean Hausdorff Distance between linA and
    linB - this is just (mhd(linA, linB) + mhd(linB, linA))/2

    """
    if linA is None or linB is None:
        return None
    else:
        return( (mhd(linA, linB) + mhd(linB, linA))/2.0 )

#}}} end smhd

def local_minima(fits, window=15): #{{{
    """
    Find the local minima within fits, and return them and their indices.

    Returns a list of indices at which the minima were found, and a list of the
    minima, sorted in order of increasing minimum.  The keyword argument window
    determines how close two local minima are allowed to be to one another.  If
    two local minima are found closer together than that, then the lowest of
    them is taken as the real minimum.  window=1 will return all locala minima.

    """
    from scipy.ndimage.filters import minimum_filter as min_filter

    minfits = min_filter(fits, size=window, mode="wrap")

    minima = [ fit for fit,minfit in zip(fits,minfits) if fit == minfit ]
    minima.sort()

    good_indices = [ fits.index(fit) for fit in minima ]
    good_fits = [ fit for fit in minima ]

    return(good_indices, good_fits)
#}}}

def generate_doppels(lins, fit_thresh=0.1, window=15): #{{{
    """
    Attempt to generate a variety of doppelganger features for a given list of
    input features.  Add the doppelgangers to the list of doppelgangers for
    each feature.

    """

    for lin in lins:
        # First we need to determine what amounts of backrotation we want to 
        # use for generation... using the good_fits method.
        good_bs, good_fits = lin.good_fits(fit_thresh=fit_thresh, window=window)

        # Wipe out the existing doppels...
        lin.doppels = []

        for b in good_bs:
            # For now, just the endpoint doppelgangers can be constructed:
            dop_E = lin.doppelganger_endpoint(propagation="east", lonshift=b)
            if dop_E is not None:
                lin.doppels.append(dop_E)
            dop_W = lin.doppelganger_endpoint(propagation="west", lonshift=b)
            if dop_W is not None:
                lin.doppels.append(dop_W)
#}}}

################################################################################
# NOT YET IMPLEMENTED:
################################################################################
def lins2shp(lins=[], shapefile=None): #{{{
    """
    Create a shapefile from a list of L{Lineament} objects, including additions
    to the attribute tables that have been calculated so that those attributes
    can be displayed easily within ArcGIS.

    """
    pass
#}}} end lins2shp

def fit_greatcircle(vertices): #{{{
    """
    Given a set of lon,lat points (vertices) on the surface of a sphere, return
    two points defining the great circle with the best least-squares fit to the
    points (see U{http://doi.acm.org/10.1145/367436.367478}).

    """
    pass
#}}} end fit_greatcircle

def lingen_random_stresscalc(stresscalc=None, init_lon="rand", init_lat="rand", max_length="rand", propagation="both", lonshift="rand", tensile_strength=0.0, time_sec=0.0, seg_len=0.05, failure="tens_frac"): #{{{
    """
    Generate a random lineament perfectly consistent with a given stress field
    defined by stresscalc.

    """
    pass
#}}}

def lingen_greatcircle(v1, v2, seg_len=0.05): #{{{
    """
    Return a L{Lineament} object closely approximating the shortest great
    circle route between v1 and v2 (two lon,lat points).

    """
    pass
#}}} end lingen_greatcircle

def lingen_power_spectrum(v1, v2, num_steps=50, pow_spect=[]): #{{{
    """
    Given two (lon,lat) points v1 and v2 defining a unique great circle, and a
    power spectrum, pow_spect, perturb the great circle using the power
    spectrum.

    Breaks the great circle arc defined by v1,v2 into num_steps segments, and
    displaces each one perpendicular to the azimuth of the great circle, by an
    amount equal to the sum of the power components at that distance from the
    end of the great circle arc.

    Returns the perturbed lineament, whose power_spectrum() method should
    return a power spectrum equivalent to the pow_spect that was passed in.

    """
    pass
#}}} end lingen_power_spectrum
