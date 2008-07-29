# Numeric Python routines, like arcsin2:
from numpy import *
from numpy import linalg
from numpy.ma.mstats import idealfourths
from pylab import *

# Open Source Geospatial libraries:
from osgeo import ogr
#from osgeo import gdal
#from osgeo import osr
#from osgeo import gdal_array
#from osgeo import gdalconst

class Lineament(object): #{{{1
    """A one dimensional feature on the surface of a spherical satellite.

    Described as a list of (longitude,latitude) points.  May be transformed as
    a whole, using translations and rotations.

    May be read in from and output to the ESRI "generate" file format.

    """
    def __init__(self, vertices, fits=[]): #{{{2
        """
        Create a lineament from a given list of (lon,lat) points.

        vertices is a list of (lon,lat) tuples in decimal degrees, representing
        the vertices making up a polyline.  East and North are taken to be
        positive.
        
        fits is a list of fits to some stress field, as the lineament is
        rotated back through pi radians of longitude.

        """

        self.vertices = radians(vertices)
        self.fits = fits

    #}}}2

    def __str__(self): #{{{2
        """
        Output the lineament as list of (lon, lat) tuples.

        """
        return("\n".join([ str(degrees(v)) for v in self.vertices ]))
    #}}}2

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

    def lonshift(self, b): #{{{2
        """
        Return the lineament shifted in longitude by 'b'.

        """
        return(Lineament([(degrees(v[0]+b), degrees(v)[1]) for v in self.vertices]))
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

    def doppelganger(self, stresscalc, time_sec=0.0, propagation="east", failuremode="tensilefracture", noise=None): #{{{2
        """
        Synthesize and return a new lineament consistent with the given
        stresscalc object, of the same length, and having one of the same
        endpoints as self (the same west endpoint if propagation is "east",
        and the same east endpoint if propagation is "west")

        """

        if propagation=="east":
            (init_lon, init_lat) = self.westend()
        else:
            assert (propagation=="west")
            (init_lon, init_lat) = self.eastend()

        return(lingen(init_lon, init_lat, stresscalc, max_length=self.length(), time_sec=time_sec, propagation=propagation, failuremode=failuremode, noise=None))

    #}}}2

    def draw(self): #{{{2
        """
        Use pylab to draw the lineament quick-and-dirty.

        """
        plot([ degrees(v[0]) for v in self.vertices ], [ degrees(v[1]) for v in self.vertices ])
        grid(True)
        axis([0,360,-90,90])
        xticks( range(0,361,30) )
        xlabel("longitude")
        yticks( range(-90,91,30) )
        ylabel("latitude")
        show()
    # }}}2

    def stresscomp(self, stresscalc, time_sec=0.0, failuremode = "tensilefracture", w_length = True, w_anisotropy = False, w_stressmag = False): # {{{2
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

        # returns a list of line segments (non-poly-lines... each defined by
        # only two endpoints on the surface of a sphere.
        segments = self.segments()

        for segment in segments:
            segment.delta = pi/2.0
            segment.w_length = w_length
            segment.w_anisotropy = w_anisotropy
            segment.w_stressmag = w_stressmag

            stress_tensor = stresscalc.tensor(theta = (pi/2.0)-segment.midpoint()[1],\
                                                phi = segment.midpoint()[0],\
                                                  t = time_sec )

            (tens_mag, tens_az, comp_mag, comp_az) = tensor2pc(stress_tensor)

            # Delta = angle between observed and expected fractures
            #   - this is the raw material for the fit-metric
            #   - set to pi/2 if stress is inapplicable to failure mode
            #   - will be scaled according to several weights below...
            #   - will be different, depending on the failure mode
            if (failuremode == "tensilefracture"):
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
        
            # W_{anisotropy}:
            #   - want to reduce the weight of segments experiencing nearly
            #     isotropic stresses, since small background perturbations can
            #     significantly re-orient the principal compoents
            #   - because the background perturbations have some finite amplitude,
            #     presumably unrelated to the magnitude of the stresses we're
            #     calculating, we can have this related linearly to the difference
            #     between the magnitudes of the two principal components.
            #   - May want to saturate at some threshold stress - possibly the
            #     strength of the ice or some fraction thereof?
            #   - for now, simply multiply by (sigma_1 - sigma_2)

            if segment.w_anisotropy:
                segment.w_anisotropy = (tens_mag - comp_mag)

            # W_{stressmagnitude}:
            #   - should be linear with the magnitude of the stress relative to the
            #     strength of the ice.
            #   - Should possibly allow saturation at some threshold above the
            #     strength of the ice... say 1x, or 2x?  Dunno.
            #   - for now, simply use ratio of stress to ice strength

            if segment.w_stressmag:
                segment.w_stressmag = abs(tens_mag/stresscalc.stresses[0].satellite.layers[-1].tensile_str)

        # Now that we have all of the error weighting functions set for each segment:
        #    - segment.delta
        #    - segment.w_length
        #    - segment.w_anisotropy
        #    - segment.w_stressmagnitude
        # we can run through the entire list of segments making up the
        # lineament, and calculate the aggregate fit-metric.  This is done
        # outside the loop that's calculating the per-segment statistics so
        # that we can play around with various ways of combining them...

        linfit = 0
        rmsfit = 0
        for segment in segments:

            segfit = segment.delta

            if segment.w_length:
                segfit *= segment.w_length
            if segment.w_anisotropy:
                segfit *= segment.w_anisotropy
            if segment.w_stressmag:
                segfit *= segment.w_stressmag

            linfit += segfit
            rmsfit += linfit**2
        
        rmsfit = sqrt(rmsfit)

        #return(rmsfit)
        return(linfit)

    # }}}2 end stresscomp

    def calc_fits(self, stress, nb, time_sec=0.0, failuremode="tensilefracture", w_length=True, w_anisotropy=False, w_stressmag=False): #{{{2
        """
        Given a stresscalc object (stress), calculate the fits between that stress
        and the lineament, at a number of backrotations (nb)

        """
        self.fits=[]
        for b in linspace(0,pi,num=nb):
            self.fits.append(self.lonshift(b).stresscomp(stress, time_sec=time_sec, failuremode = "tensilefracture", w_length=w_length, w_anisotropy = w_anisotropy, w_stressmag = w_stressmag))
    #}}}2

    def best_fit(self): #{{{2
        """
        Return the best fit listed in the fits data member, and the amount of
        backrotation at which it occured.

        """
        return(min(self.fits))
    # }}}2

    def best_b(self): #{{{2
        """
        Return the amount of backrotation at which the best fit occurs

        """
        return(linspace(0,pi,num=len(self.fits))[self.fits.index(self.best_fit())])
    # }}}2

    def fits_width(self): #{{{2
        """
        Return a metric of the width of the distribution of fits between the
        lineament and the stress field it has been compared to.  A small width
        means the lineament fits well in a narrow range of backrotations, and
        a large width means it fits equally well (or more likely, poorly) over
        a wide range of backrotations.

        For now we are using the interquartile range as the measure of fit
        distribution width.
        
        """
        return(idealfourths(self.fits))

    #}}}2

    def centroid(self): #{{{2
        """
        Return the lat/lon position of the centroid of the L{Lineament}
        object.

        CURRENTLY NOT IMPLEMENTED.
        
        """
        pass
    #}}}2 end centroid

    def mean_hausdorff_distance(self, lin2): #{{{2
        """

        """
        pass
    #}}}2 end mean_hausdorf_distance

    def symmetric_mean_hausdorff_distance(self, lin2): #{{{2
        """

        """
        pass
    #}}}2 end symmetric_mean_hausdorff_distance

#}}}1

class Segment(object): #{{{
    """
    A great-circle path line defined by two end points on the surface of a sphere.

    """

    def __init__(self, lon1, lat1, lon2, lat2, parent_lin=None):
        self.lon1 = lon1
        self.lat1 = lat1
        self.lon2 = lon2
        self.lat2 = lat2
        self.parent_lin = parent_lin

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

# Helper functions that aren't part of any object class:

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
        linlist.append(Lineament(pointlist))
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

def tensor2pc(stress_tensor): # {{{
    """
    Transform the stress tensor into principal components, and return them
    in magnitude-azimuth form.

    Always returns azimuths between 0 and pi.  The more tensile of the two
    principal components is returned first.

    """
#    (eigval_A, eigval_B), ((theta_A, phi_A), (theta_B, phi_B)) = linalg.eigh(stress_tensor)
    (eigval_A, eigval_B), ((theta_A, phi_A), (theta_B, phi_B)) = eigen2(stress_tensor)
    az_A = mod(arctan2(phi_A,-theta_A), pi)
    az_B = mod(arctan2(phi_B,-theta_B), pi)

    if (eigval_A > eigval_B):
        return(eigval_A, az_A, eigval_B, az_B)
    else:
        return(eigval_B, az_B, eigval_A, az_A)
# }}}

def eigen2(tensor): # {{{

    a = tensor[0][0]
    b = tensor[0][1]
    c = tensor[1][1]
    
    sqrt_thing = sqrt(a*a + 4.0*b*b - 2.0*a*c + c*c)

    lambda1   = (0.5)*(a + c - sqrt_thing)
    eig1theta = (-a + c + sqrt_thing)/(-2.0*b)
    eig1phi   = 1.0

    lambda2   = (0.5)*(a + c + sqrt_thing)
    eig2theta = (-a + c - sqrt_thing)/(-2.0*b)
    eig2phi   = 1.0

    mag1 = sqrt(eig1theta**2 + eig1phi**2)
    mag2 = sqrt(eig2theta**2 + eig2phi**2)

    eig1theta /= mag1
    eig1phi   /= mag1
    eig2theta /= mag2
    eig2phi   /= mag2

    return ((lambda1, lambda2), ((eig1theta, eig1phi), (eig2theta, eig2phi)))
# }}}

def lingen(stresscalc, init_lon="rand", init_lat="rand", max_length=2.0*pi, time_sec=0.0, propagation="both", failuremode="tensilefracture", noise=None, segment_length=0.01): # {{{
    """
    Generate a "perfect" lineament, given initial crack location, final length,
    and the stress field and failure mode.

    If noise is specified, add noise to the lineament, making it less than perfect.

    """
    lin_length = 0.0

    # If we're choosing a random point on the surface... let's make sure that
    # the points are evenly distributed all over the sphere:

    if init_lon == "rand":
        init_lon = 2*pi*rand()
    if init_lat == "rand":
        init_lat = arccos(2*rand()-1)-pi/2

    vertices = [(init_lon, init_lat),]

    ice_strength = stresscalc.stresses[0].satellite.layers[-1].tensile_str


    # Calculate the stress tensor at the given time and initial location
    stress_tensor = stresscalc.tensor(theta = (pi/2.0)-vertices[-1][1],\
                                        phi = vertices[-1][0],\
                                          t = time_sec )

    # Convert the stress tensor into principal components:
    (tens_mag, tens_az, comp_mag, comp_az) = tensor2pc(stress_tensor)

    # Make a note of the fact that we've got to do the other half too
    if propagation == "both":
        done = False
        propagation = "east"
    else:
        done = True

    while lin_length < max_length and tens_mag > ice_strength:
        # Convert the stress tensor into principal components:
        (tens_mag, tens_az, comp_mag, comp_az) = tensor2pc(stress_tensor)

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

        # Calculate the stress tensor at the given time and initial location
        stress_tensor = stresscalc.tensor(theta = (pi/2.0)-vertices[-1][1],\
                                            phi = vertices[-1][0],\
                                              t = time_sec )

    if len(vertices) > 1:
        if not done:
            second_half = lingen(stresscalc, init_lon=init_lon, init_lat=init_lat, max_length=max_length, time_sec=time_sec, propagation="west", failuremode=failuremode, noise=noise)
            second_half_reversed = second_half.vertices[:0:-1]
            vertices = vstack([second_half_reversed, array(vertices)])
        return Lineament(degrees(vertices))
    else:
        return None

#}}} end lingen

# Not yet implemented, or borken:
def lins2shp(lins=[], shapefile=None): #{{{
    """
    Create a shapefile from a list of L{Lineament} objects, including additions
    to the attribute tables that have been calculated so that those attributes
    can be displayed easily within ArcGIS.

    """
    # TODO
    pass
#}}}

def fixlon(lonlat, lonsign=1.0): # {{{
    """
    Removes longitude discontinuities within a single lineament.

    Currently BROKEN.

    """
    lon, lat = lonlat
    if lonsign != 0:
        while sign(lonsign) != sign(lon):
            lon += sign(lonsign)*2*pi

    return (lon, lat)

# }}} end fixlon
