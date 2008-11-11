from numpy import *
from scipy.stats.mmorestats import idealfourths
from pylab import *
from mpl_toolkits.basemap import Basemap

# Open Source Geospatial libraries:
from osgeo import ogr

class Lineament(object): #{{{1
    """
    A one dimensional feature on the surface of a spherical satellite.

    """

    def __init__(self, vertices=None, stresscalc=None, failure_mode="tens_frac", fits=None, doppels=None, is_doppel=False, proto_lin=None, proto_mhd=None, backrot=0): #{{{2
        """
        Create a lineament from a given list of (lon,lat) points.

        vertices is a list of (lon,lat) tuples in radians, representing the
        vertices making up a polyline.  East and North are taken to be
        positive.
        
        stresscalc is a L{satstress.StressCalc} object defining the field we
        are going to be comparing the lineament to.
        
        failure_mode describes the mechanism we are assuming.  Currently only
        tensile fracture is accepted.

        fits is an array or list of (b, delta_rms(b), dbar(b)) tuples defining
        how well the lineament matches its stresscalc, assuming failure_mode as
        the formation mechanism.

        doppels is a list of synthetic lineaments that are meant to approximate
        the mapped lineament.

        is_doppel indicates whether this Lineament itself is a doppelganger.

        proto_lin defines the prototype lineament this Lineament is based on,
        if it is a doppelganger.

        proto_mhd is the MHD from the prototype lineament to this lineament, in
        degrees of arc on the surface of the satellite, if it is a
        doppelganger.

        backrot is the amount of longitudinal translation that at which the
        doppelganger was created, if it is a doppelganger.

        """

        # This data is intrinsic, for the purposes of NSR comparisons
        if vertices is None:
            self.vertices = []
        else:
            self.vertices = vertices

        self.stresscalc   = stresscalc
        self.failure_mode = failure_mode

        # These are semi-private derived/cached data:
        if doppels is None:
            self.__doppelgangers = []
        else:
            self.__doppelgangers = doppels

        if fits is None:
            self.__fits = []
        else:
            self.__fits = array(fits)

        # Extra info for the doppelgangers:
        self.is_doppel = is_doppel
        self.proto_lin = proto_lin
        self.proto_mhd = proto_mhd
        self.backrot   = backrot
    #}}}2

    def __str__():
        """
        """
        pass

    def length(self): #{{{2
        """
        Return the total length of the L{Lineament} in radians of arc, along
        the surface of the sphere.

        """
        sph_len = 0

        if len(self.vertices) > 1:
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

        if len(self.vertices) > 1:
            for i in range(len(self.vertices)-1):
                lon1, lat1 = self.vertices[i]
                lon2, lat2 = self.vertices[i+1]
                seglist.append(Segment(lon1, lat1, lon2, lat2))

        return(seglist)
    #}}}2 end segments

    def midpoints(self): #{{{2
        """
        Return a list of the midpoints of each line segment in the Lineament.

        """
        mp = None

        if len(self.vertices) == 0:
            mp = None
        elif len(self.vertices) == 1:
            mp = self.vertices
        elif len(self.vertices) > 1:
            mp = [ seg.midpoint() for seg in self.segments() ]

        return(mp)
    #}}}2

    def longitudes(self): # {{{2
        """
        Return a list of the latitude values from the lineament's vertices

        """
        if len(self.vertices) == 0:
            return None
        else:
            return( [ v[0] for v in self.vertices ] )
# }}}2

    def fixed_longitudes(self): # {{{2
        """
        Return a list of the longitude values from the lineament's vertices,
        without discontinuities, and shifted to lie within -pi < lon < pi.
    
        """
        lons = self.longitudes()
        if lons is not None and len(lons) > 1:
            # Get rid of any discontinuities introduced during mapping
            for i in range(len(lons)-1):
                if abs(lons[i]-lons[i+1]) > abs(lons[i]-(lons[i+1]+2*pi)):
                    lons[i+1] += 2*pi
                elif abs(lons[i]-lons[i+1]) > abs(lons[i]-(lons[i+1]-2*pi)):
                    lons[i+1] -= 2*pi

            # Ensure that -pi < lon < pi
            while max(lons) > pi or min(lons) < -pi:
                if max(lons) > pi:
                    pi_sign = -1
                else:
                    pi_sign = +1

                lons = [ lon+(pi_sign*pi) for lon in lons ]

        return(lons)
    # }}}2

    def latitudes(self): # {{{2
        """
        Return a list of the latitude values from the lineament's vertices"

        """
        if len(self.vertices) == 0:
            return None
        else:
            return( [ v[1] for v in self.vertices ] )
    #}}}2

    def lonshift(self, b): #{{{2
        """
        Return the lineament shifted in longitude by 'b'.

        """
        shifted_vertices = [ (v[0]+b, v[1]) for v in self.vertices ]
        return(Lineament(shifted_vertices, stresscalc=self.stresscalc, failure_mode=self.failure_mode))
    #}}}2

    def poleshift(self, paleo_npole_lon, paleo_npole_lat): #{{{2
        """
        Return a Lineament object representing the location and orientation of
        self, when the point defined in modern lon,lat terms by paleo_npole_lon
        paleo_npole_lat was the north pole.

        """
        transformed_vertices = [ paleopole_transform(paleo_npole_lon, paleo_npole_lat, v[0], v[1]) for v in self.vertices ]
        return(Lineament(vertices=transformed_vertices, stresscalc=self.stresscalc, failure_mode=self.failure_mode) )
    #}}}2

    def eastend(self): #{{{2
        """
        Return a tuple representing the easternmost endpoint of the L{Lineament}.

        If only one vertex, report it.  If no vertices, return None.

        """

        # if we don't have any vertices, there's no end to report:
        if len(self.vertices) == 0:
            return None
        else:
            (lon1, lat1) = self.vertices[0]
            (lon2, lat2) = self.vertices[-1]

        # if we have more than one point, then a direction makes sense:
        if len(self.vertices) > 1 and spherical_azimuth(lon1, lat1, lon2, lat2) < pi:
            return(lon2, lat2)
        # if we only have one point, it doesn't matter which end we report
        else:
            return(lon1, lat1)
    # }}}2

    def westend(self): #{{{2
        """
        Return a tuple representing the westernmost endpoint of the L{Lineament}.

        If only one vertex, report it.  If no vertices, return None.

        """
        # if we don't have any vertices, there's no end to report:
        if len(self.vertices) == 0:
            return None
        else:
            (lon1, lat1) = self.vertices[0]
            (lon2, lat2) = self.vertices[-1]

        # if we have more than one point, then a direction makes sense:
        if len(self.vertices) > 1 and spherical_azimuth(lon1, lat1, lon2, lat2) < pi:
            return(lon1, lat1)
        # if we only have one point, it doesn't matter which end we report
        else:
            return(lon2, lat2)
    #}}}2

    def bs(self, delta_rms=None, dbar=None, delta_ismin=False, dbar_ismin=False, delta_max=None, dbar_max=None, winwidth=15): #{{{2
        """
        Return a sorted array of all the b values for which fits have been
        calculated, and which satisfy the keyword requirements.

        delta_ismin and dbar_ismin, if True, require that the b values returned
        represent local minima within the delta_rms(b) and dbar(b) curves,
        respectively.

        delta_max and dbar_max, if not None, represent maximum acceptable
        values of those variables for which to return b values.  delta_max is
        measured in radians, and dbar_max is a fraction of the overall length
        of the lineament.

        """
        # special cases:
        # find b values associated with a particular value of delta_rms or dbar
        if delta_rms is not None:
            return(self.fit()[find(self.delta_rms()==delta_rms),0])
        if dbar is not None:
            return(self.fit()[find(self.dbar()==dbar),0])

        fits = self.__fits

        # Because the local minima function still expects an evenly
        # distributed array of points, we need to find the delta minima
        # and the dbar minima separately...

        if len(fits) > 0:
            if delta_ismin:
                b_delta_mins, delta_mins = local_minima(self.bs(),self.delta_rms(), winwidth=winwidth)
            else:
                b_delta_mins, delta_mins = self.__fits[:,0], self.__fits[:,1]

            if dbar_ismin:
                b_dbar_mins, dbar_mins = local_minima(self.bs(),self.dbar(), winwidth=winwidth)
            else:
                b_dbar_mins, dbar_mins = self.__fits[:,0], self.__fits[:,2]

            fits = array([ f for f in fits if f[1] in delta_mins and f[2] in dbar_mins ])

            # Now that we've got only the fits that meet the local minima requirements,
            # we can go ahead and filter by their magnitudes:
            if delta_max is not None:
                fits = array([ f for f in fits if f[1] < delta_max ])
            if dbar_max is not None:
                fits = array([ f for f in fits if f[2] < dbar_max ])

            return(sort(fits[:,0]))

        else:
            return(array([]))
    #}}}2

    def delta_rms(self, bs=[], w_length=True, w_stress=True): #{{{2
        """
        Calculate the root mean square (RMS) angle between observed and
        predicted orientations of failure, assuming the lineament is shifted in
        longitude by b radians from it's present location, with expected
        failure orientations determined by the L{StressCalc} object, time_sec,
        and the specified failure mode.  Optionally, weight by lineament
        segment length (if w_length=True) and importance of that location
        within the stress field (if w_stress=True)

        """

        if bs == []:
            bs = self.bs()

        for x in bs:
            if x not in self.bs():
                self.fit(bs=[x,])

        deltas = []
        for b in bs:
            idx = find(b==self.__fits[:,0])
            deltas.append(self.__fits[idx,1])

        return(array(deltas))
    #}}}2

    def doppel(self, bs=[]): #{{{2
        """
        Return a list containing doppelgangers for the specified values of b.
        Calculate them if they do not already exist.  Each entry in the list
        will contain the doppelgangers for a single value of b.  If only a single
        value of b is given as an integer, return a flat list of the
        doppelgangers at that backrotation.

        """

        if bs == []:
            bs = self.bs()

        doplist = []
        for x in bs:
            # create a list of all doppelgangers at b=x:
            dops = []
            for d in self.__doppelgangers:
                if d.backrot == x:
                    dops.append(d)

            # if there were none, calculate them:
            if len(dops) == 0:
                # Generate two new endpoint doppelgangers:
                new_east_doppel = self.doppelgen_endpoint(propagation="east", lonshift=x)
                new_west_doppel = self.doppelgen_endpoint(propagation="west", lonshift=x)
                dops = [ new_east_doppel, new_west_doppel ]
                # add them to the Lineament's list of doppelganger lineaments:
                self.__doppelgangers += dops

            # Whether we had to calculate them or not, add the doppelgangers
            # to the list that we're returning:
            doplist += [ dops ]

        return(doplist)
    #}}}2

    def plot(self, map=None, fixlon=False, color='black', alpha=1.0, lw=1.0): #{{{
        """
        Plot the lineament on the provided map, or create a new one.

        """
        return(plotlinmap([self,], map=map, fixlon=fixlon, color=color, lw=lw))

    def plotdoppels(self, bs=[], backrot=True, map=None, fixlon=False, color='black', alpha=0.5, lw=1.0): #{{{
        """
        Produce a plot showing a lineament's doppelgangers, in association with
        the lineament itself.
        
        Show only those doppelgangers whose backrot is contained in the list b.
        If b is [], plot all of them.  If b is a single value, show only the
        doppelgangers having backrot==b.

        If backrot is True, shift all of the doppelgangers in longitude such
        that they share their initial starting point with one of the prototype
        lineament's endpoints.  If false, show them at the longitude where they
        were generated.

        """

        if backrot is True:
            dops2plot = [ d.lonshift(-d.backrot) for d in self.__doppelgangers if d.backrot in bs ]
        else:
            dops2plot = [ d for d in self.__doppelgangers if d.backrot in bs ]

        if len(dops2plot) > 0:
            lines_out, map_out = plotlinmap(dops2plot, map=map, fixlon=fixlon, color=color, alpha=alpha, lw=lw)
        else:
            lines_out, map_out = None, map

        return(lines_out, map_out)
    # }}}2
        
    def dbar(self, bs=[]): #{{{2
        """
        Return the mean L{mhd} from the L{Lineament} to its doppelgangers
        having the given backrotation.  If no doppelgangers have that
        backrotation, generate them.

        """

        if bs == []:
           bs = self.bs()

        for x in bs:
            if x not in self.bs():
                self.fit(bs=[x,])

        dbars=[]
        for b in bs:
            idx = find(b==self.__fits[:,0])
            dbars.append(self.__fits[idx,2])

        return(array(dbars))
    #}}}2

    def fit(self, bs=[], w_length=True, w_stress=True): #{{{2
        """
        Return a list of tuples (b, delta_rms(b), dbar(b)) describing how well
        the lineament matches the stress field in question (self.stresscalc)
        given the chosen failure mode (self.failure_mode - limited to tens_frac
        currently), at the given value, or values, of b.  If only a single
        value of b is given, return a tuple, not a list of tuples.
        
        If a fit has not been calculated at a given value of b, calculate it.

        if bs=[], return all fits.

        """

        if bs == []:
            bs = self.bs()

        for x in bs:
            # if we have no fit at this value of b, calculate it:
            if x not in self.bs():
                # first we need delta_rms
                print("  calculating delta_rms(b=%g)" % (degrees(x),))
                newdelta = self.lonshift(x).stresscomp(self.stresscalc)
                # doppelgangers are required to calculate dbar
                dops = self.doppel(bs=[x,])
                print("  calculating dbar(b=%g)" % (degrees(x),))
                newdbar = mean([ d.proto_mhd for d in dops[0] if d.backrot == x])/self.length()

                # Add the first fit by itself, instead of vstacking
                if len(self.__fits) == 0:
                    self.__fits = array([(x,newdelta,newdbar)])
                else:
                    self.__fits = vstack([self.__fits,(x,newdelta,newdbar)])

        fits = array([])
        for b in bs:
            idx = find(b==self.__fits[:,0])
            if len(fits) == 0:
                fits = self.__fits[idx]
            else:
                fits = vstack([fits,self.__fits[idx]])

        return(fits)
    #}}}2

    def doppelgen_endpoint(self, time_sec=0.0, propagation="east", lonshift=0.0): #{{{2
        """
        Synthesize and return a new lineament consistent with the given
        stresscalc object and failure mode, of the same length, and having one
        of the same endpoints as self (the same west endpoint if propagation is
        "east", and the same east endpoint if propagation is "west")

        """

        if propagation=="east":
            (init_lon, init_lat) = self.westend()
        else:
            assert (propagation=="west")
            (init_lon, init_lat) = self.eastend()

        init_lon += lonshift 

        print("  generating %s doppel(b=%g)" % (propagation, degrees(lonshift)))
        newdoppel = lingen(self.stresscalc, init_lon=init_lon, init_lat=init_lat, max_length=self.length(), time_sec=time_sec, propagation=propagation, failure=self.failure_mode)

        newdoppel.backrot = lonshift
        newdoppel.is_doppel = True
        newdoppel.proto_lin = self
        newdoppel.proto_mhd = self.mhd(newdoppel.lonshift(-newdoppel.backrot))

        return(newdoppel)
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

def plotlinmap(lins, map=None, fixlon=False, color='black', alpha=1.0, lw=1.0): #{{{
    """
    Plot a map of the lineaments listed in 'lins'.  Plot it to 'map' if given.

    """
    if map is None:
        map = Basemap()
        map.drawmapboundary(fill_color="white")
        map.drawmeridians(range(-180,181,30), labels=[1,0,0,1])
        map.drawparallels(range(-90,91,30), labels=[1,0,0,1])

    interactive(False)
    lines = []

    for lin in lins:
        if lin is not None:
            if fixlon:
                x,y = map(degrees(lin.fixed_longitudes()),degrees(lin.latitudes()))
            else:
                x,y = map(degrees(lin.longitudes()),degrees(lin.latitudes()))

        lines.append(map.plot(x, y, color=color, alpha=alpha, linewidth=lw))
    interactive(True)
    show()

    return(lines, map)
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
        pointlist = [ ogr_lin_geom.GetPoint(i)[0:2] for i in range(ogr_lin_geom.GetPointCount()) ]
        linlist.append(Lineament(radians(pointlist), stresscalc=stresscalc))
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

    if not done:
        second_half = lingen(stresscalc, init_lon=init_lon, init_lat=init_lat, max_length=max_length, time_sec=time_sec, propagation="west", failure=failure)
        second_half_reversed = second_half.vertices[:0:-1]
        vertices = vstack([second_half_reversed, array(vertices)])

    newlin = Lineament(vertices=vertices, stresscalc=stresscalc)

    if lonshift:
        if lonshift == "rand":
            lonshift = 2*pi*(rand()-0.5)

        newlin = newlin.lonshift(lonshift)

    return newlin

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
        if len(lin._Lineament__doppelgangers) > 0:
            for dop in lin._Lineament__doppelgangers:
                new_doppels.append(Lineament(dop.vertices, stresscalc=dop.stresscalc, failure_mode=dop.failure_mode, fits=dop._Lineament__fits,\
                                             is_doppel=True, proto_lin=dop.proto_lin, proto_mhd=dop.proto_mhd, backrot=dop.backrot))

        newlin = Lineament(lin.vertices, stresscalc=lin.stresscalc, failure_mode=lin.failure_mode, fits=lin._Lineament__fits, doppels=new_doppels)
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

def local_minima(x_in, y_in, winwidth=15): #{{{
    """
    Given two arrays x and y(x), find the values of x which correspond to the
    smallest local minima within the function y(x), within a window winwidth
    wide.  winwidth=0 will return all local minima, no matter how local.

    Returns two arrays, x_out, y_out, corresponding to the x and y values that
    match the given criteria.

    Returns a list of indices at which the minima were found, and a list of the
    minima, sorted in order of increasing minimum.  The keyword argument deg_win
    determines how close two local minima are allowed to be to one another.  If
    two local minima are found less than deg_win degrees apart, then the lowest
    of them is taken as the real minimum.

    For now this is hacked to assume that the x values are evenly distributed
    across the entire range of possible values, from 0-180 degrees (or 0-pi
    radians).

    """
    from scipy.ndimage.filters import minimum_filter as min_filter

    # winwidth as passed in is in degrees, but the 'size' that min_filter needs
    # is a number of neighboring points to use.  First we need to convert from
    # degrees to indices, rounding up
    deg_win = int(round(winwidth/(180.0/(len(y_in)-1))))

    # then we can calculate the minima using min_filter
    y_mins = min_filter(y_in, size=deg_win, mode="wrap")

    # but really we just want the values where y(x) = y_min
    y_out = [ y for y,y_min in zip(y_in,y_mins) if y == y_min ]

    # we need extract the x values that correspond do these y values:
    good_indices = ravel([ find(y_in==y) for y in y_out ])
    x_out = [ y_in[idx] for idx in good_indices ]

    return(x_out, y_out)
#}}}

def filter(lins, \
           min_length    = None, max_length    = None, \
           min_sinuosity = None, max_sinuosity = None, \
           min_latitude  = None, max_latitude  = None, \
           min_longitude = None, max_longitude = None, \
           min_delta_rms = None, max_delta_rms = None, \
           min_dbar      = None, max_dbar      = None): #{{{ TODO: Testing/debugging
    """
    Takes a list of lineaments, lins, and returns a list of lineaments that
    match the filtering criteria.

    Each parameter has a min_ and max_ value.  Only lineaments having at least
    the minimum, and at most the maximum, are returned.  All values are None by
    default, and if called with no keyword arguments, all lineaments in lins
    are returned.

    length: spherical length of the lineament in radians of arc on the
    satellite surface.

    sinuosity: ratio of lineament length overall to spherical distance between
    endpoints.

    latitude & longitude: requires all vertices to be greater than min_ and
    less than max_ values, stated in units of radians.  Latitude ranges from 0
    to pi.  Longitude is constrained to range from from -pi to pi.

    delta_rms: refers to the smallest value of delta_rms associated with the
    lineament, in radians.

    dbar: refers to the smallest value of dbar associated with the lineament.

    """
    if min_length is not None:
        lins = [ l for l in lins if l.length() >= min_length ]
    if max_length is not None:
        lins = [ l for l in lins if l .length() <= max_length ]

    if min_sinuosity is not None:
        lins = [ l for l in lins if l.sinuosity() >= min_sinuosity ]
    if max_sinuosity is not None:
        lins = [ l for l in lins if l.sinuosity() <= max_sinuosity ]

    if min_latitude is not None:
        lins = [ l for l in lins if min(l.latitudes()) >= min_latitude ]
    if max_latitude is not None:
        lins = [ l for l in lins if max(l.latitudes()) <= max_latitude ]

    if min_longitude is not None:
        lins = [ l for l in lins if min(l.fixed_longitudes()) >= min_longitude ]
    if max_longitude is not None:
        lins = [ l for l in lins if max(l.fixed_longitudes()) <= max_longitude ]

    if min_delta_rms is not None:
        lins = [ l for l in lins if len(l.delta_rms()) > 0 and max(l.delta_rms()) >= min_delta_rms ]
    if max_delta_rms is not None:
        lins = [ l for l in lins if len(l.delta_rms()) > 0 and min(l.delta_rms()) <= max_delta_rms ]

    if min_dbar is not None:
        lins = [ l for l in lins if len(l.dbar()) > 0 and max(l.dbar()) >= min_dbar ]
    if max_dbar is not None:
        lins = [ l for l in lins if len(l.dbar()) > 0 and min(l.dbar()) <= max_dbar ]

    return lins
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

def paleopole_transform(paleo_npole_lon, paleo_npole_lat, lon_in, lat_in): #{{{
    """
    Transforms the location of a point on the surface of a sphere, defined by
    (lon_in,lat_in) in east-positive longitude, returning the (lon,lat)
    location that point would be located at if the point defined by
    (paleo_npole_lon,paleo_npole_lat) is moved directly north until it is at
    the north pole.

    """

    # First we convert the point to be transformed into Cartesian coords.
    # Since we only care about the lon,lat position in the end, we'll treat the
    # the body as a unit sphere.  Remember that sphere2xyz needs CO-latitude:
    colat_in = pi/2 - lat_in
    xyz_in = array(sphere2xyz(1.0, colat_in, lon_in))

    # Now, remember that what we're doing is bringing a wayward pole back to
    # the top of the sphere... which means the angles we're interested in are
    # actually -colat, -lon:
    alpha = pi/2 + paleo_npole_lon
    beta  = pi/2 - paleo_npole_lat
    gamma = 0

    # The definition of the X-Z rotation matrix:
    rot_mat = array([ [ cos(alpha)*cos(gamma)-sin(alpha)*cos(beta)*sin(gamma), -cos(alpha)*sin(gamma)-sin(alpha)*cos(beta)*cos(gamma),  sin(beta)*sin(alpha) ],\
                      [ sin(alpha)*cos(gamma)+cos(alpha)*cos(beta)*sin(gamma), -sin(alpha)*sin(gamma)+cos(alpha)*cos(beta)*cos(gamma), -sin(beta)*cos(alpha) ],\
                      [                 sin(beta)*sin(gamma),                                     sin(beta)*cos(gamma),                       cos(beta)      ] ])

    # Do the matrix multiplication using dot():
    x_out, y_out, z_out = dot(xyz_in, rot_mat)

    # Transform back to spherical coordinates:
    r_out, theta_out, phi_out = xyz2sphere(x_out, y_out, z_out)

    # and back to lon/lat from colon/lat:
    lon_out = mod(phi_out + alpha,2*pi)
    lat_out = pi/2 - theta_out # co-lat to lat

    return(lon_out, lat_out)
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
