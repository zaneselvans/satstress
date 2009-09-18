#!/usr/bin/python
import networkx as nx
from . import lineament
from . import satstress
from . import nsrhist
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab
from mpl_toolkits import basemap
from osgeo import ogr
from osgeo import osr

# define the coordinate system we're working in:
# This is lame, but for the moment, let's assume that we're on Europa:
IAU2000_Europa_SRS_WKT="""GEOGCS["Europa 2000",DATUM["D_Europa_2000",SPHEROID["Europa_2000_IAU_IAG",1564130.0,488.79062499999998]],PRIMEM["Greenwich",0],UNIT["Radian",1.000000000]]"""

class GeoSupNet(nx.MultiDiGraph): #{{{
    """
    A Geological Superposition Network (GSN) encapsulates the cross-cutting
    relationships between a number of geologic features, as well as their
    geography.  It provides methods for resolving the stratigraphic
    relationships between the features, even when some of those relationships
    are ambiguous.  The class is particularly designed to deal with linear
    tectonic features of the type found on Jupiter's moon Europa and other
    geologically active icy bodies.

    """

    # In the event that we're actually generating this GSN from an ordered list
    # of features, we should store that information for later reference.
    linstack = []
    
    def lincross(self): #{{{2
        """
        Output a list of GSN edges that can be used to re-instantiate the
        connected components of the GSN.  Each edge is a 3-tuple, consisting of
        two Lineament objects (the first is the older, the second is the
        younger) and a dictionary containing the 'lon', 'lat' and 'weight'
        (confidence) values associated with the intersection.

          lincross[0] = bottom (older) feature (Lineament object)
          lincross[1] = top (younger) feature (Lineament object)
          lincross[2]['lon'] = lon of intersection (east-positive, radians)
          lincross[2]['lat'] = lat of intersection (radians)
          lincross[2]['weight'] = confidence K of intersection (0 <= K <= 1).

        """

        return self.edges(data=True)
    #}}}2

    def get_sub_GSNs(self, minsize=2): #{{{2
        """
        Return a list of all the connected components of the GSN having size
        greater than or equal to minsize.

        """

        # Convert the GSN to an undirected graph
        # and find the connected subcomponents
        U = nx.MultiGraph()
        U.add_edges_from(self.edges())
        connected_nbunches = nx.connected_components(U)
        # re-constitute those connected bunches as a new list of GSNs:
        return([ self.subgraph(nbunch) for nbunch in connected_nbunches if len(nbunch) >= minsize ])

    #}}}2

    def get_sub_linstacks(self, linstack=None): #{{{2
        """
        Given a GSN and a list of lineaments, ordered by time of formation
        (earliest first), return a list of lists of lineaments, maintaining
        the relative orderings of the linstack, but separating the features
        depending on what connected component of the GSN they are in.
        
        If a feature exists in the linstack, but not in the GSN, it is not
        included in any of the returned lists.

        if linstack is None, it is assumed that the GSN has its own embedded
        linstack (self.linstack) as a result of having been constructed from
        a known ordering, and that is used.

        """

        if linstack is None:
            assert(len(self.linstack) > 0)
            linstack = self.linstack

        sub_GSNs = self.get_sub_GSNs()
        linstacks = []

        for sub_gsn in sub_GSNs:
            sub_stack = [ lin for lin in linstack if lin in sub_gsn ]
            linstacks.append(sub_stack)

        return(linstacks)
    #}}}2

    def net_in_degree(self, with_labels=True, weighted=True, normalized=True): #{{{2
        """
        Length normalized weighted net in degree is a decent heuristic metric
        of a feature's order of formation.  If it has lots of successors in the
        graph (i.e. lots of intersections in which it is on the bottom), per
        unit length of the feature, then all else being equal, we expect it to
        be an older feature.

        The greater the lineament density (length of lineaments per unit area
        of map), and the longer the features involved, the better this metric
        will reflect reality.

        Returns either a list (if with_labels is False) or a dictionary (if
        with_labels is True) of length normalized net in degrees.  The values
        are:

            ( in_degree - out_degree )
            --------------------------
                     lennorm

        if normalized is True then lennorm = node.length (the nodes in a GSN
        are satstress.Lineament objects, so this is the length of the feature
        as measured in radians of arc)

        If weighted is True then the weighted in and out degrees are used.
        
        If with_labels is True then the dictionary returned is keyed by node.

        """

        in_degs  = self.in_degree(with_labels=with_labels, weighted=weighted)
        out_degs = self.out_degree(with_labels=with_labels, weighted=weighted)

        net_in_degs = {}
        for lin in self.nodes():
            if normalized is True:
                norm_length = lin.length
            else:
                norm_length = 1.0

            net_in_degs[lin] = (in_degs[lin]-out_degs[lin])/norm_length

        if with_labels is False:
            net_in_degs = [ net_in_degs[lin] for lin in net_in_degs.keys() ]

        return(net_in_degs)
    #}}}

    def stratigraphic_sort(self): #{{{2
        """
        Find upper and lower bounds on the stratigraphic location of each
        feature in the GSN.

        This is done by finding the sizes of the successor and predecessor
        digraphs.  Anything that isn't in one of those graphs is of
        indistinguishable age.

        """

        stratsorts = []
        for sub_gsn in self.get_sub_GSNs():
            stratsort = {}
            rev_gsn = sub_gsn.reverse()
            for lin in sub_gsn.nodes():
                ub = len(sub_gsn) - len(nx.single_source_shortest_path_length(sub_gsn, lin))
                lb = len(nx.single_source_shortest_path_length(rev_gsn, lin)) - 1
                stratsort[lin] = (lb, ub)
            stratsorts.append(stratsort)
        
        return(stratsorts)

    #}}}2 end stratigraphic_sort

    def draw_map(self, intersections=False, ax=None, cmap=pylab.cm.jet_r, minsize=2): #{{{2
        """
        Draw a map view of the GSN using Basemap and Matplotlib.

        """
        if ax is None:
            the_fig = plt.figure(figsize=(12,6))
            ax = the_fig.add_subplot(1,1,1)
            map_ax = basemap.Basemap(ax=ax)
            map_ax.drawmeridians(range(-180,181,30), labels=[1,0,0,1])
            map_ax.drawparallels(range(-90,91,30), labels=[1,0,0,1])

        sub_GSNs = self.get_sub_GSNs()
        if minsize > 2:
            sub_GSNs = [ sub_gsn for sub_gsn in sub_GSNs if len(sub_gsn) > minsize ]

        for sub_gsn,n in zip(sub_GSNs, range(len(sub_GSNs))):
            lines, gsn_map = lineament.plotlinmap(sub_gsn.nodes(), map=map_ax, linewidth=2,\
                                                  color=cmap(float(n)/len(sub_GSNs)))
            if intersections is True:
                edges = sub_gsn.edges(data=True)
                edge_lons = [ edge[2]['lon'] for edge in edges ]
                edge_lats = [ edge[2]['lat'] for edge in edges ]
                gsn_map.scatter(np.degrees(edge_lons), np.degrees(edge_lats), lw=0, color='black', s=50)
                gsn_map.scatter(np.degrees(edge_lons)+360.0, np.degrees(edge_lats), lw=0, color='black', s=50)
                gsn_map.scatter(np.degrees(edge_lons)-360.0, np.degrees(edge_lats), lw=0, color='black', s=50)

        plt.draw()

    #}}}2 end draw_map

    def draw_sort(self, orderby='mean', trueorder=None, ax=None, colorby='graph', cmap=pylab.cm.jet_r, title="", weighted=True, normalized=True, minsize=2): #{{{2
        """
        Make a plot showing the retrieved stratigraphic sort.

        The true ordering is indicated by the order of the features in the list
        trueorder.
    
        Eventually this needs to deal with things with no known ordering, by adding
        an 'orderby' parameter, which may also have the following values:

             'nid' -- net in degree (default)
            'true' -- the actual ordering, as implied by trueorder
            'mean' -- mean stratigraphic location
              'lb' -- lower bound (earliest possible formation time)
              'ub' -- upper bound (latest possible formation time)
    
        If any of weighted or normalized are set, they are passed on to the
        GeoSupNet.net_in_degree() function.

        minsize determines how many features have to be the connected subgraph
        in order for it to be displayed.  This allows you to avoid showing a
        bunch of sorts that only include 2 or 3 features if you want to.

        """

        stratsorts = self.stratigraphic_sort()

        if minsize > 2:
            stratsorts = [ stratsort for stratsort in stratsorts if len(stratsort) > minsize ]

        if ax is None:
            the_fig = plt.figure(figsize=(12,6))
            ax = the_fig.add_axes((0.05,0.1,0.9,0.8))

        N_tot = np.array([len(subsort) for subsort in stratsorts]).sum()

        x0 = 0
        for subsort,n_col in zip(stratsorts,range(len(stratsorts))):

            # Figure out what order we're going to plot the features in, and
            # generate the appropriate x, y, and (yerr_lo, yerr_hi) arrays for
            # plotting:
            if orderby == 'mean':
                sort_key = [ 0.5*(subsort[lin][0]+subsort[lin][1]) for lin in subsort.keys() ]
                lins_to_plot = subsort.keys()
            elif orderby == 'lb':
                sort_key = [ subsort[lin][0] for lin in subsort.keys() ]
                lins_to_plot = subsort.keys()
            elif orderby == 'ub':
                sort_key = [ subsort[lin][1] for lin in subsort.keys() ]
                lins_to_plot = subsort.keys()
            elif orderby == 'nid':
                lnwnids = self.net_in_degree(normalized=normalized, weighted=weighted, with_labels=True)
                sort_key = [ lnwnids[node] for node in subsort.keys() ]
                lins_to_plot = subsort.keys()
            elif orderby == 'true':
                if trueorder is None:
                    trueorder=self.linstack
                assert(len(trueorder) > 0)
                sort_key = []
                lins_to_plot = []
                for lin in trueorder:
                    if lin in subsort.keys():
                        sort_key.append(trueorder.index(lin))
                        lins_to_plot.append(lin)
            else: # Uh, we should never get here
                print("Bad orderby string found in draw_stratsort")

            # Create a structured array so we can sort by whatever key we want:
            if trueorder is not None:
                dtype = [('lin',object),('sort_key',float),('trueorder',int)]
            else:
                dtype = [('lin',object),('sort_key',float)]
            arr_to_sort = np.zeros(len(lins_to_plot), dtype=dtype)
            arr_to_sort['lin'] = lins_to_plot
            arr_to_sort['sort_key'] = sort_key

            if trueorder is not None:
                sub_trueorder=[]
                for lin in trueorder:
                    if lin in subsort.keys():
                        sub_trueorder.append(lin)

                arr_to_sort['trueorder'] = [ sub_trueorder.index(lin) for lin in lins_to_plot ]

            arr_to_sort.sort(order='sort_key')

            symb_size = 810.0/N_tot
            X = np.arange(len(arr_to_sort))
            lower_bounds = np.array([ subsort[lin][0] for lin in arr_to_sort['lin'] ])-0.5
            upper_bounds = np.array([ subsort[lin][1] for lin in arr_to_sort['lin'] ])+0.5
            ax.vlines(X+x0, lower_bounds, upper_bounds, colors=cmap(float(n_col)/len(stratsorts)), linewidth=symb_size)

            if trueorder is not None:
                ax.plot(X+x0, arr_to_sort['trueorder'], 'k_', markeredgewidth=2, markersize=0.85*symb_size, color='k', lw=0)

            x0 = x0 + len(subsort)

        ax.set_xlim(-1,N_tot)
        ax.set_xlabel('N')
        ax.set_ylim(-1,len(stratsorts[0]))
        ax.set_ylabel('N')
        sortnames = {  'ub':'upper stratigraphic bound',\
                       'lb':'lower stratigraphic bound',\
                     'mean':'mean stratigraphic location',\
                     'true':'a priori formation time',\
                      'nid':'net in degree'}
        ax.set_title("StratSort of "+title+" (ordered by "+sortnames[orderby]+")")
        ax.grid(True)

    #}}}2 end draw_sort

    def draw_graph(self, ax=None, minsize=2, cmap=pylab.cm.jet_r): #{{{2
        """
        Draw graph view of the GSN.  Plot each sub-GSN in a separate figure.
        This is janky and should be replace with something that will plot using
        graphviz instead of matplotlib.

        minsize allows you to limit which connected components are plotted, and
        avoid showing a bunch of uninteresting size 2 and 3 graphs.

        """
        sub_GSNs = self.get_sub_GSNs()
        if minsize > 2:
            sub_GSNs = [ sub_gsn for sub_gsn in sub_GSNs if len(sub_gsn) > minsize ]

        if ax is None:
            for sub_gsn,n in zip(sub_GSNs, range(len(sub_GSNs))):
                the_fig = plt.figure(figsize=(5,5))
                gax = the_fig.add_axes([0,0,1,1])

                nc_array = np.array([ cmap(float(n)/len(sub_GSNs)) for i in range(len(sub_gsn)) ])
                nx.draw_circular(sub_gsn, ax=gax, node_color=nc_array, labels={})

        else:
           nx.draw_circular(self, ax=ax, labels={})

   #}}}2 end draw_graph

    def valid_order(self, linstack): #{{{2
        """
        Given a list of features, ordered by their hypothesized times of formation
        (linstack), and stratigraphic sort (stratsort, as returned from
        GeoSupNet.stratigraphic_sort(), potentially having several members),
        determine whether the ordering of linstack is consistent with the
        constraints encoded within the stratsort.  Return True if it is, and False
        if it is not.

        """

        # Split up linstack according to which member of stratsort each feature
        # within it belongs to (if any).
        sub_GSNs = self.get_sub_GSNs()
        linstacks = self.get_sub_linstacks(linstack)

        # Go through each sub-stack and sub-GSN, and see if the orderings are
        # consistent.  If any of the implied predecessors in the stack are
        # found in a feature's *actual* successors, then the ordering is
        # invalid.  Otherwise it is valid.
        for sub_stack, sub_gsn in zip(linstacks, sub_GSNs):
            for n in range(len(sub_stack)):
                true_successors = nx.single_source_shortest_path(sub_gsn, sub_stack[n]).keys()
                supposed_predecessors = sub_stack[:n]
                for lin in supposed_predecessors:
                    if lin in true_successors:
                        return(False)
        return(True)

#}}}2 end valid_order

    def exclusion_coeff(self, iterations=100, numlins=2, minsize=3): #{{{2
        """
        Determine what fraction of possible feature orderings the GSN can
        exclude.

        iterations is the number of trials to run.

        minsize is the minimum size of the connected components you wish to
        include.  It allows you to only look at the features which actually
        have ordering information (i.e. intersections).  This value must be
        larger than numlins.

        numlins is the number of lineaments to try and validate in each
        ordering.  Those lineaments are chosen at random from the set of all
        features, and re-shuffled on each iteration.

        """

        sub_GSNs = self.get_sub_GSNs(minsize=minsize)

        if np.min([ len(sub_gsn) for sub_gsn in sub_GSNs ]) < minsize:
            assert(minsize > numlins)

        lins = []
        for sub_gsn in sub_GSNs:
            lins += sub_gsn.nodes()

        excluded = 0.0
        for i in range(iterations):
            pylab.shuffle(lins)
            if self.valid_order(lins[:numlins]) is False:
                excluded += 1.0

        return(excluded/iterations)
    #}}}

    def enumerate_cycles(self): #{{{2
        """
        Returns a list of lists of edges (ebunches, in NetworkX terminology),
        with each list representing a simple cycle present in the GSN.

        """

        pass
    #}}}2 end enumerate_cycles

    def make_dag(self): #{{{2
        """
        Uses a heuristic algorithm to find a feedback edge set, and returns
        that ebunch.

        """

        pass
    #}}}2 end enforce_dag

#}}} end GeoSupNet

def lonlat2key(lon, lat): #{{{
    """
    Take a longitude and a latitude, in radians, and turn them into a usable
    key for identifying an intersection.  Returns a 2-tuple of integers.

    """
    lon = np.mod(lon, 2*np.pi)
    lon_key = np.int(1e4*lon)
    lat_key = np.int(1e4*lat)

    return(lon_key, lat_key)
#}}}

def lincross2gsn(lincross): #{{{
    """
    Generates a GSN from a list of lineament intersections.  This is the method
    to use when creating a GSN from a map for actual stratigraphic analysis.

    Because it is possible for two lineaments to have more than one
    intersection, it's necessary to associate the location of each intersection
    with its confidence.

    lincross is a 3-tuple, such that:

       lincross[0] = bottom (older) feature (Lineament object)
       lincross[1] = top (younger) feature (Lineament object)
       lincross[2]['lon'] = lon of intersection (east-positive, radians)
       lincross[2]['lat'] = lat of intersection (radians)
       lincross[2]['weight'] = confidence K of intersection (0 <= K <= 1).

    """

    # Create an empty GSN:
    the_gsn = GeoSupNet()

    # Add each specified intersection, and thereby populate the nodes and edges
    # of the GSN
    for e in lincross:
        the_gsn.add_edge(e[0], e[1], lon=e[2]['lon'], lat=e[2]['lat'],\
                         weight=e[2]['weight'], key=lonlat2key(e[2]['lon'],e[2]['lat']) )

    return(the_gsn)

#}}} end lincross2gsn

def gsn2gsn(old_gsn): #{{{
    """
    Create a new GSN from an old one.

    """
    new_gsn = GeoSupNet()
    new_gsn.add_edges_from(old_gsn.edges(data=True, keys=True))
    new_gsn.linstack = old_gsn.linstack
    return(new_gsn)
#}}}

def linstack2gsn(linstack): #{{{
    """
    Takes a list of Lineament objects, ordered by their times of formation
    (i.e. linstack[0] is the first/oldest feature), and returns a GSN object in
    which all of the intersection confidences are 1.0.

    In order to do this all of the intersections between the features (which
    are implied by their geometry), have to be calculated.  This is done with
    the GDAL/OGR library.  The geometries of the lineaments are converted to
    OGRGeometry objects, and the OGRGeometry.intersection() method is used.

    The coordinate system associated with the lineaments is the IAU 2000
    geographic coordinate system for Europa.  If you're using another body, you
    should find the appropriate definition string here:

    http://spatialreference.org/ref/iau2000/

    """

    # See beginning of file for definition of the WKT SRS
    Europa_SRS = osr.SpatialReference(IAU2000_Europa_SRS_WKT)

    # Initialize the graph
    the_gsn = GeoSupNet()
    the_gsn.linstack = linstack

    # need to reverse the ordering of linstack if we want to keep the ordering as we said...
    backstack = linstack[::-1]
    the_gsn.add_nodes_from(backstack)

    # Because longitude is periodic, but fucking OGR doesn't understand that,
    # we have to do this a few times, over several different ranges of 2*pi:
    # Convert Lineaments to OGR Geometries:
    backstack_ogr_noshift = [ ogr.CreateGeometryFromWkt(lin.wkt(), Europa_SRS) for lin in backstack ]

    for lonshift in (-2.0*np.pi, 0.0, 2.0*np.pi):

        backstack_ogr_shifted = [ ogr.CreateGeometryFromWkt(lin.lonshift(lonshift).wkt(), Europa_SRS) for lin in backstack ]

        for n in range(len(backstack)):
            for m in range(n):
                lin_cross = backstack_ogr_noshift[n].Intersection(backstack_ogr_shifted[m])

                if lin_cross is not None:
                    for point in [ lin_cross.GetPoint(i) for i in range(lin_cross.GetPointCount()) ]:
                        lon = np.mod(point[0],2*np.pi)
                        lat = point[1]

                        the_gsn.add_edge(backstack[n], backstack[m], weight=1.0, key=lonlat2key(lon,lat), lon=lon, lat=lat)

    return(the_gsn)
#}}}

def test_sort(nlins=50, lintype='nsr', orderby='true', ncross=10, draw_map=True, draw_sort=True, draw_graph=True): #{{{
   """
   A short unit test of the sorting functionality.
   
   lintype may be one of: 'nsr', 'regular', or 'random', corresponding to the
   similarly named mapgen routines defined below.

   """

   print("Generating synthetic map")
   if lintype == 'nsr':
       test_lins = mapgen_nsr(nlins=nlins, ncross=ncross)
   elif lintype == 'regular':
       test_lins = mapgen_regular(nlins=nlins, ncross=ncross)
   else:
       assert(lintype == 'random')
       test_lins = nsrhist.random_gclins(nlins)
       for lin in test_lins:
           lin.lons = lineament.fixlons(lin.lons)

   print("Converting map to GSN")
   test_gsn = linstack2gsn(test_lins)

   print("Sorting GSN")
   test_sort = test_gsn.stratigraphic_sort()
   print("Found %d connected components" % (len(test_sort),) )

   print("Plotting results")
   if draw_sort is True:
       test_gsn.draw_sort(trueorder=test_lins, orderby=orderby, label=lintype)
   if draw_map is True:
       test_gsn.draw_map()
   if draw_graph is True:
       test_gsn.draw_graph()

   return(test_lins, test_gsn, test_sort)

#}}}

def mapgen_regular(nlins=10, ncross=None): #{{{
    """
    Create a regular map with a known order of formation for testing the GSN
    algorithms. nlins indicates how many features there are in each of the
    north-south and east-west orientations.  ncross indicates where in the
    sequence of east-west features the north-south features should appear.  If
    ncross is none, the whole ordering is randomized.

    """

    synth_lins = lineament.lingen_nsr_library(nlats=nlins)[:nlins]
    cross_lins = []

    init_lons = np.linspace(np.radians(-30), np.radians(30), nlins)
    fin_lons  = np.linspace(np.radians(-90), np.radians(90), nlins)

    for init_lon, fin_lon in zip(init_lons, fin_lons):
        cross_lins.append(lineament.lingen_greatcircle(init_lon,0,fin_lon,np.radians(89.0)))

    if ncross is None:
        outlins = synth_lins + cross_lins
        pylab.shuffle(outlins)
    else:
        outlins = synth_lins[:ncross] + cross_lins + synth_lins[ncross:]

    return(outlins)

#}}} end mapgen_regular

def mapgen_fromfile(nlins=100, linfile = 'output/lins/map_nsrfit', spin=False, tpw=False): #{{{
    """
    Draw nlins lineaments from a saved pool of features, and create a linstack
    from them.

    """

    # read in the features, and update them to make sure we've got whatever
    # interesting data they have, in the new Lineament object form.
    lins = lineament.update_lins(nsrhist.load_lins(linfile))

    if nlins == 0:
        nlins = len(lins)

    # If we're not randomizing, we can't resample with replacement, or we'll
    # get repeats that overlap themselves.
    if spin is False and tpw is False:
        # choose nlins random features from the dataset:
        # this will fail if they've asked for more features than we can get without randomizing.
        assert(nlins <= len(lins))
        pylab.shuffle(lins)
        linstack = lins[:nlins]
    else:
        linstack = nsrhist.make_crazy(nsrhist.linresample_byN(lins), tpw=tpw, spin=spin)

    return(linstack)
        
#}}}

def mapgen_nsr(nlins=50, ncross=5): #{{{
    """
    Create a regular map with a known order of formation for testing the GSN
    algorithms. nlins indicates how many features there are in each of the
    north-south and east-west orientations.  ncross indicates where in the
    sequence of east-west features the north-south features should appear.  If
    ncross is none, the whole ordering is randomized.

    """

    # read in what we need to calculate the NSR stresses:
    satfile = "input/ConvectingEuropa_GlobalDiurnalNSR.ssrun"
    europa = satstress.Satellite(open(satfile,'r'))
    nsr_stresscalc = satstress.StressCalc([satstress.NSR(europa),])

    # create some random NSR features, in the northern hemisphere

    nsr_linstack = []
    while len(nsr_linstack) < nlins:
        lon, lat = lineament.random_lonlatpoints(1)
        lon = np.mod(lon[0], np.pi)-np.pi/2.0
        lat = np.abs(lat[0])
        length = np.abs(np.random.normal())
        newlin = lineament.lingen_nsr(nsr_stresscalc, init_lon=lon, init_lat=lat, max_length=length)
        nsr_linstack.append(newlin)

    for xlon in np.linspace(-np.pi/3.0, np.pi/3.0, ncross):
        nsr_linstack.append(lineament.lingen_greatcircle(xlon,0,xlon+np.pi,np.pi/2.1))

    pylab.shuffle(nsr_linstack)

    bs = np.sort(np.random.uniform(np.pi, size=nlins+ncross))[::-1]-np.pi

    for lin, b in zip(nsr_linstack, bs):
        lin.lons = lin.lons - b

    return(nsr_linstack)

#}}} end mapgen_nsr

################################################################################
################################################################################
# ============================================================================== 
# Reading/Writing:
# ============================================================================== 
#   - Go back through the Powerpoint presentation
#   - Write up explanation of preliminary results
#   - Write up improved algorithm description

# ============================================================================== 
# GSN v1.0 (DAG only):
# ==============================================================================

# SUCCESS METRICS:
#   Given a GSN and a true order of formation, calculate the proportion
#   of the original ordering information that was recovered, i.e. the change in
#   entropy, or "entropy reduction".  Need to consult with Lada and/or Aaron on
#   this.
#
#   Exclusion Coefficient: the proportion of orderings that can be excluded
#   based on the information in the GSN.  Varies as a function of how many
#   elements there are in the ordering.

# FAILURE METRIC:
#   Given a StratSort and a true order of formation, extend the above metric
#   (proportion of ordering information recovered) to include a metric of how
#   wrong the recovered ordering is, if it is inconsistent with the true order.

# ============================================================================== 
# GSN v1.1 (non-DAG): #{{{
# ============================================================================== 

# INDUCE CYCLES ON A TEST DAG:
#   Given an acyclic GSN, generate a similar GSN, that has experienced
#   re-activation by reversing the direction of some proportion of the edges.

# ENUMERATE CYCLES:
#   Given a GSN containing cycles, return a set of "path" graphs describing all
#   of the simple cycles in the original GSN.  Give the option of removing all
#   the intersections below a particular confidence threshold first in order to
#   speed things up, and focus on the most interesting features.

# RENDER GSN ACYCLIC:
#   Given a GSN containing cycles generate an acyclic GSN, retaining as much of
#   the intersection information as is practical.  Allow the exclusion of
#   intersections with confidences below a particular threshold.  Allow setting
#   a limit on the depth of recursive calls to make.
#   
#   Three basic methods for getting rid of back-edges in pseudo-temporal sort:
#     - Remove them one-by-one until the graph is a DAG
#     - Reverse their direction if they are heavily clustered
#     - Split the feature they lie on and re-sort.

# DISPLAY GRAPH SORTING:
#   Get pyGraphViz installed and working
#
#   Given a GSN, create a figure showing the pseudo-temporal ordering of the
#   lineaments based on their net weighted (in|out)-degree.  Allow weighting of 
#   net degree by lineament length.
#   
#   Show the nodes stacked vertically, with the oldest features at the bottom,
#   and the newest ones at the top.  Put all of the forward edges on one side
#   of the stack, and all of the back edges on the other, as curved arcs.
#
#   Allow the creation of series of these stacks from left to right, showing
#   the sorting and/or DAG enforcement process.

#   When displaying the graph after sorting, show:
#     - width of edges keyed to confidence
#     - edge is solid if retained
#     - edge is dashed if removed

# IDENTIFY REACTIVATION SITES:
#   Given a GSN containing cycles, return a best guess as to which lineament
#   underwent reactivation, and which of its intersections were reversed as a
#   result (those intersections which are back-edges in the pseudo-temporal
#   sorting).
#   
#   Allow the reversal of those intersections (instead of just removing them)
#   and see if that improves the sorting results.

# LINEAMENT SPLITTING:
#   A further refinement of the above: instead of either removing or reversing
#   the back-edges, split the lineament to separate the back-edge section from
#   the rest of the feature, re-constitute the graph, and perform the sorting
#   algorithm again.

#}}}

# given an actual mapped geometry, generate a distribution of information
# retrieved, as a results of many random stratigraphic orderings.

# Create a circular or otherwise geometrically unbiased map space, and show how
# information retrieved varies with lineament density per unit area.

# Can also perform this analysis as a function of lineament length, or length
# distribution.

# This could be a 2D plot: length on one axis, density on the other, with color
# showing how much information was retrieved (sweet plot)

# Try and pitch this as why we need global high-res coverage, if it's not
# useful right now.
