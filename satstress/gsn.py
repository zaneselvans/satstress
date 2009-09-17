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

    def get_sub_GSNs(self): #{{{2
        """
        Return a list of all the disconnected components of the GSN.

        """

        # Convert the GSN to an undirected graph
        # and find the connected subcomponents
        U = nx.MultiGraph()
        U.add_edges_from(self.edges())
        connected_nbunches = nx.connected_components(U)
        # re-constitute those connected bunches as a new list of GSNs:
        return([ self.subgraph(nbunch) for nbunch in connected_nbunches ])

    #}}}2

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

    def draw_map(self, ax=None, cmap=pylab.cm.jet_r): #{{{2
        """
        Draw a map view of the GSN using Basemap and Matplotlib.

        """
        # TODO: allow coloring by things other than connected component...

        gsn_edges = self.edges(data=True)
        edge_lons = [ edge[2]['lon'] for edge in gsn_edges ]
        edge_lats = [ edge[2]['lat'] for edge in gsn_edges ]

        if ax is None:
            the_fig = plt.figure(figsize=(12,6))
            ax = the_fig.add_subplot(1,1,1)
            map_ax = basemap.Basemap(ax=ax)
            map_ax.drawmeridians(range(-180,181,30), labels=[1,0,0,1])
            map_ax.drawparallels(range(-90,91,30), labels=[1,0,0,1])

        sub_GSNs = self.get_sub_GSNs()
        for sub_gsn,n in zip(sub_GSNs, range(len(sub_GSNs))):
            lines, gsn_map = lineament.plotlinmap(sub_gsn.nodes(), map=map_ax, linewidth=2,\
                                                  color=cmap(float(n)/len(sub_GSNs)))

        gsn_map.scatter(np.degrees(edge_lons), np.degrees(edge_lats), lw=0, color='black', s=50)
        gsn_map.scatter(np.degrees(edge_lons)+360.0, np.degrees(edge_lats), lw=0, color='black', s=50)
        gsn_map.scatter(np.degrees(edge_lons)-360.0, np.degrees(edge_lats), lw=0, color='black', s=50)
    #}}}2 end draw_map

    def draw_sort(self, orderby='mean', trueorder=None, ax=None, colorby='graph', cmap=pylab.cm.jet_r, label=""): #{{{2
       """
       Draw a graph view of the GSN using Matplotlib

       """
       # TODO: Allow coloring by things other than connected components.
       # TODO: Allow ordering in ways other than a true order of formation.

       if ax is None:
           the_fig = plt.figure(figsize=(12,6))
           ax = the_fig.add_subplot(1,1,1)

       draw_stratsort(self.stratigraphic_sort(), orderby=orderby, trueorder=trueorder, ax=ax, colorby=colorby, label=label, cmap=cmap)

   #}}}2 end draw_graph

    def draw_graph(self, ax=None, cmap=pylab.cm.jet_r): #{{{2
        """
        Draw a graph view of the GSN using Matplotlib

        """
        sub_GSNs = self.get_sub_GSNs()

        if ax is None:
            for sub_gsn,n in zip(sub_GSNs, range(len(sub_GSNs))):
                the_fig = plt.figure(figsize=(5,5))
                gax = the_fig.add_axes([0,0,1,1])

                nc_array = np.array([ cmap(float(n)/len(sub_GSNs)) for i in range(len(sub_gsn)) ])
                nx.draw_circular(sub_gsn, ax=gax, node_color=nc_array, labels={})
            #nx.draw_circular(sub_gsn, ax=graph_axes[n], node_color=cmap(float(n)/len(sub_GSNs)), labels={})
            #nxs = np.ceil(np.sqrt(len(sub_GSNs)))
            #graph_axes = [ the_fig.add_subplot(nxs,nxs,n+1) for n in range(len(sub_GSNs)) ]

        else:
           nx.draw_circular(self, ax=ax, labels={})

   #}}}2 end draw_graph

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

def validate_ordering(linstack, stratsort): #{{{
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
    linstacks = []
    for sub_sort in stratsort:
        sub_stack = [ lin for lin in linstack if lin in sub_sort.keys() ]
        linstacks.append(sub_stack)

    # Go through each sub-stack and sub-sort, and see if the orderings are
    # consistent.
    for sub_stack, sub_sort in zip(linstacks, stratsort):
        for lin, n in zip(sub_stack, range(len(sub_stack))):
            pass # Still not entirely clear what to do here....

#}}} end validate_ordering

def test_sort(nlins=50, lintype='nsr', ncross=10, draw_map=True, draw_sort=True, draw_graph=True): #{{{
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

   print("Plotting results")
   if draw_sort is True:
       test_gsn.draw_sort(trueorder=test_lins, orderby='true', label=lintype)
   if draw_map is True:
       test_gsn.draw_map()
   if draw_graph is True:
       test_gsn.draw_graph()

   return(test_lins, test_gsn, test_sort)

#}}}

def draw_stratsort(stratsorts, orderby='mean', trueorder=None, ax=None, colorby='graph', cmap=pylab.cm.jet_r, label=""): #{{{
    """
    Make a plot showing the retrieved stratigraphic sort.

    The true ordering is indicated by the order of the features in the list
    trueorder.
    
    Eventually this needs to deal with things with no known ordering, by adding
    an 'orderby' parameter, which may also have the following values:

        'true' -- the actual ordering, as implied by trueorder
        'mean' -- mean stratigraphic location (default)
          'lb' -- lower bound (earliest possible formation time)
          'ub' -- upper bound (latest possible formation time)

    """

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
        elif orderby == 'true':
            assert(trueorder is not None)
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

        yerr_lo = np.array([ n-subsort[lin][0] for lin,n in zip(arr_to_sort['lin'],range(len(arr_to_sort))) ])+0.5
        yerr_hi = np.array([ subsort[lin][1]-n for lin,n in zip(arr_to_sort['lin'],range(len(arr_to_sort))) ])+0.5

        X = np.arange(len(arr_to_sort))
        ax.errorbar(X+x0, X, (yerr_lo, yerr_hi), fmt=None, capsize=0, elinewidth=(500.0/N_tot),\
                    ecolor=cmap(float(n_col)/len(stratsorts)))

        if trueorder is not None:
            ax.plot(X+x0, arr_to_sort['trueorder'], 'k_', markersize=1.2*(500.0/N_tot), markeredgewidth=3, color='k', lw=0)

        x0 = x0 + len(subsort)

    ax.set_xlim(-1,N_tot)
    ax.set_xlabel('N')
    ax.set_ylim(-1,len(stratsorts[0]))
    ax.set_ylabel('N')
    sortnames = {  'ub':'upper stratigraphic bound',\
                   'lb':'lower stratigraphic bound',\
                 'mean':'mean stratigraphic location',\
                 'true':'a priori formation time' }
    ax.set_title("StratSort of "+label+" (ordered by "+sortnames[orderby]+")")
    ax.grid(True)

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
# VISUALIZATION:
#   Allow StratSort to be ordered by any of true order, lower bound, upper
#   bound, or average stratigraphic location ((ub+lb)/2).  Always display the
#   true order of formation using the black bars.

# SORT VALIDITY:
#   Given a StratSort and a hypothesized ordering, return True or False based
#   on whether the ordering is consistent with or conflicts with the results of
#   the sort

# SUCCESS METRIC:
#   Given a StratSort and a true order of formation, calculate the proportion
#   of the original ordering information that was recovered.  Need to consult
#   with Lada and/or Aaron on this.  (Conditional entropy, mutual information).
#   
#   One possibility is to submit a large number of completely random orderings
#   for validation.  The proportion of them that are rejected would indicate
#   how much of the potential ordering space has been excluded.

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

# FAILURE METRIC:
#   Given a StratSort and a true order of formation, extend the above metric
#   (proportion of ordering information recovered) to include a metric of how
#   wrong the recovered ordering is, if it is inconsistent with the true order.

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
