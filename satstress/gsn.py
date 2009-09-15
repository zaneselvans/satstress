#!/usr/bin/python
import networkx as nx
from . import lineament
from . import satstress
import numpy as np
import pylab
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
    
    def stratigraphic_sort(self): #{{{2
        """
        Find upper and lower bounds on the stratigraphic location of each
        feature in the GSN.

        This is done by finding the sizes of the successor and predecessor
        digraphs.  Anything that isn't in one of those graphs is of
        indistinguishable age.

        """

        # Perform a separate stratsort on each one
        # Return a list of stratsorts

        # Convert the GSN to an undirected graph
        # and find the connected subcomponents
        connected_nbunches = nx.connected_components(nx.MultiGraph(self))

        # re-constitute those connected bunches as a new list of GSNs:
        sub_GSNs = [ self.subgraph(nbunch) for nbunch in connected_nbunches ]

        stratsorts = []
        for sub_gsn in sub_GSNs:
            stratsort = {}
            for lin in sub_gsn.nodes():
                # This seems to take waaaay too long for some reason
                ub = len(sub_gsn) - len(nx.dfs_tree(sub_gsn, source=lin))
                lb = len(nx.dfs_tree(sub_gsn, source=lin, reverse_graph=True))-1
                # ub = len(self) - len(self.all_successors(lin))
                # lb = len(self.all_predecessors(lin)) - 1
                stratsort[lin] = (lb, ub)
            stratsorts.append(stratsort)
        
        return(stratsorts)

    #}}}2 end stratigraphic_sort

    def make_dag(self): #{{{2
        """
        Uses a heuristic algorithm to find a feedback edge set, and returns
        that ebunch.

        """

        pass
    #}}}2 end enforce_dag

    def enumerate_cycles(self): #{{{2
        """
        Returns a list of lists of edges (ebunches, in NetworkX terminology),
        with each list representing a simple cycle present in the GSN.

        """

        pass
    #}}}2 end enumerate_cycles

    def draw_map(self): #{{{2
        """
        Draw a map view of the GSN using Basemap and Matplotlib.

        """
        # TODO: need to make this color code by connected components

        gsn_edges = self.edges(data=True)
        edge_lons = [ edge[2]['lon'] for edge in gsn_edges ]
        edge_lats = [ edge[2]['lat'] for edge in gsn_edges ]

        lines, gsn_map = lineament.plotlinmap(self.nodes())
        gsn_map.scatter(np.degrees(edge_lons), np.degrees(edge_lats), lw=0, color='r', alpha=0.75, s=30)
        gsn_map.scatter(np.degrees(edge_lons)+360.0, np.degrees(edge_lats), lw=0, color='r', alpha=0.75, s=30)
        gsn_map.scatter(np.degrees(edge_lons)-360.0, np.degrees(edge_lats), lw=0, color='r', alpha=0.75, s=30)
    #}}}2 end draw_map

    def draw_graph(self): #{{{2
       """
       Draw a graph view of the GSN using Matplotlib

       """
       the_fig = pylab.figure()
       ax = the_fig.add_subplot(1,1,1)

   #}}}2 end draw_graph

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
    
    lincross is a list of 5-tuples, such that:

        lincross[0] = bottom (older) feature (Lineament object)
        lincross[1] = top (younger) feature (Lineament object)
        lincross[2] = longitude of the intersection (east-positive, radians)
        lincross[3] = latitude of the intersection (radians)
        lincross[4] = confidence K of that intersection (0 <= K <= 1).

    """

    # Create an empty GSN:
    the_gsn = GeoSupNet()

    # Add each specified intersection, and thereby populate the nodes and edges
    # of the GSN
    for x in lincross:
        the_gsn.add_edge(x[0], x[1], lon=x[2], lat=x[3], weight=x[4], key=lonlat2key(x[2],x[3]))

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

def test_sort(nlins=50, lintype='nsr', ncross=10): #{{{
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
       test_lins = mapgen_random(nlins=nlins)

   print("Converting map to GSN")
   test_gsn  = linstack2gsn(test_lins)
   print("Sorting GSN")
   test_sort = test_gsn.stratigraphic_sort()
   print("Plotting results")
   draw_stratsort(test_sort, trueorder=test_lins)
   test_gsn.draw_map()

   return(test_lins, test_gsn, test_sort)

#}}}

def draw_stratsort(stratsorts, trueorder=None): #{{{
    """
    Make a plot showing the retrieved stratigraphic sort.

    The true ordering is indicated by the order of the features in the list
    trueorder.
    
    Eventually this needs to deal with things with no known ordering, by adding
    an 'orderby' parameter, which may also have the following values:

        'mean' -- mean stratigraphic location (default)
          'lb' -- lower bound (earliest possible formation time)
          'ub' -- upper bound (latest possible formation time)

    """

    the_fig = pylab.figure()
    ax = the_fig.add_subplot(1,1,1)

    N_tot = np.array([len(subsort) for subsort in stratsorts]).sum()

    x0 = 0
    for subsort in stratsorts:
        nlins = len(subsort)
        N = np.arange(nlins)

        suborder = [ lin for lin in trueorder if lin in subsort.keys() ]

        yerr_lo = [ n-subsort[lin][0] for lin,n in zip(suborder,N) ]
        yerr_hi = [ subsort[lin][1]-n for lin,n in zip(suborder,N) ]

        ax.errorbar(N+x0, N, (yerr_lo, yerr_hi), fmt=None, capsize=0, elinewidth=(N_tot/10.0), ecolor='k')
        ax.plot(N+x0, N, 'r_', markersize=4+(N_tot/10.0), markeredgewidth=3, color='r', lw=0)

        x0 = x0 + nlins

    ax.set_xlim(-1,N_tot)
    ax.set_xlabel('N')
    ax.set_ylim(-1,len(stratsorts[0]))
    ax.set_ylabel('N')
    ax.set_title("Stratigraphic Sort Test")
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
        nsr_linstack.append(lineament.lingen_greatcircle(xlon,0,xlon,np.pi/2.01))

    pylab.shuffle(nsr_linstack)

    bs = np.sort(np.random.uniform(np.pi/2.0, size=nlins+ncross))[::-1]-np.pi/2.0

    for lin, b in zip(nsr_linstack, bs):
        lin.lons = lin.lons - b

    return(nsr_linstack)

#}}} end mapgen_nsr

def mapgen_random(nlins): #{{{
    """
    Generate a random set of great circle segments for practice sorting.
    Should probably try and keep this within a particular geographic area just
    for simplicty.

    """
    
    pass
#}}} end mapgen_random

################################################################################
# GSN TODO:
################################################################################
#
# * Go back through the Powerpoint presentation
# * Write up explanation of preliminary results
# * Write up improved algorithm description
#
# * Get GraphViz/PyGraphViz working
# 
# * Find out how to implement a new layout algorithm
#   - linear arrangement of nodes
#     - color-code nodes by length-weighted net degree
#   - forward edges arcing over the top
#   - back-edges arcing under the bottom
#     - width (or alpha?) of edges by confidence
#     - solid if kept
#     - dotted if removed
#
# ============================================================================== 
# GSN v1.0 (DAG only):
# ==============================================================================

# SPEED UP SORTING!

# Given:
#   - A StratSort and
#   - A hypothesized ordering
# Return:
#   - True or False based on whether the ordering is consistent with the sort

# Given:
#   - A GSN and
#   - a StratSort
# Create:
#   - A map showing each feature color coded by graph component
#   - A map showing each feature color coded by average stratigraphic location.

# Given:
#   - A StratSort
# Calculate:
#   - The proportion of the ordering information that was recovered

# ============================================================================== 
# GSN v1.1 (non-DAG):
# ============================================================================== 

# Given:
#   - a GSN
# Generate:
#   - an ordered list of the nodes/edges making up each simple cycle
#     - if there are lots of cycles, may take too long
#     - allow filtering by minimum acceptable confidence 

# Given:
#   - an acyclic GSN
# Generate:
#   - a similar GSN, that has experienced re-activation
#     - reverse the z-order pair for a portion of some features

# Given:
#   - a GSN containing cycles
# Generate:
#   - an acyclic GSN
#   - retain as much information as is practical
#   - allow filtering of intersections by confidence
#   - set number of recursive calls to allow

# Given:
#   - a GSN containing cycles
# Generate:
#   - a guess as to which lineament segment reactivated
#     - consists of a set of intersections to reverse

# Given:
#   - a GSN
# Generate:
#   - A figure showing lineaments sorted by net weighted degree
#     with forward and back edges on different sides of the nodes
#     - Show the nodes stacked vertically
#     - Put oldest nodes at the bottom (high in-degree)
#     - Will allow left-to-right time series
#       - good for showing the sorting process
#       - good for showing the DAG enforcement

# Given:
#   - a StratSort
#   - a true order of formation
# Generate:
#   - a metric of how wrong the recovered ordering is, if at all.

