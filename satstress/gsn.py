#!/usr/bin/python

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

# Given:
#   - an acyclic GSN
# Generate:
#   - a StratSort of the GSN, i.e. for each lineament L:
#     - a predecessor tree (features younger than L)
#     - a successor tree (features older than L)
#     - a cohort list (features indistinguishable in age from L)
#     - nearest neighbor(s) down in the stack
#     - may have more than one of each, if there are indistinguishable
#       nodes in the graph (topologically)
#   - when removing a set of source (or sink) nodes, the next set of source
#     (or sink) nodes is their potential nearest neighbors (actual nearest
#     neighbors if they had an intersection)

# Given:
#   - a StratSort
# Generate:
#   - metric of how much ordering information was recovered
#   - a figure showing the history schematically, including the ambiguity
#   - a figure showing the history, in map form... hard.


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
#   - a GSN
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

# ============================================================================== 
# GSN Data Structures/Objects:
# ============================================================================== 
#   - StratSort (the result):
#     - for each input lineament, the nearest stratigraphic neighbors (also
#       lineaments) in the stack, both upward and downward.
#     - will need a couple of ways to visualize this:
#       - stratigraphic column style (like the schematics I made by hand)
#       - map style (like the output I made in ArcGIS before)
#
################################################################################

import sys
import string
import re

import networkx as nx
from . import lineament
import numpy as np
from osgeo import ogr
from osgeo import osr

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
        Performs a modified topological sort on the GSN, and for each node
        calculates:

          * a predecessor tree (all older features)
          * a successor tree (all younger features)
          * a peer cohort (all features indistinguishable in age)

        """

        pass
    #}}}2 end stratigraphic_sort

    def enforce_dag(self, min_conf=0.0): #{{{2
        """
        Removes any cycles in the GSN, rendering it sortable.  Only includes
        intersections which have confidences greater than or equal to min_conf.

        """

        pass
    #}}}2 end enforce_dag

    def enumerate_cycles(self, min_conf=0.0): #{{{2
        """
        Returns a list of all the simple cycles present in the GSN composed only
        of intersections having confidences greater than or equal to min_conf.

        """

        pass
    #}}}2 end enumerate_cycles

    def draw_map(self): #{{{2
        """
        Draw a map view of the GSN using Basemap and Matplotlib.

        """

        gsn_edges = self.edges(data=True)
        edge_lons = [ edge[2]['lon'] for edge in gsn_edges ]
        edge_lats = [ edge[2]['lat'] for edge in gsn_edges ]

        lines, gsn_map = lineament.plotlinmap(self.nodes())
        gsn_map.scatter(np.degrees(edge_lons), np.degrees(edge_lats))
        gsn_map.scatter(np.degrees(edge_lons)+360.0, np.degrees(edge_lats))
        gsn_map.scatter(np.degrees(edge_lons)-360.0, np.degrees(edge_lats))
    #}}}2 end draw_map

    def draw_graph(self): #{{{2
       """
       Draw a graph view of the GSN using Matplotlib

       """

       pass
   #}}}2 end draw_graph

#}}} end GeoSupNet

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

    # define the coordinate system we're working in:
    # This is lame, but for the moment, let's assume that we're on Europa:
    IAU2000_Europa_SRS_WKT="""GEOGCS["Europa 2000",DATUM["D_Europa_2000",SPHEROID["Europa_2000_IAU_IAG",1564130.0,488.79062499999998]],PRIMEM["Greenwich",0],UNIT["Decimal_Degree",0.0174532925199433]]"""
    Europa_SRS = osr.SpatialReference(IAU2000_Europa_SRS_WKT)

    # Initialize the graph
    the_gsn = GeoSupNet()
    the_gsn.add_nodes_from(linstack)

    # Because longitude is periodic, but fucking OGR doesn't understand that,
    # we have to do this a few times, over several different ranges of 2*pi:
    # Convert Lineaments to OGR Geometries:
    linstack_ogr_noshift = [ ogr.CreateGeometryFromWkt(str(lin), Europa_SRS) for lin in linstack ]

    for lonshift in (-2.0*np.pi, 0.0, 2.0*np.pi):

        linstack_ogr_shifted = [ ogr.CreateGeometryFromWkt(str(lin.lonshift(lonshift)), Europa_SRS) for lin in linstack ]

        for n in range(len(linstack)):
            for m in range(n):
                lin_cross = linstack_ogr_noshift[n].Intersection(linstack_ogr_shifted[m])

                if not lin_cross.IsEmpty():
                    for point in [ lin_cross.GetPoint(i) for i in range(lin_cross.GetPointCount()) ]:
                        lon = np.mod(np.radians(point[0]),2*np.pi)
                        lon_key = np.int(1e4*lon)

                        lat = np.radians(point[1])
                        lat_key = np.int(1e4*lat)

                        the_gsn.add_edge(linstack[n], linstack[m], weight=1.0, key=(lon_key,lat_key), lon=lon, lat=lat)

    return(the_gsn)
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
        shuffle(outlins)
    else:
        outlins = synth_lins[:ncross] + cross_lins + synth_lins[ncross:]

    return(outlins)

#}}} end mapgen_regular

def mapgen_random(nlins): #{{{
    """
    Generate a random set of great circle segments for practice sorting.
    Should probably try and keep this within a particular geographic area just
    for simplicty.

    """
    
    pass
#}}} end mapgen_random
