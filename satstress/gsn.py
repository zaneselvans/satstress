#!/usr/bin/python
"""
A Geological Superposition Network (GSN) is a data structure encoding both
geographic and stratigraphic information about features on the surface of a
planetary body.  It uses directed graphs (or "digraphs"), a construct from
discrete mathematics, to manipulate both the geography and superposition
relationships by storing each geologic feature as a "node" (or "vertex") in the
graph, and each cross-cutting relationship between features as an "edge", such
that if a feature A is seen to cross-cut another feature B, there exists a
directed edge in the graph from B to A, and in the language of graph theory, A
is said to be a successor to B.  The direction of the edge thus indicates the
arrow of time, with the implication that because A appears superimposed upon B,
B must be older than A.

Within each node of the GSN is stored the geographic information pertaining to
the geologic feature it represents.  In the applications below, the features
being analyzed are linear tectonic features, and so satstress Lineament objects
are used as nodes.

It is possible for two features to have more than one intersection, and so it
follows that in some instances there may be multiple edges between the same two
nodes.  Because of this, it is necessary for each edge to have additional data
associated with it to allow multiple edges between the same two nodes to be
distinguished from one another.  It is also useful to have the location of each
intersection stored, for display and processing purposes.  Thus, each edge also
has a latitude and longitude value.

Because in some instances the exact nature of the superposition relationship is
unclear, the ordering of each intersection is assigned an estimated probability
of its being correct.  When this value is 1, the ordering is completely clear.
When it is 0.5, no temporal information is discernible.

Because the constraints that the cross cutting relationships place on the
overall order of formation are in general fairly loose, it is usually not
practical to use a GSN to try and infer a particular chronology.  However, the
construct is still useful in answering certain kinds of questions, such as:

"What are the chances that one feature formed before or after another?"

which is really just a special case of:

"To what degree are the mapped stratigraphic relationships compatible with a
particular ordering inferred by some other method?"

A GSN may also be used to identify sets of stratigraphic relationships which
appear to be logically incompatible, and therefore indicate non-instantaneous
activity.  Such sets of intersections form loops, or "cycles" within the
directed graph.

"""

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import pylab
from mpl_toolkits import basemap
from osgeo import ogr
from osgeo import osr
from . import lineament
from . import satstress
from . import nsrhist

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

    ##########################################################################
    # Translation:
    ##########################################################################

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
          lincross[2]['weight'] = -log(conf) where 0 <= conf <= 0.5

        """

        return self.edges(data=True)
    #}}}2

    def pairs(self): #{{{2
        """
        Return a list of sets of all the pairwise stratigraphic relationships,
        as ordered tuples of lineament objects, which are determined by each
        sub graph of the GSN

        """

        pairwise_sets = []
        for sub_gsn in self.get_sub_GSNs():
            pairwise_subset = set([])
            for source in sub_gsn.nodes():
                for dest in nx.single_source_shortest_path(sub_gsn, source).keys():
                    if source != dest:
                        pairwise_subset.add((source, dest))
            pairwise_sets.append(pairwise_subset)

        return(pairwise_sets)

    #}}}2 end pairwise_orderings

    ##########################################################################
    # Display Routines:
    ##########################################################################

    def draw_map(self, intersections=False, ax=None, cmap=pylab.cm.jet_r, minsize=2, minconf=0.5): #{{{2
        """
        Draw a map view of the GSN using Basemap and Matplotlib.

        """
        if ax is None:
            the_fig = plt.figure(figsize=(12,6))
            ax = the_fig.add_subplot(1,1,1)
            map_ax = basemap.Basemap(ax=ax)
            map_ax.drawmeridians(range(-180,181,30), labels=[1,0,0,1])
            map_ax.drawparallels(range(-90,91,30), labels=[1,0,0,1])

        sub_GSNs = self.get_sub_GSNs(minsize=minsize, minconf=minconf)

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

    def draw_sort(self, orderby='mean', ax=None, colorby='graph', cmap=pylab.cm.jet_r, title="", weighted=True, normalized=True, minsize=2): #{{{2
        """
        Make a plot showing the retrieved stratigraphic sort.
        
        The features may be ordered by any of the following:
    
             'nid' -- net in degree (default)
            'true' -- the actual ordering, as implied by trueorder
            'mean' -- mean stratigraphic location
              'lb' -- lower bound (earliest possible formation time)
              'ub' -- upper bound (latest possible formation time)
    
        If the true ordering is indicated by the linstacks embedded with the
        GSN, then it will be displayed.

        If any of weighted or normalized are set, they are passed on to the
        GeoSupNet.net_in_degree() function (with orderby='nid')

        minsize determines how many features have to be the connected subgraph
        in order for it to be displayed.  This allows you to avoid showing a
        bunch of sorts that only include 2 or 3 features.

        """
        if ax is None:
            the_fig = plt.figure(figsize=(12,6))
            ax = the_fig.add_axes((0.05,0.1,0.9,0.8))

        sub_GSNs = self.get_sub_GSNs(minsize=minsize, minconf=minconf)
        N_tot = np.sum([len(sub_gsn) for sub_gsn in sub_GSNs ])

        x0 = 0
        for sub_gsn, n_col in zip(sub_GSNs, range(len(sub_GSNs))):
            sub_sort  = sub_gsn.stratigraphic_sort()
            assert (len(sub_sort) == 1)
            sub_sort = sub_sort[0] # should have only one connected component...
            sub_stack = sub_gsn.linstack

            # Figure out what order we're going to plot the features in, and
            # generate the appropriate x, y, and (yerr_lo, yerr_hi) arrays for
            # plotting:
            if orderby == 'mean':
                sort_key = [ 0.5*(sub_sort[lin][0]+sub_sort[lin][1]) for lin in sub_sort.keys() ]
                lins_to_plot = sub_sort.keys()
            elif orderby == 'lb':
                sort_key = [ sub_sort[lin][0] for lin in sub_sort.keys() ]
                lins_to_plot = sub_sort.keys()
            elif orderby == 'ub':
                sort_key = [ sub_sort[lin][1] for lin in sub_sort.keys() ]
                lins_to_plot = sub_sort.keys()
            elif orderby == 'nid':
                nids = self.net_in_degree(normalized=normalized, weighted=weighted, with_labels=True)
                sort_key = [ nids[node] for node in sub_sort.keys() ]
                lins_to_plot = sub_sort.keys()
            elif orderby == 'true':
                assert(len(sub_stack) > 0)
                sort_key = []
                lins_to_plot = []
                for lin in sub_stack:
                    if lin in sub_sort.keys():
                        sort_key.append(sub_stack.index(lin))
                        lins_to_plot.append(lin)
            else: # Uh, we should never get here
                print("Bad orderby string found in draw_stratsort")

            # Create a structured array so we can sort by whatever key we want:
            if len(sub_stack) > 0:
                dtype = [('lin',object),('sort_key',float),('trueorder',int)]
            else:
                dtype = [('lin',object),('sort_key',float)]

            arr_to_sort = np.zeros(len(lins_to_plot), dtype=dtype)
            arr_to_sort['lin'] = lins_to_plot
            arr_to_sort['sort_key'] = sort_key
            if len(sub_stack) > 0:
                arr_to_sort['trueorder'] = [ sub_stack.index(lin) for lin in lins_to_plot ]
            arr_to_sort.sort(order='sort_key')

            symb_size = 810.0/N_tot
            X = np.arange(len(arr_to_sort))
            lower_bounds = np.array([ sub_sort[lin][0] for lin in arr_to_sort['lin'] ])-0.5
            upper_bounds = np.array([ sub_sort[lin][1] for lin in arr_to_sort['lin'] ])+0.5
            ax.vlines(X+x0, lower_bounds, upper_bounds, colors=cmap(float(n_col)/len(sub_GSNs)), linewidth=symb_size)

            if len(sub_stack) > 0:
                ax.plot(X+x0, arr_to_sort['trueorder'], 'k_', markeredgewidth=2, markersize=0.85*symb_size, color='k', lw=0)

            x0 = x0 + len(sub_gsn)

        ax.set_xlim(-1,N_tot)
        ax.set_xlabel('N')
        ax.set_ylim(-1,len(sub_GSNs[0]))
        ax.set_ylabel('N')
        sortnames = {  'ub':'upper stratigraphic bound',\
                       'lb':'lower stratigraphic bound',\
                     'mean':'mean stratigraphic location',\
                     'true':'a priori formation time',\
                      'nid':'net in degree'}

        if title != "":
            title = " of "+title
        ax.set_title("StratSort"+title+" (ordered by "+sortnames[orderby]+")")
        ax.grid(True)

    #}}}2 end draw_sort

    def draw_graph(self, ax=None, minsize=2, minconf=0.5, cmap=pylab.cm.jet_r): #{{{2
        """
        Draw graph view of the GSN.  Plot each sub-GSN in a separate figure.
        This is janky and should be replace with something that will plot using
        graphviz instead of matplotlib.

        minsize allows you to limit which connected components are plotted, and
        avoid showing a bunch of uninteresting size 2 and 3 graphs.

        """
        from matplotlib.colors import rgb2hex
        from matplotlib import rcParams
        import Image


        sub_GSNs = self.get_sub_GSNs(minsize=minsize, minconf=minconf)

        rel_info = self.completeness()
        for sub_gsn,n in zip(sub_GSNs, range(len(sub_GSNs))):

            nc_array = np.array([ cmap(float(n)/len(sub_GSNs)) for i in range(len(sub_gsn)) ])
            sub_agsn = nx.to_agraph(sub_gsn)

            sub_agsn.graph_attr.update(splines  = 'true',\
                                       rankdir  = 'BT',\
                                       nslimit1 = '100',\
                                       nodesep  = '0.2',\
                                       ranksep  = '0.5 equally',\
                                       mclimit  = '10.0',\
                                       rank     = 'source')
            sub_agsn.node_attr.update(color=rgb2hex(nc_array[0][0:3]), shape='point', width='0.25')
            sub_agsn.layout(prog='dot')
            sub_agsn.draw('output/graphs/test_agsn.png')

            png_out = Image.open('output/graphs/test_agsn.png')
            dpi = rcParams['figure.dpi']
            fig_width  = png_out.size[0]/dpi
            fig_height = png_out.size[1]/dpi
            the_fig = plt.figure(figsize=(fig_width, fig_height))
            gax = the_fig.add_subplot(111, frameon=False)
            im = plt.imshow(png_out, origin='lower')
            gax.set_xticks([])
            gax.set_yticks([])
            gax.set_xlabel('I=%.3g' % (rel_info[n],))

                
    #}}}2 end draw_graph

    def draw_idist(self, iterations=100, minsize=4, cmap=pylab.cm.jet_r): #{{{2
        """
        Calculate some statistics of the possible information content within the
        GSN, depending on what order the features formed in.

        """
        sub_GSNs = self.get_sub_GSNs(minsize=minsize, minconf=minconf) 
        GIDs = geometric_information_distribution(self, iterations=iterations, minsize=minsize)
        nc_array = np.array([ cmap(float(i)/len(GIDs)) for i in range(len(GIDs)) ])
        the_fig = plt.figure(figsize=(12,6))
        hist_ax = the_fig.add_subplot(111)
        for idist,sub_gsn,nc in zip(GIDs, sub_GSNs, nc_array):
            n, bins, patches = hist_ax.hist(idist, facecolor=nc, normed=True, histtype='stepfilled', lw=2, edgecolor=nc, alpha=0.5)
            patches[0].set_label("%d features" % (len(sub_gsn),) )

        hist_ax.set_xlim((0,1))
        hist_ax.grid(True)
        hist_ax.legend(loc='upper right')
        hist_ax.set_title('%d samples' % (iterations,) )
        hist_ax.set_xlabel('I (proportion of recoverable information)')
        hist_ax.set_ylabel('N (normalized)')

        return(GIDs)
    #}}}2

    ##########################################################################
    # GSN Manipulation
    ##########################################################################

    def get_sub_GSNs(self, minconf=0.5, minsize=2): #{{{2
        """
        Return a list of all the connected components of the GSN having size
        greater than or equal to minsize, when only intersections with weights
        corresponding to confidences greater than minconf. If the original GSN
        has a linstack associated with it, add the subset of lineaments in each
        sub_GSN to its associated linstack in the same order in which they
        appear in the original GSN.

        """

        # Convert the GSN to an undirected graph
        # and find the connected subcomponents
        U = nx.MultiGraph()
        good_edges = [ e for e in self.edges(data=True) if np.exp(-e[2]['weight']) >= minconf ]
        U.add_edges_from(good_edges)
        connected_nbunches = nx.connected_components(U)
        G = GeoSupNet()
        G.add_edges_from(good_edges)

        # re-constitute those connected bunches as a new list of GSNs:
        sub_GSNs = []
        for nbunch in connected_nbunches:
            if len(nbunch) >= minsize:
                sub_gsn = G.subgraph(nbunch)
                # Add the appropraite sub-linstack if it's available:
                sub_gsn.linstack = [ lin for lin in self.linstack if lin in sub_gsn ]
                sub_GSNs.append(sub_gsn)

        return(sub_GSNs)
    #}}}2 end get_sub_GSNs()

    def find_sub_linstacks(self, linstack, minsize=2, minconf=0.5): #{{{2
        """
        Return a list of lists of lineaments containing any features in
        linstack which are also found in self, in the same order as they appear
        in linstack, but grouped according to which sub_GSN of self they're in.

        """

        linstacks = []
        for sub_gsn in self.get_sub_GSNs(minsize=minsize, minconf=minconf):
            linstack.append([ lin for lin in linstack if lin in sub_gsn ])

        return(linstacks)

    #}}}2 end get_sub_linstacks()

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

    def reorder(self, linstack): #{{{2
        """
        Given a new order of formation of the features in the GSN, defined by
        linstack, change the directions of the edges as required to reflect the
        new ordering, without needing to re-calculate the locations of all the
        intersections, which is computationally expensive.

        """
        # Make sure that linstack is a superset of nodes:
        sub_linstack = [ lin for lin in linstack if lin in self.nodes() ]
        assert(len(sub_linstack) == len(self))

        # for every edge in the GSN, check to see if it, or its reverse, is
        # consistent with the new linstack.  If its reverse is consistent, flip
        # the direction.  If neither is implied by the linstack, 
        linstack_pairs = linstack2pairs(linstack)
        old_edges=self.edges(data=True, keys=True)
        for edge in old_edges:
            if (edge[1],edge[0]) in linstack_pairs:
                self.add_edge(edge[1],edge[0],edge[2],edge[3])
                self.remove_edge(edge[0],edge[1],edge[2])
            else:
                assert((edge[0],edge[1]) in linstack_pairs)

        self.linstack = sub_linstack

    #}}}2

    ##########################################################################
    # Metrics and Analysis:
    ##########################################################################

    def conf_dist(self): #{{{2
        """
        Return an array containing the confidence values for all of the
        intersections in the GSN.

        """
        weights = np.array([e[2]['weight'] for e in self.edges(data=True)])
        return(np.exp(-weights))
    #}}}2

    def path_conf_dist(self): #{{{2
        """
        Return an array containing the confidence values for all of the
        shortest paths calculated within the GSN.

        """
        confs, paths = self.best_paths(minconf=0)
        out_confs = []
        for source in self.nodes():
            for target in confs[source].keys():
                out_confs.append(confs[source][target])

        return(out_confs)

    #}}}2

    def completeness(self): #{{{2 TODO: log-weights
        """
        Return a list of the proportion of pairwise orderings which are
        specified for each sub-GSN.  Indicates how information rich the
        GSNs are.

        """
        
        comps = []
        sub_gsns   = self.get_sub_GSNs()
        sub_npairs = [ len(pairs) for pairs in self.pairs() ]
        for sub_gsn, sub_npair in zip(sub_gsns, sub_npairs):
            L = len(sub_gsn)
            comps.append(sub_npair/((L*(L-1))/2.0))
        return(comps)
    #}}}2

    def best_paths(self, minconf=0.5): #2{{{
        """
        For each ordered pair of nodes (A,B) in the GSN, find the shortest
        weighted path from A to B.  Return both the paths, and their
        confidences, as two dictionaries of dictionaries such that for example:

        conf[A][B]

        would be the confidence of the path from A to B, which are part of the
        largest sub-GSN, and:

        path[A][B]

        would be the corresponding path, represented as an ordered list of
        nodes.

        Paths with confidences less than minconf are excluded.

        """

        path = {}
        conf = {}
        nodelist = self.nodes()

        for source in nodelist:
            conf[source], path[source] = nx.single_source_dijkstra(self, source)
            # restrict the output to only sufficiently confident paths
            for target in conf[source].keys():
                conf[source][target] = np.exp(-conf[source][target])
                if conf[source][target] < minconf:
                    del path[source][target]
                    del conf[source][target]

        return(conf, path)
    #2}}}

    def enumerate_cycles(self, minconf=0.5): #{{{2
        """
        Find all pairs of nodes (A,B) for which there exist paths both from A
        to B and B to A with confidences greater than minconf.  Many of these
        paths will be isomorphic with each other (i.e. passing through all of
        the same nodes... but with different starting and ending points A and
        B).  Compare all of the paths, and return one from each isomorphic
        group, as a list of edge tuples, (lin, lin, dict).

        """

        cycles = []
        cycle_confs = []
        conf, path = self.best_paths(minconf=np.sqrt(minconf))

        for n1 in self.nodes():
            for n2 in self.nodes():
                try:
                    if conf[n1][n2] and conf[n2][n1] and n1 != n2:
                        cycles.append(path[n1][n2][:-1]+path[n2][n1][:-1])
                        cycle_confs.append(conf[n1][n2]*conf[n2][n1])

                except(KeyError):
                    continue

        # only include one copy of each isomorphic cycle
        unique_cycles = []
        unique_cycle_confs = []
        cycles_as_sets = []
        for i in range(len(cycles)):
            cyc_set = set(cycles[i])
            if cyc_set not in cycles_as_sets:
                cycles_as_sets.append(cyc_set)
                unique_cycles.append(cycles[i])
                unique_cycle_confs.append(cycle_confs[i])

        # TODO: sort the returned cycles first by confidence, and then by number of nodes
        dtype = [('cycle',object),('path_length',int),('conf',float)]
        cycle_sort = np.zeros(len(unique_cycles), dtype=dtype)
        cycle_sort['cycle'] = unique_cycles
        cycle_sort['path_length'] = np.array([ len(c) for c in unique_cycles ])
        cycle_sort['conf'] = -np.array(unique_cycle_confs)
        cycle_sort.sort(order=['conf','path_length'])

        return(cycle_sort['cycle'], -cycle_sort['conf'])

    #}}}2 end enumerate_cycles

    def weighted_pairs(self, minsize=2, minconf=0.5): #{{{ TODO
        """
        Find the shortest weighted path between all pairs of nodes in the GSN.

        """

    #}}}


    ##########################################################################
    # Obsolete:
    ##########################################################################

    def valid_order(self, linstack): #{{{2
        """
        Given a list of features, ordered by their hypothesized times of formation
        (linstack), and stratigraphic sort (stratsort, as returned from
        GeoSupNet.stratigraphic_sort(), potentially having several members),
        determine whether the ordering of linstack is consistent with the
        constraints encoded within the stratsort.  Return True if it is, and False
        if it is not.

        """

        sub_GSNs = self.get_sub_GSNs()
        sub_linstacks = self.find_sub_linstacks(linstack)

        # Go through each sub-stack and sub-GSN, and see if the orderings are
        # consistent.  If any of the implied predecessors in the stack are
        # found in a feature's *actual* successors, then the ordering is
        # invalid.  Otherwise it is valid.
        for sub_stack, sub_gsn in zip(sub_linstacks, sub_GSNs):
            # create a dictionary of the successors for each lineament in the
            # GSN, so we're not doing the same work over and over again
            # TODO: Actually, for exclusion_coeff... this only needs to be done once!
            true_successors = {}
            for lin in sub_stack:
                true_successors[lin] = nx.single_source_shortest_path(sub_gsn, lin).keys()

            for n in range(len(sub_gsn.linstack)):
                # run through all the lineaments in the stack below n
                for lin in sub_stack[:n]:
                    # if any of those features are a successor of the nth lin, the ordering fails
                    if lin in true_successors[sub_stack[n+1]]:
                        return(False)
        return(True)

#}}}2 end valid_order

    def exclusion_coeff(self, numlins=2, iterations=100, minsize=3): #{{{2
        """
        Determine what fraction of possible feature orderings the GSN can
        exclude.

        numlins is the number of lineaments to try and validate in each
        ordering.  Those lineaments are chosen at random from the set of all
        features, and re-shuffled on each iteration.

        iterations is the number of trials to run.

        minsize is the minimum size of the connected components you wish to
        include.  It allows you to only look at the features which actually
        have ordering information (i.e. intersections).  This value must be
        larger than numlins.

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

    def net_in_degree(self, with_labels=True, weighted=True, normalized=False): #{{{2
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

#}}} end GeoSupNet

################################################################################
# GSN Generation and Translation:
################################################################################

def shp2gsn(lin_shp, cross_shp): #{{{
    the_gsn = lincross2gsn(shp2lincross(lin_shp, cross_shp))
    return(the_gsn)
#}}}

def shp2lincross(lin_shp, cross_shp): #{{{
    """
    Given the paths to two shapefiles, one defining a set of mapped lineaments
    (lin_shp) and one describing the superposition relationships at their
    intersections (cross_shp) return a list of lincross tuples, from which a GSN
    may be constructed.

    Each feature within lin_shp must have a field in its attribute table named
    'gid', which contains an integer uniquely identifying it within the
    layer.

    Each feature within cross_shp must have three fields in its attribute
    table, named 'top_lin', 'bottom_lin', and 'conf', indicating the gid of the
    features within lin_shp which participate in the intersection, which one
    appears to be on top, and which on the bottom, and the confidence with
    which that assignment is made, as a subjective probability, whose value
    ranges between 0.5 and 1.0, with conf=0.5 indicating that the imaging
    from which the features was digitized did not allow any determination of
    the relative ages of the features at that point.  Probabilities less than
    0.5 are prohibited, because they would in effect be indicating that the
    ordering of features ought to be reversed.

    The set of lineament crossings returned is compatible with lincross2gsn,
    defined below, and takes the form of a list of 3-tuples containing two
    lineament objects, and a dictionary of data with which to construct an edge
    in the GSN:

    (bottom_lin, top_lin, {'weight':-log(conf),'lon':LON,'lat':LAT})

    The confidence is transformed to -log(conf) in the graph weighting to allow
    us to use the sums of weights to infer the products of probabilities.

    """

    # Create a dictionary of all the lineaments defined in lin_shp, keyed by
    # their LIN_ID values:
    lins, lin_ids = lineament.shp2lins(lin_shp, lin_ids=True)

    # check and make sure that all of the lin_id values are unique
    # How do I do that?  Need to look it up.

    lindict = {}
    for lin, lin_id in zip(lins, lin_ids):
        lindict[lin_id] = lin

    # Iterate over the set of points described in cross_shp, and for each one
    # look up the lineament objects involved in the intersection in the
    # dictionary we just defined, and associate those Lineament objects with the
    # edge data we can read from the points
    cross_pt_data = ogr.Open(cross_shp, update=0)
    cross_pt_layer = cross_pt_data.GetLayer(0)

    lincross = []
    cross_pt = cross_pt_layer.GetNextFeature()
    while cross_pt is not None:
        bottom_lin_idx = cross_pt.GetFieldIndex('bottom_lin')
        bottom_lin_id  = cross_pt.GetField(bottom_lin_idx)
        bottom_lin     = lindict[bottom_lin_id]

        top_lin_idx = cross_pt.GetFieldIndex('top_lin')
        top_lin_id  = cross_pt.GetField(top_lin_idx)
        top_lin     = lindict[top_lin_id]

        # The confidence is the probability (0.5 < P < 1.0) that the 
        # lineaments have the observed temporal ordering.  P < 0.5 is
        # prohibited because it would be equivalent to P=(1-P), with
        # the observed temporal ordering reversed.
        conf_idx = cross_pt.GetFieldIndex('conf')
        conf     = cross_pt.GetField(conf_idx)
        # Because most weighted graph algorithms sum the weights of the
        # edges, and what we're really interested in is their product
        # (the probability of a particular path being true) we need to
        # transform the confidence into a weight, logarithmically:
        pt_wt = -np.log(conf)

        cross_pt_geom = cross_pt.GetGeometryRef()
        pt_lon, pt_lat, no_z_coord = cross_pt_geom.GetPoint()
        pt_lon = np.radians(pt_lon)
        pt_lat = np.radians(pt_lat)

        lincross.append( (bottom_lin, top_lin, {'weight':pt_wt, 'lon':pt_lon, 'lat':pt_lat}) )

        # grab the next intersection:
        cross_pt = cross_pt_layer.GetNextFeature()

    return(lincross)

#}}} end shp2lincross

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
    new_gsn.linstack = [ lin for lin in old_gsn.linstack ]
    return(new_gsn)
#}}}

def linstack2pairs(linstack): #{{{
    """
    Takes a list of Lineament objects, ordered by their times of formation, and
    returns a set of tuples describing all of the pairwise orderings implied by
    the stack, such that (A,B) means that A formed before B.
    """

    pairs = set([])
    for n in range(len(linstack)):
        for m in range(n+1,len(linstack)):
            pairs.add((linstack[n],linstack[m]))
    return(pairs)
#}}} end linstack2pairs

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

def linstack_regular(nlins=10, overpole=False): #{{{
    """
    Create a regular map with a known order of formation for testing the GSN
    algorithms. nlins indicates how many features there are in each of the
    north-south and east-west orientations.  ncross indicates where in the
    sequence of east-west features the north-south features should appear.  If
    ncross is none, the whole ordering is randomized.

    """

    northsouth_init_lons = np.linspace(-np.pi/3.1, np.pi/3.1, nlins)
    northsouth_fin_lons  = np.linspace(-np.pi/3.1, np.pi/3.1, nlins)
    northsouth_init_lats = np.array([-np.pi/2.1,]*nlins)
    northsouth_fin_lats  = np.array([ np.pi/2.1,]*nlins)

    eastwest_init_lons = np.array([-np.pi/3.0]*nlins)
    eastwest_fin_lons  = np.array([ np.pi/3.0]*nlins)
    eastwest_lats      = np.linspace(-np.pi/3.0, np.pi/3.0, nlins)

    if overpole is True:
        northsouth_fin_lons += np.pi
        northsouth_fin_lats = np.array([ np.pi/2.01,]*nlins)

    init_lons = np.concatenate([northsouth_init_lons, eastwest_init_lons])
    init_lats = np.concatenate([northsouth_init_lats, eastwest_lats])
    fin_lons  = np.concatenate([northsouth_fin_lons,  eastwest_fin_lons])
    fin_lats  = np.concatenate([northsouth_fin_lats,  eastwest_lats])

    outlins = []
    for lon1, lat1, lon2, lat2 in zip(init_lons, init_lats, fin_lons, fin_lats):
        outlins.append(lineament.lingen_greatcircle(lon1, lat1, lon2, lat2))

    pylab.shuffle(outlins)

    return(outlins)

#}}} end linstack_regular 

def linstack_file(nlins=0, linfile='output/lins/map_nsrfit', spin=False, tpw=False): #{{{
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
    else:
        linstack = nsrhist.make_crazy(nsrhist.linresample_byN(lins), tpw=tpw, spin=spin)

    linstack = lins[:nlins]

    return(linstack)
        
#}}}

def linstack_nsr(nlins=50, ncross=5): #{{{
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

#}}} end linstack_nsr

def linstack2lincross(linstack): #{{{
    """
    Take an ordered list of lineaments, and return it as a list of GSN edges.

    """

    new_gsn = linstack2gsn(linstack)
    lincross = new_gsn.edges(data=True)
    return(lincross)

#}}}

def linstack_random(nlins=100): #{{{
   """
   Create a random collection of great circle segments.

   """

   test_lins = nsrhist.random_gclins(nlins)
   for lin in test_lins:
       lin.lons = lineament.fixlons(lin.lons)
   return(test_lins)
#}}}

def random_gsn(nlins=100): #{{{
    """
    Create a random collection of great circle segments, with intersections
    having ambiguous superposition relationships.

    """
    lincross = linstack2lincross(linstack_random(nlins=nlins))

    for cross in lincross:
        cross[2]['weight'] = -np.log(np.random.uniform(low=0.5, high=1.0))

    return(lincross2gsn(lincross))

#}}}

################################################################################
# GSN Metrics
################################################################################

def pairs_relevance(A, B): #{{{
    """
    Within two sets of statements (about the ordering of features), A and B,
    each statement has to be in one of three classes: either they agree (the
    statement is in both A and B), they disagree (the statement is in A, and
    its reverse is in B) or they are unrelated (the statement appears in A or
    B, but not in the other)

    This function returns the number of statements which the two sets share
    either in agreement or disagreement) divided by the number of statements
    within A, and so indicates how relevant B as a whole is to A.

    """
    
    agree    = A.intersection(B)
    disagree = A.intersection(set([ pair[::-1] for pair in B ]))

    return(float(len(agree)+len(disagree))/len(A))
#}}}

def pairs_agreement(A, B): #{{{
    """
    Within two sets of statements (about the ordering of features), A and B,
    each statement has to be in one of three classes: either they agree (the
    statement is in both A and B), they disagree (the statement is in A, and
    its reverse is in B) or they are unrelated (the statement appears in A or
    B, but not in the other)

    This function returns the number of statements on which the two sets agree,
    divided by the number of statements which they share at all (either in
    agreement or disagreement)

    """

    agree    = A.intersection(B)
    disagree = A.intersection(set([ pair[::-1] for pair in B ]))

    return(len(agree)/(len(agree)+len(disgree)))
#}}}

def pairs_disagreement(A, B): #{{{
    """
    Within two sets of statements (about the ordering of features), A and B,
    each statement has to be in one of three classes: either they agree (the
    statement is in both A and B), they disagree (the statement is in A, and
    its reverse is in B) or they are unrelated (the statement appears in A or
    B, but not in the other)

    This function returns the number of statements on which the two sets
    disagree, divided by the number of statements which they share at all
    (either in agreement or disagreement).

    """

    return(1.0-pairs_disagreement(A,B))
#}}}

def geometric_information_distribution(gsn, iterations=100, minsize=2): #{{{
    """
    Calculate the relative information content of a GSN repeatedly, re-ordering
    its true order of formation (linstack) repeatedly, while holding its
    geometry constant.  Do this for each connected component separately, and
    return 2D array of the results, with one row for each sub GSN.

    """

    # Make a copy of the GSN so we don't alter the original:
    test_gsn = gsn2gsn(gsn)

    sub_GSNs = test_gsn.get_sub_GSNs(minsize=minsize, minconf=minconf)
    sub_comps = np.zeros((len(sub_GSNs),iterations))
    for sub_gsn,n in zip(sub_GSNs,range(len(sub_GSNs))):
        for i in range(iterations):
            pylab.shuffle(sub_gsn.linstack)
            sub_gsn.reorder(sub_gsn.linstack)
            sub_comps[n,i] = sub_gsn.completeness()[0]

    return(sub_comps)
#}}}

################################################################################
# Helper Functions:
################################################################################
def wkt2proj(wkt_str=None): #{{{
    """
    Take a "Well Known Text" (WKT) spatial reference system (SRS) and convert
    it to a set of Proj.4 command line arguments, for consumption by qGIS.

    """
    if wkt_str is None:
        wkt_str = """PROJCS["Europa_Mercator_AUTO",GEOGCS["Europa 2000",DATUM["D_Europa_2000",SPHEROID["Europa_2000_IAU_IAG",1564130.0,488.79062499999998]],PRIMEM["Greenwich",0],UNIT["Decimal_Degree",0.0174532925199433]],PROJECTION["Mercator_1SP"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",180],PARAMETER["Standard_Parallel_1",0],UNIT["Meter",1]]"""
    coord_ref_sys = osr.SpatialReference(wkt_str)
    return(coord_ref_sys.ExportToProj4())
#}}}

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

################################################################################
# Testing Functions:
################################################################################

def test(nlins=100, lintype='real', orderby='mean', ncross=None, spin=False, tpw=False, minsize=4, draw_map=True, draw_sort=True, draw_graph=True, draw_idist=True): #{{{
   """
   A short unit test of the sorting functionality.
   
   lintype may be one of: 'nsr', 'regular', or 'random', corresponding to the
   similarly named linstack routines defined below.

   """

   print("Generating linstack")
   if lintype == 'nsr':
       test_lins = linstack_nsr(nlins=nlins, ncross=ncross)
   elif lintype == 'regular':
       test_lins = linstack_regular(nlins=nlins)
   elif lintype == 'real':
       test_lins = linstack_file(nlins=nlins, spin=spin, tpw=spin)
   else:
       assert(lintype == 'random')
       test_lins = linstack_random(nlins=nlins)

   print("Converting map to GSN")
   test_gsn = linstack2gsn(test_lins)

   print("Plotting results")
   if draw_sort is True:
       test_gsn.draw_sort(orderby=orderby, title=lintype, minsize=minsize)
   if draw_map is True:
       test_gsn.draw_map(minsize=minsize)
   if draw_graph is True:
       test_gsn.draw_graph(minsize=minsize)
   if draw_idist is True:
       test_gsn.draw_idist(minsize=minsize, iterations=iterations)

   return(test_gsn)

#}}} end test()

################################################################################
# TODO:
################################################################################
# Mapping:
# --------
# add descriptive attributes to lineaments:
#   - width/form:
#     - fracture (hairline, no apparent topography)
#     - ridge (single apparent topographic feature, or significant width)
#     - multi-ridge (2-3 apparent topographic features)
#     - band (many sub-parallel features taking up significant area)
#   - color:
#     - light
#     - mid
#     - dark

# ============================================================================== 
# GSN v1.1 (intersection confidences):
# ==============================================================================

# Introduce weights on the edges which are less than unity, to represent the
# uncertainty in the relationships they depict.

# Use bi-directed graph with P and 1-P as the weights for any given
# intersection, where P (>= 0.5) is the confidence (probability) that the edge
# goes in the more probable direction.

# Find the shortest path between all ordered pairs of nodes in the graph using
# Dijkstra's algorithm.  As weights for the edges, use -log(P), such that the
# sums of the weights correspond to the product of the probabilities.

# With bi-directed graph, the interesting question, regarding cycles, becomes
# for a given threshold (maximum path length) do there exist paths both from A
# to B *and* from B to A such that the sum of the weights of the edges involved
# is below the threshold.  The right threshold is probably 0.5

# The distribution of shortest path lengths (instead of proportion of
# relationships which are defined) becomes the interesting thing to look at,
# because all relationships within a connected network of features are
# "defined", though some of them have miniscule (or zero) probability.

# We want to be able ask questions like:
#   - given three classes of features (X,Y,Z) how consistent are the intersections
#     we mapped with X forming first, then Y, then Z?
#   - Are there any clear instances of lineament re-activation?
#   - How consistent are the cross cutting relationships with the order of formation
#     one would infer from NSR?  Is it much more consistent with NSR than with a
#     random order of formation?

# A binary ordering of two features is well defined if P(A,B) >> P(B,A).  Can
# define a metric of clarity which is the magnitude of the difference between
# the two: C = |P(A,B) - P(B,A)|.
#   - If one is large, and the other small, then C is large.
#   - If both are small, then C is small.
#   - If both are large, things get weird.  What if instead, we use:
#     max(P(A,B),P(B,A))
#     - |1.0-0.0| - |0.5-0.5| =  1.0
#     - |0.9-0.1| - |0.6-0.4| =  0.6
#     - |0.8-0.2| - |0.7-0.3| =  0.2
#     - |0.7-0.3| - |0.8-0.2| = -0.2
#     - |0.6-0.4| - |0.9-0.1| = -0.6
#     - |0.5-0.5| - |1.0-0.0| = -1.0

# How does retrieveable information vary as a function of geometry
# Create a circular or otherwise geometrically unbiased map space, and show how
# information retrieved varies with:
#   - lineament density per unit area.
#   - lineament length (as a fraction of the scale of the area)
#
# This could be a 2D plot: length on one axis, density on the other, with color
# showing how much information was retrieved (sweet plot)

################################################################################

# ============================================================================== 
# GSN v1.2 (allow cycles): #{{{
# ============================================================================== 

# INDUCE CYCLES ON A TEST DAG:
# ----------------------------
#   Given an acyclic GSN, generate a similar GSN, that has experienced
#   re-activation by reversing the direction of some proportion of the edges.

# ENUMERATE CYCLES:
# ----------------------------
#   Given a GSN containing cycles, return a set of "path" graphs describing all
#   of the simple cycles in the original GSN.  Give the option of removing all
#   the intersections below a particular confidence threshold first in order to
#   speed things up, and focus on the most interesting features.
#
#   To first order, this is probably just finding all the strongly connected
#   components of the digraph.

# IDENTIFY REACTIVATION SITES:
# ----------------------------
#   Figure out exactly what algorithm the 'dot' program uses to draw
#   hierarchical graphs.  It's probably about as good as we're going to get at
#   identifying the back-edges.

#   Test whether multiple back edges clustered along a given lineament are a
#   good indication of reactivation having reversed them, by using test
#   datasets.

#}}}

# ============================================================================== 
# GSN v1.3 (cycle resolution): #{{{ FUTURE WORK
# ============================================================================== 

# LINEAMENT SPLITTING:
#   Given a best guess at which intersections have been affected by reactivation,
#   determine where to cut the lineament most involved in order to separate out
#   the reactivated portion, and allow it to sort independently.

# RENDER GSN ACYCLIC:
#   - Either remove back edges one-by-one until the graph is acyclic or
#   - Split the feature they lie on and attempt to re-sort.

#}}}

