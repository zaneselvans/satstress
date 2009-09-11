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
# * Come up with a set of relatively simple sorting test cases just to make
#   sure that the old code does what it's supposed to.
#
# Things GSNs need to do (v1.0):
# ===========================================
# Given:
#   - a list of Lineament objects, ordered by time of formation
# Generate:
#   - ordered lists of intersections along each lineament
#     (This is actually the primary dataset)
#     - location (lon,lat)
#     - z-order (bottomlin_id, toplin_id)
#     - confidence of intersection (0 < K < 1)
#   - unweighted GSN corresponding to that map
# Why:
#   - allows easy creation of test datasets
#
# Given a GSN:
#   - enumerate all the simple cycles
#     - may need to use a confidence cutoff to keep it computable...
#   - render the GSN acyclic
#   - display the GSN graphically (see above)
#   - perform a stratigraphic sort on the GSN
#     - bottom-up topologic sort defines upper bounds
#     - top-down topologic sort defines lower bounds
#
# Given a stratigraphic sort:
#   - Create a graphical representation of the sort
#   - Calculate a metric of how much of the ordering information was retrieved
#   - Calculate mean value of best_b for each set of equivalent features
#
################################################################################

"""linsort: a stratigraphic classification tool by zane.crawford@colorado.edu

   linsort [-h|--help]
     * prints this help message

   linsort -f <infile> -t <conf_thresh> -i <outfile1> -l <outfile2>
     * take the features and intersections described in <infile>
     * remove all the intersections less confident than <conf_thresh>
     * perform a stratigraphic sort of the features
     * output the fates of the intersections to outfile1
     * output the sort of the features to outfile2
"""
import sys
import getopt
import string
import re
import networkx as nx

class XDG_Node:
    """XDG_Node is an object class which contains data and methods associated with
    the mapped features being sorted.  They are represented by nodes in the graph.
    """
    def __init__(self, feature_id):
        """This function initializes the XDG_Node, and takes a single argument - the
        unique ID associated with the feature.  It also sets the stratigraphic 
        equivalence class to -1, a value which indicates the feature has not yet
        been placed in the sort.
        """
        self.feature_id = feature_id
        self.strat_eq_class = -1

    def __hash__(self):
        """The NetworkX module requires that node objects be hashable, so we must
        define the __hash__ method for our XDG_Node object.
        """
        return hash(self.feature_id)

    def weighted_connectedness(self, in_xdg):
        """The "weighted connectedness" of a node in the graph is simply the sum
        of all of the weights (in our case confidences) of the edges which are
        incident upon it.  It is a measure of how connected to the rest of the
        graph the node is - in this application, it gives us some idea of how
        well defined the node's place in the stratigraphic sort is.  Because one
        node object could be included in multiple graphs, and the connectedness
        depends not just on the node, but on the structure of the graph in which
        it is embedded, this method must be given the graph in which the node is
        embedded as an argument.
        """
        wc = 0
        for edge in in_xdg.edges(self):
            wc = wc + edge[2].confidence

        return wc

    def weighted_degree(self, in_xdg):
        """
        The weighted degree of a node is the sum of the weights of the incoming
        edges, minus the sum of the weights of the outgoing edges.  It gives a
        rough estimate of where in the final sort we should expect the node
        to appear, and so is useful in ordering the nodes prior to sorting, in
        a way which allows us to remove a small number of edges to make the graph
        acyclic.
        """
        wt_in_deg  = 0
        wt_out_deg = 0

        for node in in_xdg.successors(self):
            out_edge = in_xdg.get_edge(self, node)
            wt_out_deg = wt_out_deg + out_edge.confidence

        for node in in_xdg.predecessors(self):
            in_edge = in_xdg.get_edge(node, self)
            wt_in_deg = wt_in_deg + in_edge.confidence

        return(wt_in_deg - wt_out_deg)
        
    def lower_bound(self, in_xdg):
        """
        This method returns the lowest stratigraphic class (which unintuitively
        means a higher number - the deeper in the stratigraphy you are, the 
        higher the number associated with the stratigraphic class) of which the
        node could conceivably be a member, based on the stratigraphic class of
        its nearest descendant node.  If the node has no descendants, then the
        lower bound is the deepest stratigraphic class.
        """
        if (self.strat_eq_class == -1):
            raise("UnsortedGraph")

        successors = in_xdg.successors(self)

        if(len(successors) > 0):
            least = successors[0].strat_eq_class
            for node in successors:
                if(node.strat_eq_class < least):
                    least = node.strat_eq_class
            return(least - 1)

        else:
            least = 0
            for node in in_xdg:
                if node.strat_eq_class > least:
                    least = node.strat_eq_class
            return(least)

class XDG_Edge:
    def __init__(self, intersection_id, confidence):
        self.intersection_id = intersection_id
        self.confidence = float(confidence)
        self.status = "TBD"

    def __cmp__(self, other):
        return cmp(self.confidence, other.confidence)

    def __hash__(self):
        return hash(self.intersection_id)


def csv_to_xdg(infile):
# Create an empty graph representing the crosscutting relationships
# which are clear enough to be used in the sorting algorithm
    INFILE = open(infile)
    new_xdg = nx.XDiGraph(multiedges=False)

    for line in INFILE.readlines():
        # strip out the comments
        if(re.match("^#",line)):
            continue
        else:
            # get rid of newlines and whitespace
            line = re.sub("(\n|\s)","",line)
            # split the input line by commas
            fields = string.splitfields(line, ',')

            top_lin = XDG_Node(fields[1])
            bot_lin = XDG_Node(fields[2])

            for node in new_xdg:
                if node.feature_id == fields[1]:
                    top_lin = node
                if node.feature_id == fields[2]:
                    bot_lin = node

            new_edge = XDG_Edge(fields[0], fields[3])
            new_xdg.add_edge(top_lin, bot_lin, new_edge)

    return (new_xdg)

def filter_by_confidence(in_xdg, conf_thresh):
    filtered_xdg = nx.XDiGraph(multiedges=False)
    ignored_xdg  = nx.XDiGraph(multiedges=False)

    for edge in in_xdg.edges():
        if(edge[2].confidence < conf_thresh):
            edge[2].status = "ignored"
            ignored_xdg.add_edge(edge)
        else:
            filtered_xdg.add_edge(edge)

    return (filtered_xdg, ignored_xdg)

def is_DAG(in_xdg):
    if nx.topological_sort_recursive(in_xdg):
        return True
    else:
        return False

def enforce_DAG(in_xdg):

    acyclic_xdg = nx.XDiGraph(multiedges=False)
    acyclic_xdg.add_edges_from(in_xdg.edges())

    removed_xdg = nx.XDiGraph(multiedges=False)

# create a list of all edges in the graph
    sorted_edges = acyclic_xdg.edges()

# this somewhat obscure line sorts the edges in the list using the XDG_Edge
# object as the key, and there is a __cmp__ method defined for that object
# which compares the confidences associated with the intersection - so this
# sort is sorting the intersections by confidence, in ascending order.
    sorted_edges.sort(lambda x,y: cmp(x[2],y[2]))

# So long as the graph is not acyclic
    while not is_DAG(acyclic_xdg):
# iterate up through the list of edges sorted by confidence
# find the lowest confidence edge which is pointing the "wrong way".
# remove it, and escape the for loop to test and see if the graph is now
# acyclic.  The way this is written right now it is adjusting the weighted
# degree of all of the nodes each time an edge is removed, which I suspect
# will improve the algorithm (vs. using the same old graph we started with)
        for edge in sorted_edges:
            if ((edge[0].weighted_degree(acyclic_xdg) > edge[1].weighted_degree(acyclic_xdg))):
                edge[2].status = "removed"
                removed_xdg.add_edge(edge)
                acyclic_xdg.delete_edge(edge)
                continue

    for edge in acyclic_xdg.edges():
        edge[2].status = "included"

    return (acyclic_xdg, removed_xdg)

# - create an empty graph (sorted_xdg) which will ultimate be filled in with
#   the nodes that have their strat_eq_class set.
# - create a temporary graph (tmp_xdg) which we will destroy during the sort
# - then, so long as we have any nodes left in the graph
# - create an empty graph (condemned).
# - look for all of the "top" nodes in the temporary graph
# - when we find one, set its strat_eq_class, then add it to the output
#   graph and also add it to the list of nodes we'll delete this time through
#   We can't remove the "top" nodes one by one because that would change who
#   was on top next time around!
# - once we've found all the "top" nodes, remove them from the temporary graph
#   and increment the stratigraphic equivalence class
# - when we're really finished, we need to add all of the original edges to
#   our output graph...  since we couldn't add them during the sort.
def stratigraphic_sort(in_xdg):

    sorted_xdg = nx.XDiGraph(multiedges=False)

    tmp_xdg = nx.XDiGraph(multiedges=False)
    tmp_xdg.add_edges_from(in_xdg.edges())

    strat_eq_class = 0

    while(tmp_xdg.order() > 0):

        condemned = nx.XDiGraph(multiedges=False)

        for node in tmp_xdg.nodes():
            if(tmp_xdg.in_degree(node) == 0):
                node.strat_eq_class = strat_eq_class
                sorted_xdg.add_node(node)
                condemned.add_node(node)

        tmp_xdg.delete_nodes_from(condemned)
        strat_eq_class = strat_eq_class + 1

    sorted_xdg.add_edges_from(in_xdg.edges())
    return(sorted_xdg)

def output_sorted_lineaments(in_xdg, outfile):

    OUTFILE = open(outfile, 'w')
    OUTFILE.write("# LIN_ID, UPPER_BOUND, LOWER_BOUND, WEIGHTED_CONNECTEDNESS\n")

    node_list = in_xdg.nodes()
    node_list.sort(lambda x,y: cmp(x.strat_eq_class,y.strat_eq_class))

    for node in node_list:
        out = []
        out.append(str(node.feature_id))
        out.append(str(node.strat_eq_class))
        out.append(str(node.lower_bound(in_xdg)))
        out.append(str(node.weighted_connectedness(in_xdg)))
        out = ", ".join(out)
        out = out + "\n"
        OUTFILE.write(out)

    OUTFILE.close()

def output_intersection_fates(in_xdg, outfile):

    OUTFILE = open(outfile, 'w')
    OUTFILE.write("# INTERSECTION_ID, TOP_LIN, BOT_LIN, CONFIDENCE, STATUS\n")

    for edge in in_xdg.edges():
        out = []
        out.append(str(edge[2].intersection_id))
        out.append(str(edge[0].feature_id))
        out.append(str(edge[1].feature_id))
        out.append(str(edge[2].confidence))
        out.append(str(edge[2].status))
        out = ", ".join(out)
        out = out + "\n"
        OUTFILE.write(out)

    OUTFILE.close()

# The main() function controls the overall flow of the program
def main(argv=None):
    if argv is None:
        argv = sys.argv
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:i:l:t:", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        return(2)

    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            return(0)
        if o in ('-f'):
            infile = a
        if o in ('-i'):
            outfile_intersections = a
        if o in ('-l'):
            outfile_lineaments = a
        if o in ('-t'):
            conf_thresh = float(a)
    # process arguments
    #for arg in args:
    #    process(arg) # process() is defined elsewhere

    # create an all-inclusive graph with the data in our input file
    all_input_xdg = csv_to_xdg(infile)
    print("read in %d lineaments..." % (all_input_xdg.order()) )
    print("having %d intersections..." % (len(all_input_xdg.edges())) )

    # filter the graph based on the minimum acceptable confidence which
    # the user has specified, removing all edges with less than the threshold
    (filtered_xdg, ignored_xdg) = filter_by_confidence(all_input_xdg, conf_thresh)
    print("ignored %d intersections (with confidence < %.2f)..." % (len(ignored_xdg.edges()), conf_thresh))

    (acyclic_xdg, removed_xdg) = enforce_DAG(filtered_xdg)
    print("removed %d intersections to enforce an acyclic graph..." % (len(removed_xdg.edges())))

    # perform a stratigraphic sort on the graph, setting strat_eq_class
    # for each of the lineaments.  Removed edges get dumped into 
    # the removed_xdg graph so that we can see what happened to them.
    sorted_xdg = stratigraphic_sort(acyclic_xdg)
    print("sorted %d lineaments..." % (sorted_xdg.order()))
    print("having %d intersections..." % (len(sorted_xdg.edges())))

    # spit out the results of the sort
    output_sorted_lineaments(sorted_xdg, outfile_lineaments)
    print("Sending lineament sort results to %s" % (outfile_lineaments))

    # build the output graph, containing all of the edges
    all_output_xdg = nx.XDiGraph(multiedges=False)
    all_output_xdg.add_edges_from(ignored_xdg.edges())
    all_output_xdg.add_edges_from(removed_xdg.edges())
    all_output_xdg.add_edges_from(sorted_xdg.edges())

    # print the fates of all the intersections
    output_intersection_fates(all_output_xdg, outfile_intersections)
    print("Sending intersection fates to %s" % (outfile_intersections))

    # how many lineaments are in each stratigraphic class
    strat_classes = [ node.strat_eq_class for node in sorted_xdg ]
    strat_hist = {}
    for n in strat_classes:
        if strat_hist.has_key(n):
            strat_hist[n] += 1
        else:
            strat_hist[n] = 1

    strat_levels = strat_hist.keys()
    strat_levels.sort()

    for level in strat_levels:
        print("Stratigraphic class %d has %d members" % (level, strat_hist[level]))
    
    removed_confidences = [ edge[2].confidence for edge in removed_xdg.edges() ]

    # average confidence of intersections removed
    print("Average confidence of removed intersections: %.2f" % (sum(removed_confidences)/float(len(removed_confidences))))

    # highest confidence of any removed intersection
    print("Highest confidence of a removed intersection: %.2f" % (max(removed_confidences)))

    return(0)

# this allows the program to be run as a non-interactive, stand-alone, tool
if __name__ == "__main__":
    sys.exit(main())


class CrossCutGraph: #{{{
    """A graph (network) describing the superposition relationships between features.

    At each point where two features (probably Lineaments) intersect it is possible
    to infer a cross-cutting relationship.  Each edge in the CrossCutGraph defines
    one of these relationships.  Each node in the graph represents a feature
    (Lineament).  An edge from A to B means A is on top of B.
    Each relationship so defined can have a confidence associated with it.

    This is a stub class, to be done after the NSR Fit Histogram paper...

    """
    pass
#}}} end class CrossCutGraph

