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
