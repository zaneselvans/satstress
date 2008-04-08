class Lineament: #{{{
    """A one dimensional feature on the surface of a spherical satellite.

    Described as a list of (co-lat,lon) points.  May be transformed as a whole,
    using translations and rotations.

    May be read in from and output to the ESRI "generate" file format.

    """
#}}}

class FailureGeom: #{{{
    """Describes the geometry relating a stress to the tectonic feature it results in.

    Assumes that a fault propagates at some rate across the surface, which may be
    fast enough relative to the forcing timescale that it can be treated as 
    instantaneous (as generally is the case with NSR).
    """
#}}} end class FailureGeom

class TensileFracture(FailureGeom): #{{{
    """Failure that occurs in tension, parallel to least tensile stress component."""
#}}} end class TensileFracture

class StrikeSlip(FailureGeom): #{{{
    """Failure that occurs when principal stress components have opposite sign.

    This is a stub class.
    """
    pass
#}}} end class StrikeSlip

class ThrustFault(FailureGeom): #{{{
    """Failure that occurs when stress principal components are both compressive.

    This is a stub class.
    """
    pass
#}}} end class ThurstFault
