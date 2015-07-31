import numpy as np

def a2p(o,z,x):
    y = np.cross(z, x)
    r = np.eye(4)
    r[:-1,:-1] = x,y,z
    r[-1,:-1] = o
    return r.T

def axis2placement(plc):
    z = np.array(plc.Axis.DirectionRatios if plc.Axis else (0,0,1))
    x = np.array(plc.RefDirection.DirectionRatios if plc.RefDirection else (1,0,0))
    o = plc.Location.Coordinates
    return a2p(o,z,x)
    
def local_placement(plc):
    if plc.PlacementRelTo is None:
        parent = np.eye(4)
    else:
        parent = local_placement(plc.PlacementRelTo)
    return np.dot(axis2placement(plc.RelativePlacement), parent)
