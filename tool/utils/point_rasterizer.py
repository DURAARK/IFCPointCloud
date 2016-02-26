import math
import pickle
import numpy as np
import scipy.interpolate    

from collections import namedtuple

from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from scipy.spatial import cKDTree

grid = namedtuple('grid', 'placement unum vnum uspacing vspacing'.split())
grid_values = namedtuple('grid_values', 'grid values'.split())

class rasterizer(object):
    def __init__(self, face_bounds):
        self.face_bounds = face_bounds

    def rasterize(self, eid, fid, uvws, grid_spacing=0.1):
        
        (u1, v1), (u2, v2) = self.face_bounds[(eid, fid)]
        minu, minv, maxu, maxv = min(uvws[:,0]), min(uvws[:,1]), max(uvws[:,0]), max(uvws[:,1])
        numu, numv = (maxu - minu) * (u2 - u1) / grid_spacing, (maxv - minv) * (v2 - v1) / grid_spacing
        numu, numv = map(lambda f: int(math.ceil(f)), (numu, numv))
        numu, numv = numu + 1, numv + 1
        
        grid_sample_points = np.meshgrid(
            np.linspace(minu, maxu, numu),
            np.linspace(minv, maxv, numv)
        )
        
        uspacing = (maxu - minu) / (numu - 1.)
        vspacing = (maxv - minv) / (numv - 1.)
        placement = minu, minv
        
        try:
            values = scipy.interpolate.griddata(uvws[:,0:2], uvws[:,2], tuple(grid_sample_points), 'linear')
        except:
            values = scipy.interpolate.griddata(uvws[:,0:2], uvws[:,2], tuple(grid_sample_points), 'nearest')
        
        area_in_wcs = ((v2 - v1) * (u2 - u1))
        threshold = area_in_wcs / math.sqrt(len(uvws)) / 2. # ** 2
        
        # Keep these UV's are normalized here, so for distances
        # a non-uniformly scaled copy needs to be made.
        uvs = uvws[:,0:2].copy()
        uvs *= (u2 - u1), (v2 - v1)
        grid_sample_points[0] *= (u2 - u1)
        grid_sample_points[1] *= (v2 - v1)
        tree = cKDTree(uvs)
        gridpoints = _ndim_coords_from_arrays(tuple(grid_sample_points), ndim=uvs.shape[1])
        dists = tree.query(gridpoints)[0]

        values[dists > threshold] = np.nan
            
        return grid_values(grid(placement, numu, numv, uspacing, vspacing), values)

if __name__ == '__main__':
    import point_reader
    import matplotlib.pyplot as plt
    
    POINT_PATTERN = "../intermediate_files/Haus30/associated_points/points-subset-%d.bin"
    FACE_BOUNDS_FILE = "../intermediate_files/Haus30/face_bounds.dat"
    
    r = rasterizer(FACE_BOUNDS_FILE)
    
    associated, unassociated = point_reader.read(POINT_PATTERN, range(2))
    for k, ds in associated.iteritems():
        (u1, v1), (u2, v2) = faces[k].uv_bounds
        if ((v2 - v1) * (u2 - u1)) < 5.: continue
        a = ds.merge()
        if len(a) < 128: continue
        result = r.rasterize(k[0], k[1], a).values
        plt.figimage(result)
        plt.show()
