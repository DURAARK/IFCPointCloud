import sys
import json

from . import geometry_cache
from . import face_parameterization

tupelize_occ = lambda occ_xyz: tuple((occ_xyz.X(), occ_xyz.Y(), occ_xyz.Z()))

def obtain(fn):

    cache = geometry_cache.open(fn); f = cache.file

    di = {}

    for p in f.by_type("IfcProduct"):
        if p.Representation is not None:
            eid = p.id()
            try:
                s = face_parameterization.shape(cache.create_shape(p, force_from_cache=True))
            except: continue
            for fid, fa in enumerate(s.faces):
                try:
                    (u1, v1), (u2, v2) = fa.uv_bounds
                    di[(fid, eid)] = tuple(map(tupelize_occ, (fa.origin, fa.normal, fa.axis))) + (0., 0., u2 - u1, v2 - v1)
                except:
                    # fails for non-planar surfs... let's hope they do not have associated points
                    pass

    return json.dumps(list(di.items()))
