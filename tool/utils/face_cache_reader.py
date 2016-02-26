from . import face_parameterization
from . import geometry_cache

def read(fn):
    c = geometry_cache.open(fn); f = c.file
    for p in f.by_type("IfcProduct"):
        if p.Representation is not None:
            eid = p.id()
            try:
                s = face_parameterization.shape(c.create_shape(p, force_from_cache=True))
            except: continue
            for fid, fa in enumerate(s.faces): yield (eid, fid), fa

class lazy(object):
    def __init__(self, fn):
        self.c = geometry_cache.open(fn); 
        self.f = self.c.file
        self.d = {}
    def __getitem__(self, i):
        if i not in self.d:
            def obtain():
                eid, rid = i
                try: inst = self.f[eid]
                except: return None
                if not inst.is_a("IfcProduct"): return None
                if not inst.Representation: return None
                try:
                    s = face_parameterization.shape(self.c.create_shape(inst, force_from_cache=True))
                except: return None
                for fid, fa in enumerate(s.faces):
                    self.d[(eid, fid)] = fa
            obtain()
        return self.d[i]
        