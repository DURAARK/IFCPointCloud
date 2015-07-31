import os, gzip
import ifcopenshell, ifcopenshell.geom
import OCC.BRepTools

# Specify to return pythonOCC shapes from ifcopenshell.geom.create_shape()
settings = ifcopenshell.geom.settings()
settings.set(settings.USE_PYTHON_OPENCASCADE, True)

class c(object):
    def __init__(self, fn):
        p = os.path.join("intermediate_files", "%s.__cache__" % os.path.basename(fn))
        if not os.path.exists(p):
            os.mkdir(p)
        self.dirname = p
        self.file = ifcopenshell.open(fn)
    def create_shape(self, p, force_from_cache=False):
        ss = OCC.BRepTools.BRepTools_ShapeSet()
        fn = os.path.join(self.dirname, "%d.brep.gz" % p.id())
        if os.path.exists(fn):
            ss.ReadFromString(gzip.open(fn).read())
            return ss.Shape(ss.NbShapes())
        elif not force_from_cache:
            shp = ifcopenshell.geom.create_shape(settings, p).geometry
            ss.Add(shp)
            gzip.open(fn, 'w').write(ss.WriteToString())
            return shp
        else:
            raise ValueError("Cache entry not found")
    def is_empty(self):
        return len(os.listdir(self.dirname)) == 0
    def __getattr__(self, k):
        return getattr(self.file, k)

def open(fn):
    return c(fn)
