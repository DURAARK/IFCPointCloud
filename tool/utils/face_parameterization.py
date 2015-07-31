import sys

import OCC.gp
import OCC.BRep
import OCC.TopoDS
import OCC.TopAbs
import OCC.TopExp
import OCC.ShapeAnalysis
import OCC.BRepTopAdaptor

# Helper class to help with point projection on face surface 
# and cache classification info to speed up point projection.
class face(object):

    # Attributes that can be dynamically accessed by means
    # of __getattr__(), they are computed once by means of
    # the getter functions named according to the attribute.
    attrs = {'surface', 'parametrizer', 'uv_bounds', 'classifier', 'origin', 'normal', 'axis', 'reversed'}
    
    def __init__(self, face):
        self.face = face
        
    def get_reversed(self):
        return self.face.Orientation() == OCC.TopAbs.TopAbs_REVERSED
        
    def get_surface(self):
        return OCC.BRep.BRep_Tool.Surface(self.face)
        
    def get_parametrizer(self):
        return OCC.ShapeAnalysis.ShapeAnalysis_Surface(self.surface)
        
    def get_uv_bounds(self):
        exp = OCC.TopExp.TopExp_Explorer(self.face, OCC.TopAbs.TopAbs_VERTEX)
        umin, umax, vmin, vmax = 1e9, -1e9, 1e9, -1e9
        while exp.More():
            vertex = OCC.TopoDS.topods.Vertex(exp.Current())
            xyz = OCC.BRep.BRep_Tool.Pnt(vertex)
            u, v = self.uvd(xyz)[:2]
            if u < umin: umin = u
            if u > umax: umax = u
            if v < vmin: vmin = v
            if v > vmax: vmax = v
            exp.Next()
        return ((umin, vmin), (umax, vmax))
        
    def get_classifier(self):
        return OCC.BRepTopAdaptor.BRepTopAdaptor_FClass2d(self.face, 1e-9)
        
    def get_origin(self):
        (a,b), (c,d) = self.uv_bounds
        return self.parametrizer.Value(OCC.gp.gp_Pnt2d(a,b))
        
    def get_axis(self):
        
        (a,b), (c,d) = self.uv_bounds
        o = self.parametrizer.Value(OCC.gp.gp_Pnt2d(a,b))
        p = self.parametrizer.Value(OCC.gp.gp_Pnt2d(c,b))
        ax = OCC.gp.gp_Dir(p.XYZ() - o.XYZ())
        
        surf = self.surface
        assert surf.GetObject().DynamicType().GetObject().Name() == "Geom_Plane"
        plane = OCC.Geom.Handle_Geom_Plane.DownCast(surf).GetObject()
        bx = plane.Pln().XAxis().Direction()

        return ax
        
    def get_normal(self):
        surf = self.surface
        assert surf.GetObject().DynamicType().GetObject().Name() == "Geom_Plane"
        plane = OCC.Geom.Handle_Geom_Plane.DownCast(surf).GetObject()
        ax = plane.Axis().Direction()
        # if self.reversed: ax.Reverse()
        # if not plane.Pln().Direct(): 
        #     ax.Reverse()
        return ax
        
    def uvd(self, xyz):
        uv = self.parametrizer.ValueOfUV(xyz, 1e-9)
        u,v,d = uv.X(), uv.Y(), None
        
        st = self.classifier.Perform(uv)
        if st == OCC.TopAbs.TopAbs_IN or st == OCC.TopAbs.TopAbs_ON:
            projected = self.parametrizer.Value(uv)
            d = projected.Distance(xyz)
            
            try:
                dir = OCC.gp.gp_Vec(xyz.XYZ() - projected.XYZ())
                if dir.XYZ().Dot(self.normal.XYZ()) < 0: d *= -1
            except:
                # TODO: For non-planar surfaces only the unsigned distance is computed
                pass

        return u, v, d
        
    def __getattr__(self, k):
        if k in face.attrs:
            x = getattr(self, "get_"+k)()
            setattr(self, k, x)
            return x
    
# A ordered collections of faces
class shape(object):
    def __init__(self, s):
        self.s, self.faces = s, []
        exp = OCC.TopExp.TopExp_Explorer(self.s, OCC.TopAbs.TopAbs_FACE)
        while exp.More():
            self.faces.append(face(OCC.TopoDS.topods.Face(exp.Current())))
            exp.Next()
