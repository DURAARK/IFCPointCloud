import math
import argparse
import operator
import itertools
import functools

import h5py

import numpy as np

import OCC
import OCC.BRep
import OCC.Display.SimpleGui
import OCC.TopExp
import OCC.BRepGProp
import OCC.GProp

from PyQt4 import QtGui

from utils import geometry_cache, ifc_math

argp = argparse.ArgumentParser()
argp.add_argument('hdf_file')
argp.add_argument('--show-model', required=False, action='store_true')
argp.add_argument('--show-unassociated', required=False, action='store_true')
argp.add_argument('coloring', choices=['color_deviation', 'color_parameterization'])
argp.add_argument('lod', type=int)
ns = argp.parse_args()

# Configurable options
SHOW_MODEL = ns.show_model
SHOW_UNASSOCIATED = ns.show_unassociated
COLOR_D = ns.coloring = 'color_deviation'
LOD = ns.lod
GRID_DOWNSAMPLE = int(math.sqrt(LOD))

# Hard-coded options
MIN_AREA = 2.
MIN_POINTS = 256
NUM_COLOR_SAMPLES = 128
NUM_COLOR_SAMPLES_SQRT = int(math.sqrt(NUM_COLOR_SAMPLES))

# Initialize display
handle, main_loop, add_menu, add_function_to_menu = OCC.Display.SimpleGui.init_display()
handle.set_bg_gradient_color(*(255,)*6)

GRAY = OCC.Quantity.Quantity_Color(0.5, 0.5, 0.5, OCC.Quantity.Quantity_TOC_RGB)

def render_depthcolor(c, xyzs):
    def color_from_w(c):
        r = lambda x: max(-x, 0.)
        g = lambda x: 1. - abs(x)
        b = lambda x: max( x, 0.)
        x = float(c) / NUM_COLOR_SAMPLES
        return r(x), g(x), b(x)
    bu = OCC.BRep.BRep_Builder()
    sh = OCC.TopoDS.TopoDS_Compound()
    bu.MakeCompound(sh)
    r, g, b = color_from_w(c)
    clr = OCC.Quantity.Quantity_Color(r, g, b, OCC.Quantity.Quantity_TOC_RGB)
    for xyz in xyzs:
        vertex = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeVertex(OCC.gp.gp_Pnt(*map(float, tuple(xyz)))).Vertex()
        bu.Add(sh, vertex)
    handle.DisplayShape(sh, color=clr)
    
def render_uvcolor(uvc, xyzs):
    bu = OCC.BRep.BRep_Builder()
    sh = OCC.TopoDS.TopoDS_Compound()
    bu.MakeCompound(sh)
    clr = OCC.Quantity.Quantity_Color(uvc[0], uvc[1], 0., OCC.Quantity.Quantity_TOC_RGB)
    for xyz in xyzs:
        vertex = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeVertex(OCC.gp.gp_Pnt(*map(float, tuple(xyz)))).Vertex()
        bu.Add(sh, vertex)
    handle.DisplayShape(sh, color=clr)
    
def render_constcolor(clr, xyzs):
    bu = OCC.BRep.BRep_Builder()
    sh = OCC.TopoDS.TopoDS_Compound()
    bu.MakeCompound(sh)
    for xyz in xyzs:
        vertex = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeVertex(OCC.gp.gp_Pnt(*map(float, tuple(xyz)))).Vertex()
        bu.Add(sh, vertex)
    handle.DisplayShape(sh, color=clr)
    
class hdf5_population(object):

    INST_DTYPE = np.dtype([('_HDF5_dataset_index_', '<i2'), ('_HDF5_instance_index_', '<i4')])

    class hdf5_entity_instance(object):
        def __init__(self, file, data):
            self.file = file
            self.data = data
        def __getattr__(self, attr):
            attr_value = self.data[attr]
            if isinstance(attr_value, h5py.h5r.Reference):
                attr_value = self.file.file[attr_value]
            if attr_value.dtype == hdf5_population.INST_DTYPE:
                if attr_value.shape == ():
                    return self.file.resolve(attr_value)
                else:
                    return map(self.file.resolve, attr_value)
            elif getattr(attr_value, 'dtype', None) == "<f4":
                # This is necessary in case of single precision HDF
                # files, because they do not have a common ancestor
                # with Python's float type, which is a double.
                dims = len(attr_value.shape)
                if dims == 0: return np.array([attr_value], dtype="f8")[0]
                else: return np.array(attr_value, dtype="f8")                    
            else:
                return attr_value

    def __init__(self, file, population_group_name=None):
        self.file = file
        if population_group_name is None:
            population_group_name = list(set(map(operator.itemgetter(0), file.items())) - {u'IFC2X3_encoding'})[0]
        self.root_group = self.file[population_group_name]
        self.dataset_names = self.root_group.attrs["iso_10303_26_data_set_names"]
        
    def __getattr__(self, entity):
        try:
            return map(
                functools.partial(hdf5_population.hdf5_entity_instance, self), 
                self.root_group["%s_objects" % entity]["%s_instances" % entity]
            )
        except: return []
        
    def resolve(self, idx):
        ds_idx, inst_idx = idx
        if ds_idx > 0:
            entity = self.dataset_names[ds_idx]
            return hdf5_population.hdf5_entity_instance(
                self, 
                self.root_group["%s_objects" % entity]["%s_instances" % entity][inst_idx]
            )
        
        
if SHOW_MODEL:
    all_faces = []
    ifc_file = cache.open("Plan3D_Haus30_PREVIEW_NEW.ifc")
    for p in ifc_file.by_type("IfcProduct"):
        if p.is_a("IfcSpace") or p.is_a("IfcOpeningElement"): continue
        try: shape = ifc_file.create_shape(p, force_from_cache=True)
        except: continue
        all_faces.append(shape)
        exp = OCC.TopExp.TopExp_Explorer(shape, OCC.TopAbs.TopAbs_FACE)
        while exp.More():
            face = OCC.TopoDS.topods.Face(exp.Current())
            props = OCC.GProp.GProp_GProps()                    
            OCC.BRepGProp.brepgprop_SurfaceProperties(face, props)
            exp.Next()

    for i, f in enumerate(all_faces):
        if i % 64 == 0:
            handle.FitAll()
            QtGui.QApplication.processEvents()
        handle.Context.SetTransparency(handle.DisplayShape(f, color=GRAY), 0.8)
        
# Open the file
f = hdf5_population(h5py.File(ns.hdf_file, "r"))

if SHOW_UNASSOCIATED:
    try:
        render_constcolor(GRAY, f.IfcCartesianPointList3D[0].CoordList[LOD])
    except: pass

for inst in f.IfcRelAssignsToProduct:
    pc_elem = inst.RelatedObjects[0]

    # TODO:
    # product = inst.RelatingProduct
    # placement = pc_elem.ObjectPlacement
    # matrix = local_placement(placement)
    
    definition = pc_elem.Representation
    shaperep = definition.Representations[0]
    
    point_cloud = shaperep.Items[0]
    coords = point_cloud.Coordinates
    if not hasattr(coords, 'Surface'): # GridOffsetList
        if len(coords.Offsets) < MIN_POINTS: continue
        
        grid = coords.Grid
        
        if hasattr(grid, "Representation"): # IfcGrid
            trimmed_surf = grid.Representation.Representations[0].Items[0]
            u0 = grid.UAxes[0].Pnt.Coordinates[0]
            v0 = grid.VAxes[0].Pnt.Coordinates[1]
            distances = lambda arr: np.array(map(lambda a: a.data["Distance"], arr[1:]))
            us = np.concatenate(([0], distances(grid.UAxes))) + u0
            vs = np.concatenate(([0], distances(grid.VAxes))) + v0
        else: # IfcSurfaceGrid
            trimmed_surf = grid.Location.BasisSurface
            u0 = grid.Location.PointParameterU
            v0 = grid.Location.PointParameterV
            us = np.arange(grid.NumUAxes) * grid.USpacing
            vs = np.arange(grid.NumVAxes) * grid.VSpacing
            
        du, dv = trimmed_surf.U2 - trimmed_surf.U1, trimmed_surf.V2 - trimmed_surf.V1
        if du*dv < MIN_AREA: continue      
        
        offsets = coords.Offsets[:].flatten()
        
        is_int = str(offsets.dtype).startswith('uint')
        if is_int:
            non_nans = offsets != 255 # always uint8
            offsets = np.array(offsets, dtype='f')
            offsets = offsets / coords.SampleWidth * (coords.MaxOrthogonalOffset - coords.MinOrthogonalOffset) + coords.MinOrthogonalOffset
        else:        
            non_nans = offsets != 1.
        
        uvs = np.array(np.meshgrid(vs,us)[::-1]).T.flatten().reshape(-1,2)
        uvs *= du, dv
        uvds = np.c_[uvs, offsets]
        
        if GRID_DOWNSAMPLE > 1:
            maskx, masky = np.meshgrid(np.arange(len(us)),np.arange(len(vs)))
            mask = ((maskx % GRID_DOWNSAMPLE == 0) & (masky % GRID_DOWNSAMPLE == 0)).flatten()
            uvds = uvds[mask & non_nans]
        else:
            uvds = uvds[non_nans]
            
        max_offset, min_offset = coords.MaxOrthogonalOffset, coords.MinOrthogonalOffset
    else:
        if coords.Values.shape[0] < MIN_POINTS: continue
        
        trimmed_surf = coords.Surface
        du, dv = trimmed_surf.U2 - trimmed_surf.U1, trimmed_surf.V2 - trimmed_surf.V1
        if du*dv < MIN_AREA: continue
        
        uvds = coords.Values[::LOD,:]
        
        is_int = str(uvds.dtype).startswith('uint')
        if is_int:
            uvds = np.array(uvds, dtype='f') / 2.**16
            uvds[:,:-1] *= du, dv
        
        max_offset, min_offset = coords.WMaxOffset, coords.WMinOffset
        
    plane_matrix = ifc_math.axis2placement(trimmed_surf.BasisSurface.Position)
    
    if not SHOW_MODEL:
        pln = OCC.gp.gp_Pln(OCC.gp.gp_Ax3(
            OCC.gp.gp_Pnt(*trimmed_surf.BasisSurface.Position.Location.Coordinates),
            OCC.gp.gp_Dir(*trimmed_surf.BasisSurface.Position.Axis.DirectionRatios),
            OCC.gp.gp_Dir(*trimmed_surf.BasisSurface.Position.RefDirection.DirectionRatios)))
            
        f = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeFace(pln, 0., trimmed_surf.U2, 0., trimmed_surf.V2).Face()
        handle.Context.SetTransparency(handle.DisplayShape(f, color=GRAY), 0.8)
        
        exp = OCC.TopExp.TopExp_Explorer(f, OCC.TopAbs.TopAbs_WIRE)
        while exp.More():
            handle.Context.SetTransparency(handle.DisplayShape(exp.Current(), color=GRAY), 0.2)
            exp.Next()
            
    uvds = np.c_[uvds, np.ones(len(uvds))]
    uvds[:,2] *= max_offset
    
    if COLOR_D:
        xyzs = np.array((np.asmatrix(plane_matrix) * np.asmatrix(uvds).T).T)
        xyzs = np.c_[xyzs[:,:-1], np.array(uvds[:,2] * (NUM_COLOR_SAMPLES / max_offset), dtype='i')]
        del uvds
        for i in range(-NUM_COLOR_SAMPLES, NUM_COLOR_SAMPLES+1):
            mask = np.repeat(xyzs[:,-1:] == i, 4, axis=1)
            masked = xyzs[mask].flatten().reshape((-1, 4))
            render_depthcolor(i, masked[:,:-1])
    else:
        xyzs = np.array((np.asmatrix(plane_matrix) * np.asmatrix(uvds).T).T)
        uvds[:,0:2] /= trimmed_surf.U2, trimmed_surf.V2
        xyzs = np.c_[xyzs[:,:-1], np.array(uvds[:,0:2] * NUM_COLOR_SAMPLES_SQRT, dtype='i')]
        del uvds
        for ij in itertools.product(range(NUM_COLOR_SAMPLES_SQRT+1), repeat=2):
            mask = np.apply_along_axis(lambda a: (a[-2:] == ij).all(), 1, xyzs).flatten().reshape((-1, 1))
            mask = np.repeat(mask, 5, axis=1)
            masked = xyzs[mask].flatten().reshape((-1, 5))
            render_uvcolor(map(float, np.array(ij, dtype='f') / NUM_COLOR_SAMPLES_SQRT), masked[:,:-2])
        
    QtGui.QApplication.processEvents()
    
    handle.FitAll()

main_loop()
