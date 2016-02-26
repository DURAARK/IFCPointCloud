import os
import sys
import json
import time
import struct
import random
import platform
import argparse
import itertools

if platform.architecture()[0] == '64bit':
    import ifcopenshell_x64
    import ifcopenshell_x64.geom
    ifcopenshell = ifcopenshell_x64
else:
    import ifcopenshell
    import ifcopenshell.geom

import OCC
import OCC.gp
import OCC.Bnd
import OCC.TopExp
import OCC.TopAbs
import OCC.BRepBndLib
import OCC.ShapeAnalysis
import OCC.BRepTopAdaptor
import OCC.BRepGProp
import OCC.GProp
import OCC.BRepBuilderAPI

import numpy as np
import numpy.linalg as linalg

from collections import defaultdict, namedtuple

from rtree import index

from utils import pcd_file, geometry_cache, face_parameterization, face_information

POINT_PASSES = 100
MAX_DISTANCE = (0.01, 0.05, 0.25)
MAX_MAX_DISTANCE = max(MAX_DISTANCE)
ENTITIES = {"IfcWall", "IfcSlab", "IfcStairFlight", "IfcStair"}

argp = argparse.ArgumentParser()
argp.add_argument('--use-display', required=False, action='store_true')
argp.add_argument('ifc_file')
argp.add_argument('pcd_files')
argp.add_argument('alignment_matrix')
ns = argp.parse_args()

USE_DISPLAY, IFC_FILE, PCD_FILES, ALIGNMENT_MATRIX_FILE = \
    ns.use_display, ns.ifc_file, ns.pcd_files, ns.alignment_matrix
    
if USE_DISPLAY:
    from PyQt4 import QtGui
    
ALIGNMENT_MATRIX = np.array(json.load(open(ALIGNMENT_MATRIX_FILE)))

IFC_FILE_BASE = os.path.splitext(os.path.basename(IFC_FILE))[0]

INTERMEDIATE_DIR = 'intermediate_files'

OUTPUT_DIR = os.path.join(INTERMEDIATE_DIR, IFC_FILE_BASE, "associated_points")

if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)

INDEX_FILE = os.path.join(INTERMEDIATE_DIR, '%s_%d_%d_3d_index' % (IFC_FILE_BASE, time.time(), os.getpid()))

index_properties = index.Property()
index_properties.dimension = 3
spatial_index = index.Index(INDEX_FILE, properties=index_properties)
index_counter = 1

if USE_DISPLAY:
    # Initialize a graphical display window
    occ_display = ifcopenshell.geom.utils.initialize_display()

    default_color = (OCC.Quantity.Quantity_Color(.5, .5, .5, OCC.Quantity.Quantity_TOC_RGB), 0.9)
    RED = OCC.Quantity.Quantity_Color(.8, .1, .1, OCC.Quantity.Quantity_TOC_RGB)
    BLUE = OCC.Quantity.Quantity_Color(.1, .1, .8, OCC.Quantity.Quantity_TOC_RGB)
    colors = {
        'IfcWindow':           (OCC.Quantity.Quantity_Color(.5, .6, .7, OCC.Quantity.Quantity_TOC_RGB), 0.95),
        'IfcSlab':             (OCC.Quantity.Quantity_Color(.4, .4, .4, OCC.Quantity.Quantity_TOC_RGB), 0.9),
        'IfcDoor':             (OCC.Quantity.Quantity_Color(.6, .5, .4, OCC.Quantity.Quantity_TOC_RGB), 0.95)
    }

    # Render a transparent shape for the product into the display
    # window, based on entity type and colors defined above
    def process_shape(product, shape):
        color, transparency = colors.get(product.is_a(), default_color)
        disp = occ_display.DisplayShape(shape, color=color, material=OCC.Graphic3d.Graphic3d_NOM_SATIN, update=True)
        if transparency:
            ifcopenshell.geom.utils.set_shape_transparency(disp, transparency)
else:
    # An identity function is created for the function that
    # would otherwise render the shape to the display window
    process_shape = lambda *args: None

face_dict = {}
assoc_dict = {}

# Add axis-aligned bounding box information of
# the face to the spatial indexing structure.
def process_face(instance_id, face_id, bbox, face_shape):
    global index_counter
    object_key = (instance_id, face_id)
    spatial_index.insert(index_counter, bbox, object_key)
    
    face_dict[object_key] = face_parameterization.face(face_shape)
    index_counter += 1

ifc_file = None

# Compute a precise bounding box, disregard tolerances in the BRep model
def precise_bbox(s):
    exp = OCC.TopExp.TopExp_Explorer(s, OCC.TopAbs.TopAbs_VERTEX)
    x1,y1,z1,x2,y2,z2 = 1e9,1e9,1e9,-1e9,-1e9,-1e9
    while exp.More():
        p = OCC.BRep.BRep_Tool.Pnt(OCC.TopoDS.topods.Vertex(exp.Current()))
        x,y,z = p.X(), p.Y(), p.Z()
        if x < x1: x1 = x
        if y < y1: y1 = y
        if z < z1: z1 = z
        if x > x2: x2 = x
        if y > y2: y2 = y
        if z > z2: z2 = z
        exp.Next()
    return x1,y1,z1,x2,y2,z2
    
def load_ifc_file_per_face():
    global ifc_file
    # Open the IFC file using IfcOpenShell
    ifc_file = geometry_cache.open(IFC_FILE)
    cache_initialized = not ifc_file.is_empty()
    # Display the geometrical contents of the file using Python OpenCascade
    products = ifc_file.by_type("IfcProduct")
    for i, product in enumerate(products):
        if USE_DISPLAY:
            QtGui.QApplication.processEvents()
        if any(product.is_a(e) for e in ENTITIES) and product.Representation:
            try:
                shape = ifc_file.create_shape(product, force_from_cache=cache_initialized)
            except: continue
            process_shape(product, shape)
            # every face is stored separately in the indexing structure
            exp = OCC.TopExp.TopExp_Explorer(shape, OCC.TopAbs.TopAbs_FACE)
            face_index = 0
            
            while exp.More():
                bbox = OCC.Bnd.Bnd_Box()
                face = OCC.TopoDS.topods.Face(exp.Current())
                OCC.BRepBndLib.brepbndlib_Add(face, bbox)
                process_face(product.id(), face_index, precise_bbox(face), face)
                face_index += 1
                exp.Next()

                
ASSOCIATED_OUTPUT_FN = os.path.join(OUTPUT_DIR, "points-subset-%d.bin")

def load_and_associate_point_cloud(scan_number, point_pass):
    global associated_point_file
    pcdf = pcd_file.open(PCD_FILES % scan_number, mod=(POINT_PASSES, point_pass))
    sys.stderr.write("\r[pass %d] processing scan %d  " % (j, i))
    for x,y,z in pcdf:
        a = np.dot(ALIGNMENT_MATRIX, np.array((x,y,z,1.)))
        x,y,z = a[:3]
        
        associated_point_file.write(struct.pack("@ddd", x, y ,z))
        p = OCC.gp.gp_Pnt(x,y,z)
        
        # Progressively increase the bounding box size
        for bbs in MAX_DISTANCE:
        
            # bbs represents half of the bounding cube size used
            # for intersection testing with our indexing structure
        
            xs = list(spatial_index.intersection((x-bbs, y-bbs, z-bbs, x+bbs, y+bbs, z+bbs), objects=True))
            
            U, V, D, O, F = None, None, 1e9, None, None
            umin, vmin, umax, vmax = (None,) * 4
            for x_ in xs:
                f = face_dict[x_.object]
                # project point on surface and get u,v and distance from surface
                try:
                    # TODO: fails for non-planar surfs because the normal is not evaluated at UV
                    u,v,d = f.uvd(p)
                except:
                    continue
                
                if d is not None and abs(d) < bbs and d < D:
                    (umin, vmin), (umax, vmax) = f.uv_bounds
                    if (umax-umin) > 1e-9 and (vmax - vmin) > 1e-9:
                        U,V,D,O,F = u,v,d,x_.object,f
                    
                    
            if U and V:
                un, vn, dn = (U-umin) / (umax-umin), (V-vmin) / (vmax-vmin), D / MAX_MAX_DISTANCE
                associated_point_file.write(struct.pack("@qqqddd", 1, x_.object[0], x_.object[1], un, vn, dn))           
                
            if USE_DISPLAY:
                if U and V:
                    # color vertex based on U (red), V (green), Distance (blue)
                    un, vn, dn = (U-umin) / (umax-umin), (V-vmin) / (vmax-vmin), D / (bbs * 2.) + 0.5
                    clr = OCC.Quantity.Quantity_Color(un, vn, dn, OCC.Quantity.Quantity_TOC_RGB)
                else:
                    clr = default_color[0]
                vertex = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeVertex(p).Vertex()            
                occ_display.DisplayShape(vertex, color=clr, material=OCC.Graphic3d.Graphic3d_NOM_SATIN, update=random.random() < 0.01)
                QtGui.QApplication.processEvents()
                
            if U and V: break # break in order for the else block not to kick in
        else:
            # Normally not a big fan of for...else but it does the job here
            associated_point_file.write(struct.pack("@q", 0))
            
print "Loading IFC file"
load_ifc_file_per_face()

assoc_files = []

try:
    for j in range(POINT_PASSES):
        # This will likely go out of memory at some point, but the script is re-entrant.
        fn, lock_fn = (ASSOCIATED_OUTPUT_FN % j), (ASSOCIATED_OUTPUT_FN + ".lock") % j
        if not os.path.exists(fn) and not os.path.exists(lock_fn):
            with open(lock_fn, "wb"): pass
            associated_point_file = open(fn, "wb")
            for i in itertools.count():
                try: load_and_associate_point_cloud(i, j)
                except IOError: break
                associated_point_file.flush()
            os.unlink(lock_fn)
            associated_point_file.close()
            assoc_files.append(fn)
except KeyboardInterrupt:
    print "\n[Ctrl]+[C] / Assocation interrupted by user"
        
sys.stderr.write("\nDone :)\n")

print "Writing face information"
with open(os.path.join(INTERMEDIATE_DIR, IFC_FILE_BASE, "face_information.json"), "w") as f:
    f.write(face_information.obtain(IFC_FILE))

del spatial_index
if os.path.exists(INDEX_FILE + ".dat"):
    os.unlink(INDEX_FILE + ".dat")
    os.unlink(INDEX_FILE + ".idx")
    
if USE_DISPLAY:
    ifcopenshell.geom.utils.main_loop()
