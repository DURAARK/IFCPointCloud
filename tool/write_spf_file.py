import os
import sys
import uuid
import time
import json
import struct
import platform
import argparse
import functools
import itertools
import subprocess

import numpy as np

if platform.architecture()[0] == '64bit':
    import ifcopenshell_x64
    import ifcopenshell_x64.geom
    ifcopenshell = ifcopenshell_x64
else:
    import ifcopenshell
    import ifcopenshell.geom

from collections import namedtuple, defaultdict

from utils import point_reader

# from ifc_math import *

ALL_POINTS_AS_UNASSOCIATED = '--all_points_as_unassociated' in sys.argv
ONLY_INCLUDE_UNASSOCIATED = '--only_include_unassociated' in sys.argv
ONLY_INCLUDE_UNASSOCIATED_SEGMENTED = '--only_include_unassociated_segmented' in sys.argv
SINGLE_ASSIGNMENT = False # due to h5py vlen limitations

INTERMEDIATE_DIR = "intermediate_files"

assert sum([ALL_POINTS_AS_UNASSOCIATED, ONLY_INCLUDE_UNASSOCIATED]) <= 1

def get_argparser():
    argp = argparse.ArgumentParser()
    argp.add_argument('--all_points_as_unassociated', required=False, action='store_true')
    argp.add_argument('--only_include_unassociated', required=False, action='store_true')
    argp.add_argument('ifc_file')
    return argp

if not ALL_POINTS_AS_UNASSOCIATED and not ONLY_INCLUDE_UNASSOCIATED and not ONLY_INCLUDE_UNASSOCIATED_SEGMENTED:
    argp = get_argparser()
    argp.add_argument('grid_spacing', type=float)
    argp.add_argument('store_unassociated', choices=['without_unassociated', 'unassociated_cartesian', 'unassociated_segmented'])
    argp.add_argument('param_format',       choices=['continuous_param', 'discrete_param']       )
    argp.add_argument('grid_format',        choices=['legacy_grid', 'raster_grid']               )
    ns = argp.parse_args()
    
    GRID_SPACING = ns.grid_spacing
    STORE_UNASSOCIATED = ns.store_unassociated != 'without_unassociated'
    STORE_UNASSOCIATED_SEGMENTED = ns.store_unassociated == 'unassociated_segmented'
    USE_DISCRETE_PARAM = ns.param_format == 'discrete_param'
    USE_LEGACY_GRID = ns.grid_format == 'legacy_grid'
    
    flags = [
        '',
        ('grid_%0.3f' % GRID_SPACING) if GRID_SPACING > 1e-5 else "no_grid",
        ns.store_unassociated,
        ns.grid_format,
        ns.param_format
    ]
else: 
    argp = get_argparser()
    ns = argp.parse_args()
    
    if ALL_POINTS_AS_UNASSOCIATED: 
        flag_name = 'without-association'
        STORE_UNASSOCIATED = True
        # This has to be set to false, because segmentation has not been
        # performed for the unassociated points.
        STORE_UNASSOCIATED_SEGMENTED = False
        GRID_SPACING = 1e-9
    elif ONLY_INCLUDE_UNASSOCIATED: 
        flag_name = 'only-unassociated'
        STORE_UNASSOCIATED = True
        STORE_UNASSOCIATED_SEGMENTED = False
        GRID_SPACING = 1e-9
    elif ONLY_INCLUDE_UNASSOCIATED_SEGMENTED: 
        flag_name = 'only-unassociated-segmented'
        STORE_UNASSOCIATED = True
        STORE_UNASSOCIATED_SEGMENTED = True
        GRID_SPACING = 1e-9
    else: raise Exception()
    
    flags = [
        '', 
        flag_name
    ]
    
IFC_FN        = ns.ifc_file
IFC_FILE_BASE = os.path.splitext(os.path.basename(IFC_FN))[0]
ASSOC_FN      = "intermediate_files/%s/associated_points/points-subset-%%d.bin" % IFC_FILE_BASE
SEGMENTED_FN  = "intermediate_files/%s/associated_points/segmented/plane-%%d.bin" % IFC_FILE_BASE
REMAINING_FN  = "intermediate_files/%s/associated_points/segmented/remaining-unassociated-points.bin" % IFC_FILE_BASE
IFC_FN_FLAGS  = "generated_files/spf/" + IFC_FILE_BASE + "-".join(flags)
IFC_PC_FN     = IFC_FN_FLAGS + ".ifc"
IFC_FACE_INFO = os.path.join(INTERMEDIATE_DIR, IFC_FILE_BASE, "face_information.json")
IFC_PRODUCTS  = os.path.join(INTERMEDIATE_DIR, IFC_FILE_BASE, "face_information.json")

if GRID_SPACING > 1e-5:
    from utils import point_rasterizer, face_bounds
    rasterizer = point_rasterizer.rasterizer(face_bounds.get_mapping(IFC_FN))

if not os.path.exists(os.path.dirname(IFC_FN_FLAGS)): os.makedirs(os.path.dirname(IFC_FN_FLAGS))

decomp_relations = defaultdict(list)
t0 = int(time.time())

# add range=[0] to speed things up
datasets, global_ds = point_reader.read(ASSOC_FN, associate=not ALL_POINTS_AS_UNASSOCIATED) 

print "\nReading face parametrization data"
shapes = dict((tuple(k),v) for k, v in json.load(open(IFC_FACE_INFO)))

print "Reading original SPF file"
f = ifcopenshell.open(IFC_FN)

print "Creating new IFC representations"    

new_guid = lambda: ifcopenshell.guid.compress(uuid.uuid4().hex)

tupelize = lambda np_array: tuple(map(lambda a: tuple(float(x) for x in a), np_array))
tupelize1d = lambda np_array: tuple(map(float, np_array))
tupelize2d = lambda np_array: tuple(map(tupelize1d, np_array))

tupelize1d_map = lambda fn, np_array: tuple(map(fn, np_array))
tupelize2d_map = lambda fn, np_array: tuple(map(functools.partial(tupelize1d_map, fn), np_array))

discretisize0d = lambda p: lambda f: int(f*(p-1)) if f == f else 0
discretisize1d = lambda p: lambda a: tuple(map(discretisize0d(p), a))
discretisize2d = lambda p: lambda a: tuple(map(discretisize1d(p), a))


context = f.createIfcGeometricRepresentationContext("PointCloud", "PointCloud", 3, 1e-3, f.createIfcAxis2Placement3D())
f.by_type("IfcProject")[0].RepresentationContexts += (context,)
    
person = f.createIfcPerson(GivenName="Thomas", FamilyName="Krijnen")
organisation = f.createIfcOrganization(Name="Eindhoven University of Technology")
pando = f.createIfcPersonAndOrganization(person, organisation)
app = f.createIfcApplication(organisation, "0.2", "IFC+PointCloud prototypical implementation", "IfcPointCloud app v0.2")
owner_hist = f.createIfcOwnerHistory(pando, app, None, "ADDED", t0, pando, app, t0)

root_pc_elem = f.createIfcPointCloudElement(new_guid(), owner_hist, "Container for point clouds")

if STORE_UNASSOCIATED:
    if STORE_UNASSOCIATED_SEGMENTED:
        reps = []
        # idxs=range(64) to speed things up
        for m, points in point_reader.read_param_points(SEGMENTED_FN):
            if m is not None:
                p = f.createIfcCartesianPoint(tupelize1d( m[0:3,3]))
                z = f.createIfcCartesianPoint(tupelize1d(-m[0:3,2]))
                x = f.createIfcCartesianPoint(tupelize1d( m[0:3,0]))
                ax = f.createIfcAxis2Placement3d(p, z, x)
                pln = f.createIfcPlane(ax)
            
            points = np.array(list(points))
            
            us, vs, ws = points.T            
            u1, u2 = np.min(us), np.max(us)
            v1, v2 = np.min(vs), np.max(vs)
            w1, w2 = np.min(ws), np.max(ws)
            
            surf = f.createIfcRectangularTrimmedSurface(pln, u1, v1, u2, v2)
            
            # Move to trimmed surface coords
            us -= u1
            vs -= v1
            
            if USE_DISCRETE_PARAM:
                ws -= w1
                
                # Clamp to [0,1] range
                us /= (u2 - u1)
                vs /= (v2 - v1)
                ws /= (w2 - w1)
                
                us = discretisize1d(2**16)(us)
                vs = discretisize1d(2**16)(vs)
                ws = discretisize1d(2**16)(ws)
                points = f.createIfcDiscreteParameterValueList(w1, w2, surf, tuple(zip(us, vs, ws)), 2**16, 2**16, 2**16)
            else:
                points = f.createIfcContinuousParameterValueList(w1, w2, surf, tupelize(points))
            
            cloud = f.createIfcPointCloud(Coordinates=points, Attributes=())
            reps.append(cloud)
        
        # Add the remaining unsegmented points
        points = np.fromfile(REMAINING_FN).reshape((-1,3))
        points = points[::32]
        points = f.createIfcCartesianPointList3D(tupelize(points))
        cloud = f.createIfcPointCloud(Coordinates=points, Attributes=())        
        reps.append(cloud)
        
        elem = f.createIfcPointCloudElement(new_guid(), owner_hist, "Unassociated point clouds")
        rep = f.createIfcShapeRepresentation(context, "PointCloud", "PointCloud", reps)
        elem.Representation = f.createIfcProductDefinitionShape(Representations=(rep,))
    else:
        points = f.createIfcCartesianPointList3D(tupelize(global_ds.merge()))
        cloud = f.createIfcPointCloud(Coordinates=points, Attributes=())
        elem = f.createIfcPointCloudElement(new_guid(), owner_hist, "Unassociated point clouds")
        rep = f.createIfcShapeRepresentation(context, "PointCloud", "PointCloud", (cloud,))
        elem.Representation = f.createIfcProductDefinitionShape(Representations=(rep,))
        
        decomp_relations[root_pc_elem.id()].append(elem)

del global_ds

if ONLY_INCLUDE_UNASSOCIATED or ONLY_INCLUDE_UNASSOCIATED_SEGMENTED: datasets = {}

# Carefully iterate over a copy of the keys so that elements can be freed as IFC representations of them are created in order not to consume more memory
# Additionally: <eid, fid> pairs are sorted so that IfcRelAssignsToProduct can be grouped per RelatingProduct

previous_eid = None
related_objects = []

def write_assignment_relationship(objects, product):
    f.createIfcRelAssignsToProduct(new_guid(), owner_hist, RelatedObjects=objects, RelatingProduct=product)

for dataset_key in sorted(datasets.keys()):

    # oops.. reversed key in get_face_info.py:20
    eid, fid = dataset_key
    shape_key = fid, eid
    if shape_key not in shapes:
        print "No surface matched for %dth face of:\n%r" % (fid, f[eid])
        continue
    o,z,x,u1,v1,u2,v2 = shapes[shape_key]
    
    if not SINGLE_ASSIGNMENT:
        if eid != previous_eid and len(related_objects):
            write_assignment_relationship(related_objects, f[previous_eid])
            related_objects[:] = []
    previous_eid = eid
    
    # TODO
    # plane_matrix = a2p(o,z,x)
    # product_matrix = obtain inverse product placement matrix
    # rel_plane_matrix = np.dot(product_matrix, plane_matrix)
    
    # Initially I was under the impression that IfcRectangularTrimmedSurface would normalize the 
    # parametric range of the base surface to [(0,0);(1,1)]. However, this is not the case, but the
    # association code writes parameter values in this range. Hence they are multiplied.
    
    ax = f.createIfcAxis2Placement3D(
            f.createIfcCartesianPoint(o),
                 f.createIfcDirection(z),
                 f.createIfcDirection(x))
    pln = f.createIfcPlane(ax)
    surf = f.createIfcRectangularTrimmedSurface(pln, u1, v1, u2, v2)

    f2d = lambda f: lambda a: f(map(f, a))
    min2d = f2d(np.nanmin)
    max2d = f2d(np.nanmax)
            
    numpy_array = datasets[dataset_key].merge()
    
    du = max(numpy_array[:,0]) - min(numpy_array[:,0])
    dv = max(numpy_array[:,1]) - min(numpy_array[:,1])
    density = float(len(numpy_array)) / (du * dv)
    
    try_rasterization = GRID_SPACING > 1e-5 and not np.isinf(density) and len(numpy_array) > 4 # and density > (0.1 / (GRID_SPACING * GRID_SPACING)) 
    
    is_smaller = big_enough = percentage_nan_low = False
    if try_rasterization:
        rasterized = rasterizer.rasterize(eid, fid, numpy_array, grid_spacing=GRID_SPACING)
        is_smaller = (rasterized.values.nbytes < numpy_array.nbytes // 2)
        num_non_nans = sum(~np.isnan(rasterized.values.flatten()))
        num_nans = sum( np.isnan(rasterized.values.flatten()))
        big_enough = num_non_nans > 16  # rasterized.values.size > 16
        percentage_nan_low = num_non_nans > num_nans 
        
    if try_rasterization and is_smaller and big_enough and percentage_nan_low:
        # Some debugging information:
        # print 'Rasterized %d / %d' % (num_non_nans, num_non_nans + num_nans)
        # print '%d -> %d bytes' % (numpy_array.nbytes, rasterized.values.nbytes)
        
        if USE_LEGACY_GRID:
            grid_origin = f.createIfcCartesianPoint(tupelize1d(rasterized.grid.placement))
            first_u_axis = f.createIfcLine(grid_origin, f.createIfcVector(f.createIfcDirection((1., 0.)), 1.))
            first_v_axis = f.createIfcLine(grid_origin, f.createIfcVector(f.createIfcDirection((0., 1.)), 1.))
            uaxes, vaxes = [first_u_axis], [first_v_axis]
            for i in range(1, rasterized.grid.unum):
                uaxes.append(f.createIfcOffsetCurve2D(first_u_axis, i * rasterized.grid.uspacing))
            for i in range(1, rasterized.grid.vnum):
                vaxes.append(f.createIfcOffsetCurve2D(first_v_axis, i * rasterized.grid.vspacing))
            # TODO: Needs to be IfcPCurve, but is IFC4 entity not in pc schema
            # Therefore 'Hackhackhack' workaround is introduced to add the surface as a representation item on the grid
            rep = f.createIfcShapeRepresentation(context, "Hackhackhack", "Hackhackhack", (surf,))
            grid = f.createIfcGrid(new_guid(), owner_hist, Representation=rep, UAxes=uaxes, VAxes=vaxes)
            grid.Representation = f.createIfcProductDefinitionShape(Representations=(rep,))
        else:
            pu, pv = tupelize1d(rasterized.grid.placement)
            grid_origin = f.createIfcPointOnSurface(surf, pu, pv)
            g = rasterized.grid
            grid = f.createIfcSurfaceGrid(grid_origin, g.uspacing, g.vspacing, g.unum, g.vnum)
            
        # The rasterizer returns uncovered cells as NaN, in 
        # the schema these are written as a separate attribute
        mask = None
        values = rasterized.values
        nans = np.isnan(values)
        num_nans = np.sum(nans)
        min_value, max_value = map(float, (min2d(values), max2d(values)))
        if num_nans > 0:
            values[nans] = (min_value + max_value) / 2.
            mask = tupelize2d_map(bool, nans)
                
        if USE_DISCRETE_PARAM:
            values = discretisize2d(256)((values - min_value) / (max_value - min_value))
            points = f.createIfcDiscreteGridOffsetList(
                Grid                = grid, 
                MinOrthogonalOffset = min_value, 
                MaxOrthogonalOffset = max_value,
                CoverageMask        = mask,
                Offsets             = values, 
                SampleWidth         = 256
            )
        else:
            values = tupelize2d(values)
            # TODO: TYPO ALERT
            points = f.createIfcContiousGridOffsetList(
                Grid                = grid, 
                MinOrthogonalOffset = min_value, 
                MaxOrthogonalOffset = max_value, 
                CoverageMask        = mask,
                Offsets             = values
            )
    else:
        if USE_DISCRETE_PARAM:
            us, vs, ws = numpy_array.T
            # us *= (u2-u1) NB: DO NOT MULTIPLY.. should be in
            # vs *= (v2-v1)     range [0,1) for discretization
            us = discretisize1d(2**16)(us)
            vs = discretisize1d(2**16)(vs)
            wmin, wmax = map(float, (min(ws), max(ws)))
            ws = discretisize1d(2**16)((ws - wmin) / (wmax - wmin))
            points = f.createIfcDiscreteParameterValueList(wmin, wmax, surf, tuple(zip(us, vs, ws)), 2**16, 2**16, 2**16)
        else:
            numpy_array[:,:-1] *= (u2-u1), (v2-v1)
            wmin, wmax = map(float, (min(numpy_array[:,2]), max(numpy_array[:,2])))
            points = f.createIfcContinuousParameterValueList(wmin, wmax, surf, tupelize(numpy_array))
    
    cloud = f.createIfcPointCloud(Coordinates=points, Attributes=())
    elem = f.createIfcPointCloudElement(new_guid(), owner_hist, "Associated points")
    rep = f.createIfcShapeRepresentation(context, "PointCloud", "PointCloud", (cloud,))
    elem.Representation = f.createIfcProductDefinitionShape(Representations=(rep,))
    # TODO
    # elem.ObjectPlacement = f[eid].ObjectPlacement
    if SINGLE_ASSIGNMENT:
        write_assignment_relationship((elem,), f[eid])
    else:
        related_objects.append(elem)
    decomp_relations[root_pc_elem.id()].append(elem)
    
    del datasets[dataset_key]

if len(related_objects):
    write_assignment_relationship(related_objects, f[eid])
    
for roid, objects in decomp_relations.iteritems():
    f.createIfcRelDecomposes(new_guid(), owner_hist, RelatingObject=f[roid], RelatedObjects=objects)
            
print "Writing SPF file"
f.write(IFC_PC_FN)

print "Done :)"

del f