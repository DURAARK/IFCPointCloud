import os
import sys
import uuid
import time
import json
import struct
import argparse
import itertools
import subprocess

import numpy as np

import ifcopenshell

import h5py

from collections import namedtuple, defaultdict

from schemas import ifc2x3_pc_cached as ifc
from utils import hdf5_mapping

argp = argparse.ArgumentParser()
argp.add_argument('ifc_file')
argp.add_argument('compressed', choices=['compressed', 'uncompressed']          )
argp.add_argument('precision',  choices=['single_precision', 'double_precision'])

ns = argp.parse_args()

flags = [
    '', 
    ns.precision,
    ns.compressed,
]

if not os.path.exists("generated_files/hdf"):
    os.makedirs("generated_files/hdf")

IFC_FN        = ns.ifc_file
IFC_FILE_BASE = os.path.splitext(os.path.basename(IFC_FN))[0]
HDF_FN        = "generated_files/hdf/" + IFC_FILE_BASE + "-".join(flags) + ".hdf"

SINGLE_PRECISION = ns.precision == 'single_precision'
USE_COMPRESSION  = ns.compressed == 'compressed'

print "Creating HDF5 schema mapping"

if os.path.exists(HDF_FN):
    os.unlink(HDF_FN)
f = h5py.File(HDF_FN, "w")

# Some attributes are ignored in order to prevent vlen attributes, which would render entities to have to be emulated
ignored = [('IfcPointCloud', 'Attributes'), ('IfcGrid', 'WAxes')]
# Some attributes are artificially set to have a fixed length, in order to prevent the aforementioned
fixed_length = {('IfcCartesianPoint', 'Coordinates'): 3,
                ('IfcDirection', 'DirectionRatios'): 3,
                ('IfcProductDefinitionShape', 'Representations'): 1,
                ('IfcShapeRepresentation', 'Items'): 1,
                ('IfcRelAssignsToProduct', 'RelatedObjects'): 1}
xref = [('IfcGridOffsetList', 'Offsets'),
        ('IfcGrid', 'UAxes'),
        ('IfcGrid', 'VAxes')]
precision = 32 if SINGLE_PRECISION else 64
mapping = hdf5_mapping.hdf5_mapping(ifc.schema, f, ignore=ignored, fix_length=fixed_length, force_xref=xref, precision=precision)

print "Reading SPF data"

ifcf = ifcopenshell.open(IFC_FN)

print "Creating HDF5 instance population"

mapping.create_population("Haus30_population", ifcf, compressed=USE_COMPRESSION)

f.close()

print "Done :)"
