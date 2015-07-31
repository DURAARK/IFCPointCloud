import os
import sys
import logging
import operator
import itertools
import collections

from collections import defaultdict

import numpy as np

import h5py

# Module packages reference each other cyclically,
# so this absolute import mechanism has to be used
def relative_import(path):
    abs_path = os.path.dirname(os.path.abspath(os.path.join(os.path.dirname(__file__), path)))
    name = os.path.basename(path)
    if name in sys.modules:
        globals()[name] = sys.modules[name]
    else:
        sys.path.insert(0, abs_path)
        globals()[name] = __import__(name)
        sys.path[0:1] = []

relative_import("../schemas/ifc2x3_pc_cached")
relative_import("../schemas/nodes")
relative_import("../dependencies/toposort")
relative_import("../ifcopenshell")

# isinstance() does not work well due to importing hack
isAggregationType = lambda n: n.__class__.__name__ == 'AggregationType'
isBinaryType = lambda n: n.__class__.__name__ == 'BinaryType'

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.ERROR)

# Some remaining questions
# ------------------------

# 6.6 Mappings for EXPRESS entity data types
# [...] The name of each member of the HDF5 Compound Type representing an EXPRESS attribute shall be the name of the EXPRESS explicit attribute in upper case.
# -- This does not seem to be the case in any of the examples, nor the appendix source code

# 6.8.3 EXPRESS array datatype values as HDF5
# Table 2 - Summary of mapping of EXPRESS array to HDF5
# ENUMERATION | same as for INTEGER
# 6.9.2 Mappings for EXPRESS enumeration data types as HDF5
# [...] Because of EXPRESS Edition 2 Extensible Enumeration Types, there is no ordering that can be related to an EXPRESS enumeration literal. Therefore, references must be made using the symbolic name.
# -- What is the symbolic name, how is it written as an INTEGER?

# 6.4 Mappings for EXPRESS simple data types
# Table 1 - Summary of mapping of EXPRESS data types to HDF5
# BINARY | HDF5 OPAQUE if FIXED or a width is specified, otherwise HDF5 VARYING OPAQUE
# -- HDF5 VARYING OPAQUE does not exist, some vlen attribute I suppose?

# 6.9.3.4 EXPRESS select of mixed types as HDF5
# -- name_of_typed_aggregate is not specified

class hdf5_mapping:

    compound_attribute = collections.namedtuple('compound_attribute', ('name', 'type'))
    
    # Emulated entity denotes a singleton which specifies a entity has a vlen attribute and needs to be emulated.
    # Emulation entails writing the entity in numpy's object serialization into a root vlen field of uint8
    class emulated_entity: pass
    emulated_entity = emulated_entity()
    
    # This is a quick annotation to denote a vlen type. The original resides in the 'type' attribute
    class vlen_attribute:
        def __init__(self, vlen_type): self.type = vlen_type

    def attributes(self, entity):
        supertype_attributes = self.attributes(self.schema.entities[entity.supertypes[0]]) if len(entity.supertypes) == 1 else []
        return supertype_attributes + entity.attributes
        
    def derived_attributes(self, entity):
        supertype_derived = self.derived_attributes(self.schema.entities[entity.supertypes[0]]) if len(entity.supertypes) == 1 else []
        own_derived = [str(s) for s in entity.derive.elements] if entity.derive else []
        return supertype_derived + own_derived
        
    def create_compound(self, name, attrs):
        size = sum(map(lambda a: a.type.get_size(), attrs))
        compound = h5py.h5t.create(h5py.h5t.COMPOUND, size)
        offset = 0
        for attr in attrs:
            compound.insert(attr.name, offset, attr.type)
            offset += attr.type.get_size()
        if name is not None:
            compound.commit(self.group.id, name)
        return compound
        
    # This part of the standard is not enforced. Iff a multidimensional aggregate (LIST OF LIST OF ...) an object reference is written 
    def create_aggregate_reference(self, aggregate_type):
        return self.create_compound(None, (
            hdf5_mapping.compound_attribute('obj_ref_or_vlen', h5py.h5t.NATIVE_INT8),
            hdf5_mapping.compound_attribute('object_reference', h5py.h5t.STD_REF_OBJ),
            hdf5_mapping.compound_attribute('vlen_array', h5py.h5t.vlen_create(aggregate_type))
        ))
        
    def flatten_aggregate_name(self, t):
        name = ""
        while isAggregationType(t):
            name += "%s-of-" % t.aggregate_type
            t = t.type
        return name + t
        
    def follow_simpletype(self, simpletype):
        while self.schema.is_simpletype(simpletype):
            simpletype = self.schema.simpletypes[simpletype]
        return simpletype
        
    # 6.5 Mappings for EXPRESS schema declarations and interface specifications
    def __init__(self, schema, file, ignore=None, fix_length=None, force_xref=None, precision=64):
        self.types = {}
        # vlen flags for simpletypes, checked in the mapping. big problem for IfcRoot
        self.vlen_flags = {}
        # count the number of select possibilities (which result into compound attribute fields) for default initialization
        self.compound_attribute_count = {}
        self.compound_attributes = defaultdict(list)
        self.schema = schema
        self.file = file
        self.group = file.create_group("%s_encoding" % schema.name)
        self.group.attrs['iso_10303_26_schema'] = schema.name
        self.xref_attributes = defaultdict(dict)
        if force_xref is not None:
            for entity, attr_name in force_xref:
                self.xref_attributes[entity][attr_name] = True
        self.ignored_attributes = set() if ignore is None else set(ignore)
        self.fixed_length_attributes = {} if fix_length is None else fix_length
        
        # Create the instance reference type according to:
        # 6.10.4 EXPRESS entity instance references
        self.instance_reference = self.create_compound('_HDF_INSTANCE_REFERENCE_HANDLE_', (
            hdf5_mapping.compound_attribute('_HDF5_dataset_index_', h5py.h5t.NATIVE_INT16),
            hdf5_mapping.compound_attribute('_HDF5_instance_index_', h5py.h5t.NATIVE_INT32)
        ))
        
        self.boolean_type = h5py.h5t.enum_create(h5py.h5t.NATIVE_INT8)
        for i, v in enumerate(("BOOLEAN-FALSE", "BOOLEAN-TRUE")):
            self.boolean_type.enum_insert(v, i)
        
        self.logical_type = h5py.h5t.enum_create(h5py.h5t.NATIVE_INT8)
        for i, v in enumerate(("LOGICAL-UNKNOWN", "LOGICAL-FALSE", "LOGICAL-TRUE")):
            self.logical_type.enum_insert(v, i - 1)
            
        # h5py.h5t.special_dtype(vlen=...)?
        self.binary_type = h5py.h5t.vlen_create(h5py.h5t.create(h5py.h5t.OPAQUE, 1)) # 1?
        
        self.float_type = h5py.h5t.NATIVE_DOUBLE if precision==64 else h5py.h5t.NATIVE_FLOAT
        
        self.byte       = np.dtype("int8")
        self.byte_array = h5py.special_dtype(vlen=self.byte)
        
        self.object_reference_type = h5py.h5t.py_create(h5py.special_dtype(ref=h5py.Reference), logical=True)
        
        self.string_type = h5py.h5t.C_S1.copy()
        self.string_type.set_size(32)
        
        for enum in schema.enumerations.items(): 
            self.types[enum[0]] = self.map_enumeration(*enum)
            
        type_dependencies = {}
        for type_name, simpletype in schema.simpletypes.items():
            while isAggregationType(simpletype):
                simpletype = simpletype.type
            if isBinaryType(simpletype):
                simpletype = 'binary'
            if isinstance(simpletype, str): type_dependencies[type_name] = {simpletype}
            else:
                raise ValueError("Unexpected %r of type %r" % (simpletype, type(simpletype)))
                
        for type_name in toposort.toposort_flatten(type_dependencies):
            simpletype = schema.simpletypes.get(type_name)
            if simpletype:
                self.types[type_name] = self.map_simpletype(type_name, simpletype)
            else:
                # one of the express data types
                pass
                
        for type_name, select in schema.selects.items():
            dtype = self.map_select_type(type_name, select)
            if dtype is not None:
                self.types[type_name] = dtype
            else:
                # Express select types that resolve solely to entities (category b) are transparent in HDF5 in the sence that instances are populated with just a _HDF_INSTANCE_REFERENCE_HANDLE_
                self.types[type_name] = self.instance_reference
                pass
                
        for entity in schema.entities.values():
            self.types[entity.name] = self.map_entity(entity)
                
        # self.file.create_dataset("booleans", (100,), dtype=self.types['IfcWall'])
        
    def map_express_type_allow_vlen(self, *args):
        attr = self.map_express_type(*args)
        if isinstance(attr, hdf5_mapping.vlen_attribute):
            attr = attr.type
        return attr
    
    # 6.9.3 Mappings for EXPRESS select data types as HDF5
    def map_select_type(self, typename, select_type):
        
        self.compound_attribute_count[typename] = 2
        
        def follow_dependent_types(s):
            # print "follow_dependent_types(%r)" % s
            deps = set()
            for v in s.values:
                sl = self.schema.selects.get(v)
                if sl is None:
                    deps |= {v}
                else:
                    deps |= follow_dependent_types(sl)
            return deps
        def select_category(select_type):
            # a) all resolve to same underlying express type [TK: isn't this a weird corner case?]
            # b) only entity instance refs
            # c) mixed
            deps = follow_dependent_types(select_type)
            if all(map(self.schema.is_entity, deps)):
                return 'b'
            else:
                return 'c'
        cat = select_category(select_type)
        if cat == 'a':
            raise ValueError("Unimplemented selection type returned <%s>" % cat)
        elif cat == 'b':
            return None
        elif cat == 'c':
            # Due to type_path always needs to be emulated?
            # TODO: Check whether this is the case, because vlen string is special case
            # It is. But strings are fixed width now
            # return hdf5_mapping.emulated_entity
            
            # TODO: A type_path is not always necessary
            attrs = [
                hdf5_mapping.compound_attribute('select_bitmap', h5py.h5t.NATIVE_INT8),
                hdf5_mapping.compound_attribute('type_path', self.string_type) # type path should be LIST OF STRING, but will be delimited instead
            ]
            attr_names = set()
            deps = follow_dependent_types(select_type)
            
            def insert_attribute(name, value):
                if name not in attr_names:
                    attrs.append(hdf5_mapping.compound_attribute(name, value))
                    attr_names.add(name)
                    self.compound_attribute_count[typename] += 1
                    return True
                return False
                    
            def create_dependent_basetype(t):
                if self.schema.is_entity(t):
                    if insert_attribute('instance-value', self.instance_reference):
                        self.compound_attributes[typename].append(lambda: (0,0))
                elif self.schema.is_enumeration(t):
                    if insert_attribute('%s-value' % t, self.types[t]):
                        self.compound_attributes[typename].append(int)
                elif self.schema.is_simpletype(t):
                    t = self.follow_simpletype(t)
                    create_dependent_basetype(t)
                elif isAggregationType(t):
                    # do not allow vlen either?
                    # insert_attribute(self.flatten_aggregate_name(t), self.map_express_type_allow_vlen(t))
                    if insert_attribute(self.flatten_aggregate_name(t), self.map_express_type(t)):
                        self.compound_attributes[typename].append(str)
                else:
                    # do not allow vlen either?
                    # insert_attribute('%s-value' % t, self.map_express_type_allow_vlen(t))
                    if insert_attribute('%s-value' % t, self.map_express_type(t)):
                        self.compound_attributes[typename].append({
                            'number': float,
                            'real': float,
                            'integer': int,
                            'string': str
                        }.get(t))
                # else:
                #     raise ValueError("Unexpected dependent type encountered <%r>" % t)
            for dep in deps: create_dependent_basetype(dep)
            
            # If the list of attributes contains a variable length attribute, the select has to be emulated
            attr_types = [a.type.__class__ for a in attrs]
            if hdf5_mapping.vlen_attribute in attr_types:
                return hdf5_mapping.emulated_entity                
            
            return self.create_compound(typename, attrs)
        else:
            raise ValueError("Undefined selection type returned <%s>" % cat)
        

    def map_express_type(self, ty, levels=0, attr_name=None, entity_name=None):
        def f(ty):
            if isBinaryType(ty): ty = 'binary'
            if isinstance(ty, str):
                # Entity instance references are not statically typed in HDF5
                if ty in self.schema.entities: 
                    return self.instance_reference                
                
                mapped_type = self.types.get(ty)
                if mapped_type: 
                    if self.vlen_flags.get(ty, False):
                        # it is still necessary to flag vlen attributes for simple types.
                        # unfortunately this rules out any IfcRoot descendent [because of the vlen GlobalId (which is not really vlen), but also Name and Description] to have a proper h5py implementation
                        return hdf5_mapping.vlen_attribute(mapped_type)
                    else:
                        return mapped_type
                
                express_type = {
                    'boolean' : self.boolean_type,
                    'logical' : self.logical_type,
                    'integer' : h5py.h5t.NATIVE_INT32,
                    'real'    : self.float_type,
                    'number'  : self.float_type,
                   #'string'  : hdf5_mapping.vlen_attribute(self.string_type),
                   #'binary'  : hdf5_mapping.vlen_attribute(self.string_type)
                    'string'  : self.string_type,
                    'binary'  : self.string_type
                   # TODO: Binaries just treated as strings for now
                   #'binary'  : self.binary_type                   
                }.get(ty)
                # TK: Take care not to evaluate thruthiness, the string_type evaluates to false
                if express_type is not None: return express_type
                
            elif isAggregationType(ty):
                # why does this return a <numpy.dtype>?
                # TODO: Arrays should *not* be mapped as vlen attributes
                # TODO: Array unset elements, contains possible null references
                def fixed_width(bs):
                    l, u = map(int, (bs.lower, bs.upper))
                    return u == l
                
                artificially_fixed_array_bound = self.fixed_length_attributes.get((entity_name, attr_name))
                if artificially_fixed_array_bound is not None:
                    """" 
                    print >>sys.stderr, "[NOTICE] On entity:"
                    print >>sys.stderr, "", self.schema.entities[entity_name]
                    print >>sys.stderr, "remapped:"
                    print >>sys.stderr, " %s: %s" % (attr_name, ty)
                    print >>sys.stderr, "to:\n LIST[%d:%d]" % (artificially_fixed_array_bound, artificially_fixed_array_bound)
                    print >>sys.stderr, ""
                    """
                    return h5py.h5t.array_create(self.map_express_type(ty.type), (artificially_fixed_array_bound,))
                
                if ty.aggregate_type == "list":
                    if isAggregationType(ty.type) and ty.type.aggregate_type == "list":
                        #NB multi dim aggregates always x-reffed
                        return self.object_reference_type
                        # return h5py.h5t.STD_REF_OBJ
                        
                        # obsolete:
                        if fixed_width(ty.bounds) and fixed_width(ty.type.bounds):
                            dims = tuple(map(int, (ty.bounds.upper, ty.type.bounds.upper)))
                            return h5py.h5t.array_create(self.map_express_type(ty.type.type), dims)
                    else:
                        pass
                        # FIXME: Implement, but does not exist in schema
                
                # Not mapped as a fixed length attributes.
                # vlen instead, the horror.
                np_dtype = h5py.h5t.special_dtype(vlen=self.map_express_type(ty.type, levels+1))
                return hdf5_mapping.vlen_attribute(h5py.h5t.py_create(np_dtype))

            raise ValueError("Unexpected type <%r>" % ty)
        d = f(ty)
        ### print "%s<%s.%s> -> <%s>" % ("    " * levels, self.schema.name, ty, type(d).__name__)
        return d
        
    # 6.9.2 Mappings for EXPRESS enumeration data types as HDF5
    def map_enumeration(self, enum_name, enum):
        dtype = h5py.h5t.enum_create(h5py.h5t.NATIVE_INT8) # 16?
        for i, v in enumerate(sorted(enum.values)):
            dtype.enum_insert("/".join((enum_name, v)), i)
        dtype.commit(self.group.id, enum_name)
        return dtype
        
    # 6.9.4 Mappings for EXPRESS simple defined types
    def map_simpletype(self, type_name, simpletype):
        # header_print("%s.%s" % (self.schema.name, type_name))
        dtype = self.map_express_type_allow_vlen(simpletype)
        self.vlen_flags[type_name] = isinstance(self.map_express_type(simpletype), hdf5_mapping.vlen_attribute)
        # print type(dtype), dtype
        dtype = dtype.copy()
        dtype.commit(self.group.id, type_name)
        return dtype
        
    # 6.6 Mappings for EXPRESS entity data types
    def map_entity(self, entity):
        logging.debug("Start mapping %s.%s" % (self.schema.name, entity.name))
        attributes = self.attributes(entity)
        derived = set(self.derived_attributes(entity))
        
        attrs = [
            hdf5_mapping.compound_attribute('set_unset_bitmap', h5py.h5t.NATIVE_INT16),
            hdf5_mapping.compound_attribute('Entity-Instance-Identifier', h5py.h5t.NATIVE_INT32)
        ]
        
        for attr in attributes:
            # print "%s - %s" % (entity.name, attr.name)
            # print attr
            if attr.name in derived: 
                # print "*derived*"
                continue
            if (entity.name, attr.name) in self.ignored_attributes:
                """
                print >>sys.stderr, "[NOTICE] On entity:"
                print >>sys.stderr, "", entity
                print >>sys.stderr, "ignored:"
                print >>sys.stderr, "", attr
                print >>sys.stderr, ""
                """
                continue
            
            if self.xref_attributes[entity.name].get(attr.name, False):
                mapped_type = h5py.h5t.STD_REF_OBJ 
            else:
                mapped_type = self.map_express_type(attr.type, entity_name=entity.name, attr_name=attr.name)
            
                if mapped_type is h5py.h5t.STD_REF_OBJ:
                    self.xref_attributes[entity.name][attr.name] = True
                
            if isinstance(mapped_type, hdf5_mapping.vlen_attribute) or mapped_type is hdf5_mapping.emulated_entity:
                return hdf5_mapping.emulated_entity
                
            attrs.append(hdf5_mapping.compound_attribute(attr.name, mapped_type))

        ### print ""
        return self.create_compound(entity.name, attrs)
            
    # 6.3.3 EXPRESS-driven data populations as HDF5
    def create_population(self, population_name, instances, compressed=False):
        dataset_names = sorted(set(i.is_a() for i in instances))
        dataset_name_mapping = dict((j,i) for i,j in enumerate(dataset_names))
        dataset_instance_indices = dict((i, dict((r,q) for q, r in enumerate(sorted(map(operator.itemgetter(1), j))))) for i,j in itertools.groupby(sorted((i.is_a(), i.id()) for i in instances), operator.itemgetter(0)))
        dataset_instance_index = lambda inst: dataset_instance_indices.get(inst.is_a()).get(inst.id())
        entity_instance_to_hdf5 = lambda inst: (dataset_name_mapping[inst.is_a()], dataset_instance_index(inst))
        
        group = self.file.create_group(population_name)
        
        if False:
            included_entities = set(['IfcProductDefinitionShape', 'IfcPlane', 'IfcDirection', 'IfcAxis2Placement3D', 'IfcPointCloudElement', 'IfcRectangularTrimmedSurface',
                'IfcParameterValueList', 'IfcPointCloud', 'IfcShapeRepresentation', 'IfcCartesianPoint', 'IfcRelAssignsToProduct', 'IfcRelDecomposes', 'IfcGrid',
                'IfcGridOffsetList', 'IfcLine', 'IfcOffsetCurve2D', 'IfcContinuousParameterValueList', 'IfcDiscreteParameterValueList', 'IfcSurfaceGrid',
                'IfcContiousGridOffsetList', 'IfcDiscreteGridOffsetList', 'IfcCartesianPointList3D'])
        else:
            included_entities = dataset_names
        
        for entity in included_entities:
            entity_group = group.create_group("%s_objects" % entity)
            
            sys.stdout.write("\r%s                    " % entity)
            sys.stdout.flush()
            def f():
                for inst in sorted(instances.by_type(entity), key=lambda i: i.id()):
                    # additional check to filter out subtypes
                    if inst.is_a() != entity: continue
                    logging.debug("> %r" % inst)
                    def g():
                        attr_values = [inst[i] for i in range(len(inst))]
                        set_unset = map(lambda x: x is not None, attr_values)
                        set_unset = sum(map(lambda a: a[1] << a[0], enumerate([inst[i] is not None for i in range(len(inst))])))
                        yield set_unset
                        yield inst.id()
                        def m(schema_attr, attr, length_fixed=False):
                            flattened_type = self.follow_simpletype(schema_attr.type)
                            result = None
                            if (entity, schema_attr.name) in self.fixed_length_attributes and not length_fixed:
                                assert isinstance(attr, tuple)
                                # TK: NB this has only been implemented for floats and instances
                                flattened_aggre_elem_type = self.follow_simpletype(schema_attr.type.type)
                                fixed_length = self.fixed_length_attributes[(entity, schema_attr.name)]
                                if flattened_aggre_elem_type == 'real':
                                    return attr[:fixed_length] + (float('nan'),) * (fixed_length - len(attr))
                                else:
                                    attr = tuple(map(entity_instance_to_hdf5, attr))
                                    attr = attr[:fixed_length] + ((-1,-1),) * (fixed_length - len(attr))
                                    # TK: This is a bit peculiar, with fixed length aggregates of length 1 numpy does not expect a sequence:
                                    #           >>> np.ndarray(1, np.dtype([('a', [('b', 'i', 1)])]))[0]
                                    #           ((5832796,),)                                 ^
                                    #           >>> np.ndarray(1, np.dtype([('a', [('b', 'i', 2)])]))[0]
                                    #           (([5832792, 47120432],),)                     ^
                                    # ie. the level of nested-ness is different depending on the aggregate size
                                    if len(attr) == 1:
                                        return attr[0]
                                    else: return attr
                            elif entity in self.xref_attributes and self.xref_attributes[entity].get(schema_attr.name) is not None:
                                if attr is None:
                                    return group.ref
                                else:
                                    # Hackhackahackahakchack
                                    if isinstance(attr, tuple):
                                        if set(map(type, attr)) == {ifcopenshell.entity_instance}:
                                            attr = np.array(map(entity_instance_to_hdf5, attr), dtype=self.instance_reference.dtype)
                                        elif set(map(type, attr)) == {int} or set(map(type, attr)) == {long}:
                                            nbytes = int(np.ceil(np.log2(max(attr)))) // 8
                                            attr = np.array(attr, dtype="u%d" % nbytes)
                                        elif set(map(type, attr)) == {tuple} and (set(map(type, attr[0])) == {int} or set(map(type, attr[0])) == {long}):
                                            max_value = max(map(max, attr))
                                            nbytes = int(np.ceil(np.log2(max_value))) // 8
                                            attr = np.array(attr, dtype="u%d" % nbytes)
                                        elif set(map(type, attr)) == {float}:
                                            attr = np.array(attr, dtype=self.float_type.dtype)
                                        elif set(map(type, attr)) == {tuple} and set(map(type, attr[0])) == {float}:
                                            attr = np.array(attr, dtype=self.float_type.dtype)
                                    
                                    kwargs = {}
                                    if compressed:
                                        kwargs["shuffle"] = True
                                        kwargs["compression"] = "gzip"
                                        kwargs["compression_opts"] = 9
                                    ds = entity_group.create_dataset("Aggr_%s_%d" % (schema_attr.name, inst.id()), data=attr, **kwargs)
                                    return ds.ref
                            elif isinstance(attr, ifcopenshell.entity_instance):
                                if attr.id() == 0:
                                    # These are simple types:
                                    # See whether this select is ambiguous
                                    if isinstance(self.types[flattened_type], h5py.h5t.TypeCompoundID):
                                        # NB: Not implemented
                                        return (0,'') + tuple(a() for a in self.compound_attributes[flattened_type])
                                    else:                                        
                                        return m(schema_attr, attr.wrapped_data)
                                else:
                                    result = entity_instance_to_hdf5(attr)
                            elif flattened_type in {'boolean', 'logical'}:
                                # booleans and logicals are also enumerations in HDF5
                                # int representation of Python bools corresponds to HDF5 enumeration values
                                # logical not properly implemented in IfcOpenShell, only {True, False}
                                if attr is None: return 0
                                else: return int(attr)
                            elif self.schema.is_enumeration(flattened_type):
                                if attr is None: return 0
                                result = self.schema.enumerations[flattened_type].values.index(attr)
                            elif isinstance(attr, list):
                                # TK: This is obsolete probably?
                                result = tuple(attr)
                            elif isinstance(attr, str):
                                result = attr.encode()
                            elif attr is None:
                                # empty fields do not exist. hence values are default initialized
                                if self.schema.is_entity(flattened_type):
                                    result = (-1,-1)
                                elif self.schema.is_select(flattened_type):
                                    if self.types[flattened_type] == self.instance_reference:
                                        result = (-1,-1)
                                    elif isinstance(self.types[flattened_type], h5py.h5t.TypeCompoundID):
                                        # evaluate the default constructors:
                                        result = (0,'') + tuple(a() for a in self.compound_attributes[flattened_type])
                                    else:
                                        result = None
                                else:
                                    result = {
                                        'string' : '',
                                        'real'   : float('NaN'),
                                        'number' : float('NaN'),
                                        'integer': 0
                                    }.get(flattened_type)
                            else:
                                result = attr
                            logging.debug("| %r (%s) -> %r" % (schema_attr, flattened_type, result))
                            return result
                        schema_entity = self.schema.entities[inst.is_a()]
                        derived = set(self.derived_attributes(schema_entity))
                        for schema_attr in zip(self.attributes(schema_entity), attr_values):
                            attr_name = schema_attr[0].name
                            if attr_name not in derived:
                                if (entity, attr_name) not in self.ignored_attributes:
                                    yield m(*schema_attr)
                    if self.types[inst.is_a()] is hdf5_mapping.emulated_entity:
                        yield np.array(map(ord, np.array(list(g()), dtype=object).tostring()))
                    else:
                        try:
                            yield np.array(tuple(g()), dtype=self.types[inst.is_a()])
                        except Exception as e:
                            print >>sys.stderr, "\n[ERROR] Failed to encode:"
                            print >>sys.stderr, inst
                            print >>sys.stderr, "Reason: '%s'" % str(e)
                            print >>sys.stderr, "Storing:\n%r" % (tuple(g()),)
                            print >>sys.stderr, "In:\n%s" % str(self.types[inst.is_a()].dtype)
                        
            def create_instance_dataset(**kwargs):
                if compressed:
                    kwargs["shuffle"] = True
                    kwargs["compression"] = "gzip"
                    kwargs["compression_opts"] = 9
                return entity_group.create_dataset("%s_instances" % entity, **kwargs)
            
            if self.types[entity] is hdf5_mapping.emulated_entity:
                d = np.array(list(f()))
                ds = create_instance_dataset(shape=(d.shape[0], 1), dtype=self.byte_array)
                for i, elem in enumerate(d):
                    ds[i] = elem
            elif entity in self.xref_attributes:
                values = list(f())
                # TK: Immediate initialization throws a numpy exception:
                #     ValueError: Setting void-array with object members using buffer.
                #     More info here? https://github.com/h5py/h5py/issues/463
                #     Assignment in a loop seems to work
                if len(values):
                    ds = create_instance_dataset(dtype=values[0].dtype, shape=(len(values),))
                    for i, v in enumerate(values):
                        ds[i] = v
            else:
                d = np.array(list(f()))
                create_instance_dataset(data=d)
                
        # For newline
        print ""
    
        group.attrs["iso_10303_26_data"] = self.schema.name
        group.attrs["iso_10303_26_data_set_names"] = dataset_names
        # TODO attributes:
        # iso_10303_26_data           = <schema_id> : HDF5 STRING
        # iso_10303_26_data_set_names = ? 
        # optional:
        # iso_10303_26_description: which has a data value of the <user_defined_description>
        # iso_10303_26_timestamp: which has a data value corresponding to the extended format for the complete calendar date as specified in ISO 8601 concatenated to the extended format for the time of the day also as specified in ISO 8601. The date and time shall be separated by the capital letter "T" as specified in ISO 8601, which also defines alternate formats that permit the optional inclusion of a time zone specifier.
        # iso_10303_26_author: which has a data value of the <user>
        # iso_10303_26_organization: which has a data value of the <user_organization>
        # iso_10303_26_originating_system: which has a data value of the <software_system_name>
        # iso_10303_26_preprocessor_version: which has a data value of the <software_application_and_version>
        # iso_10303_26_context: which has a data value of the <context within which data is applicable>
        # iso_10303_26_language: which has a data value of the <default language for string values> where the name of the language shall be encoded using the Alpha-3 bibliographic code specified in ISO 639-2.
        pass
