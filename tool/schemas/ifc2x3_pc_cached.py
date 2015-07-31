import os, pickle

cache_file = "intermediate_files/ifc2x3_pc_schema.cache"

if os.path.isfile(cache_file):
    with open(cache_file, "rb") as f:
        schema = pickle.load(f)
else:
    import ifc2x3_pc
    schema = ifc2x3_pc.schema
    with open(cache_file, "wb") as f:
        pickle.dump(schema, f, protocol=0)
