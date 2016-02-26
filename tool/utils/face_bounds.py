from . import face_cache_reader

class proxy(object):
    def __init__(self, fcr):
        self.fcr = fcr
    def __getitem__(self, i):
        return self.fcr[i].uv_bounds

def get_mapping(fn):
    return proxy(face_cache_reader.lazy(fn))
    return dict((k, fa.uv_bounds) for k, fa in face_cache_reader.read(fn))
