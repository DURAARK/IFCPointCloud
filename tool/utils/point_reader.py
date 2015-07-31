import os
import sys
import struct
import itertools
import numpy as np
from collections import namedtuple

associated_point = namedtuple("associated_point", list("xyzefuvd"))
free_point = namedtuple("associated_point", list("xyz"))

def obtain_points(filename_pattern, range=None):
    if range is None: range = itertools.count()

    for i in range:
        if os.path.exists(filename_pattern % i) and not os.path.exists((filename_pattern % i) + ".lock"):
            with open(filename_pattern % i, "rb") as f:
                eof = os.fstat(f.fileno()).st_size
                while f.tell() < eof:
                    if (f.tell() // 1000) % 1000 == 0:
                        sys.stderr.write("\r[pass %d] Reading association data: %d kb   " % (i, f.tell() // 1024))
                    r = lambda fmt: struct.unpack(fmt, f.read(struct.calcsize(fmt)))
                    xyz = r("@ddd")
                    has_assoc = r("@q")[0]
                    if has_assoc:
                        efuvd = r("@qqddd")
                        yield associated_point(*(xyz+efuvd))
                    else:
                        yield free_point(*xyz)
        else: break
    print ""
        
class buffered_dataset:
    def __init__(self, xshape=3):
        self.buffer = []
        self.ds = np.ndarray((0,3))
    def merge(self):
        if len(self.buffer):
            b = self.buffer
            while len(b) > 1:
                nb = []
                def add(*args): 
                    nb.append(np.vstack(args))
                for i, current in enumerate(b):
                    if i % 2: previous = add(previous, current)
                    else: previous = current
                if previous is not None: nb.append(previous)
                b = nb
            self.ds = np.vstack((self.ds, b[0]))
            self.buffer = []
        return self.ds
    def add(self, *args):
        self.buffer.append(np.array(args))
        if len(self.buffer) > 1000: self.merge()
        
def read(filename_pattern, range=None):
    datasets = {}
    global_ds = buffered_dataset()
    
    for p in obtain_points(filename_pattern, range):
        if isinstance(p, associated_point):
            k = p.e, p.f
            ds = datasets.get(k)
            if ds is None: ds = datasets[k] = buffered_dataset()
            ds.add(p.u, p.v, p.d)
        else:
            global_ds.add(p.x, p.y, p.z)
            
    return datasets, global_ds
