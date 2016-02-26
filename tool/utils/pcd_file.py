import gzip

def open(fn, mod=None):
    f = gzip.open(fn, 'r')
    for l in ['', 'VERSION', 'FIELDS', 'SIZE', 'TYPE', 'COUNT', 'WIDTH', 'HEIGHT', 'VIEWPOINT', 'POINTS', 'DATA']:
        f.readline()

    if mod is None:
        for l in f:
            xyz = tuple(map(float, l.split()[0:3]))
            yield xyz
    else:
        for i, l in enumerate(f):
            if i % mod == 0: 
                xyz = tuple(map(float, l.split()[0:3]))
                yield xyz
        