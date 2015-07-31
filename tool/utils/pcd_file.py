import gzip

def open(fn):
    f = gzip.open(fn, 'r')
    for l in ['', 'VERSION', 'FIELDS', 'SIZE', 'TYPE', 'COUNT', 'WIDTH', 'HEIGHT', 'VIEWPOINT', 'POINTS', 'DATA']:
        f.readline()

    for l in f:
        xyz = tuple(map(float, l.split()[0:3]))
        yield xyz
