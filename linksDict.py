# linksDict.py encoding: utf-8
import numpy as np

class LinksDict(dict):
    def __init__(self, nodes):
        self.nodes = nodes
        self.geo_hash = {}
        super(dict, self).__init__()

    def __setitem__(self, key, value):
        # optional processing here

        # Update geo cache
        n0, n1 = key
        ix0, iy0 = [int(z/100.0) for z in self.nodes[n0].p]
        ix1, iy1 = [int(z/100.0) for z in self.nodes[n1].p]
        for ix in range(min(ix0, ix1), max(ix0, ix1) + 1):
            for iy in range(min(iy0, iy1), max(iy0, iy1) + 1):
                self.geo_hash[(ix,iy)] = self.geo_hash.get((ix, iy), [])
                self.geo_hash[(ix,iy)].append((n0, n1, value))

        super(LinksDict, self).__setitem__(key, value)

    def __delitem__(self, key):
        # Update geo cache
        value = self[key]
        n0, n1 = key
        ix0, iy0 = [int(z/100.0) for z in self.nodes[n0].p]
        ix1, iy1 = [int(z/100.0) for z in self.nodes[n1].p]
        for ix in range(min(ix0, ix1), max(ix0, ix1) + 1):
            for iy in range(min(iy0, iy1), max(iy0, iy1) + 1):
                self.geo_hash[(ix,iy)] = self.geo_hash.get((ix, iy), [])
                self.geo_hash[(ix,iy)].remove((n0, n1, value))

        return super(LinksDict, self).__delitem__(key)

    def getLinksInside(self, p0, p1):
        ix0, iy0 = [int(z/100.0) for z in p0]
        ix1, iy1 = [int(z/100.0) for z in p1]
        links = dict()
        for ix in range(min(ix0, ix1), max(ix0, ix1) + 1):
            for iy in range(min(iy0, iy1), max(iy0, iy1) + 1):
                if (ix, iy) in self.geo_hash:
                    for n0, n1, link in self.geo_hash[(ix, iy)]:
                        links[(n0, n1)] = link
        return links

    def getLinksNear(self, p0, distance):
        return self.getLinksInside(np.add(p0, -distance), np.add(p0, distance))

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError("update expected at most 1 arguments, "
                                "got %d" % len(args))
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]
