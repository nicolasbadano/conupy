# GIS functions for ArcGis 10.1
import math
import pickle
from itertools import tee, izip
import numpy as np

def saveOnFile( data, name ):
    oF = open(name + ".pkl", 'wb')
    pickle.dump(data, oF)
    oF.close()


def readFromFile( name ):
    iF = open(name + ".pkl", 'rb')
    data = pickle.load(iF)
    iF.close()
    return data


def dist(p0, p1):
    return ((p0[0]-p1[0])**2+(p0[1]-p1[1])**2)**0.5


def distToSegment(p0, p10, p11):
    r = p0 - p10
    n = p11 - p10

    n_norm = np.linalg.norm(n)
    nn = n / n_norm
    alpha = np.dot(nn, r)

    if alpha < 0:
        return np.linalg.norm(r)
    elif alpha > n_norm:
        return np.linalg.norm(p0 - p11)

    return np.linalg.norm(r - nn * alpha)


def distSq(p0, p1):
    return ((p0[0]-p1[0])**2+(p0[1]-p1[1])**2)


def interpolate(p0, p1, alpha):
    return [p0[0] + (p1[0]-p0[0])*alpha, p0[1] + (p1[1]-p0[1])*alpha]


def insertPoints(points, maxLength = 50):
    i = 0
    while (i<len(points)-1):
        p0 = points[i]
        p1 = points[i+1]
        length = dist(p0, p1)
        if (length > maxLength):
            numSpans = int(math.floor(length / maxLength + 1))
            newLength = length / numSpans
            for j in xrange(1,numSpans):
                points.insert(i+j, interpolate(p0, p1, float(j)/float(numSpans)))
            i += numSpans
        else:
            i += 1


def removePoints(points, maxLength = 25):
    i = 1
    while (i<len(points)-1):
        p2 = points[i-1]
        p0 = points[i]
        p1 = points[i+1]
        lengthAnterior = dist(p2, p0)
        lengthSiguiente = dist(p0, p1)
        if lengthAnterior < maxLength and lengthSiguiente < maxLength and p0[-1] != "snapped":
            points.pop(i)
        else:
            i += 1


def atraviesaArroyo(p0, p1, nodos, links):
    for (n0, n1), link in links.iteritems():
        if link["type"] == "channel":
            if (intersect(p0, p1, nodos[n0], nodos[n1])):
                print n0, n1, link
                return n0, n1, link

    return None


def ccw(A,B,C):
    return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])


def intersection(P1,P2,P3,P4):
    return [
        ((P1[0]*P2[1] - P1[1]*P2[0]) * (P3[0] - P4[0]) - (P1[0] - P2[0]) * (P3[0]*P4[1] - P3[1]*P4[0])) / ((P1[0] - P2[0]) * (P3[1] - P4[1]) - (P1[1] - P2[1]) * (P3[0] - P4[0])),
        ((P1[0]*P2[1] - P1[1]*P2[0]) * (P3[1] - P4[1]) - (P1[1] - P2[1]) * (P3[0]*P4[1] - P3[1]*P4[0])) / ((P1[0] - P2[0]) * (P3[1] - P4[1]) - (P1[1] - P2[1]) * (P3[0] - P4[0]))]


def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


def addNode(nodos, punto, tipo, geo_hash, tolSq = 30):
    if not any(geo_hash):
        # The geo_hash is empty, populate it
        for nodo in nodos:
            ix, iy = int(nodo[0]/1000.0), int(nodo[1]/1000.0)
            geo_hash[(ix,iy)] = geo_hash.get((ix, iy), [])
            geo_hash[(ix,iy)].append(len(nodos)-1)

    ix, iy = int(punto[0]/1000.0), int(punto[1]/1000.0)
    nodosCercanos = []
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            a = geo_hash.get((ix+i, iy+j), [])
            nodosCercanos.extend(a)


    for i in nodosCercanos:
        nodo = nodos[i]
        if distSq(nodo, punto) < tolSq:
            if nodo[2] == "conduit":
                # Si llega al menos un arroyo a un nodo conducto, pasa a ser nodo arroyo
                if tipo == "channel":
                    nodos[i][2] = "channel"
                    return i
                elif tipo == "conduit":
                    return i
                # Si se quiere crear un nodo esquina cerca de uno conducto, no los une
                elif tipo == "corner":
                    continue

            elif nodo[2] == "corner":
                if tipo == "channel":
                    #Esto no debería ocurrir, porque los arroyos se crean antes que las esquinas
                    return i
                elif tipo == "conduit":
                    #Esto no debería ocurrir, porque los arroyos se crean antes que las esquinas
                    continue
                elif tipo == "corner":
                    return i

            elif nodo[2] == "channel":
                if tipo == "channel":
                    return i
                elif tipo == "conduit":
                    return i
                elif tipo == "corner":
                    return i

    # Si no hay nodos cercanos, agregar uno
    punto.append(tipo)
    nodos.append(punto)
    geo_hash[(ix,iy)] = geo_hash.get((ix, iy), [])
    geo_hash[(ix,iy)].append(len(nodos)-1)
    return len(nodos)-1;


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)
