# GIS functions for ArcGis 10.1
import math
import pickle

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


def distSq(p0, p1):
    return ((p0[0]-p1[0])**2+(p0[1]-p1[1])**2)


def interpolate(p0, p1, alpha):
    return [p0[0] + (p1[0]-p0[0])*alpha, p0[1] + (p1[1]-p0[1])*alpha]


def insertPoints(puntos, distMax = 50):
    i = 0
    while (i<len(puntos)-1):
        p0 = puntos[i]
        p1 = puntos[i+1]
        largo = dist(p0, p1)
        if (largo > distMax):
            numTramos = int(math.floor(largo / distMax + 1))
            largoNuevo = largo / numTramos
            for j in xrange(1,numTramos):
                puntos.insert(i+j, interpolate(p0, p1, float(j)/float(numTramos)))
            i += numTramos
        else:
            i += 1


def removePoints(puntos, distMax = 25):
    i = 1
    while (i<len(puntos)-1):
        p2 = puntos[i-1]
        p0 = puntos[i]
        p1 = puntos[i+1]
        largoAnterior = dist(p2, p0)
        largoSiguiente = dist(p0, p1)
        if largoAnterior < distMax and largoSiguiente < distMax and p0[-1] != "snapped":
            puntos.pop(i)
        else:
            i += 1


def atraviesaArroyo(p0, p1, nodos, lineas):
    for linea in lineas:
        if (linea[2] == "arroyo"):
            if (intersect(p0, p1, nodos[linea[0]], nodos[linea[1]])):
                return linea[0]
        if (linea[2] == "calle"):
            break

    return -1


def ccw(A,B,C):
    return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])


def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


def addNode(nodos, punto, tipo = "esquina", tolSq = 30):
    for (i, nodo) in enumerate(nodos):
        if distSq(nodo, punto) < tolSq:
            if nodo[2] == "conducto":
                # Si llega al menos un arroyo a un nodo conducto, pasa a ser nodo arroyo
                if tipo == "arroyo":
                    nodos[i][2] = "arroyo"
                    return i
                elif tipo == "conducto":
                    return i
                # Si se quiere crear un nodo esquina cerca de uno conducto, no los une
                elif tipo == "esquina":
                    continue

            elif nodo[2] == "esquina":
                if tipo == "arroyo":
                    #Esto no debería ocurrir, porque los arroyos se crean antes que las esquinas
                    return i
                elif tipo == "conducto":
                    #Esto no debería ocurrir, porque los arroyos se crean antes que las esquinas
                    continue
                elif tipo == "esquina":
                    return i

            elif nodo[2] == "arroyo":
                if tipo == "arroyo":
                    return i
                elif tipo == "conducto":
                    return i
                elif tipo == "esquina":
                    return i

    # Si no hay nodos cercanos, agregar uno
    punto.append(tipo)
    nodos.append(punto)
    return len(nodos)-1;
