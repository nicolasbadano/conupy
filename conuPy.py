# conuPy.py encoding: utf-8

import os
import swmmout
import numpy as np
from collections import OrderedDict
from munch import Munch as Bunch

engine = "qgis" #"arcgis10" | "qgis_standalone"
if engine == "arcgis10":
    import engine_arcgis10 as gis
elif engine == "qgis":
    import engine_qgis as gis

from funciones import (dist, pairwise, insertPoints, removePoints, addNode,
    saveOnFile, readFromFile, atraviesaArroyo, intersection, distToSegment,
    intersect)
from nodesList import NodesList
from linksDict import LinksDict
from timing import timing

# Output fles
shpFileNodesDrainageNetwork = "nodesDrainageNetwork.shp"
shpFileNodesAlongNetwork    = "nodesAlongNetwork.shp"
shpFileNodesNetworkDepth    = "nodesNetworkDepth.shp"
shpFileNodosSample          = "nodosSample.shp"
shpFileNodos                = "nodos.shp"
shpFileCentros              = "centros.shp"
shpFileLineas               = "lineas.shp"
shpFileSumideros            = "sumideros.shp"
shpFileVertederos           = "vertederos.shp"
subcuencasShpFile           = "subcuencas.shp"
subcuencasClipShpFile       = "subcuencasClip.shp"
subcuencasAreaShpFile       = "subcuencasArea.shp"

# Parameters
params = {}
# Max dist for which channel and stream nodes are snapped together
params["maxDistSnapStreamNodes"] = 20.0
params["targetDx"] = 50.0
# Max length stream spans can have without been subdivided
params["maxLengthForStreamSpanDivide"] = params["targetDx"] * 4.0/3
# Min length stream spans can have without been joined
params["minLengthForStreamSpanJoin"] = params["targetDx"] * 2.0/3
# Max length street spans can have without been subdivided
params["maxLengthForStreetSpanDivide"] = params["targetDx"] * 4.0/3
# Min length street spans can have without been joined
params["minLengthForStreetSpanJoin"] = params["targetDx"] * 2.0/3
# Default coverage for conduits when no depths or levels were defined
params["minCoverage"] = 0.3
# Max distance required for gutters to be created from a corner to a conduit node
params["maxDistGutter"] = params["targetDx"] * 1.1
# Max distance required for weirs to be created from a corner to a channel segment
params["maxDistWeir"] = params["targetDx"] * 0.5
# Max dist for which outfall nodes are connected to regular nodes
params["maxDistConnectOutfallNodes"] = params["targetDx"]
# % of the basin area in which water is stored in a node
params["juApondPer"] = 0.75
# Water depth at start of simulation (ft or m) (default is 0).
params["juY0"] = 0
params["evap"] = 4.0 # mm/dia
# Manning's n for overland flow over the impervious sub-area.
params["saNImp"] = 0.025
# Manning's n for overland flow over the pervious sub-area.
params["saNPerv"] = 0.05
# depression storage for impervious sub-area (inches or mm).
params["saSImp"] = 0 #mm
# depression storage for pervious sub-area (inches or mm).
params["saSPerv"] = 0 #mm
# percent of impervious area with no depression storage.
params["saZero"] = 100 #%
# Infiltration Horton
params["inF0"] = 50.0 #mm/hr
params["inFf"] = 5.0 #mm/hr
params["inCoefDecaim"] = 2.0
# Time it takes for fully saturated soil to dry  (days).
params["inDryTime"] = 5.0
# Maximum infiltration volume possible (0 if not applicable) (in or mm)
params["inMaxInf"] = 0
# Value of n (i.e., roughness parameter) in Manning's equation.
params["coN"] = 0.04
# Weirs parameters
params["weAlturaCordon"] = 0.05
params["weCd"] = 3.0
# Orifice parameters
params["orCd"] = 0.65
# XSection parameters
params["outfallXsWidth"] = 50
params["xsG1"], params["xsG3"], params["xsG4"] = 10, 0, 0
params["xsSumideroH"], params["xsSumideroW"] = 0.15, 2.0
params["xsVertederoH"], params["xsVertederoW"] = 10.0, 20.0
# Transect parameters
params["traNConducto"] = 0.03
params["traNArroyoPlanicie"] = 0.05
params["traNArroyoCauce"] = 0.03
params["traNCalle"] = 0.04
params["traAnchoMargenArroyo"] = 20
params["traAnchoVereda"] = 2
params["traAltoCordon"] = 0.2
# 2d zones parameters
params["cellSize2D"] = 10
params["maxDist2DConnection"] = 1.5 * params["cellSize2D"]
# Number of discrete gages
params["numDiscreteGages"] = 21

def print_decorate(func):
    def wrapper(*args, **kwargs):
        print "INFO: ENTERED %s" % func.__name__
        ret = func(*args, **kwargs)
        print "INFO: EXITED %s" % func.__name__
        return ret
    return wrapper

@print_decorate
def mainCleanWorkspace(workspace):
    print "Limpiando el directorio de trabajo..."

    file_list = os.listdir(workspace)

    for f in file_list:
        try:
            os.remove(f)
        except:
            print "ERROR: " + f + " no pudo ser borrado."

    print "Finalizado el limpiado del directorio de trabajo."

@print_decorate
def mainPrepareDrainageNetwork(shpFileDrainageOriginal,
    shpFileDrainagePrepared, rasterFileDEM):
    # Read the original drainage network
    streams = gis.leer_shp_polilineas(shpFileDrainageOriginal,
        ['Ancho', 'Alto', 'Tipo', 'depthIni', 'depthFin', 'levelIni', 'levelFin'])
    spatial_ref = gis.leer_spatial_reference(shpFileDrainageOriginal)

    # Write a shape file with the nodes of the network
    nodesDrainageNetwork = []
    for stream in streams:
        nodesDrainageNetwork.append(stream[0][0])
        nodesDrainageNetwork.append(stream[0][-1])
    gis.escribir_shp_puntos(shpFileNodesDrainageNetwork, nodesDrainageNetwork, {}, spatial_ref)

    # Sample the terrain on the nodes
    nodesTerrainLevels = gis.sample_raster_on_nodes(shpFileNodesDrainageNetwork, rasterFileDEM)

    # Prepare the drainage network elevations base on available data for each node
    inode = 0
    for stream in streams:
        points, w, h, typ, depthIni, depthFin, levelIni, levelFin = stream
        if w is None:
            print "ERROR: Missing w value on stream"
            return
        if h is None:
            print "ERROR: Missing h value on stream"
            return
        if typ is None:
            print "ERROR: Missing type value on stream"
            return
        if str(typ).lower() in ["entubado", "conducto", "conduit"]:
            typ = "conduit"
        else:
            typ = "channel"
        if levelIni is None:
            if depthIni is None:
                if typ == "conduit":
                    depthIni = h + params["minCoverage"]
                else:
                    depthIni = h
            levelIni = nodesTerrainLevels[inode] - depthIni
        else:
            depthIni = nodesTerrainLevels[inode] - levelIni

        if levelFin is None:
            if depthFin is None:
                if typ == "conduit":
                    depthFin = h + params["minCoverage"]
                else:
                    depthFin = h
            levelFin = nodesTerrainLevels[inode+1] - depthFin
        else:
            depthFin = nodesTerrainLevels[inode+1] - levelFin

        stream[3], stream[4], stream[5], stream[6], stream[7] = typ, float(depthIni), float(depthFin), float(levelIni), float(levelFin)
        inode += 2

    # Calculate length of each stream
    for stream in streams:
        length = sum(dist(p0, p1) for p0, p1 in pairwise(stream[0]))
        stream.append(length)

    # Write the prepared drainage network
    polylines = [stream[0] for stream in streams]
    fields = OrderedDict()
    fields["w"]        = [stream[1] for stream in streams]
    fields["h"]        = [stream[2] for stream in streams]
    fields["type"]     = [str(stream[3]) for stream in streams]
    fields["depthIni"] = [stream[4] for stream in streams]
    fields["depthFin"] = [stream[5] for stream in streams]
    fields["levelIni"] = [stream[6] for stream in streams]
    fields["levelFin"] = [stream[7] for stream in streams]
    fields["slope"]    = [float((stream[6] - stream[7]) / stream[8]) for stream in streams]
    gis.escribir_shp_polilineas(shpFileDrainagePrepared, polylines, fields, spatial_ref)

    # Create points along streams
    pointsAlongNetwork = []
    levelsAlongNetwork = []
    for stream in streams:
        points, _, _, _, _, _, levelIni, levelFin, length = stream
        points = points[:]
        insertPoints(points, 50)

        prog = 0
        levels = [float(levelIni)]
        for p0, p1 in pairwise(points):
            prog += dist(p0,p1)
            levels.append(float(levelIni + prog / length * (levelFin - levelIni)))

        pointsAlongNetwork.extend(points)
        levelsAlongNetwork.extend(levels)

    # Write points for sampling
    gis.escribir_shp_puntos(shpFileNodesAlongNetwork, pointsAlongNetwork, {}, spatial_ref)
    # Sample the terrain on the nodes
    pointsTerrainLevels = gis.sample_raster_on_nodes(shpFileNodesAlongNetwork, rasterFileDEM)
    # Write depths of points along the network
    gis.escribir_shp_puntos(shpFileNodesNetworkDepth, pointsAlongNetwork, {
        "depth": [pointsTerrainLevels[i] - level for i, level in enumerate(levelsAlongNetwork)],
        "levelBot": levelsAlongNetwork,
        "levelTer": pointsTerrainLevels
    }, spatial_ref)


@print_decorate
def mainCreateSWMMModel(shpFileDrainagePrepared, shpFileCalles, shpFileCuenca,
    rasterFileDEM, rasterFileSlope, rasterFileImpermeabilidad,
    shpFileNodosBorde, gageMethod, gageFileName, rasterFileCoeficiente,
    gagesFileName, stationsFileName, swmmInputFileName):

    # Create nodes and links for the model
    nodos = NodesList()
    links = LinksDict(nodos)

    # Read spatial reference for the project
    spatial_ref = gis.leer_spatial_reference(shpFileCalles)

    # Read the drainage network and create nodes and links
    readDrainageNetwork(nodos, links, shpFileDrainagePrepared)
    # Read the street network and create nodes and links
    readStreets(nodos, links, shpFileCalles)
    # Create gutters and weirs
    createGutters(nodos, links)
    createWeirs(nodos, links)
    # Read 2d zones and create nodes and links
    # generate2dZones(nodos, links, "F:/Trabajo/Dropbox/Federico/ConuPy_version_27-4-2016/Ejemplo/Zona2D/zona2D.shp")
    # Calculate elevations
    calculateElevations(nodos, links, shpFileNodosSample, rasterFileDEM, rasterFileSlope, rasterFileImpermeabilidad, spatial_ref)

    # Create outfalls
    nodosOutfall, lineasOutfall = createOutfallNodes(nodos, shpFileNodosBorde)

    # Calculate dead depths
    calculateDeadDepths(nodos, links, lineasOutfall)

    # Write the network
    writeNetworkShapes(nodos, links, shpFileNodos, shpFileLineas, spatial_ref)

    # Create subcatchments
    centros, subcuencas = createSubcatchments(nodos, shpFileCuenca, spatial_ref)


    # Create raingages and map them to each subcatchment
    if gageMethod == "createRainGagesMethod0":
        gages, subcatchmentGages = createRainGagesMethod0(centros, gageFileName, rasterFileCoeficiente, gagesFileName)
    else:
        gages, subcatchmentGages = createRainGagesMethod1(centros, stationsFileName)

    saveOnFile(nodos, 'nodos')
    saveOnFile(links, 'links')
    saveOnFile(centros, 'centros')
    saveOnFile(lineasOutfall, 'lineasOutfall')

    # Write the model file
    writeSWMMFile(nodos, links, centros, subcuencas, nodosOutfall, lineasOutfall, gages, subcatchmentGages, swmmInputFileName)

@print_decorate
def readDrainageNetwork(nodos, links, shpFileDrainagePrepared):
    # Read the prepared drainage network
    streams = gis.leer_shp_polilineas(shpFileDrainagePrepared, ['w', 'h', 'type', 'depthIni', 'depthFin', 'levelIni', 'levelFin'])
    streams = [Bunch(points = stream[0],
                     w = float(stream[1]),
                     h = float(stream[2]),
                     typ = str(stream[3]).lower(),
                     levelIni = float(stream[6]),
                     levelFin = float(stream[7])) for stream in streams]

    # Subdivide the streams in spans shorter than maxLengthForStreamSpanDivide
    for stream in streams:
        insertPoints(stream.points, params["maxLengthForStreamSpanDivide"])

    snappedPoints = []
    # Snap the end nodes of each stream to nearby nodes of other streams
    for stream in streams:
        def snap(p, stream, streams):
            for targetStream in streams:
                # Doing stream == targetStream triggers problems when comparing numpy arrays
                if (stream.__dict__ == targetStream.__dict__):
                    continue
                mindist, minj = params["maxDistSnapStreamNodes"], -1
                for (j, p2) in enumerate(targetStream.points):
                    d = dist(p, p2)
                    if d < mindist:
                        mindist, minj = d, j
                if minj == -1:
                    continue
                # Convert coordinates to tuple because numpy arrays don't
                # work properly with the "in" operand
                snappedPoints.append(tuple(targetStream.points[minj]))
                return targetStream.points[minj]
            return p

        stream.points[0]  = snap(stream.points[0],  stream, streams)
        stream.points[-1] = snap(stream.points[-1], stream, streams)


    # Join streams spans shorter than minLengthForStreamSpanJoin unless they've been "snapped"
    for stream in streams:
        removePoints(stream.points, params["minLengthForStreamSpanJoin"], snappedPoints)

    for stream in streams:
        tipoTramo = "conduit" if stream.typ in ["entubado", "conducto", "conduit"] else "channel"

        nodesn = [addNode(nodos, p, tipoTramo) for p in stream.points]
        length = sum(dist(nodos[n0].p, nodos[n1].p) for n0, n1 in pairwise(nodesn))
        progFin = 0
        for n0, n1 in pairwise(nodesn):
            progIni = progFin
            progFin = progIni + dist(nodos[n0].p, nodos[n1].p)
            if n0 == n1:
                continue
            # Create a new link
            links[(n0, n1)] = {"type": tipoTramo,
                               "w": stream.w,
                               "h": stream.h,
                               "levelIni": stream.levelIni + (stream.levelFin - stream.levelIni) * progIni / length,
                               "levelFin": stream.levelIni + (stream.levelFin - stream.levelIni) * progFin / length}

    print "\tNumber of nodes for the drainage network: %i" % len(nodos)
    print "\tNumber of links for the drainage network: %i" % len(links)
    print "FINISHED: Drainage network construction"

@print_decorate
def readStreets(nodos, links, shpFileCalles):
    streets = gis.leer_shp_polilineas(shpFileCalles, ['ANCHO'])
    streets = [Bunch(points = street[0],
                     w = float(street[1])) for street in streets]

    for street in streets:
        if len(street.points) < 2:
            continue

        # Subdivide the street in spans shorter than maxLengthForStreetSpanDivide
        insertPoints(street.points, params["maxLengthForStreetSpanDivide"])
        # Join streams spans shorter than minLengthForStreetSpanJoin
        removePoints(street.points, params["minLengthForStreetSpanJoin"], [])

    for (i, street) in enumerate(streets):
        if i % 100 == 0:
            print "Procesando calle " + str(i) + " de " + str(len(streets))

        # Si el ancho es nulo
        if (street.w == 0):
            continue

        if len(street.points) < 2:
            continue

        for p0, p1 in pairwise(street.points):

            # Calculate
            # Crear u obtener los nodos extremos
            with timing.getTimer("addNode"):
                n0 = addNode(nodos, p0, "corner")
            with timing.getTimer("addNode"):
                n1 = addNode(nodos, p1, "corner")

            # Si la polilinea es demasiado corta
            if n0 == n1:
                continue

            # Verificar si la calle atraviesa un arroyo
            with timing.getTimer("atraviesaArroyo"):
                channel_data = atraviesaArroyo(nodos[n0].p, nodos[n1].p, nodos, links)

            if channel_data is not None:
                nch0, nch1, channel_link = channel_data
                p4 = intersection(nodos[n0].p, nodos[n1].p, nodos[nch0].p, nodos[nch1].p)

                nchannel = nch0 if dist(nodos[nch0].p, p4) <= dist(nodos[nch1].p, p4) else nch1

                # Dividir primero el arroyo si vale la pena (distancia al nodo más
                # cercano > a 15% maxLengthForStreamSpanDivide)
                if (min(dist(nodos[nch0].p, p4), dist(nodos[nch1].p, p4)) >
                    params["maxLengthForStreamSpanDivide"] * 0.15):

                    with timing.getTimer("addNode"):
                        nchannel = addNode(nodos, p4, "channel", 0)

                    alpha = dist(nodos[nch0].p, nodos[nchannel].p) / dist(nodos[nch0].p, nodos[nch1].p)
                    levelMid = (1 - alpha) * channel_link["levelIni"] + alpha * channel_link["levelFin"]

                    links[(nch0, nchannel)] = {"type": channel_link["type"],
                                               "w": channel_link["w"],
                                               "h": channel_link["h"],
                                               "levelIni": channel_link["levelIni"],
                                               "levelFin": levelMid}
                    links[(nchannel, nch1)] = {"type": channel_link["type"],
                                               "w": channel_link["w"],
                                               "h": channel_link["h"],
                                               "levelIni": levelMid,
                                               "levelFin": channel_link["levelFin"]}
                    del links[(nch0, nch1)]

                # Atraviesa un arroyo --> crear conexión entre el las dos esquinas y el arroyo
                if n0 != nchannel:
                    # Crear un nuevo vertedero
                    links[(n0, nchannel)] = {"type":"weir", "w":street.w}
                if n1 != nchannel:
                    # Crear un nuevo vertedero
                    links[(n1, nchannel)] = {"type":"weir", "w":street.w}
            else:
                # No atraviesa --> crear una calle comun
                links[(n0, n1)] = {"type":"street", "w":street.w}

@print_decorate
def createGutters(nodos, links):
    # Crear sumideros
    print "Creando sumideros"
    for (i, nodo) in enumerate(nodos):
        if nodo.type != "corner":
            continue

        with timing.getTimer("crearSumidero"):
            # Buscar el nodo conducto mas cercano
            nearINodes = nodos.getINodesNear(nodo.p, params["maxDistGutter"])
            mindist, minj = params["maxDistGutter"], -1
            for j in nearINodes:
                nodo2 = nodos[j]
                if nodo2.type == "conduit":
                    d = dist(nodo.p, nodo2.p)
                    if d < mindist:
                        mindist = d
                        minj = j
            if minj == -1:
                continue
            # Existe un nodo conducto cerca (< 80m)
            n0, n1 = i, minj
            # Si ya existe una conexión entre los nodos
            if (n0, n1) in links or (n1, n0) in links:
                continue
            # Crear un sumidero
            links[(n0, n1)] = {"type":"gutter", "w":params["xsSumideroW"]}

@print_decorate
def createWeirs(nodos, links):
    # Crear vertederos
    print "Creando vertederos"
    for (i, nodo) in enumerate(nodos):
        if nodo.type != "corner":
            continue

        with timing.getTimer("crearVertedero"):
            # Buscar el tramo de arroyo mas cercano
            nearLinks = links.getLinksNear(nodo.p, params["maxDistWeir"])
            mindist, minj = params["maxDistWeir"], -1
            for (n10, n11), link  in nearLinks.items():
                if link["type"] != "channel":
                    continue
                p10, p11 = nodos[n10].p, nodos[n11].p
                d = distToSegment(nodo.p, p10, p11)
                if d < mindist:
                    mindist = d
                    minj = n10 if dist(nodo.p, p10) <= dist(nodo.p, p11) else n11

            if minj == -1:
                continue

            # Existe un nodo arroyo cerca, conectar
            n0, n1 = i, minj
            # Si ya existe una conexión entre los nodos
            if (n0, n1) in links or (n1, n0) in links:
                continue
            # Crear un vertedero
            links[(n0, n1)] = {"type":"weir", "w":params["xsVertederoW"]}

    timing.dump()


@print_decorate
def generate2dZones(nodos, links, shpFile2dZones):
    minX, maxX, minY, maxY = [int(x) for x in gis.get_extent(shpFile2dZones)]
    junctions = []
    for x in range(minX, maxX, params["cellSize2D"]):
        for y in range(minY, maxY, params["cellSize2D"]):
            junctions.append(Bunch(p = np.array([x, y]), inside = False))
    gis.mark_points_inside_features(junctions, "F:/Trabajo/Dropbox/Federico/ConuPy_version_27-4-2016/Ejemplo/Zona2D/zona2D.shp")

    for junction in junctions:
        if junction.inside:
            addNode(nodos, junction.p, "2d")

    def atraviesaLinksTipo(p0, p1, tipos, nodos, links):
        for (n0, n1), link in links.getLinksInside(p0, p1).iteritems():
            if link["type"] in tipos:
                if (intersect(p0, p1, nodos[n0].p, nodos[n1].p)):
                    return True
        return False

    # Create 2D connections
    for y in range(minY, maxY, params["cellSize2D"]):
        for x0, x1 in pairwise(range(minX, maxX, params["cellSize2D"])):
            in0 = nodos.getINode(np.array([float(x0),float(y)]))
            in1 = nodos.getINode(np.array([float(x1),float(y)]))
            if in0 is None or in1 is None:
                continue
            if (in0, in1) in links or (in1, in0) in links:
                continue
            if atraviesaLinksTipo(nodos[in0].p,
                                    nodos[in1].p,
                                    ["channel", "street"],
                                    nodos,
                                    links):
                continue
            links[(in0, in1)] = {"type":"2d", "w":params["cellSize2D"]}

    for x in range(minX, maxX, params["cellSize2D"]):
        for y0, y1 in pairwise(range(minY, maxY, params["cellSize2D"])):
            in0 = nodos.getINode(np.array([float(x),float(y0)]))
            in1 = nodos.getINode(np.array([float(x),float(y1)]))
            if in0 is None or in1 is None:
                continue
            if (in0, in1) in links or (in1, in0) in links:
                continue
            if atraviesaLinksTipo(nodos[in0].p,
                                    nodos[in1].p,
                                    ["channel", "street"],
                                    nodos,
                                    links):
                continue
            links[(in0, in1)] = {"type":"2d", "w":params["cellSize2D"]}

    for (i, nodo) in enumerate(nodos):
        if nodo.type != "2d":
            continue

        # Buscar el nodo esquina o calle más cercano
        nearINodes = nodos.getINodesNear(nodo.p, params["maxDist2DConnection"])
        mindist, minj = params["maxDist2DConnection"], -1
        for j in nearINodes:
            nodo2 = nodos[j]
            if nodo2.type in ["corner", "channel"]:
                d = dist(nodo.p, nodo2.p)
                if d < mindist:
                    mindist = d
                    minj = j
        if minj == -1:
            continue
        # Existe un nodo conducto cerca (< 80m)
        n0, n1 = i, minj
        # Si ya existe una conexión entre los nodos
        if (n0, n1) in links or (n1, n0) in links:
            continue

        if atraviesaLinksTipo(nodo.p + 0.1 * (nodos[minj].p - nodo.p),
                              nodo.p + 0.9 * (nodos[minj].p - nodo.p),
                              ["channel", "street", "2d"],
                              nodos,
                              links):
            continue
        # Crear un sumidero
        links[(n0, n1)] = {"type":"2d", "w":params["cellSize2D"]}

@print_decorate
def calculateElevations(nodos, links, shpFileNodos, rasterFileDEM, rasterFileSlope, rasterFileImpermeabilidad, spatial_ref):
    gis.escribir_shp_puntos(shpFileNodos, [nodo.p for nodo in nodos], {}, spatial_ref)
    nodosElev = gis.sample_raster_on_nodes(shpFileNodos, rasterFileDEM)
    nodosSlope = gis.sample_raster_on_nodes(shpFileNodos, rasterFileSlope)
    nodosImpermeabilidad = gis.sample_raster_on_nodes(shpFileNodosSample, rasterFileImpermeabilidad)
    for i, nodo in enumerate(nodos):
        nodo.elev = float(nodosElev[i])
        nodo.slope = float(nodosSlope[i])
        nodo.imper = float(nodosImpermeabilidad[i])
        nodo.offset = 0
        nodo.length = 0

    for (n0, n1), link in links.iteritems():
        # Acumulate the length of segments incoming to the node
        if link["type"] in ["street", "channel", "2d"]:
            length = dist(nodos[n0].p, nodos[n1].p)
            nodos[n0].length += length
            nodos[n1].length += length

        # Update invert offset if appropiate
        if link["type"] in ["conduit", "channel"]:
            offset0 = link["levelIni"] - nodos[n0].elev
            offset1 = link["levelFin"] - nodos[n1].elev
            nodos[n0].offset = min(nodos[n0].offset, offset0)
            nodos[n1].offset = min(nodos[n1].offset, offset1)

    # Set elevations for streets, gutters and weirs
    for (n0, n1), link in links.iteritems():
        if link["type"] == "street":
            link["levelIni"] = nodos[n0].elev
            link["levelFin"] = nodos[n1].elev
        elif link["type"] in ["weir", "gutter"]:
            link["levelIni"] = nodos[n0].elev
            link["levelFin"] = nodos[n1].elev + nodos[n1].offset

            if link["type"] in "weir":
                if link["levelFin"] > link["levelIni"]:
                    print "WARNING: Weir %s - channel node invert is higher than the weir crest. This creates instabilties in SWMM."

@print_decorate
def writeNetworkShapes(nodos, links, shpFileNodos, shpFileLineas, spatial_ref):
    print u"Número de nodos: %d" % len(nodos)
    print u"Número de links: %d" % len(links)

    # Escribir shape con la posicion de los nodos
    campos = OrderedDict()
    campos["name"] = ['NODO%d' % i for i, _ in enumerate(nodos)]
    campos["type"] = [nodo.type for nodo in nodos]
    campos["elev"] = [float(nodo.elev) for nodo in nodos]
    campos["offs"] = [float(nodo.offset) for nodo in nodos]
    campos["inve"] = [float(nodo.elev + nodo.offset) for nodo in nodos]
    campos["dead"] = [float(nodo.tirante) for nodo in nodos]
    gis.escribir_shp_puntos(shpFileNodos, [nodo.p for nodo in nodos], campos, spatial_ref)
    # Escribir shape con los links
    polilineas = []
    campos = OrderedDict()
    campos["name"] = []
    campos["n0"] = []
    campos["n1"] = []
    campos["type"]  = []
    campos["w"] = []
    campos["elev0"] = []
    campos["elev1"] = []
    for i, (n0, n1) in enumerate(links):
        link = links[(n0, n1)]
        polilineas.append([nodos[n0].p, nodos[n1].p])
        campos["name"].append(str(link["type"]) + str(i))
        campos["n0"].append(int(n0))
        campos["n1"].append(int(n1))
        campos["type"].append(str(link["type"]))
        try:
            campos["w"].append(float(link.get("w", -1.0)))
            campos["elev0"].append(float(link["levelIni"]))
            campos["elev1"].append(float(link["levelFin"]))
        except:
            print "ERROR: Link (%d,%d) with type %s is missing w, levelIni or levelFin data." % (n0, n1, link["type"])
            print link
    gis.escribir_shp_polilineas(shpFileLineas, polilineas, campos, spatial_ref)

@print_decorate
def createSubcatchments(nodos, shpFileCuenca, spatial_ref):
    # Crear centros de cuencas
    centros = []
    for (i, nodo) in enumerate(nodos):
        if nodo.type != "conduit":
            centros.append([nodo.p[0], nodo.p[1], i])
    print "Numero de cuencas:\t%d" % len(centros)
    # Escribir shape con la posicion de los baricentros de subcuencas
    gis.escribir_shp_puntos(shpFileCentros, centros, {}, spatial_ref)

    # Escribir shape con la posicion de los baricentros de subcuencas
    gis.create_thiessen_polygons(shpFileCentros, subcuencasShpFile)
    gis.clip_feature(subcuencasShpFile, shpFileCuenca, subcuencasClipShpFile)

    subcuencas = []
    if engine == "arcgis10":
        gis.calculate_areas(subcuencasClipShpFile, subcuencasAreaShpFile)
        subcuencas = gis.leer_shp_poligonos(subcuencasAreaShpFile, ["Input_FID", "F_AREA"])

    elif engine == "qgis":
        areas = gis.read_areas(subcuencasClipShpFile)
        subcuencas = gis.leer_shp_poligonos(subcuencasClipShpFile, ["FID"])

        for i, subcuenca in enumerate(subcuencas):
            subcuenca.append(areas[i])

    subcuencasDict = {}
    for subcuenca in subcuencas:
        poligono, fid, area = subcuenca[0], subcuenca[1], subcuenca[2]
        subcuencasDict[fid] = [poligono, area]

    # Completar si falta alguna y eliminar duplicadas si existieran
    subcuencasCompletas = []
    for i in range(0,len(centros)):
        subcuencasCompletas.append(subcuencasDict.get(i, [[], 0]))

    return centros, subcuencasCompletas

@print_decorate
def createOutfallNodes(nodos, shpFileNodosBorde):
    posicionesNodosOutfall = gis.leer_shp_puntos(shpFileNodosBorde)

    nodosOutfall = []
    lineasOutfall = []
    for (i, nodo) in enumerate(nodos):
        # Buscar la posicion outfall mas cercana
        mindist, minj = params["maxDistConnectOutfallNodes"], -1
        for (j, posNO) in enumerate(posicionesNodosOutfall):
            if dist(nodo.p, posNO) < mindist:
                mindist = dist(nodo.p, posNO)
                minj = j
        if minj == -1:
            continue

        nodoOutfall = Bunch(p = np.array(posicionesNodosOutfall[minj]),
                            elev = nodo.elev,
                            offset = nodo.offset)
        nodosOutfall.append(nodoOutfall)

        lineasOutfall.append( [i, len(nodosOutfall)-1, params["outfallXsWidth"]] )

    return nodosOutfall, lineasOutfall

@print_decorate
def createRainGagesMethod0(centros, gageFileName, rasterFileCoeficiente, gagesFileName):
    print "Leyendo pluviometro maestro..."

    with open(gageFileName, "r") as iFile:
        lineas = iFile.readlines()


    print "Creando estaciones..."
    gages = []
    with open(gagesFileName, "w") as tF:
        for i in xrange(0, params["numDiscreteGages"]):
            coef = i * (100/(params["numDiscreteGages"]-1))

            gage = {}
            gage["name"] = 'GAGE'+str(coef)
            _, gage["file"] = os.path.split(gagesFileName)
            gage["interval"] = '0:10'
            gages.append(gage)

            for linea in lineas:
                datos = linea.split()
                if len(datos) == 0:
                    continue
                datos[0]  = gage["name"]
                datos[-1] = float(datos[-1])*(float(coef)/100.0)
                tF.write(("").join([ str(x).ljust(15, ' ') for x in datos]))
                tF.write("\n")

    print "Leyendo mapa de decaimiento"
    nodosDecaimiento = gis.sample_raster_on_nodes(shpFileNodos, rasterFileCoeficiente)

    print "Seleccionando pluviómetro para cada subcuenca"
    subcatchmentGages = []
    for (i, centro) in enumerate(centros):
        coef = int(nodosDecaimiento[centro[2]])
        gageName = 'GAGE' + str(coef - (coef%(100/(params["numDiscreteGages"]-1))))
        subcatchmentGages.append(gageName)

    return gages, subcatchmentGages

@print_decorate
def createRainGagesMethod1(centros, stationsFileName):
    print "Leyendo lista de pluviometros..."
    gages = []
    with open(stationsFileName, "r") as f:
        for i, line in enumerate(f):
            data = line.split()
            gage = Bunch()
            gage.coord = np.array([float(data[0]), float(data[1])])
            gage.name = data[2]
            gage.file = data[3]
            gage.interval = '0:05'
            gages.append(gage)

    print "Seleccionando pluviómetro para cada subcuenca..."
    subcatchmentGages = []
    for (i, centro) in enumerate(centros):
        minDist, minGage = 1e10, None

        for gage in gages:
            gDist = dist(gage.coord, np.array(centro[0:2]))

            if gDist < minDist:
                minGage, minDist = gage, gDist

        subcatchmentGages.append(minGage.name)

    return gages, subcatchmentGages

@print_decorate
def writeSWMMFile(nodos, links, centros, subcuencas, nodosOutfall, lineasOutfall, gages, subcatchmentGages, swmmInputFileName):
    for nodo in nodos:
        nodo.area = 1.167
    for (i, centro) in enumerate(centros):
        nodos[centro[2]].area = params["juApondPer"] * subcuencas[i][1]

    with open(swmmInputFileName, "w") as tF:
        tF.write("[TITLE]\n")
        tF.write("Conurbano\n")

        tF.write("\n")
        tF.write("[OPTIONS]\n")
        tF.write("FLOW_UNITS           CMS\n")
        tF.write("INFILTRATION         HORTON\n") #CURVE_NUMBER
        tF.write("FLOW_ROUTING         DYNWAVE\n")
        tF.write("LINK_OFFSETS         ELEVATION\n")
        tF.write("START_DATE           01/01/2017\n")
        tF.write("START_TIME           00:00\n")
        tF.write("END_DATE             01/02/2017\n")
        tF.write("END_TIME             06:00\n")
        tF.write("WET_STEP             00:00:10\n")
        tF.write("DRY_STEP             00:00:10\n")
        tF.write("ROUTING_STEP         00:00:01\n")
        tF.write("ALLOW_PONDING        YES\n")
        tF.write("INERTIAL_DAMPING     FULL\n")
        tF.write("MAX_TRIALS           100\n")
        tF.write("HEAD_TOLERANCE       0.001\n")
        tF.write("MINIMUM_STEP         1\n")
        tF.write("THREADS              4\n")
        # Minimum Surface Area - This is a minimum surface area used at nodes when computing changes in water depth. If 0 is entered, then the default value of 12.566 ft2 (1.167 m2) is used. This is the area of a 4-ft diameter manhole. The value entered should be in square feet for US units or square meters for SI units.
        tF.write("MIN_SURFAREA         1.167\n") #m2, equivalente a 10 m de diametro
        tF.write("\n")

        tF.write("[FILES]\n")
        tF.write("SAVE RAINFALL        rainfall.rff\n")
        tF.write("SAVE RUNOFF          runoff.rof\n")
        tF.write("SAVE OUTFLOWS        outflows.txt\n")


        tF.write("\n")
        tF.write("[RAINGAGES]\n")
        tF.write(";;Name         Format         Interval       SCF            DataSrc        SourceName     Sta            Units\n")
        tF.write(";;============================================================================================================\n")
        for gage in gages:
            list = [gage["name"], 'INTENSITY', gage["interval"], 1.0, 'FILE', gage["file"], gage["name"], "MM"]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")

        tF.write("\n")
        tF.write("[EVAPORATION]\n")
        tF.write("CONSTANT  %.3f\n" % params["evap"])

        tF.write("\n")
        tF.write("[SUBCATCHMENTS]\n")
        tF.write(";;Name         Raingage       Outlet         Area           %ImperV        Width          Slope          Curve Length\n")
        tF.write(";;======================================================================================================================\n")
        for (i, centro) in enumerate(centros):
            n1 = centro[2]
            nodo = nodos[n1]
            gageName = subcatchmentGages[i]
            list = ['CUENCA'+str(i), gageName, 'NODO'+str(n1), "%.3f" % (float(subcuencas[i][1])/10000.0), "%.3f" % nodo.imper, "%.3f" % (nodo.length/2), "%.3f" % (nodo.slope), "%.3f" % (subcuencas[i][1]**0.5)]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")

        tF.write("\n")
        tF.write("[SUBAREAS]\n")
        tF.write(";;Name         N_Imp          N_Perv         S_Imp          S_Perv         %ZER           RouteTo\n")
        tF.write(";;=======================================================================================================\n")
        for (i, centro) in enumerate(centros):
            list = ['CUENCA'+str(i), params["saNImp"], params["saNPerv"], params["saSImp"], params["saSPerv"], params["saZero"], 'OUTLET']
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")

        tF.write("\n")
        tF.write("[INFILTRATION]\n")
        tF.write(";;Subcat       MaxRate        MinRate        Decay          DryTime        Max Inf\n")
        tF.write(";;========================================================================================\n")
        for (i, centro) in enumerate(centros):
            list = ['CUENCA'+str(i), params["inF0"], params["inFf"], params["inCoefDecaim"], params["inDryTime"], params["inMaxInf"]]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")


        # Depth from ground to invert elevation (ft or m) (default is 0).
        # juYmax = 1.0
        # maximum additional head above ground elevation that manhole junction can sustain under surcharge conditions (ft or m) (default is 0).
        # juYsur = 0
        # % del area de la subcuenca en que se almacena el agua en el nodo
        params["juApondPer"] = 0.75
        tF.write("\n")
        tF.write("[JUNCTIONS]\n")
        tF.write(";;Name         Elev           Ymax           Y0             Ysur           Apond) \n")
        tF.write(";;========================================================================================\n")
        for (i, nodo) in enumerate(nodos):
            if not nodo.type in ["conduit", "2d"]:
                continue
            list = ['NODO%d' % i,
                    "%.3f" % (nodo.elev + nodo.offset),
                    "%.3f" % (-nodo.offset + 20),
                    "%.3f" % params["juY0"],
                    "%.3f" % 0,
                    "%.3f" % nodo.area]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")


        tF.write("\n")
        tF.write("[STORAGE]\n")
        tF.write(";;Name         Elev           Ymax           Y0             TABULAR        Apond          ) \n")
        tF.write(";;========================================================================================\n")
        for (i, nodo) in enumerate(nodos):
            if nodo.type in ["conduit", "2d"]:
                continue
            tF.write(("").join([ str(x).ljust(15, ' ') for x in [
                'NODO%d' % i,
                "%.3f" % (nodo.elev + nodo.offset),
                "%.3f" % (-nodo.offset+20),
                "%.3f" % params["juY0"],
                'TABULAR',
                'STORAGE%d' % i,
                "%.3f" % nodo.area]]))
            tF.write("\n")


        tF.write("\n")
        tF.write("[CURVES]\n")
        tF.write(";;Name         Type           x-value        y-value\n")
        tF.write(";;========================================================================================\n")
        for (i, nodo) in enumerate(nodos):
            if (nodo.offset < 0):
                list = [
                    'STORAGE%d' % i,
                    'STORAGE',
                    "%.3f" % 0, 1.167,
                    "%.3f" % (-nodo.offset),      "%.3f" % 1.167,
                    "%.3f" % (-nodo.offset + 1),  "%.3f" % (nodo.area/10),
                    "%.3f" % (-nodo.offset + 2),  "%.3f" % nodo.area,
                    "%.3f" % (-nodo.offset + 20), "%.3f" % nodo.area]
            else:
                list = [
                    'STORAGE%d' % i,
                    'STORAGE',
                    "%.3f" % 0,                   "%.3f" % 1.167,
                    "%.3f" % 1,                   "%.3f" % (nodo.area/10),
                    "%.3f" % 2,                   "%.3f" % nodo.area,
                    "%.3f" % 20,                  "%.3f" % nodo.area]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")


        tF.write("\n")
        tF.write("[OUTFALLS]\n")
        tF.write(";;Name         Elev           Type           Gate\n")
        tF.write(";;==========================================================\n")
        for (i, nodo) in enumerate(nodosOutfall):
            list = [
                'NODOOUT'+str(i),
                nodo.elev + nodo.offset,
                "FREE",
                "NO"]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")

        tF.write("\n")
        tF.write("[CONDUITS]\n")
        tF.write(";;Name         Node1          Node2          Length         N              Z1             Z2             Q0\n")
        tF.write(";;======================================================================================================================\n")
        for i, ((in0, in1), link) in enumerate(links.iteritems()):
            name = link["type"] + str(i)
            length = dist(nodos[in0].p, nodos[in1].p)
            if link["type"] == "street":
                list = [
                    name,
                    'NODO'+str(in0),
                    'NODO'+str(in1),
                    "%.3f" % length,
                    "%.3f" % params["coN"],
                    "%.3f" % nodos[in0].elev,
                    "%.3f" % nodos[in1].elev,
                    0]
            elif link["type"] == "2d":
                list = [
                    name,
                    'NODO'+str(in0),
                    'NODO'+str(in1),
                    "%.3f" % length,
                    "%.3f" % params["coN"],
                    "%.3f" % nodos[in0].elev,
                    "%.3f" % nodos[in1].elev,
                    0]
            elif link["type"] in ["channel", "conduit"]:
                list = [
                    name,
                    'NODO'+str(in0),
                    'NODO'+str(in1),
                    "%.3f" % length,
                    "%.3f" % params["coN"],
                    "%.3f" % link["levelIni"],
                    "%.3f" % link["levelFin"],
                    0]
            else:
                continue
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")
        for (i, linea) in enumerate(lineasOutfall):
            in0, in1, ancho = linea
            name = 'SALIDA' + str(i)
            length = dist(nodos[in0].p, nodosOutfall[in1].p)
            list = [
                name,
                'NODO'+str(in0),
                'NODOOUT'+str(in1),
                "%.3f" % length,
                "%.3f" % params["coN"],
                "%.3f" % (nodos[in0].elev + nodos[in0].offset),
                "%.3f" % (nodos[in0].elev + nodos[in0].offset),
                0]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")


        tF.write("\n")
        tF.write("[WEIRS]\n")
        tF.write(";;Name           From Node        To Node          Type         CrestHt    Qcoeff     Gated    EndCon   EndCoeff   Surcharge \n")
        tF.write(";;-------------- ---------------- ---------------- ------------ ---------- ---------- -------- -------- ---------- ----------\n")
        for i, ((in0, in1), link) in enumerate(links.iteritems()):
            if link["type"] != "weir":
                continue
            name = "weir" + str(i)
            list = [
                name,
                'NODO'+str(in0),
                'NODO'+str(in1),
                "SIDEFLOW",
                "%.3f" % (nodos[in0].elev + params["weAlturaCordon"]),
                params["weCd"],
                "NO",
                0,
                0,
                "NO"]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")


        tF.write("\n")
        tF.write("[ORIFICES]\n")
        tF.write(";;Name         Node1          Node2          Type           Offset         Cd             Flap           Orate\n")
        tF.write(";;======================================================================================================================\n")
        for i, ((in0, in1), link) in enumerate(links.iteritems()):
            if link["type"] != "gutter":
                continue
            name = "gutter" + str(i)
            list = [
                name,
                'NODO'+str(in0),
                'NODO'+str(in1),
                "SIDE",
                "%.3f" % nodos[in0].elev,
                params["weCd"],
                "NO",
                0]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")


        tF.write("\n")
        tF.write("[XSECTIONS]\n")
        tF.write(";;Link         Type           G1             G2             G3             G4\n")
        tF.write(";;========================================================================================\n")
        transectas = {}
        for i, ((in0, in1), link) in enumerate(links.iteritems()):
            name = link["type"] + str(i)

            if link["type"] == "street":
                tname = link["type"] + str(int(link["w"]))
                transectas[tname] = [link["type"], tname, link["w"], link.get("h",0)]
                list = [name, 'IRREGULAR', tname, 0, 0, 0]
            elif link["type"] == "2d":
                list = [name, 'RECT_OPEN', 100, 20.0, 2, 0, 1]
            elif link["type"] == "channel":
                tname = link["type"] + str(int(link["w"])) + "x" + str(int(link["h"]))
                transectas[tname] = [link["type"], tname, link["w"], link.get("h",0)]
                list = [name, 'IRREGULAR', tname, 0, 0, 0]
            elif link["type"] == "conduit":
                list = [name, 'RECT_CLOSED', link["h"], link["w"], 0, 0]
            elif link["type"] == "weir":
                list = [link["type"]+str(i), 'RECT_OPEN', params["xsVertederoH"], params["xsVertederoW"], 0, 0]
            elif link["type"] == "gutter":
                list = [link["type"]+str(i), 'RECT_CLOSED', params["xsSumideroH"], params["xsSumideroW"], 0, 0]
            else:
                continue
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")
        for (i, linea) in enumerate(lineasOutfall):
            in0, in1, ancho = linea
            list = ['SALIDA'+str(i), 'RECT_OPEN', params["xsG1"], ancho, params["xsG3"], params["xsG4"]]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")


        tF.write("\n")
        tF.write("[TRANSECTS]\n")
        tF.write(";;NC           Nleft          Nright         Nchanl         x              x \n")
        tF.write(";;X1           Name           Nsta           Xleft          Xright         0 0 0 0 0 0 \n")
        tF.write(";;GR           Elev           GR             Elev           GR             Elev....\n")
        tF.write(";;========================================================================================\n")
        for key, value in transectas.items():
            tipo, tname, ancho, alto = value
            if (tipo == "channel"):
                traTiranteArroyo = alto
                list = ['NC', params["traNArroyoPlanicie"], params["traNArroyoPlanicie"], params["traNArroyoCauce"]]
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")
                list = ['X1', tname, 6, -ancho * 0.5, ancho * 0.5, 0, 0, 0, 0, 0, 0]
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")
                list = ['GR',
                        traTiranteArroyo + 20, -params["traAnchoMargenArroyo"] - ancho*0.5,
                        traTiranteArroyo, -params["traAnchoMargenArroyo"] - ancho*0.5,
                        traTiranteArroyo, -ancho*0.5,
                        0, -ancho*0.5 + 0.25,
                        0, ancho*0.5 - 0.25,
                        traTiranteArroyo, ancho*0.5,
                        traTiranteArroyo, params["traAnchoMargenArroyo"] + ancho*0.5,
                        traTiranteArroyo + 20, params["traAnchoMargenArroyo"] + ancho*0.5]
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")
                tF.write(";;-------------------------------------------\n")
            # elif (tipo == "conduit"):
                # list = ['NC', params["traNConducto"], params["traNConducto"], params["traNConducto"]]
                # tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                # tF.write("\n")
                # list = ['X1', tname, 6, -ancho * 0.5, ancho * 0.5, 0, 0, 0, 0, 0, 0]
                # tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                # tF.write("\n")
                # list = ['GR', -ancho*0.5, alto, -ancho*0.5, 0, ancho*0.5, 0, ancho*0.5, alto]
                # tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                # tF.write("\n")
                # tF.write(";;-------------------------------------------\n")
            elif (tipo == "street"):
                list = ['NC', params["traNCalle"], params["traNCalle"], params["traNCalle"]]
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")
                list = ['X1', tname, 7, -ancho * 0.5, ancho * 0.5, 0, 0, 0, 0, 0, 0]
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")
                list = ['GR',
                        params["traAltoCordon"] + 20, -params["traAnchoVereda"] - ancho*0.5,
                        params["traAltoCordon"], -params["traAnchoVereda"] - ancho*0.5,
                        params["traAltoCordon"], -ancho*0.5,
                        0, -ancho*0.5,
                        params["traAltoCordon"], 0,
                        0, ancho*0.5,
                        params["traAltoCordon"], ancho*0.5,
                        params["traAltoCordon"], params["traAnchoVereda"] + ancho*0.5,
                        params["traAltoCordon"] + 20, params["traAnchoVereda"] + ancho*0.5]
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")
                tF.write(";;-------------------------------------------\n")


        tF.write("\n")
        tF.write("[REPORT]\n")
        tF.write("CONTINUITY               YES\n")
        tF.write("FLOWSTATS                YES\n")
        tF.write("CONTROLS                 NO\n")
        tF.write("SUBCATCHMENTS            ALL\n")
        tF.write("NODES                    ALL\n")
        tF.write("LINKS                    ALL\n")


        tF.write("\n\n\n\n")
        tF.write("[MAP]\n")
        minx = min(nodo.p[0] for nodo in nodos)
        maxx = max(nodo.p[0] for nodo in nodos)
        miny = min(nodo.p[1] for nodo in nodos)
        maxy = max(nodo.p[1] for nodo in nodos)
        distx, disty = maxx - minx, maxy - miny
        tF.write(";;             minx           miny           maxx           maxy\n")
        tF.write(";;=========================================================================\n")
        list = ['DIMENSIONS', minx - 0.1 * distx, miny - 0.1 * disty, maxx + 0.1 * distx, maxy + 0.1 * disty]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\nUNITS METERS\n")

        tF.write("\n")
        tF.write("[COORDINATES]\n")
        tF.write(";;Node         xcoord         ycoord         \n")
        tF.write(";;===========================================\n")
        for (i, nodo) in enumerate(nodos):
            list = ['NODO%d' % i, nodo.p[0], nodo.p[1]]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")
        for (i, nodo) in enumerate(nodosOutfall):
            list = ['NODOOUT%d' % i, nodo.p[0], nodo.p[1]]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")

        tF.write("\n")
        tF.write("[VERTICES]\n")
        tF.write(";;LInk         xcoord         ycoord         \n")
        tF.write(";;===========================================\n")

        tF.write("\n")
        tF.write("[POLYGONS]\n")
        tF.write(";;Subcat       xcoord         ycoord         \n")
        tF.write(";;===========================================\n")
        for (i, subcuenca) in enumerate(subcuencas):
            for punto in subcuenca[0]:
                list = ['CUENCA'+str(i), punto[0], punto[1]]
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")

        tF.write("\n")
        tF.write("[SYMBOLS]\n")
        tF.write(";;Gage         xcoord         ycoord         \n")
        tF.write(";;===========================================\n")
        list = ['GAGE1', 0.5*(minx+maxx), 0.5*(miny+maxy)]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")

@print_decorate
def mainReadSWMMResultsDepths(swmmOuputFileName):
    nodos = readFromFile('nodos')

    print "Numero de nodos: ", len(nodos)

    outfile = swmmout.open(swmmOuputFileName)

    query_vars = ['head']
    query_nodes = ['NODO%d' % i for i, nodo in enumerate(nodos)]
    data = outfile.get_values('nodes', query_nodes, query_vars)

    # Conseguir la referencia geografica
    spatial_ref = gis.leer_spatial_reference(shpFileNodos)

    # Escribir shape con los niveles y profundidades en cada paso de tiempo
    for i, dataline in enumerate(data):
        campos = OrderedDict()
        campos["elev"]  = [dataline[i+1][1] for i, nodo in enumerate(nodos)]
        campos["depth"] = [max(dataline[i+1][1] - (nodo.elev + nodo.offset), 0) for i, nodo in enumerate(nodos)]
        gis.escribir_shp_puntos("nodeDepth%04d.shp" % i, [nodo.p for nodo in nodos], campos, spatial_ref)

    for i, nodo in enumerate(nodos):
        nodo.maxHead = max([dataline[i+1][1] for dataline in data])

    # Escribir shape con los niveles y profundidades máximas
    campos = OrderedDict()
    campos["elev"]  = [nodo.maxHead for nodo in nodos]
    campos["depth"] = [max(nodo.maxHead - (nodo.elev + nodo.offset), 0) for nodo in nodos]
    gis.escribir_shp_puntos(workspace + "/" + "nodeDepthMax.shp", [nodo.p for nodo in nodos], campos, spatial_ref)

@print_decorate
def calculateDeadDepths(nodos, links, lineasOutfall):
    for nodo in nodos:
        nodo.tirante = 100
    for linea in lineasOutfall:
        nodo0 = linea[0]
        nodos[nodo0].tirante = 0

    i = 0
    maxbajada = 100
    while maxbajada > 0.01:
        maxbajada = 0

        def igualar(nodo0, nodo1, link):
            linklevel = max(link["levelIni"], link["levelFin"])

            elev0 = nodos[nodo0].elev + nodos[nodo0].offset + nodos[nodo0].tirante
            elev1 = nodos[nodo1].elev + nodos[nodo1].offset + nodos[nodo1].tirante

            newEle0 = min(elev0, max(linklevel, elev1))
            newEle1 = min(elev1, max(linklevel, elev0))

            bajada = max(elev0 - newEle0, elev1 - newEle1)

            nodos[nodo0].tirante = newEle0 - (nodos[nodo0].elev + nodos[nodo0].offset)
            nodos[nodo1].tirante = newEle1 - (nodos[nodo1].elev + nodos[nodo1].offset)

            return bajada

        for n0, n1 in links:
            link = links[(n0, n1)]
            bajada = igualar(n0, n1, link)
            maxbajada = max(maxbajada, bajada)

        i += 1
        print "Iteracion %i - Max Bajada %f" % (i, maxbajada)


if __name__ == '__main__':
    # Origen de datos
    dataFolder = "F:/Desarrollo/Utilidades/conuPy/Ejemplos/CONUPY_cuenca completa1/"
    # Archivos de entrada
    shpFileDrainageOriginal     = dataFolder + "Arroyos/Arroyos2.shp"
    shpFileDrainagePrepared     = dataFolder + "Arroyos/Arroyos2Preparados.shp"
    shpFileCalles               = dataFolder + "Calles/Calles_Recortado2.shp"
    shpFileCuenca               = dataFolder + "Cuenca/Cuencas_Laferrere.shp"
    shpFileNodosBorde           = dataFolder + "NodosBorde/NodosBorde.shp"
    rasterFileDEM               = dataFolder + "MDT/MDT"
    rasterFileSlope             = dataFolder + "Pendientes/Pendientes"
    rasterFileCoeficiente       = dataFolder + "Lluvias/lluvias.img"
    rasterFileImpermeabilidad   = dataFolder + "Impermeabilidad/Imp"
    # workspace temporal
    workspace                   = dataFolder + "ws"
    # Directorio del modelo
    modelFolder                 = dataFolder + "Modelo/"
    swmmInputFileName           = modelFolder + "conurbano.inp"
    defaultGageFileName         = modelFolder + "lluvia.in"
    defaultGagesFileName        = modelFolder + "pluviom.dat"
    swmmOuputFileName           = modelFolder + "conurbano.out"
    gageMethod                  = "createRainGagesMethod0"
    gageFileName                = ""
    rasterFileCoeficiente       = ""
    gagesFileName               = ""
    stationsFileName            = ""
    swmmInputFileName           = ""

    os.chdir(workspace)
    gis.gis_init()

    # Opciones de corrida
    print "Que desea hacer?"
    print " 0 - Preparar red de drenaje"
    print " 1 - Crear modelo SWMM"
    print " 2 - Leer resultados y escribir shp con profundidad y elevacion en nodos"
    x = input("Opcion:")
    if (x == 0):
        mainPrepareDrainageNetwork(shpFileDrainageOriginal, shpFileDrainagePrepared, rasterFileDEM)
    elif (x == 1):
        mainCreateSWMMModel(shpFileDrainagePrepared, shpFileCalles, shpFileCuenca,
            rasterFileDEM, rasterFileSlope, rasterFileImpermeabilidad,
            shpFileNodosBorde, gageMethod, gageFileName, rasterFileCoeficiente,
            gagesFileName, stationsFileName, swmmInputFileName)
    elif (x == 2):
        mainReadSWMMResultsDepths(swmmOuputFileName)
