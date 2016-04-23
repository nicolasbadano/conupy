# conuPy.py encoding: utf-8

engine = "qgis" #"arcgis10" | "qgis_standalone"

import sys, os
from collections import OrderedDict

if engine == "arcgis10":
    from engine_arcgis10 import *
elif engine == "qgis":
    from engine_qgis import *
elif engine == "qgis_standalone":
    from engine_qgis import *
    import inspect
    print inspect.getfile(inspect.currentframe()) # script filename (usually with path)
    print os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
    print os.path.dirname(inspect.getfile(inspect.currentframe()))
    sys.path.append(os.path.dirname(inspect.getfile(inspect.currentframe())))

from funciones import *
import swmmout

# Origen de datos
dataFolder = "F:/Desarrollo/Utilidades/conuPy/Ejemplos/CONUPY_cuenca completa1/"

# workspace temporal
workspace                        = dataFolder + "ws"
# Directorio del modelo
modelFolder                      = dataFolder + "Modelo/"

# Archivos de entrada
defaultShpFileArroyos            = dataFolder + "Arroyos/Arroyos2.shp"
defaultShpFileCalles             = dataFolder + "Calles/Calles_Recortado2.shp"
defaultShpFileCuenca             = dataFolder + "Cuenca/Cuencas_Laferrere.shp"
defaultShpFileNodosBorde         = dataFolder + "NodosBorde/NodosBorde.shp"

defaultRasterFileDEM             = dataFolder + "MDT/MDT"
defaultRasterFileSlope           = dataFolder + "Pendientes/Pendientes"
defaultRasterFileCoeficiente     = dataFolder + "Lluvias/lluvias.img"
defaultRasterFileImpermeabilidad = dataFolder + "Impermeabilidad/Imp"

swmmCorregirArroyosFileName  = modelFolder + "conurbano_3.inp"

# Output fles
shpFileNodos                = "nodos.shp"
shpFileCentros              = "centros.shp"
shpFileLineas               = "lineas.shp"
shpFileSumideros            = "sumideros.shp"
shpFileVertederos           = "vertederos.shp"
subcuencasShpFile           = "subcuencas2.shp"
subcuencasClipShpFile       = "subcuencasClip.shp"
subcuencasAreaShpFile       = "subcuencasArea.shp"

defaultSwmmInputFileName    = modelFolder + "conurbano.inp"
defaultGageFileName         = modelFolder + "lluvia.in"
defaultGagesFileName        = modelFolder + "pluviom.dat"
defaultSwmmOuputFileName    = modelFolder + "conurbano.out"

# Numero de Gages
numGages = 21


def mainCleanWorkspace(workspace):
    print "Limpiando el directorio de trabajo..."

    file_list = os.listdir(workspace)

    for f in file_list:
        try:
            os.remove(f)
        except Exception, err:
            print "\t" + f + " no pudo ser borrado."

    print "Finalizado el limpiado del directorio de trabajo."


def mainReadRivers(shpFileArroyos):
    print "Proceso de creado de arroyos..."

    # Leer arroyos
    arroyos = leer_shp_polilineas(shpFileArroyos, ['Ancho', 'Alto', 'Tipo'])

    # Dividir los subtramos en longitudes menores a 100m
    for arroyo in arroyos:
        insertPoints(arroyo[0], 100)

    # Hacer coincidir los nodos extremos de los arroyos con nodos intermedios cercanos de otros arroyos
    for arroyo in arroyos:
        def snap(p, arroyo, arroyos):
            for arroyo2 in arroyos:
                if (arroyo != arroyo2):
                    mindistSq, minj = 50*50, -1
                    for (j, p2) in enumerate(arroyo2[0]):
                        dSq = distSq(p, p2)
                        if dSq < mindistSq:
                            mindistSq = dSq
                            minj = j
                    if minj != -1:
                        arroyo2[0][minj].append("snapped")
                        return arroyo2[0][minj]
            return p

        arroyo[0][0]  = snap(arroyo[0][0],  arroyo, arroyos)
        arroyo[0][-1] = snap(arroyo[0][-1], arroyo, arroyos)


    # Unir subtramos consecutivos de longitudes menores a 75m que cuyo nodo central no es "snapped"
    for arroyo in arroyos:
        removePoints(arroyo[0], 75)

    for arroyo in arroyos:
        puntos, ancho, alto, tipo = arroyo[0], arroyo[1], arroyo[2], arroyo[3]
        for punto in puntos:
            if punto[-1] == "snapped":
                del punto[-1]

    # Crear nodos arroyos y lineas arroyo
    nodos = []
    lineas = []
    for arroyo in arroyos:
        puntos, ancho, alto, tipo = arroyo[0], arroyo[1], arroyo[2], arroyo[3]
        tipoTramo = "arroyo" if (tipo == "Canal") else "conducto"

        punto = puntos[0]
        puntoInicial = addNode(nodos, punto, tipoTramo)
        for i in xrange(1, len(puntos)):
            punto = puntos[i]
            puntoFinal = addNode(nodos, punto, tipoTramo)
            if ( puntoInicial != puntoFinal ):
                 lineas.append( [puntoInicial, puntoFinal, tipoTramo, ancho, alto] )
            puntoInicial = puntoFinal

    print "Numero de nodos: ", len(nodos)
    print "Numero de lineas: ", len(lineas)

    # Write list files
    saveOnFile(nodos, "nodosArroyo")
    saveOnFile(lineas, "lineasArroyo")

    print "Finalizado proceso de creado de arroyos"


def mainReadStreets(shpFileCalles):
    print "Proceso de creado de calles"
    tF = open("log", "w")
    tF.write("Proceso de creado de calles\n")

    spatial_ref = leer_spatial_reference(shpFileCalles)
    calles = leer_shp_polilineas(shpFileCalles, ['ANCHO'])

    # Crear nodos esquina y lineas calle
    nodos = readFromFile('nodosArroyo')
    lineas = readFromFile('lineasArroyo')
    for (i,calle) in enumerate(calles):
        if i % 100 == 0:
            print "Procesando calle " + str(i) + " de " + str(len(calles))
            tF.write("Procesando calle " + str(i) + " de " + str(len(calles)) + "\n")

        puntos, ancho = calle[0], calle[1]
        if (ancho == 0):
            continue

        if len(puntos) < 2:
            print "Error"
            print puntos

        p = atraviesaArroyo(puntos[0], puntos[-1], nodos, lineas)
        if (p != -1):
            # Atraviesa un arroyo --> crear dos semicalles hacia el nodo arroyo

            #puntoInicial = addNode(nodos, puntos[0], "esquina")
            puntoFinal   = p
            #if ( puntoInicial != puntoFinal ):
            #    lineas.append( [puntoInicial, puntoFinal, "calle", ancho] )
            #
            #puntoInicial = addNode(nodos, puntos[-1], "esquina")
            #puntoFinal   = p
            #if ( puntoInicial != puntoFinal ):
            #    lineas.append( [puntoInicial, puntoFinal, "calle", ancho] )
        else:
            # No atraviesa --> crear una calle comun
            puntoInicial = addNode(nodos, puntos[0], "esquina")
            puntoFinal   = addNode(nodos, puntos[-1], "esquina")
            if ( puntoInicial != puntoFinal ):
                lineas.append( [puntoInicial, puntoFinal, "calle", ancho] )
                puntoInicial = puntoFinal


    # Crear sumideros
    print "Creando sumideros"
    tF.write("Creando sumideros\n")
    sumideros = []
    for (i, nodo) in enumerate(nodos):
        if nodo[2] == "esquina":
            # Buscar el nodo conducto mas cercano
            mindist, minj = 80, -1
            for (j, nodo2) in enumerate(nodos):
                if nodo2[2] == "conducto":
                    d = dist(nodo, nodo2)
                    if d < mindist:
                        mindist = d
                        minj = j

                if nodo2[2] == "esquina":
                    break

            if minj != -1:
                # Existe un nodo conducto cerca (< 80m)
                sumideros.append( [i, minj] )


    # Crear vertederos
    print "Creando vertederos"
    tF.write("Creando vertederos\n")
    vertederos = []
    for (i, nodo) in enumerate(nodos):
        if nodo[2] == "esquina":
            # Buscar el nodo arroyo mas cercano
            mindist, minj = 50, -1
            for (j, nodo2) in enumerate(nodos):
                if nodo2[2] == "arroyo":
                    d = dist(nodo, nodo2)
                    if d < mindist:
                        mindist = d
                        minj = j

                if nodo2[2] == "esquina":
                    break

            if minj != -1:
                # Existe un nodo conducto cerca (< 80m)
                vertederos.append( [i, minj] )

    #Crear centros de cuencas
    centros = []
    for (i, nodo) in enumerate(nodos):
        if nodo[2] != "conducto":
            centros.append([nodo[0], nodo[1], i])

    print "Numero de nodos:   ", len(nodos)
    print "Numero de cuencas: ", len(centros)
    print "Numero de lineas:  ", len(lineas)
    print "Numero de sumideros:  ", len(sumideros)
    print "Numero de vertederos:  ", len(vertederos)

    # Escribir shape con la posicion de los nodos
    escribir_shp_puntos(shpFileNodos, nodos, {}, spatial_ref)
    # Escribir shape con la posicion de los baricentros de subcuencas
    escribir_shp_puntos(shpFileCentros, centros, {}, spatial_ref)

    # Escribir shape con las conducciones
    polilineas = [[nodos[linea[0]], nodos[linea[1]]] for linea in lineas]
    campos = OrderedDict()
    campos["nodo0"] = [int(linea[0]) for linea in lineas]
    campos["nodo1"] = [int(linea[1]) for linea in lineas]
    campos["tipo"]  = [str(linea[2]) for linea in lineas]
    campos["ancho"] = [float(linea[3]) for linea in lineas]
    escribir_shp_polilineas(shpFileLineas, polilineas, campos, spatial_ref)

    # Escribir shape de sumideros
    polilineas = [[nodos[linea[0]], nodos[linea[1]]] for linea in sumideros]
    escribir_shp_polilineas(shpFileSumideros, polilineas, {}, spatial_ref)

    # Escribir shape de vertederos
    polilineas = [[nodos[linea[0]], nodos[linea[1]]] for linea in vertederos]
    escribir_shp_polilineas(shpFileVertederos, polilineas, {}, spatial_ref)

    # Write list files
    saveOnFile(nodos, "nodos")
    saveOnFile(centros, "centros")
    saveOnFile(lineas, "lineas")
    saveOnFile(sumideros, "sumideros")
    saveOnFile(vertederos, "vertederos")

    print "Finalizado proceso de creado de calles"
    tF.write("Finalizado proceso de creado de calles\n")
    tF.close()


def mainGetSubcatchments(shpFileCuenca):
    print "Proceso de creado de cuencas"

    create_thiessen_polygons(shpFileCentros, subcuencasShpFile)

    clip_feature(subcuencasShpFile, shpFileCuenca, subcuencasClipShpFile)

    print "Finalizado proceso de creado de cuencas"

    subcuencas = []
    if engine == "arcgis10":
        calculate_areas(subcuencasClipShpFile, subcuencasAreaShpFile)
        subcuencas = leer_shp_poligonos(subcuencasAreaShpFile, ["Input_FID", "F_AREA"])


    elif engine == "qgis":
        areas = read_areas(subcuencasClipShpFile)
        subcuencas = leer_shp_poligonos(subcuencasClipShpFile, ["FID"])

        for i, subcuenca in enumerate(subcuencas):
            subcuenca.append(i)
            subcuenca.append(areas[i])

    subcuencas.sort(key=lambda x: float(x[1]))

    # Completar si falta alguna y eliminar duplicadas si existieran
    centros = readFromFile('centros')
    i = 0
    while (i<len(centros)):
        fid = subcuencas[i][1]
        if (fid > i):
            subcuencas.insert(i,[[], i, 0])
            i += 1
        elif (fid < i):
            subcuencas.remove(i)
        else:
            i += 1

    saveOnFile(subcuencas, "subcuencas")

    print "Finalizado proceso de creado de cuencas"


def mainSampleNodeData(rasterFileDEM, rasterFileSlope, rasterFileCoeficiente, rasterFileImpermeabilidad):
    print "Proceso de muestreo de datos..."

    nodosElev = sample_raster_on_nodes(shpFileNodos, rasterFileDEM)
    nodosSlope = sample_raster_on_nodes(shpFileNodos, rasterFileSlope)
    nodosCoeficiente = sample_raster_on_nodes(shpFileNodos, rasterFileCoeficiente)
    nodosImpermeabilidad = sample_raster_on_nodes(shpFileNodos, rasterFileImpermeabilidad)

    # Bajar nodos arroyo
    #nodos = readFromFile('nodos')
    #for i in xrange(0,len(nodosElev)):
    #    if (nodos[i][2] == "arroyo"):
    #        nodosElev[i] -= 1.5

    saveOnFile(nodosElev, "nodosElev")
    saveOnFile(nodosSlope, "nodosSlope")
    saveOnFile(nodosCoeficiente, "nodosCoeficiente")
    saveOnFile(nodosImpermeabilidad, "nodosImpermeabilidad")
    print "Finalizado proceso de generacion de muestreo"


def mainCreateOutfallNodes(shpFileNodosBorde):
    print "Proceso de generacion de nodos outfall..."

    posicionesNodosOutfall = leer_shp_puntos(shpFileNodosBorde)

    nodos = readFromFile('nodos')
    nodosElev = readFromFile('nodosElev')

    nodosOutfall = []
    nodosOutfallElev = []
    lineasOutfall = []
    print posicionesNodosOutfall

    for (i, nodo) in enumerate(nodos):
        # Buscar la posicion outfall mas cercana
        mindist, minj = 50, -1
        for (j, posNO) in enumerate(posicionesNodosOutfall):
            if dist(nodo, posNO) < mindist:
                mindisdt = dist(nodo, posNO)
                minj = j
        if minj != -1:
            nodosOutfall.append(list(posicionesNodosOutfall[minj]))
            nodosOutfallElev.append(nodosElev[i])
            lineasOutfall.append( [i, len(nodosOutfall)-1, 50] )

    print nodosOutfall, nodosOutfallElev, lineasOutfall

    saveOnFile(nodosOutfall, "nodosOutfall")
    saveOnFile(nodosOutfallElev, "nodosOutfallElev")
    saveOnFile(lineasOutfall, "lineasOutfall")
    print "Finalizado proceso de generacion de nodos outfall"


def mainCalculateInvertOffsets():
    print "Proceso de calculo de inverts..."

    coTapada = 0.3

    nodos = readFromFile('nodos')
    lineas = readFromFile('lineas')

    nodosInvElevOffset = [0] * len(nodos)
    for (i, linea) in enumerate(lineas):
        in0, in1, tipo, ancho, alto = linea[0], linea[1], linea[2], linea[3], (0 if len(linea) < 5 else linea[4])
        offset = 0 if tipo == "calle" else -alto - coTapada
        if (offset < nodosInvElevOffset[in0]):
            nodosInvElevOffset[in0] = offset;
        if (offset < nodosInvElevOffset[in1]):
            nodosInvElevOffset[in1] = offset;
    saveOnFile(nodosInvElevOffset, "nodosInvElevOffset")

    print "Finalizado el calculo de inverts"


def mainCorregirArroyos():
    nodos = readFromFile('nodos')
    nodosElev = readFromFile('nodosElev')
    nodosInvElevOffset = readFromFile('nodosInvElevOffset')

    # Leer archivo imp anterior
    tF = open(swmmCorregirArroyosFileName, "r")
    lineas = tF.readlines()
    tF.close()
    # Leer nodos e invElevsAnteriores
    invElevsAnteriores = dict()
    for i in xrange(0,len(lineas)):
        if lineas[i].startswith("[STORAGE]") :
            i = i + 1
            while (lineas[i].startswith("[") == False):
                if (lineas[i].startswith(";") == False):
                    partes = lineas[i].split()
                    if len(partes) >= 2:
                        invElevsAnteriores[partes[0]] = float(partes[1])
                i += 1
            break
    # Leer coordenadas anteriores
    coordenadasAnteriores = []
    for i in xrange(0,len(lineas)):
        if lineas[i].startswith("[COORDINATES]") :
            i = i + 1
            while (lineas[i].startswith("[") == False):
                if (lineas[i].startswith(";") == False):
                    partes = lineas[i].split()
                    if len(partes) == 3:
                        coordenadasAnteriores.append( [float(partes[1]), float(partes[2]), partes[0]] )
                i += 1
            break

    # Bajar los inverts actuales de los nodos arroyo y conducto si hace falta
    for i, nodo in enumerate(nodos):
        if nodo[2] != "esquina":
            # Buscar nodo cercano
            for coordenada in coordenadasAnteriores:
                #Buscar el mismo nodo en el archivo anterior
                if (distSq(nodo, coordenada) < 1.0):
                    invElevAnterior = invElevsAnteriores[coordenada[2]]
                    #Bajar el terreno si es necesario
                    if (invElevAnterior < nodosElev[i] + nodosInvElevOffset[i] - 0.02 ) :
                        print "Nodo ", i, " bajado", nodosElev[i] + nodosInvElevOffset[i] - invElevAnterior
                        nodosInvElevOffset[i] = invElevAnterior - nodosElev[i]
                    break

    # Grabar invertOffsets modificados
    saveOnFile(nodosInvElevOffset, "nodosInvElevOffset")


def mainCreateSWMM(swmmInputFileName):
    print "Proceso de escritura de archivo SWMM..."

    evap = 4 # mm/dia
    # percent imperviousness of subcatchment.
    #scImperv = 25 #%
    # subcatchment slope (percent).
    #scSlope = 0.2 #%

    # Manning's n for overland flow over the impervious sub-area.
    saNImp = 0.025
    # Manning's n for overland flow over the pervious sub-area.
    saNPerv = 0.05
    # depression storage for impervious sub-area (inches or mm).
    saSImp = 0 #mm
    # depression storage for pervious sub-area (inches or mm).
    saSPerv = 0 #mm
    # percent of impervious area with no depression storage.
    saZero = 100 #%

    # Infiltration Horton
    inF0 = 50 #mm/hr
    inFf = 5 #mm/hr
    inCoefDecaim = 2
    # Time it takes for fully saturated soil to dry  (days).
    inDryTime = 5
    # Maximum infiltration volume possible (0 if not applicable) (in or mm)
    inMaxInf = 0

    # Depth from ground to invert elevation (ft or m) (default is 0).
    juYmax = 1.0
    # Water depth at start of simulation (ft or m) (default is 0).
    juY0 = 0
    # maximum additional head above ground elevation that manhole junction can sustain under surcharge conditions (ft or m) (default is 0).
    juYsur = 0
    # % del area de la subcuenca en que se almacena el agua en el nodo
    juApondPer = 0.75

    # value of n (i.e., roughness parameter) in Manning's equation.
    coN = 0.04
    coTapada = 0.3

    # Weirs parameters
    weAlturaCordon = 0.05
    weCd = 3.0
    # Orifices parameters
    orCd = 0.6

    # XSection parameters
    xsG1, xsG3, xsG4 = 10, 0, 0
    xsSumideroH, xsSumideroW = 0.15, 2.0
    xsVertederoH, xsVertederoW = 10.0, 20.0

    # Transect parameters
    traNConducto = 0.03
    traNArroyoPlanicie = 0.05
    traNArroyoCauce = 0.03
    traNCalle = 0.04
    traAnchoMargenArroyo = 20
    #traTiranteArroyo = 1.5
    traAnchoVereda = 2
    traAltoCordon = 0.2

    nodos = readFromFile('nodos')
    centros = readFromFile('centros')
    lineas = readFromFile('lineas')
    sumideros = readFromFile('sumideros')
    vertederos = readFromFile('vertederos')
    subcuencas = readFromFile('subcuencas')
    nodosOutfall = readFromFile('nodosOutfall')
    nodosOutfallElev = readFromFile('nodosOutfallElev')
    lineasOutfall = readFromFile('lineasOutfall')

    nodosInvElevOffset = readFromFile('nodosInvElevOffset')

    nodosElev = readFromFile('nodosElev')
    nodosSlope = readFromFile('nodosSlope')
    nodosCoeficiente = readFromFile('nodosCoeficiente')
    nodosImpermeabilidad = readFromFile('nodosImpermeabilidad')

    nodosLongitudLineas = [0] * len(nodos)
    for (i, nodo) in enumerate(nodos):
        for linea in lineas:
            if linea[0] == i or linea[1] == i:
                nodosLongitudLineas[i] = nodosLongitudLineas[i] + dist(nodos[linea[0]], nodos[linea[1]])


    areasNodo = [1.167] * len(nodos)
    for (i, centro) in enumerate(centros):
        areasNodo[centro[2]] = juApondPer * subcuencas[i][2]

    tF = open(swmmInputFileName, "w")
    tF.write("[TITLE]\n")
    tF.write("Conurbano\n")

    tF.write("\n")
    tF.write("[OPTIONS]\n")
    tF.write("FLOW_UNITS         CMS\n")
    tF.write("INFILTRATION       HORTON\n") #CURVE_NUMBER
    tF.write("FLOW_ROUTING       DYNWAVE\n")
    tF.write("LINK_OFFSETS       ELEVATION\n")
    tF.write("START_DATE         7/8/2013\n")
    tF.write("START_TIME         00:00\n")
    tF.write("END_TIME           06:00\n")
    tF.write("WET_STEP           00:00:30\n")
    tF.write("DRY_STEP           00:01:00\n")
    tF.write("ROUTING_STEP       00:00:30\n")
    # Minimum Surface Area - This is a minimum surface area used at nodes when computing changes in water depth. If 0 is entered, then the default value of 12.566 ft2 (1.167 m2) is used. This is the area of a 4-ft diameter manhole. The value entered should be in square feet for US units or square meters for SI units.
    #tF.write("MIN_SURFAREA       75.54\n") #m2, equivalente a 10 m de diametro
    tF.write("MIN_SURFAREA       1.167\n") #m2, equivalente a 10 m de diametro
    tF.write("\n")

    tF.write("[FILES]\n")
    tF.write("SAVE RAINFALL      rainfall.rff\n")
    tF.write("SAVE RUNOFF        runoff.rof\n")
    tF.write("SAVE OUTFLOWS      outflows.txt\n")


    tF.write("\n")
    tF.write("[RAINGAGES]\n")
    tF.write(";;Name         Format         Interval       SCF            DataSrc        SourceName     Sta            Units\n")
    tF.write(";;============================================================================================================\n")
    for i in xrange(0,numGages):
        gageName = 'GAGE'+str(i * (100/(numGages-1)))
        list = [gageName, 'INTENSITY', '0:05', 1.0, 'FILE', "pluviom.dat", gageName, "MM"]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")

    tF.write("\n")
    tF.write("[EVAPORATION]\n")
    tF.write("CONSTANT  "+str(evap)+"\n")

    tF.write("\n")
    tF.write("[SUBCATCHMENTS]\n")
    tF.write(";;Name         Raingage       Outlet         Area           %ImperV        Width          Slope          Curve Length\n")
    tF.write(";;======================================================================================================================\n")
    for (i, centro) in enumerate(centros):
        numNodo = centro[2]
        coef = int(nodosCoeficiente[numNodo])
        gageName = 'GAGE' + str(coef - (coef%(100/(numGages-1))))
        list = ['CUENCA'+str(i), gageName, 'NODO'+str(numNodo), "%.3f" % (float(subcuencas[i][2])/10000.0), "%.3f" % nodosImpermeabilidad[numNodo], "%.3f" % (nodosLongitudLineas[numNodo]/2), "%.3f" % (nodosSlope[numNodo]), "%.3f" % (subcuencas[i][2]**0.5)]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")

    tF.write("\n")
    tF.write("[SUBAREAS]\n")
    tF.write(";;Name         N_Imp          N_Perv         S_Imp          S_Perv         %ZER           RouteTo\n")
    tF.write(";;=======================================================================================================\n")
    for (i, centro) in enumerate(centros):
        list = ['CUENCA'+str(i), saNImp, saNPerv, saSImp, saSPerv, saZero, 'OUTLET']
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")

    tF.write("\n")
    tF.write("[INFILTRATION]\n")
    tF.write(";;Subcat       MaxRate        MinRate        Decay          DryTime        Max Inf\n")
    tF.write(";;========================================================================================\n")
    for (i, centro) in enumerate(centros):
        list = ['CUENCA'+str(i), inF0, inFf, inCoefDecaim, inDryTime, inMaxInf]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")


    # tF.write("\n")
    # tF.write("[JUNCTIONS]\n")
    # tF.write(";;Name         Elev           Ymax           Y0             Ysur           Apond) \n")
    # tF.write(";;========================================================================================\n")
    # areasNodo = [75.54] * len(nodos)
    # for (i, centro) in enumerate(centros):
        # areasNodo[centro[2]] = juApondPer * subcuencas[i][2]
    # for (i, nodo) in enumerate(nodos):
        # list = ['NODO'+str(i), nodosElev[i]+nodosInvElevOffset[i], juYmax-nodosInvElevOffset[i], juY0, juYsur, "%.3f" % areasNodo[i]]
        # tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        # tF.write("\n")


    tF.write("\n")
    tF.write("[STORAGE]\n")
    tF.write(";;Name         Elev           Ymax           Y0             TABULAR        Apond          ) \n")
    tF.write(";;========================================================================================\n")
    for (i, nodo) in enumerate(nodos):
        list = ['NODO'+str(i), "%.3f" % (nodosElev[i]+nodosInvElevOffset[i]), "%.3f" % (-nodosInvElevOffset[i]+20), juY0, 'TABULAR', 'STORAGE'+str(i), "%.3f" % (areasNodo[i])]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")


    tF.write("\n")
    tF.write("[CURVES]\n")
    tF.write(";;Name         Type           x-value        y-value\n")
    tF.write(";;========================================================================================\n")
    for (i, nodo) in enumerate(nodos):
        if (nodosInvElevOffset[i] < 0):
            list = ['STORAGE'+str(i), 'STORAGE', 0, 1.167, -nodosInvElevOffset[i], 1.167, -nodosInvElevOffset[i]+1, "%.3f" % (areasNodo[i]/10),-nodosInvElevOffset[i]+2, "%.3f" % (areasNodo[i]),-nodosInvElevOffset[i]+20, "%.3f" % (areasNodo[i])]
        else:
            list = ['STORAGE'+str(i), 'STORAGE', 0, 1.167, 1, "%.3f" % (areasNodo[i]/10), 2, "%.3f" % (areasNodo[i]), 20, "%.3f" % (areasNodo[i])]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")


    tF.write("\n")
    tF.write("[OUTFALLS]\n")
    tF.write(";;Name         Elev           Type           Gate\n")
    tF.write(";;==========================================================\n")
    for (i, nodo) in enumerate(nodosOutfall):
        list = ['NODOOUT'+str(i), nodosOutfallElev[i], "FREE", "NO"]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")

    tF.write("\n")
    tF.write("[CONDUITS]\n")
    tF.write(";;Name         Node1          Node2          Length         N              Z1             Z2             Q0\n")
    tF.write(";;======================================================================================================================\n")
    for (i, linea) in enumerate(lineas):
        in0, in1, tipo, ancho, alto = linea[0], linea[1], linea[2], linea[3], (0 if len(linea) < 5 else linea[4])
        offset = 0 if tipo == "calle" else -alto - coTapada
        name = tipo + str(i)
        if (tipo == "calle"):
            list = [name, 'NODO'+str(in0), 'NODO'+str(in1), "%.3f" % dist(nodos[in0], nodos[in1]), "%.3f" % coN, "%.3f" % (nodosElev[in0]+offset), "%.3f" % (nodosElev[in1]+offset), 0]
        else:
            list = [name, 'NODO'+str(in0), 'NODO'+str(in1), "%.3f" % dist(nodos[in0], nodos[in1]), "%.3f" % coN, "*", "*", 0]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")
    for (i, linea) in enumerate(lineasOutfall):
        in0, in1, ancho = linea
        list = ['SALIDA'+str(i), 'NODO'+str(in0), 'NODOOUT'+str(in1), "%.3f" % dist(nodos[in0], nodosOutfall[in1]), "%.3f" % coN, "%.3f" % (nodosElev[in0]+nodosInvElevOffset[in0]), "%.3f" % (nodosElev[in0]+nodosInvElevOffset[in0]), 0]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")
    for (i, vertedero) in enumerate(vertederos):
        in0, in1 = vertedero[0], vertedero[1]
        list = ["vertedero" + str(i), 'NODO'+str(in0), 'NODO'+str(in1), 20, "%.3f" % coN, "%.3f" % (nodosElev[in0]+weAlturaCordon), "*", 0]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")

    # tF.write("\n")
    # tF.write("[WEIRS]\n")
    # tF.write(";;Name         Node1          Node2          Type           Offset         Cd             Flap  (EC Cd2)          \n")
    # tF.write(";;======================================================================================================================\n")
    # for (i, vertedero) in enumerate(vertederos):
        # in0, in1 = vertedero[0], vertedero[1]

        # name = "vertedero" + str(i)
        # list = [name, 'NODO'+str(in0), 'NODO'+str(in1), "SIDEFLOW", "%.3f" % (nodosElev[in0]+weAlturaCordon), weCd, "NO"]
        # tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        # tF.write("\n")


    tF.write("\n")
    tF.write("[ORIFICES]\n")
    tF.write(";;Name         Node1          Node2          Type           Offset         Cd             Flap           Orate\n")
    tF.write(";;======================================================================================================================\n")
    for (i, sumidero) in enumerate(sumideros):
        in0, in1 = sumidero[0], sumidero[1]

        name = "sumidero" + str(i)
        list = [name, 'NODO'+str(in0), 'NODO'+str(in1), "SIDE", 0, "%.3f" % (nodosElev[in0]-traAltoCordon), "NO", 0]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")


    tF.write("\n")
    tF.write("[XSECTIONS]\n")
    tF.write(";;Link         Type           G1             G2             G3             G4\n")
    tF.write(";;========================================================================================\n")
    transectas = {}
    for (i, linea) in enumerate(lineas):
        in0, in1, tipo, ancho, alto = linea[0], linea[1], linea[2], linea[3], (0 if len(linea) < 5 else linea[4])
        name = tipo + str(i)
        if (tipo == "conducto"):
            list = [name, 'RECT_CLOSED', alto, ancho, 0, 0]
        else:
            tname = tipo + str(int(ancho)) + "x" + str(int(alto))
            transectas[tname] = [tipo, tname, ancho, alto]
            list = [name, 'IRREGULAR', tname, 0, 0, 0]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")
    for (i, linea) in enumerate(lineasOutfall):
        in0, in1, ancho = linea
        list = ['SALIDA'+str(i), 'RECT_OPEN', xsG1, ancho, xsG3, xsG4]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")
    for (i, vertedero) in enumerate(vertederos):
        list = ['vertedero'+str(i), 'RECT_OPEN', xsVertederoH, xsVertederoW, 0, 0]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")
    for (i, sumidero) in enumerate(sumideros):
        list = ['sumidero'+str(i), 'RECT_CLOSED', xsSumideroH, xsSumideroW, 0, 0]
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
        if (tipo == "arroyo"):
            traTiranteArroyo = alto + coTapada
            list = ['NC', traNArroyoPlanicie, traNArroyoPlanicie, traNArroyoCauce]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")
            list = ['X1', tname, 6, -ancho * 0.5, ancho * 0.5, 0, 0, 0, 0, 0, 0]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")
            # list = ['GR', -traAnchoMargenArroyo-ancho*0.5, traTiranteArroyo+20, -traAnchoMargenArroyo-ancho*0.5, traTiranteArroyo, -ancho*0.5, traTiranteArroyo, -ancho*0.5+1, 0, ancho*0.5-1, 0, ancho*0.5,traTiranteArroyo, traAnchoMargenArroyo+ancho*0.5, traTiranteArroyo, traAnchoMargenArroyo+ancho*0.5, traTiranteArroyo+20]
            list = ['GR', traTiranteArroyo+20, -traAnchoMargenArroyo-ancho*0.5, traTiranteArroyo, -traAnchoMargenArroyo-ancho*0.5, traTiranteArroyo, -ancho*0.5, 0, -ancho*0.5+0.25, 0, ancho*0.5-0.25, traTiranteArroyo, ancho*0.5, traTiranteArroyo, traAnchoMargenArroyo+ancho*0.5, traTiranteArroyo+20, traAnchoMargenArroyo+ancho*0.5]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")
            tF.write(";;-------------------------------------------\n")
        # elif (tipo == "conducto"):
            # list = ['NC', traNConducto, traNConducto, traNConducto]
            # tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            # tF.write("\n")
            # list = ['X1', tname, 6, -ancho * 0.5, ancho * 0.5, 0, 0, 0, 0, 0, 0]
            # tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            # tF.write("\n")
            # list = ['GR', -ancho*0.5, alto, -ancho*0.5, 0, ancho*0.5, 0, ancho*0.5, alto]
            # tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            # tF.write("\n")
            # tF.write(";;-------------------------------------------\n")
        elif (tipo == "calle"):
            list = ['NC', traNCalle, traNCalle, traNCalle]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")
            list = ['X1', tname, 7, -ancho * 0.5, ancho * 0.5, 0, 0, 0, 0, 0, 0]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")
            # list = ['GR', -traAnchoVereda-ancho*0.5, traAltoCordon+20, -traAnchoVereda-ancho*0.5, traAltoCordon, -ancho*0.5, traAltoCordon, -ancho*0.5, 0, 0, traAltoCordon, ancho*0.5, 0, ancho*0.5, traAltoCordon, traAnchoVereda+ancho*0.5, traAltoCordon, traAnchoVereda+ancho*0.5, traAltoCordon+20]
            list = ['GR', traAltoCordon+20, -traAnchoVereda-ancho*0.5, traAltoCordon, -traAnchoVereda-ancho*0.5, traAltoCordon, -ancho*0.5, 0, -ancho*0.5, traAltoCordon, 0, 0, ancho*0.5, traAltoCordon, ancho*0.5, traAltoCordon, traAnchoVereda+ancho*0.5, traAltoCordon+20, traAnchoVereda+ancho*0.5]
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
    minx, miny, maxx, maxy = 1e9, 1e9, -1e9,-1e9
    for n in nodos:
        if (n[0] < minx):
            minx = n[0]
        if (n[0] > maxx):
            maxx = n[0]
        if (n[1] < miny):
            miny = n[1]
        if (n[1] > maxy):
            maxy = n[1]
    distx = maxx - minx;
    disty = maxy - miny;
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
        list = ['NODO'+str(i), nodo[0], nodo[1]]
        tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
        tF.write("\n")
    for (i, nodo) in enumerate(nodosOutfall):
        list = ['NODOOUT'+str(i), nodo[0], nodo[1]]
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

    tF.close()

    print "Finalizado el proceso de escritura de archivo SWMM"


def mainReadSWMMResultsDepths(swmmOuputFileName):
    nodos = readFromFile('nodos')
    nodosElev = readFromFile('nodosElev')
    nodosInvElevOffset = readFromFile('nodosInvElevOffset')

    print "Numero de nodos: ", len(nodos)

    outfile = swmmout.open(swmmOuputFileName)

    query_vars = ['depth']
    query_nodes = []
    for i in xrange(0,len(nodos)):
        query_nodes.append( 'NODO'+str(i) )
    data = outfile.get_values('nodes', query_nodes, query_vars)

    nodeMaxDepths = []
    for i in xrange(0,len(nodos)):
        nodeMaxDepths.append(max([d[i+1][1] for d  in data]))
        nodeMaxDepths[i] = nodeMaxDepths[i] + nodosInvElevOffset[i]
        if (nodeMaxDepths[i] < 0):
            nodeMaxDepths[i] = 0

    # Conseguir la referencia geografica
    spatial_ref = leer_spatial_reference(shpFileNodos)

    # Escribir shape con las profundidades maximas
    campos = OrderedDict()
    campos["depth"] = [maxDepth for maxDepth in nodeMaxDepths]
    campos["elev"] = [maxDepth + nodosElev[j] for j, maxDepth in enumerate(nodeMaxDepths)]
    escribir_shp_puntos(workspace + "/" + "nodeDepthMax.shp", nodos, campos, spatial_ref)

    # Escribir shape con las profundidades en cada paso de tiempo
    for i, dataline in enumerate(data):
        campos = OrderedDict()
        campos["depth"] = [max(dataline[j][1] + nodosInvElevOffset[j-1], 0) for j in range(1,len(dataline))]
        campos["elev"] = [max(dataline[j][1] + nodosInvElevOffset[j-1], 0) + nodosElev[j-1] for j in range(1,len(dataline))]
        escribir_shp_puntos("nodeDepth%04d.shp" % i, nodos, campos, spatial_ref)


def mainCreateGenerateRain(gageFileName, gagesFileName):
    print "Leyendo pluviometro..."
    iFile = open(gageFileName, "r")
    lineas = iFile.readlines()
    iFile.close()

    print "Escribiendo estaciones..."
    tF = open(gagesFileName, "w")
    for i in xrange(0,numGages):
        coef = i * (100/(numGages-1))
        gageName = 'GAGE'+str(coef)
        for linea in lineas:
            datos = linea.split()
            if len(datos) == 0:
                continue
            datos[0]  = gageName
            datos[-1] = float(datos[-1])*(float(coef)/100.0)
            tF.write(("").join([ str(x).ljust(15, ' ') for x in datos]))
            tF.write("\n")

    print "Estaciones listas"


def mainCalculateDeadDepths():
    nodos = readFromFile('nodos')
    nodosElev = readFromFile('nodosElev')
    print len(nodos)
    print len(nodosElev)
    nodosInvElevOffset = readFromFile('nodosInvElevOffset')
    lineas = readFromFile('lineas')
    sumideros = readFromFile('sumideros')
    vertederos = readFromFile('vertederos')
    lineasOutfall = readFromFile('lineasOutfall')


    tirantes = [100 for nodo in nodos]
    print lineasOutfall
    i = 0
    maxbajada = 100
    while maxbajada > 0.01:
        maxbajada = 0
        for linea in lineasOutfall:
            nodo0 = linea[0]
            tirantes[nodo0] = 0

        def igualar(nodo0, nodo1):
            elev0 = nodosElev[nodo0] + tirantes[nodo0]
            elev1 = nodosElev[nodo1] + tirantes[nodo1]
            if elev0 <= elev1:
                elev11 = max(nodosElev[nodo1], elev0)
                tirantes[nodo1] = elev11 - nodosElev[nodo1]
                bajada = elev1 - elev11
            else:
                elev01 = max(nodosElev[nodo0], elev1)
                tirantes[nodo0] = elev01 - nodosElev[nodo0]
                bajada = elev0 - elev01

            return bajada

        for linea in lineas:
            nodo0, nodo1 = linea[0], linea[1]
            bajada = igualar(nodo0, nodo1)
            maxbajada = max(maxbajada, bajada)

        for linea in sumideros:
            nodo0, nodo1 = linea[0], linea[1]
            bajada = igualar(nodo0, nodo1)
            maxbajada = max(maxbajada, bajada)

        for linea in vertederos:
            nodo0, nodo1 = linea[0], linea[1]
            bajada = igualar(nodo0, nodo1)
            maxbajada = max(maxbajada, bajada)

        i += 1
        print "Iteracion %i - Max Bajada %f" % (i, maxbajada)

    # Conseguir la referencia geografica
    spatial_ref = leer_spatial_reference(shpFileNodos)

    # Escribir shape con las profundidades maximas
    campos = {
                "depth" : tirantes,
                "nelev" : nodosElev
             }
    escribir_shp_puntos("nodesProfMuerta.shp", nodos, campos, spatial_ref)





if __name__ == '__console__' :
    os.chdir(workspace)
    mainReadRivers()

if __name__ == '__main__':
    os.chdir(workspace)
    gis_init()

    # Opciones de corrida
    print "Que desea hacer?"
    print " 0 - Leer arroyos"
    print " 1 - Leer las calles"
    print " 2 - Crear y leer subcuencas a partir de nodos"
    print " 3 - Leer los rasters en cada nodo"
    print " 4 - Leer Nodos de borde y generar outfalls"
    print " 5 - Calcular los inverts"
    print " 6 - (Opcional)Corregir los perfiles de arroyos y conductos usando otro archivo .inp"
    print " 7 - Generar archivos de SWMM"
    print " 8 - Todo"
    print " 9 - Leer resultados y escribir shp con profundidad y elevacion en nodos"
    print "11 - Crear pluviometros"
    print "12 - Analisis de tirantes muertos"
    x = input("Opcion:")
    if (x == 0):
        mainReadRivers(defaultShpFileArroyos)
    elif (x == 1):
        mainReadStreets(defaultShpFileCalles)
    elif (x == 2):
        mainGetSubcatchments(defaultShpFileCuenca)
    elif (x == 3):
        mainSampleNodeData(defaultRasterFileDEM, defaultRasterFileSlope, defaultRasterFileCoeficiente, defaultRasterFileImpermeabilidad)
    elif (x == 4):
        mainCreateOutfallNodes(defaultShpFileNodosBorde)
    elif (x == 5):
        mainCalculateInvertOffsets()
    elif (x == 6):
        mainCorregirArroyos()
    elif (x == 7):
        mainCreateSWMM(defaultSwmmInputFileName)
    elif (x == 8):
        mainReadRivers(defaultShpFileArroyos)
        mainReadStreets(defaultShpFileCalles)
        mainGetSubcatchments(defaultShpFileCuenca)
        mainSampleNodeData(defaultRasterFileDEM, defaultRasterFileSlope, defaultRasterFileCoeficiente, defaultRasterFileImpermeabilidad)
        mainCreateOutfallNodes(defaultShpFileNodosBorde)
        mainCalculateInvertOffsets()
        mainCreateSWMM(defaultSwmmInputFileName)
    elif (x == 9):
        mainReadSWMMResultsDepths(defaultSwmmOuputFileName)
    elif (x == 11):
        mainCreateGenerateRain(defaultGageFileName, defaultGagesFileName)
    elif (x == 12):
        mainCalculateDeadDepths()

