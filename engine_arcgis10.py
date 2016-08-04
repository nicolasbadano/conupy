# GIS functions for ArcGis 10.1

import arcpy
from arcpy import env
from arcpy.sa import *
import os
import numpy as np

default_spatial_reference = None

def gis_init():
    return


def crear_datos_temporales():
    return arcpy.CreateScratchName(workspace=arcpy.env.scratchGDB)


def leer_spatial_reference(fileName):
    desc = arcpy.Describe(fileName)
    return desc.SpatialReference


def leer_shp_puntos(shp_file, lista_campos = []):
    resultados = []
    print "Leyendo shape de puntos ", shp_file
    try:
        sRow, sCur, sFeat = None, None, None
        shapeName = arcpy.Describe( shp_file ).ShapeFieldName
        sCur = arcpy.SearchCursor( shp_file )
        sRow = sCur.next()
        pnt = arcpy.CreateObject("Point")

        while (sRow):
            sFeat = sRow.getValue(shapeName)

            pnt = sFeat.getPart(0)
            punto = np.array([pnt.X, pnt.Y])
            for campo in lista_campos:
                punto.append( sRow.getValue(campo) )

            resultados.append(punto)
            sRow = sCur.next()

    except Exception, err:
        print err.message
        arcpy.AddError(err.message)

    finally:
        if sRow:
            del sRow
        if sCur:
            del sCur
        if sFeat:
            del sFeat
        print "Finalizado de leer shape de puntos."
        return resultados


def leer_shp_polilineas(shp_file, lista_campos = []):
    resultados = []
    print "Leyendo shape de polilineas ", shp_file
    try:
        sRow, sCur, sFeat = None, None, None
        shapeName = arcpy.Describe( shp_file ).ShapeFieldName
        sCur = arcpy.SearchCursor( shp_file )
        sRow = sCur.next()
        pnt = arcpy.CreateObject("Point")

        while (sRow):
            sFeat = sRow.getValue(shapeName)

            partcount = sFeat.partCount
            if (partcount == 1):
                part = sFeat.getPart(0)

                poly = []
                puntos = []
                pnt = part.next()
                while pnt:
                    punto = np.array([pnt.X, pnt.Y])
                    puntos.append(punto)
                    pnt = part.next()

                poly.append(puntos)
                for campo in lista_campos:
                    poly.append( sRow.getValue(campo) )

                resultados.append(poly)

            sRow = sCur.next()

    except Exception, err:
        print err.message
        arcpy.AddError(err.message)

    finally:
        if sRow:
            del sRow
        if sCur:
            del sCur
        if sFeat:
            del sFeat
        print "Finalizado de leer shape de polilineas."
        return resultados


def leer_shp_poligonos(shp_file, lista_campos = []):
    resultados = []
    print "Leyendo shape de poligonos ", shp_file
    try:
        sRow, sCur, sFeat = None, None, None
        shapeName = arcpy.Describe( shp_file ).ShapeFieldName
        sCur = arcpy.SearchCursor( shp_file )
        sRow = sCur.next()
        pnt = arcpy.CreateObject("Point")

        while (sRow):
            sFeat = sRow.getValue(shapeName)

            partcount = sFeat.partCount
            if (partcount == 1):
                part = sFeat.getPart(0)

                poly = []
                puntos = []
                pnt = part.next()
                while pnt:
                    punto = np.array([pnt.X, pnt.Y])
                    puntos.append(punto)
                    pnt = part.next()

                poly.append(puntos)
                for campo in lista_campos:
                    poly.append( sRow.getValue(campo) )

                resultados.append(poly)

            sRow = sCur.next()

    except Exception, err:
        print err.message
        arcpy.AddError(err.message)

    finally:
        if sRow:
            del sRow
        if sCur:
            del sCur
        if sFeat:
            del sFeat
        print "Finalizado de leer shape de polígonos."
        return resultados


def escribir_shp_puntos(fileName, puntos, campos, spatial_reference):
    print "Escribiendo shape de puntos: " + fileName
    try:
        # Assign empty values to cursor and row objects
        iCur, iFeat = None, None

        outPath, outFC = os.path.split(fileName)
        arcpy.env.workspace = outPath
        arcpy.env.OverWriteOutput = True

        if not arcpy.Exists(fileName):
            arcpy.CreateFeatureclass_management(outPath, outFC, "Point", None, "", "", spatial_reference)

            for nombre in campos:
                dato = campos[nombre][0]
                type = ""
                if isinstance(dato, int):
                    type = "SHORT"
                    arcpy.AddField_management(fileName, nombre, type, 9, "", "", "", "NON_NULLABLE")
                elif isinstance(dato, float):
                    type = "FLOAT"
                    arcpy.AddField_management(fileName, nombre, type, 9, 6, "", "", "NON_NULLABLE")
                elif isinstance(dato, str):
                    type = "TEXT"
                    arcpy.AddField_management(fileName, nombre, type, 9, "", "", "", "NON_NULLABLE")
        else:
            arcpy.DeleteFeatures_management(fileName)

        iCur = arcpy.InsertCursor(fileName)
        pnt = arcpy.CreateObject("Point")

        for (i, punto) in enumerate(puntos):
            (pnt.X, pnt.Y) = (punto[0], punto[1])
            iFeat = iCur.newRow()
            iFeat.setValue("Shape", pnt)
            for nombre in campos:
                iFeat.setValue(nombre, campos[nombre][i])

            iCur.insertRow(iFeat)

    except Exception, err:
        print err.message
        arcpy.AddError(err.message)

    finally:
        print "Finalizado escribir shape de puntos."
        if iCur:
            del iCur
        if iFeat:
            del iFeat


def escribir_shp_polilineas(fileName, polilineas, campos, spatial_reference):
    print "Escribiendo shape de polilineas: " + fileName
    try:
        # Assign empty values to cursor and row objects
        iCur, iFeat = None, None

        outPath, outFC = os.path.split(fileName)
        arcpy.env.workspace = outPath
        arcpy.env.OverWriteOutput = True

        if not arcpy.Exists(fileName):
            arcpy.CreateFeatureclass_management(outPath, outFC, "Polyline", None, "", "", spatial_reference)

            for nombre in campos:
                dato = campos[nombre][0]
                type = ""
                if isinstance(dato, int):
                    type = "SHORT"
                    arcpy.AddField_management(fileName, nombre, type, 9, "", "", "", "NON_NULLABLE")
                elif isinstance(dato, float):
                    type = "FLOAT"
                    arcpy.AddField_management(fileName, nombre, type, 9, 6, "", "", "NON_NULLABLE")
                elif isinstance(dato, str):
                    type = "TEXT"
                    arcpy.AddField_management(fileName, nombre, type, 9, "", "", "", "NON_NULLABLE")
        else:
            arcpy.DeleteFeatures_management(fileName)

        iCur = arcpy.InsertCursor(fileName)
        pnt = arcpy.Point()
        arrayObj = arcpy.Array()

        for (i, polilinea) in enumerate(polilineas):
            for (j, punto) in enumerate(polilinea):
                (pnt.X, pnt.Y) = (punto[0], punto[1])
                arrayObj.add(pnt)

            iFeat = iCur.newRow()
            iFeat.setValue("Shape", arrayObj)
            for nombre in campos:
                iFeat.setValue(nombre, campos[nombre][i])

            iCur.insertRow(iFeat)
            arrayObj.removeAll()

    except Exception, err:
        print err.message
        arcpy.AddError(err.message)

    finally:
        print "Finalizado escribir shape de polilineas."
        if iCur:
            del iCur
        if iFeat:
            del iFeat


def escribir_shp_poligonos(fileName, poligonos, campos, spatial_reference):
    print "Escribiendo shape de polígonos: " + fileName
    try:
        # Assign empty values to cursor and row objects
        iCur, iFeat = None, None

        outPath, outFC = os.path.split(fileName)
        arcpy.env.workspace = outPath
        arcpy.env.OverWriteOutput = True

        if not arcpy.Exists(fileName):
            arcpy.CreateFeatureclass_management(outPath, outFC, "Polygon", None, "", "", spatial_reference)

            for nombre in campos:
                dato = campos[nombre][0]
                type = ""
                if isinstance(dato, int):
                    type = "SHORT"
                    arcpy.AddField_management(fileName, nombre, type, 9, "", "", "", "NON_NULLABLE")
                elif isinstance(dato, float):
                    type = "FLOAT"
                    arcpy.AddField_management(fileName, nombre, type, 9, 6, "", "", "NON_NULLABLE")
                elif isinstance(dato, str):
                    type = "TEXT"
                    arcpy.AddField_management(fileName, nombre, type, 9, "", "", "", "NON_NULLABLE")
        else:
            arcpy.DeleteFeatures_management(fileName)

        iCur = arcpy.InsertCursor(fileName)
        pnt = arcpy.Point()
        arrayObj = arcpy.Array()

        for (i, poligono) in enumerate(poligonos):
            for (j, punto) in enumerate(poligono[0]):
                (pnt.X, pnt.Y) = (punto[0], punto[1])
                arrayObj.add(pnt)

            iFeat = iCur.newRow()
            iFeat.setValue("Shape", arrayObj)
            for nombre in campos:
                iFeat.setValue(nombre, campos[nombre][i])

            iCur.insertRow(iFeat)
            arrayObj.removeAll()

    except Exception, err:
        print err.message
        arcpy.AddError(err.message)

    finally:
        print "Finalizado escribir shape de polígonos."
        if iCur:
            del iCur
        if iFeat:
            del iFeat


def sample_raster_on_nodes(nodesFile, rasterFile):
    print "Muestreando el raster..."

    inRasters = [rasterFile]
    inLocations = nodesFile
    outTable = arcpy.CreateScratchName(workspace=arcpy.env.scratchGDB)

    arcpy.CheckOutExtension("Spatial")
    Sample(inRasters, inLocations, outTable, "NEAREST")

    field_list = arcpy.ListFields(outTable)
    print [field.name for field in field_list]

    print "Leyendo los valores muestreandos..."
    value_list = []

    sRow, sCur = None, None
    sCur = arcpy.SearchCursor(outTable)

    sRow = sCur.next()
    i=0
    while (sRow):
        valor = [sRow.getValue(field.name) for field in field_list]
        value_list.append( valor )
        sRow = sCur.next()
        i += 1


    # Quedarse solo con el valor de la elevación
    values = [None for v in value_list]
    pairs = {};
    for row in value_list:
        pairs[row[1]] = float(row[-1])

    values = [pairs[i] for i in range(0, len(values))]

    print "Finalizado el muestreado de valores."
    return values


def create_thiessen_polygons(nodesFile, fileName):
    print "Creando polígonos de Thiessen. Destino: " + fileName
    try:
        outPath, outFC = os.path.split(fileName)
        arcpy.env.workspace = outPath
        arcpy.env.OverWriteOutput = True

        arcpy.CreateThiessenPolygons_analysis(nodesFile, outFC, "ONLY_FID")

    except Exception, err:
        print err.message
        arcpy.AddError(err.message)

    finally:
        print "Finalizada creación de polígonos."


def clip_feature(originalFile, clipPolygonFile, fileName):
    print "Recortando feature. Destino: " + fileName
    try:
        outPath, outFC = os.path.split(fileName)
        arcpy.env.workspace = outPath
        arcpy.env.OverWriteOutput = True

        arcpy.Clip_analysis(originalFile, clipPolygonFile, outFC, "")

    except Exception, err:
        print err.message
        arcpy.AddError(err.message)

    finally:
        print "Finalizado el recortado."


def calculate_areas(polygonFile, fileName):
    print "Calculando areas. Destino: " + fileName
    try:
        outPath, outFC = os.path.split(fileName)
        arcpy.env.workspace = outPath
        arcpy.env.OverWriteOutput = True

        arcpy.CalculateAreas_stats(polygonFile, outFC)

    except Exception, err:
        print err.message
        arcpy.AddError(err.message)

    finally:
        print "Finalizado el cálculo de areas."
