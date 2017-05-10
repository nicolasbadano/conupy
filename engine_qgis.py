# GIS functions for QGIS 2.10.1

import qgis.core
from PyQt4 import QtCore
import os
import processing
import numpy as np

default_spatial_reference = None
arcpy = []

def gis_init():
    # supply path to where is your qgis installed
    app = qgis.core.QgsApplication([],True)
    print app
    print qgis.core.QgsApplication.setPrefixPath("C:/Program Files/QGIS Pisa", True)

    # load providers
    print qgis.core.QgsApplication.initQgis()


def crear_datos_temporales():
    return arcpy.CreateScratchName(workspace=arcpy.env.scratchGDB)


def leer_spatial_reference(fileName):
    layer = qgis.core.QgsVectorLayer(fileName, "capa1", "ogr")
    spatial_reference = layer.crs()
    print spatial_reference
    return spatial_reference


def leer_shp_puntos(shp_file, lista_campos = []):
    layer = qgis.core.QgsVectorLayer(shp_file, "capa1", "ogr")

    resultados = []
    print "Leyendo shape de puntos ", shp_file

    features = layer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        x = geom.asPoint()

        poly = []
        punto = np.array([x[0], x[1]])
        poly.append(punto)

        # fetch attributes
        attrs = feature.attributes()
        # attrs is a list. It contains all the attribute values of this feature
        for campo in lista_campos:
            idx = layer.fieldNameIndex(campo)
            val = attrs[idx] if idx != -1 else None
            if isinstance(val, QtCore.QPyNullVariant):
                val = None
            poly.append(val)

        resultados.append(poly)

    print "Finalizado de leer shape de puntos."
    return resultados


def leer_shp_polilineas(shp_file, lista_campos = []):
    layer = qgis.core.QgsVectorLayer(shp_file, "capa1", "ogr")

    resultados = []
    print "Leyendo shape de polilineas ", shp_file

    features = layer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        x = geom.asPolyline()

        poly = []
        puntos = []
        for p in x:
            punto = np.array([p[0], p[1]])
            puntos.append(punto)

        if len(puntos) < 2:
            continue

        poly.append(puntos)

        # fetch attributes
        attrs = feature.attributes()
        # attrs is a list. It contains all the attribute values of this feature
        for campo in lista_campos:
            idx = layer.fieldNameIndex(campo)
            val = attrs[idx] if idx != -1 else None
            if isinstance(val, QtCore.QPyNullVariant):
                val = None
            poly.append(val)

        resultados.append(poly)

    print "Finalizado de leer shape de polilineas."
    return resultados


def leer_shp_poligonos(shp_file, lista_campos = []):
    layer = qgis.core.QgsVectorLayer(shp_file, "capa1", "ogr")

    resultados = []
    print "Leyendo shape de poligonos ", shp_file

    features = layer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        x = geom.asPolygon()

        poly = []
        puntos = []

        # Keep only the first ring (e.g. no holes)
        if len(x) < 1:
            continue

        for p in x[0]:
            punto = np.array([p[0], p[1]])
            puntos.append(punto)

        if len(puntos) < 2:
            continue

        poly.append(puntos)

        # fetch attributes
        attrs = feature.attributes()
        # attrs is a list. It contains all the attribute values of this feature
        for campo in lista_campos:
            idx = layer.fieldNameIndex(campo)
            val = attrs[idx] if idx != -1 else None
            if isinstance(val, QtCore.QPyNullVariant):
                val = None
            poly.append(val)

        resultados.append(poly)

    print "Finalizado de leer shape de polígonos."
    return resultados


def escribir_shp_puntos(fileName, puntos, campos, spatial_reference):
    print "Escribiendo shape de puntos: " + fileName
    outPath, outFC = os.path.split(fileName)

    fields = qgis.core.QgsFields()
    for nombre in campos:
        dato = campos[nombre][0]
        if isinstance(dato, int):
            fields.append(qgis.core.QgsField(nombre, QtCore.QVariant.Int))
        elif isinstance(dato, float):
            fields.append(qgis.core.QgsField(nombre, QtCore.QVariant.Double))
        elif isinstance(dato, str):
            fields.append(qgis.core.QgsField(nombre, QtCore.QVariant.String))

    writer = qgis.core.QgsVectorFileWriter(fileName, "CP1250", fields,
        qgis.core.QGis.WKBPoint, spatial_reference, "ESRI Shapefile")

    if writer.hasError() != qgis.core.QgsVectorFileWriter.NoError:
        print "Error when creating shapefile: ",  writer.errorMessage()

    # Add features
    for (i, punto) in enumerate(puntos):
        fet = qgis.core.QgsFeature()

        fet.setGeometry(qgis.core.QgsGeometry.fromPoint(qgis.core.QgsPoint(punto[0], punto[1])))
        fet.setAttributes([campos[nombre][i] for nombre in campos])

        writer.addFeature(fet)
        del fet

    del writer
    print "Finalizado escribir shape de puntos."


def escribir_shp_polilineas(fileName, polilineas, campos, spatial_reference):
    print "Escribiendo shape de polilineas: " + fileName

    fields = qgis.core.QgsFields()
    for nombre in campos:
        dato = campos[nombre][0]
        if isinstance(dato, int):
            fields.append(qgis.core.QgsField(nombre, QtCore.QVariant.Int))
        elif isinstance(dato, float):
            fields.append(qgis.core.QgsField(nombre, QtCore.QVariant.Double))
        elif isinstance(dato, str):
            fields.append(qgis.core.QgsField(nombre, QtCore.QVariant.String))

    writer = qgis.core.QgsVectorFileWriter(fileName, "CP1250", fields,
        qgis.core.QGis.WKBLineString, spatial_reference, "ESRI Shapefile")

    if writer.hasError() != qgis.core.QgsVectorFileWriter.NoError:
        print "Error when creating shapefile: ",  writer.errorMessage()

    # Add features
    for (i, polilinea) in enumerate(polilineas):
        fet = qgis.core.QgsFeature()
        puntos = []
        for (j, punto) in enumerate(polilinea):
            puntos.append(qgis.core.QgsPoint(float(punto[0]), float(punto[1])))

        fet.setGeometry(qgis.core.QgsGeometry.fromPolyline(puntos))
        fet.setAttributes([campos[nombre][i] for nombre in campos])

        writer.addFeature(fet)
        del fet

    del writer
    print "Finalizado escribir shape de polilineas."


def escribir_shp_poligonos(fileName, poligonos, campos, spatial_reference):
    print "Escribiendo shape de polígonos: " + fileName

    fields = qgis.core.QgsFields()
    for nombre in campos:
        dato = campos[nombre][0]
        if isinstance(dato, int):
            fields.append(qgis.core.QgsField(nombre, QtCore.QVariant.Int))
        elif isinstance(dato, float):
            fields.append(qgis.core.QgsField(nombre, QtCore.QVariant.Double))
        elif isinstance(dato, str):
            fields.append(qgis.core.QgsField(nombre, QtCore.QVariant.String))

    writer = qgis.core.QgsVectorFileWriter(fileName, "CP1250", fields,
        qgis.core.QGis.WKBPolygon, spatial_reference, "ESRI Shapefile")

    if writer.hasError() != qgis.core.QgsVectorFileWriter.NoError:
        print "Error when creating shapefile: ",  writer.errorMessage()

    # Add features
    for (i, polilinea) in enumerate(poligonos):
        fet = qgis.core.QgsFeature()
        puntos = []
        for (j, punto) in enumerate(polilinea):
            puntos.append(qgis.core.QgsPoint(punto[0], punto[1]))

        fet.setGeometry(qgis.core.QgsGeometry.fromPolygon(puntos))
        fet.setAttributes([campos[nombre][i] for nombre in campos])

        writer.addFeature(fet)
        del fet

    del writer
    print "Finalizado escribir shape de polígonos."


def sample_raster_on_nodes(nodesFile, rasterFile):
    print "Muestreando el raster:" + rasterFile

    nodesLayer = qgis.core.QgsVectorLayer(nodesFile, "nodeLayer", "ogr")
    rasterLayer = qgis.core.QgsRasterLayer(rasterFile, "rasterLayer")

    values = []
    features = nodesLayer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        point = geom.asPoint()

        res = rasterLayer.dataProvider().identify(point, qgis.core.QgsRaster.IdentifyFormatValue)
        if res.isValid():
            values.append(float(res.results().values()[0]))
        else:
            values.append(None)

    print "Finalizado el muestreado de valores."
    return values



def create_thiessen_polygons(nodesFile, fileName):
    print "Creando polígonos de Thiessen. Destino: " + fileName

    pointLayer = qgis.core.QgsVectorLayer(nodesFile, "capa1", "ogr")
    processing.runalg("qgis:voronoipolygons", pointLayer, 1, fileName)

    print "Finalizada creación de polígonos."


def clip_feature(originalFile, clipPolygonFile, fileName):
    print "Recortando feature. Destino: " + fileName

    originalLayer = qgis.core.QgsVectorLayer(originalFile, "capa3", "ogr")
    clipLayer = qgis.core.QgsVectorLayer(clipPolygonFile, "capa4", "ogr")
    processing.runalg("qgis:clip", originalLayer, clipLayer, fileName)

    print "Finalizado el recortado."


def read_areas(polygonFile):
    print "Leyendo areas: " + polygonFile

    layer = qgis.core.QgsVectorLayer(polygonFile, "capa1", "ogr")

    resultados = []

    features = layer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        resultados.append(geom.area())

    return resultados


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

def get_extent(shp_file):
    layer = qgis.core.QgsVectorLayer(shp_file, "capa1", "ogr")
    extent = layer.extent()
    return (extent.xMinimum(),
            extent.xMaximum(),
            extent.yMinimum(),
            extent.yMaximum())

def read_points_inside_features(shp_file_points, shp_file_polygons, lista_campos = []):
    layer = qgis.core.QgsVectorLayer(shp_file_points, "capa1", "ogr")
    layer_polygons = qgis.core.QgsVectorLayer(shp_file_polygons, "capa1", "ogr")

    resultados = []
    print "Leyendo shape de puntos ", shp_file_points

    features = layer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        x = geom.asPoint()
        punto = np.array([x[0], x[1]])

        for polyfeature in layer_polygons.getFeatures():
            if polyfeature.geometry().contains(x):
                # fetch attributes
                attrs = feature.attributes()
                # attrs is a list. It contains all the attribute values of this feature
                for campo in lista_campos:
                    idx = layer.fieldNameIndex(campo)
                    val = attrs[idx] if idx != -1 else None
                    qgis.core.poly.append(val if val is not None else None)

                resultados.append(punto)

    print "Finalizado de leer shape de puntos."
    return resultados

def mark_points_inside_features(puntos, shp_file_polygons):
    layer_polygons = qgis.core.QgsVectorLayer(shp_file_polygons, "capa1", "ogr")

    for punto in puntos:
        p = qgis.core.QgsPoint(punto.p[0], punto.p[1])
        punto.inside = False
        for polyfeature in layer_polygons.getFeatures():
            if polyfeature.geometry().contains(p):
                punto.inside = True
