# GIS functions for QGIS 2.10.1

from qgis.core import *
from PyQt4 import QtCore
import os
import processing

default_spatial_reference = None
arcpy = []

def gis_init():
    # supply path to where is your qgis installed
    app = QgsApplication([],True)
    print app
    print QgsApplication.setPrefixPath("C:/Program Files/QGIS Pisa", True)

    # load providers
    print QgsApplication.initQgis()


def crear_datos_temporales():
    return arcpy.CreateScratchName(workspace=arcpy.env.scratchGDB)


def leer_spatial_reference(fileName):
    layer = QgsVectorLayer(fileName, "capa1", "ogr")
    spatial_reference = layer.crs()
    print spatial_reference
    return spatial_reference


def leer_shp_puntos(shp_file, lista_campos = []):
    layer = QgsVectorLayer(shp_file, "capa1", "ogr")

    resultados = []
    print "Leyendo shape de puntos ", shp_file

    features = layer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        x = geom.asPoint()
        punto = [x[0], x[1]]

        # fetch attributes
        attrs = feature.attributes()
        # attrs is a list. It contains all the attribute values of this feature
        for campo in lista_campos:
            idx = layer.fieldNameIndex(campo)
            val = attrs[idx] if idx != -1 else None
            poly.append(val if val else None)

        resultados.append(punto)

    print "Finalizado de leer shape de puntos."
    return resultados


def leer_shp_polilineas(shp_file, lista_campos = []):
    layer = QgsVectorLayer(shp_file, "capa1", "ogr")

    resultados = []
    print "Leyendo shape de polilineas ", shp_file

    features = layer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        x = geom.asPolyline()

        poly = []
        puntos = []
        for p in x:
            punto = [p[0], p[1]]
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
            poly.append(val if val else None)

        resultados.append(poly)

    print "Finalizado de leer shape de polilineas."
    return resultados


def leer_shp_poligonos(shp_file, lista_campos = []):
    layer = QgsVectorLayer(shp_file, "capa1", "ogr")

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
            punto = [p[0], p[1]]
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
            poly.append(val if val else None)

        resultados.append(poly)

    print "Finalizado de leer shape de polígonos."
    return resultados


def escribir_shp_puntos(fileName, puntos, campos, spatial_reference):
    print "Escribiendo shape de puntos: " + fileName
    outPath, outFC = os.path.split(fileName)

    fields = QgsFields()
    for nombre in campos:
        dato = campos[nombre][0]
        type = ""
        if isinstance(dato, int):
            fields.append(QgsField(nombre, QtCore.QVariant.Int))
        elif isinstance(dato, float):
            fields.append(QgsField(nombre, QtCore.QVariant.Double))
        elif isinstance(dato, str):
            fields.append(QgsField(nombre, QtCore.QVariant.String))

    writer = QgsVectorFileWriter(fileName, "CP1250", fields, QGis.WKBPoint, spatial_reference, "ESRI Shapefile")

    if writer.hasError() != QgsVectorFileWriter.NoError:
        print "Error when creating shapefile: ",  writer.errorMessage()

    # Add features
    for (i, punto) in enumerate(puntos):
        fet = QgsFeature()

        fet.setGeometry(QgsGeometry.fromPoint(QgsPoint(punto[0], punto[1])))
        fet.setAttributes([campos[nombre][i] for nombre in campos])

        writer.addFeature(fet)
        del fet

    del writer
    print "Finalizado escribir shape de puntos."


def escribir_shp_polilineas(fileName, polilineas, campos, spatial_reference):
    print "Escribiendo shape de polilineas: " + fileName

    fields = QgsFields()
    for nombre in campos:
        dato = campos[nombre][0]
        type = ""
        if isinstance(dato, int):
            fields.append(QgsField(nombre, QtCore.QVariant.Int))
        elif isinstance(dato, float):
            fields.append(QgsField(nombre, QtCore.QVariant.Double))
        elif isinstance(dato, str):
            fields.append(QgsField(nombre, QtCore.QVariant.String))

    writer = QgsVectorFileWriter(fileName, "CP1250", fields, QGis.WKBLineString, spatial_reference, "ESRI Shapefile")

    if writer.hasError() != QgsVectorFileWriter.NoError:
        print "Error when creating shapefile: ",  writer.errorMessage()

    # Add features
    for (i, polilinea) in enumerate(polilineas):
        fet = QgsFeature()
        puntos = []
        for (j, punto) in enumerate(polilinea):
            puntos.append(QgsPoint(punto[0], punto[1]))

        fet.setGeometry(QgsGeometry.fromPolyline(puntos))
        fet.setAttributes([campos[nombre][i] for nombre in campos])

        writer.addFeature(fet)
        del fet

    del writer
    print "Finalizado escribir shape de polilineas."


def escribir_shp_poligonos(fileName, poligonos, campos, spatial_reference):
    print "Escribiendo shape de polígonos: " + fileName

    fields = QgsFields()
    for nombre in campos:
        dato = campos[nombre][0]
        type = ""
        if isinstance(dato, int):
            fields.append(QgsField(nombre, QtCore.QVariant.Int))
        elif isinstance(dato, float):
            fields.append(QgsField(nombre, QtCore.QVariant.Double))
        elif isinstance(dato, str):
            fields.append(QgsField(nombre, QtCore.QVariant.String))

    writer = QgsVectorFileWriter(fileName, "CP1250", fields, QGis.WKBPolygon, spatial_reference, "ESRI Shapefile")

    if writer.hasError() != QgsVectorFileWriter.NoError:
        print "Error when creating shapefile: ",  writer.errorMessage()

    # Add features
    for (i, polilinea) in enumerate(poligonos):
        fet = QgsFeature()
        puntos = []
        for (j, punto) in enumerate(polilinea):
            puntos.append(QgsPoint(punto[0], punto[1]))

        fet.setGeometry(QgsGeometry.fromPolygon(puntos))
        fet.setAttributes([campos[nombre][i] for nombre in campos])

        writer.addFeature(fet)
        del fet

    del writer
    print "Finalizado escribir shape de polígonos."


def sample_raster_on_nodes(nodesFile, rasterFile):
    print "Muestreando el raster:" + rasterFile

    nodesLayer = QgsVectorLayer(nodesFile, "nodeLayer", "ogr")
    rasterLayer = QgsRasterLayer(rasterFile, "rasterLayer")

    bandcount = rasterLayer.bandCount()

    values = []
    features = nodesLayer.getFeatures()
    for feature in features:
        geom = feature.geometry()
        point = geom.asPoint()

        res = rasterLayer.dataProvider().identify(point, QgsRaster.IdentifyFormatValue)
        if res.isValid():
            values.append(res.results().values()[0])
        else:
            values.append(None)

    print "Finalizado el muestreado de valores."
    return values



def create_thiessen_polygons(nodesFile, fileName):
    print "Creando polígonos de Thiessen. Destino: " + fileName

    pointLayer = QgsVectorLayer(nodesFile, "capa1", "ogr")
    processing.runalg("qgis:voronoipolygons", pointLayer, 1, fileName)

    print "Finalizada creación de polígonos."


def clip_feature(originalFile, clipPolygonFile, fileName):
    print "Recortando feature. Destino: " + fileName

    originalLayer = QgsVectorLayer(originalFile, "capa1", "ogr")
    clipLayer = QgsVectorLayer(clipPolygonFile, "capa2", "ogr")
    processing.runalg("qgis:clip", originalLayer, clipLayer, fileName)

    print "Finalizado el recortado."


def read_areas(polygonFile):
    print "Leyendo areas: " + polygonFile

    layer = QgsVectorLayer(polygonFile, "capa1", "ogr")

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
