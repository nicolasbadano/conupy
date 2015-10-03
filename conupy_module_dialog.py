# -*- coding: utf-8 -*-
"""
/***************************************************************************
 ConuPyDialog
                                 A QGIS plugin
 ConuPy for QGIS
                             -------------------
        begin                : 2015-10-02
        git sha              : $Format:%H$
        copyright            : (C) 2015 by Nicolás Diego Badano
        email                : nicolas.d.badano@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os
import sys
import pickle
from PyQt4 import QtCore, QtGui, uic
import conuPy

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'conupy_module_dialog_base.ui'))

path = os.path.dirname(__file__)

class WriteStreamInmediate(object):
    def __init__(self, textBox):
        self.textBox = textBox

    def write(self, text):
        self.textBox.moveCursor(QtGui.QTextCursor.End)
        self.textBox.insertPlainText(text)
        QtGui.QApplication.processEvents()


class ConuPyDialog(QtGui.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(ConuPyDialog, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)

        self.last_path = path

        # Enlazar los Event Handlers de los botones
        self.actionNuevo.clicked.connect(self.nuevo)
        self.actionAbrir.clicked.connect(self.abrir)
        self.actionGuardar.clicked.connect(self.guardar)
        self.actionGuardarComo.clicked.connect(self.guardar_como)

        self.btnArroyos.clicked.connect(self.btnArroyos_click)
        self.btnCalles.clicked.connect(self.btnCalles_click)
        self.btnCuenca.clicked.connect(self.btnCuenca_click)
        self.btnNodosBorde.clicked.connect(self.btnNodosBorde_click)
        self.btnMDE.clicked.connect(self.btnMDE_click)
        self.btnSlope.clicked.connect(self.btnSlope_click)
        self.btnCoeficiente.clicked.connect(self.btnCoeficiente_click)
        self.btnImpermeabilidad.clicked.connect(self.btnImpermeabilidad_click)
        self.btnWorkspace.clicked.connect(self.btnWorkspace_click)
        self.btnModelDirectory.clicked.connect(self.btnModelDirectory_click)
        self.btnRain.clicked.connect(self.btnRain_click)

        self.btnActionLimpiarWorkspace.clicked.connect(self.btnActionLimpiarWorkspace_click)
        self.btnActionArroyosCalles.clicked.connect(self.btnActionArroyosCalles_click)
        self.btnActionSubcuencas.clicked.connect(self.btnActionSubcuencas_click)
        self.btnActionSample.clicked.connect(self.btnActionSample_click)
        self.btnActionBordes.clicked.connect(self.btnActionBordes_click)
        self.btnActionInverts.clicked.connect(self.btnActionInverts_click)
        self.btnActionSWMM.clicked.connect(self.btnActionSWMM_click)
        self.btnActionGenerateRain.clicked.connect(self.btnActionGenerateRain_click)
        self.btnActionProfMuertas.clicked.connect(self.btnActionProfMuertas_click)
        self.btnActionExtraerProf.clicked.connect(self.btnActionExtraerProf_click)

        # Enlazar los Event Handlers de las cajas de texto
        self.connect(self.textShpFileArroyos,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textShpFileCalles,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textShpFileCuenca,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textShpFileNodosBorde,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textRasterFileDEM,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textRasterFileSlope,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textRasterFileCoeficiente,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textRasterFileImpermeabilidad,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textWorkspace,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textModelDirectory,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textModelFileName,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.connect(self.textRainFileName,
                        QtCore.SIGNAL("textEdited(QString)"),
                        self.leer_datos_generales)
        self.nuevo()

    def nuevo(self):
        self.nombreArchivo = None

        self.shpFileArroyos = ""
        self.shpFileCalles = ""
        self.shpFileCuenca = ""
        self.shpFileNodosBorde = ""
        self.rasterFileDEM = ""
        self.rasterFileSlope = ""
        self.rasterFileCoeficiente = ""
        self.rasterFileImpermeabilidad = ""
        self.workspace = ""
        self.modelFolder = ""
        self.modelFileName = "model"
        self.rainFileName = ""

        self.actualizar_datos_generales()

    def devolver_diccionario(self):
        variables = [ \
            "shpFileArroyos", \
            "shpFileCalles", \
            "shpFileCuenca", \
            "shpFileNodosBorde", \
            "rasterFileDEM", \
            "rasterFileSlope", \
            "rasterFileCoeficiente", \
            "rasterFileImpermeabilidad", \
            "workspace", \
            "modelFolder", \
            "modelFileName", \
            "rainFileName"]

        return dict([ (name, getattr(self, name)) for name in variables])

    def abrir(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar proyecto de ConuPy',
            self.last_path, u"Proyecto ConuPy (*.conuPy)")
        print text
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.nuevo()

            try:
                iF = open(text, 'rb')
                header = pickle.load(iF)
                md = pickle.load(iF)
                iF.close()
            except:
                print "Error al leer el Proyecto ConuPy"

            # Cargar variables principales
            if not md is None:
                for name, value in md.iteritems():
                    setattr(self, name, value)

            self.nombreArchivo = text
            self.actualizar_datos_generales()

    def guardar(self):
        if self.nombreArchivo is None:
            self.guardar_como()
        else:
            md = self.devolver_diccionario()

            try:
                oF = open(self.nombreArchivo, 'wb')
                pickle.dump("CONUPY2", oF)
                pickle.dump(md, oF)
                oF.close()
            except:
                print "Error al escribir el Proyecto ConuPy"

    def guardar_como(self):
        text = QtGui.QFileDialog.getSaveFileName(self, u'Seleccionar proyecto de ConuPy',
            self.last_path, u"Proyecto ConuPy (*.conuPy)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.nombreArchivo = text
            self.guardar()

    def actualizar_datos_generales(self):
        try:
            self.textShpFileArroyos.setText(str(self.shpFileArroyos))
            self.textShpFileCalles.setText(str(self.shpFileCalles))
            self.textShpFileCuenca.setText(str(self.shpFileCuenca))
            self.textShpFileNodosBorde.setText(str(self.shpFileNodosBorde))
            self.textRasterFileDEM.setText(str(self.rasterFileDEM))
            self.textRasterFileSlope.setText(str(self.rasterFileSlope))
            self.textRasterFileCoeficiente.setText(str(self.rasterFileCoeficiente))
            self.textRasterFileImpermeabilidad.setText(str(self.rasterFileImpermeabilidad))
            self.textWorkspace.setText(str(self.workspace))
            self.textModelDirectory.setText(str(self.modelFolder))
            self.textModelFileName.setText(str(self.modelFileName))
            self.textRainFileName.setText(str(self.rainFileName))
        except: pass

    def leer_datos_generales(self):
        try:
            pass
            self.shpFileArroyos = self.textShpFileArroyos.text()
            self.shpFileCalles = self.textShpFileCalles.text()
            self.shpFileCuenca = self.textShpFileCuenca.text()
            self.shpFileNodosBorde = self.textShpFileNodosBorde.text()
            self.rasterFileDEM = self.textRasterFileDEM.text()
            self.rasterFileSlope = self.textRasterFileSlope.text()
            self.rasterFileCoeficiente = self.textRasterFileCoeficiente.text()
            self.rasterFileImpermeabilidad = self.textRasterFileImpermeabilidad.text()
            self.workspace = self.textWorkspace.text()
            self.modelFolder = self.textModelDirectory.text()
            self.modelFileName = self.textModelFileName.text()
            self.rainFileName = self.textRainFileName.text()
        except: pass

    # Event handlers de btnArroyos
    def btnArroyos_click(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar shape de arroyos',
            self.last_path, u"Shape (*.shp)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.textShpFileArroyos.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnCalles
    def btnCalles_click(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar shape de calles',
            self.last_path, u"Shape (*.shp)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.textShpFileCalles.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnCuenca
    def btnCuenca_click(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar shape de cuenca',
            self.last_path, u"Shape (*.shp)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.textShpFileCuenca.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnNodosBorde
    def btnNodosBorde_click(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar shape de nodos de borde',
            self.last_path, u"Shape (*.shp)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.textShpFileNodosBorde.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnMDE
    def btnMDE_click(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar raster de mde',
            self.last_path, u"Raster (*.adf *.img)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.textRasterFileDEM.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnSlope
    def btnSlope_click(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar raster de pendiente',
            self.last_path, u"Raster (*.adf *.img)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.textRasterFileSlope.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnCoeficiente
    def btnCoeficiente_click(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar raster de coeficiente de precipitación',
            self.last_path, u"Raster (*.adf *.img)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.textRasterFileCoeficiente.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnImpermeabilidad
    def btnImpermeabilidad_click(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar raster de Impermeabilidad de precipitación',
            self.last_path, u"Raster (*.adf *.img)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.textRasterFileImpermeabilidad.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnWorkspace
    def btnWorkspace_click(self):
        dialog = QtGui.QFileDialog()
        dialog.setFileMode(QtGui.QFileDialog.Directory)
        dialog.setOption(QtGui.QFileDialog.ShowDirsOnly)
        dialog.setDirectory(self.last_path)
        if dialog.exec_():
            text = dialog.selectedFiles()[0]
            self.last_path = text

            self.textWorkspace.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnModelDirectory
    def btnModelDirectory_click(self):
        dialog = QtGui.QFileDialog()
        dialog.setFileMode(QtGui.QFileDialog.Directory)
        dialog.setOption(QtGui.QFileDialog.ShowDirsOnly)
        dialog.setDirectory(self.last_path)
        if dialog.exec_():
            text = dialog.selectedFiles()[0]
            self.last_path = text

            self.textModelDirectory.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnRain
    def btnRain_click(self):
        text = QtGui.QFileDialog.getOpenFileName(self, u'Seleccionar archivo precipitación',
            self.last_path, u"Precipitación (*.in)")
        if text is not None and not str(text) == "":
            self.last_path = os.path.dirname(text)
            self.textRainFileName.setText(text)
            self.leer_datos_generales()

    # Event handlers de btnActionLimpiarWorkspace
    def btnActionLimpiarWorkspace_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainCleanWorkspace(self.workspace)

    # Event handlers de btnActionArroyosCalles
    def btnActionArroyosCalles_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainReadRivers(self.shpFileArroyos)
        conuPy.mainReadStreets(self.shpFileCalles)

    # Event handlers de btnActionSubcuencas
    def btnActionSubcuencas_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainGetSubcatchments(self.shpFileCuenca)

    # Event handlers de btnActionSample
    def btnActionSample_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainSampleNodeData(self.rasterFileDEM,
                                  self.rasterFileSlope,
                                  self.rasterFileCoeficiente,
                                  self.rasterFileImpermeabilidad)

    # Event handlers de btnActionBordes
    def btnActionBordes_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainCreateOutfallNodes(self.shpFileNodosBorde)

    # Event handlers de btnActionInverts
    def btnActionInverts_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainCalculateInvertOffsets()

    # Event handlers de btnActionSWMM
    def btnActionSWMM_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainCreateSWMM(self.modelFolder + "/" + self.modelFileName + ".inp")

    # Event handlers de btnActionGenerateRain
    def btnActionGenerateRain_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainCreateGenerateRain(self.rainFileName, self.modelFolder + "/" + "pluviom.dat")

    # Event handlers de btnActionProfMuertas
    def btnActionProfMuertas_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainCalculateDeadDepths()

    # Event handlers de btnActionExtraerProf
    def btnActionExtraerProf_click(self):
        # Redireccionar stdout a la caja de texto
        sys.stdout = WriteStreamInmediate(self.textInfo)

        os.chdir(self.workspace)
        conuPy.mainReadSWMMResultsDepths(self.modelFolder + "/" + self.modelFileName + ".out")
