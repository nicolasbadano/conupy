# ConuPyWindow

import sys
import os
from PyQt4 import QtCore, QtGui, uic

from math import sqrt, pow, exp
import pickle

# Load the UI
path = os.getcwd() + "/"
print path
path = "F:/Desarrollo/Utilidades/conuPy/ConuPy/"

form_class = uic.loadUiType(path + "conuPy.ui")[0]

class ConuPyWindow(QtGui.QMainWindow, form_class):

    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.setupUi(self)

        # # Enlazar los Event Handlers de los botones
        # self.botonNuevaQuebrada.clicked.connect(self.botonNuevaQuebrada_click)
        # self.botonBorrarQuebrada.clicked.connect(self.botonBorrarQuebrada_click)
        # self.botonNuevoTalud.clicked.connect(self.botonNuevoTalud_click)
        # self.botonBorrarTalud.clicked.connect(self.botonBorrarTalud_click)
        # self.botonCalcular.clicked.connect(self.botonCalcular_click)
        # self.botonSeleccionarShapeTraza.clicked.connect(self.botonSeleccionarShapeTraza_click)
        # self.botonSeleccionarMDE.clicked.connect(self.botonSeleccionarMDE_click)
        # self.botonCopiarSeccion.clicked.connect(self.botonCopiarSeccion_click)
        # self.botonPegarSeccion.clicked.connect(self.botonPegarSeccion_click)
        # self.botonCopiarHietograma.clicked.connect(self.botonCopiarHietograma_click)
        # self.botonPegarHietograma.clicked.connect(self.botonPegarHietograma_click)
        # # Enlazar los Event Handlers de las cajas de texto
        # self.connect(self.textNombreQuebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.textProgresivaQuebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.textNivelConductoQuebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.textAreaQuebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.textCoefEscorrentiaQuebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.textPendLongQuebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.textNManningQuebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.textCoefBetaQuebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.textD50Quebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.textGammaSecoQuebrada,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_quebrada)
        # self.connect(self.comboSedimQuebrada,
        #                 QtCore.SIGNAL("currentIndexChanged(int)"),
        #                 self.leer_datos_quebrada)

        # self.connect(self.textNombreTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)
        # self.connect(self.textProgresivaTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)
        # self.connect(self.textHTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)
        # self.connect(self.textAlfaTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)
        # self.connect(self.textRhoGranoTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)
        # self.connect(self.textFiTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)
        # self.connect(self.textCTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)
        # self.connect(self.textPorosidadTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)
        # self.connect(self.textHumedadAntesTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)
        # self.connect(self.textTasaInfiltracionTalud,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_talud)

        # self.connect(self.textArchivoTraza,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_generales)
        # self.connect(self.textArchivoMDE,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_generales)
        # self.connect(self.textK,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_generales)
        # self.connect(self.textC,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_generales)
        # self.connect(self.textP,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_generales)
        # self.connect(self.textRho,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_generales)
        # self.connect(self.textAnchoPista,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_generales)
        # self.connect(self.textAnchoZanja,
        #                 QtCore.SIGNAL("textEdited(QString)"),
        #                 self.leer_datos_generales)

        # # Enlazar los Event Handlers
        # self.connect(self.tableSeccion,
        #                 QtCore.SIGNAL("itemChanged(QTableWidgetItem*)"),
        #                 self.actualizar_puntos_seccion)
        # self.actualizando_tabla = False
        # self.connect(self.tableHietograma,
        #                 QtCore.SIGNAL("itemChanged(QTableWidgetItem*)"),
        #                 self.actualizar_puntos_hietograma)
        # self.actualizando_tabla_hietograma = False

        # self.connect(self.actionNuevo,
        #                 QtCore.SIGNAL("triggered()"),
        #                 self.nuevo)
        # self.connect(self.actionAbrir,
        #                 QtCore.SIGNAL("triggered()"),
        #                 self.abrir)
        # self.connect(self.actionGuardar,
        #                 QtCore.SIGNAL("triggered()"),
        #                 self.guardar)
        # self.connect(self.actionGuardarComo,
        #                 QtCore.SIGNAL("triggered()"),
        #                 self.guardar_como)
        # self.connect(self.actionSalir,
        #                 QtCore.SIGNAL("triggered()"),
        #                 self.close)

        # # Agregar validadores a los campos
        # self.textProgresivaQuebrada.setValidator(QtGui.QDoubleValidator(0,1e6, 6, self.textProgresivaQuebrada));
        # self.textAreaQuebrada.setValidator(QtGui.QDoubleValidator(0,1e6, 6, self.textAreaQuebrada));
        # self.textNivelConductoQuebrada.setValidator(QtGui.QDoubleValidator(-1000,8000, 3, self.textNivelConductoQuebrada));
        # self.textCoefEscorrentiaQuebrada.setValidator(QtGui.QDoubleValidator(0,1e6, 6, self.textCoefEscorrentiaQuebrada));
        # self.textPendLongQuebrada.setValidator(QtGui.QDoubleValidator(0,1e6, 6, self.textPendLongQuebrada));
        # self.textNManningQuebrada.setValidator(QtGui.QDoubleValidator(0,1e6, 6, self.textNManningQuebrada));
        # self.textCoefBetaQuebrada.setValidator(QtGui.QDoubleValidator(0.5,1.5, 6, self.textCoefBetaQuebrada));
        # self.textD50Quebrada.setValidator(QtGui.QDoubleValidator(1e-4,1000, 6, self.textD50Quebrada));
        # self.textGammaSecoQuebrada.setValidator(QtGui.QDoubleValidator(0,10, 6, self.textGammaSecoQuebrada));

        # self.textProgresivaTalud.setValidator(QtGui.QDoubleValidator(0,1e6, 6, self.textProgresivaTalud));
        # self.textHTalud.setValidator(QtGui.QDoubleValidator(0,100, 6, self.textHTalud));
        # self.textAlfaTalud.setValidator(QtGui.QDoubleValidator(0,90, 6, self.textAlfaTalud));
        # self.textFiTalud.setValidator(QtGui.QDoubleValidator(0,90, 6, self.textFiTalud));
        # self.textRhoGranoTalud.setValidator(QtGui.QDoubleValidator(0,10000, 2, self.textRhoGranoTalud));
        # self.textCTalud.setValidator(QtGui.QDoubleValidator(0,1000, 6, self.textCTalud));
        # self.textPorosidadTalud.setValidator(QtGui.QDoubleValidator(0.01,0.99, 6, self.textPorosidadTalud));
        # self.textHumedadAntesTalud.setValidator(QtGui.QDoubleValidator(0.0,2.0, 6, self.textHumedadAntesTalud));
        # self.textTasaInfiltracionTalud.setValidator(QtGui.QDoubleValidator(0.001,1000, 6, self.textTasaInfiltracionTalud));

        # self.textK.setValidator(QtGui.QDoubleValidator(0,100000, 3, self.textK));
        # self.textC.setValidator(QtGui.QDoubleValidator(0,100000, 3, self.textC));
        # self.textP.setValidator(QtGui.QDoubleValidator(0,100000, 3, self.textP));
        # self.textRho.setValidator(QtGui.QDoubleValidator(0,10000, 2, self.textRho));
        # self.textAnchoPista.setValidator(QtGui.QDoubleValidator(0,100, 2, self.textAnchoPista));
        # self.textAnchoZanja.setValidator(QtGui.QDoubleValidator(0,100, 2, self.textAnchoZanja));

        # self.quebradas = []
        # self.taludes = []

        # self.quebradasModelo = ModeloLista(None, self.quebradas)
        # self.listViewQuebradas.setModel(self.quebradasModelo)
        # self.listViewQuebradas.selectionModel().selectionChanged.connect(self.listViewQuebradas_selectionChanged)

        # self.taludesModelo = ModeloLista(None, self.taludes)
        # self.listViewTaludes.setModel(self.taludesModelo)
        # self.listViewTaludes.selectionModel().selectionChanged.connect(self.listViewTaludes_selectionChanged)

        self.nuevo()

    def nuevo(self):
        # self.quebrada_seleccionada = None
        # self.talud_seleccionado = None
        # self.nombreArchivo = None
        # self.archivo_traza = ""
        # self.archivo_MDE = ""
        # self.k = 0.003
        # self.c = 0.45
        # self.p = 1.0
        # self.rho = 1900.0
        # self.ancho_pista = 10.0
        # self.ancho_zanja = 1.0
        # self.actualizar_datos_generales()
        pass

    def devolver_diccionario(self):
        variables = [ \
            "archivo_traza", \
            "archivo_MDE", \
            "k", \
            "c", \
            "p", \
            "rho", \
            "ancho_pista", \
            "ancho_zanja", \
            "hietograma"]

        return dict([ (name, getattr(self, name)) for name in variables])

    def abrir(self):
        text = QtGui.QFileDialog.getOpenFileName(self, 'Seleccionar proyecto de estrata',
            path, "Proyecto ConuPy (*.conuPy)")
        print text
        if text is not None and not str(text) == "":
            self.nuevo()

            try:
                iF = open(text, 'rb')
                header = pickle.load(iF)
                qds = pickle.load(iF)
                tds = pickle.load(iF)
                md = pickle.load(iF)
                iF.close()
            except:
                print "Error al leer el Proyecto Estrata"

            # Cargar variables principales
            if not md is None:
                for name, value in md.iteritems():
                    setattr(self, name, value)
            # Crear quebradas a partir de los diccionarios
            if not qds is None:
                for d in qds:
                    self.quebradas.append(Quebrada(d))
            # Crear taludes a partir de los diccionarios
            if not tds is None:
                for d in tds:
                    self.taludes.append(Talud(d))

            self.quebrada_seleccionada = None
            self.talud_seleccionado = None
            self.nombreArchivo = text
            self.actualizar_datos_generales()

    def guardar(self):
        if self.nombreArchivo is None:
            self.guardar_como()
        else:
            md = self.devolver_diccionario()
            qds = [q.devolver_diccionario() for q in self.quebradas]
            tds = [t.devolver_diccionario() for t in self.taludes]

            try:
                oF = open(self.nombreArchivo, 'wb')
                pickle.dump("ESTRATA2", oF)
                pickle.dump(qds, oF)
                pickle.dump(tds, oF)
                pickle.dump(md, oF)
                oF.close()
            except:
                print "Error al escribir el Proyecto Estrata"

    def guardar_como(self):
        text = QtGui.QFileDialog.getSaveFileName(self, 'Seleccionar proyecto de estrata',
            path, "Proyecto ConuPy (*.conuPy)")
        if text is not None and not str(text) == "":
            self.nombreArchivo = text
            self.guardar()


    # Event handlers de botonSeleccionarShapeTraza
    def botonSeleccionarShapeTraza_click(self):
        # text = QtGui.QFileDialog.getOpenFileName(self, 'Seleccionar shape de traza',
        #     path, "Shape (*.shp)")
        # if text is not None and not str(text) == "":
        #     self.textArchivoTraza.setText(text)
        #     self.leer_datos_generales()

    # Event handlers de botonCalcular
    # def botonCalcular_click(self):
    #     # Crear traza

    #     print "Calculo finalizado"
        pass

    def actualizar_datos_generales(self):
        try:
            pass
            # self.textK.setText(str(self.k))
            # self.textC.setText(str(self.c))
            # self.textP.setText(str(self.p))
            # self.textRho.setText(str(self.rho))
            # self.textAnchoPista.setText(str(self.ancho_pista))
            # self.textAnchoZanja.setText(str(self.ancho_zanja))
            # self.textArchivoTraza.setText(self.archivo_traza)
            # self.textArchivoMDE.setText(self.archivo_MDE)
        except: pass

    def leer_datos_generales(self):
        try:
            pass
            # self.k = float(self.textK.text())
            # self.c = float(self.textC.text())
            # self.p = float(self.textP.text())
            # self.rho = float(self.textRho.text())
            # self.ancho_pista = float(self.textAnchoPista.text())
            # self.ancho_zanja = float(self.textAnchoZanja.text())
            # self.archivo_traza = str(self.textArchivoTraza.text())
            # self.archivo_MDE = str(self.textArchivoMDE.text())
        except: pass

if __name__ == '__main__':
    # Procedimiento principal
    aplicacion = QtGui.QApplication(sys.argv)
    ventana_principal = VentanaPrincipal(None)
    ventana_principal.show()

    sys.exit(aplicacion.exec_())
