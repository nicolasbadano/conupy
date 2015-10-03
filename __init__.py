# -*- coding: utf-8 -*-
"""
/***************************************************************************
 ConuPy
                                 A QGIS plugin
 ConuPy for QGIS
                             -------------------
        begin                : 2015-10-02
        copyright            : (C) 2015 by Nicolás Diego Badano
        email                : nicolas.d.badano@gmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load ConuPy class from file ConuPy.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .conupy_module import ConuPy
    return ConuPy(iface)
