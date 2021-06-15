#!/usr/bin/env python

import sys
import CronosWidgetDesign
import CronosWidgeteer
from PyQt4 import QtCore, QtGui


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    Form = QtGui.QWidget()
    WidgetTier = CronosWidgeteer.CronosWidgeteer()
    WidgetTier.setupWidget(Form)
    WidgetTier.show()
    sys.exit(app.exec_())
    

