from PySide2.QtWidgets import *
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import *
from PySide2.QtGui import QPixmap

from main import MainWindow
import homeui

# pyside2-uic -x main.ui -o mainui.py


class MainWidget (QWidget, homeui.Ui_Form):

    def __init__(self, parent=None):
        super(MainWidget, self).__init__(parent)

        #f = QFile("home.ui")
        # f.open(QFile.ReadOnly)

        #loader = QUiLoader()
        #self.ui = loader.load(f,  self)

        # f.close()
        #logo = self.ui.findChild(QLabel, 'logoLabel')
        # logo.setPixmap(QPixmap('qt.png'))

        #version = self.ui.findChild(QLabel, 'versionLabel')
        # version.setText()

        #importButton = self.ui.findChild(QPushButton, 'createProjectButton')
        # importButton.clicked.connect(self.openMainWindow)

        self.setupUi(self)
        self.logoLabel.setPixmap(QPixmap('qt.png'))
        self.createProjectButton.clicked.connect(self.openMainWindow)

    def openMainWindow(self):
        self.mainUIWindow = MainWindow()
        self.hide()
        print('A new window is opened!')


if __name__ == '__main__':
    QCoreApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
    app = QApplication([])
    window = MainWidget()
    window.show()
    app.exec_()
