import operator
import sys
import csv
import numpy as np
import pandas as pd

from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.preprocessing import PolynomialFeatures


from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from PySide2.QtUiTools import QUiLoader

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvas,  NavigationToolbar2QT as NavigationToolbar

# pyside2-uic - x main.ui - o mainui.py
# export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"

from matplotlib.figure import Figure
from matplotlib import style
style.use('ggplot')


class MainWindow(QMainWindow, QObject):

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.center()

        f = QFile('main.ui')
        f.open(QFile.ReadOnly)

        loader = QUiLoader()
        self.ui = loader.load(f)
        self.ui.show()
        f.close()
        
        #---------------------------------Home---------------------------------
        self.optiDrill = self.ui.findChild(QLabel, 'titleLabel')
        logo = self.ui.findChild(QLabel, 'logoLabel')
        logo.setPixmap(QPixmap('logo.png').scaled(600, 400))

        version = self.ui.findChild(QLabel, 'versionLabel')
        version.setText('Version 0.1.1 \n All Right Reserved')

        # 1: Length conversion:
        #lengthUnit = ['--select unit--', 'mm', 'cm', 'm', 'km', 'in', '1/64in', 'ft', 'yard']
        self.lengthUnit = {'m': 1, 'mm': 1000, 'cm': 100, 'km': 0.001}
        units = sorted(self.lengthUnit.keys())

        self.lengthFromCB = self.ui.findChild(QComboBox, 'lengthFromComboBox')
        self.lengthFromCB.addItems(units)
        self.lengthToCB = self.ui.findChild(QComboBox, 'lengthToComboBox')
        self.lengthToCB.addItems(units)
        self.lengthFromLE = self.ui.findChild(QLineEdit, 'lengthFromLineEdit')
        self.lengthToLE = self.ui.findChild(QLineEdit, 'lengthToLineEdit')

        self.connect(self.lengthFromCB, SIGNAL('currentIndexChanged(int)'), self.convertLength)
        self.connect(self.lengthToCB, SIGNAL('currentIndexChanged(int)'), self.convertLength)
        self.connect(self.lengthFromLE, SIGNAL('editingFinished()'), self.convertLength)

        # 2: Area Conversion:
        # areaUnit = [u'ft\u00b2', u'cm\u00b2', 'acre', u'm\u00b2', 'ha', u'in\u00b2', u'cm\u00b2', u'mm\u00b2']
        self.areaUnits = {'sq.cm': 10000, 'sq. meter': 1, 'sq. feet': 10.7639, 'ha': 1e-4}
        #sortedAreaUnit = sorted(self.areaUnits.keys())
        areaUnits = list(self.areaUnits.keys())

        self.areaFromCB = self.ui.findChild(QComboBox, 'areaFromComboBox')
        self.areaFromCB.addItems(areaUnits)
        self.areaToCB = self.ui.findChild(QComboBox, 'areaToComboBox')
        self.areaToCB.addItems(areaUnits)
        self.areaFromLE = self.ui.findChild(QLineEdit, 'areaFromLineEdit')
        self.areaToLE = self.ui.findChild(QLineEdit, 'areaToLineEdit')

        self.connect(self.areaFromCB, SIGNAL('currentIndexChanged(int)'), self.convertArea)
        self.connect(self.areaToCB, SIGNAL('currentIndexChanged(int)'), self.convertArea)
        self.connect(self.areaFromLE, SIGNAL('editingFinished()'), self.convertArea)

        # 3: Volume Calculation
        # volumeUnit = [Bcf, Tcf, acre.ft, MMcf, ft3, yard3, galUK, cm3, L, dm3, bbl, m3, quart, in3, MCF, galUS]
        self.volumeDict = {'ml': 1, 'm3': 1, 'ft3': 35.3147, 'litre': 1000}
        #sortedVolumeUnits = sorted(self.areaUnits.keys())
        volumeUnits = list(self.volumeDict.keys())

        self.volumeFromCB = self.ui.findChild(QComboBox, 'volumeFromComboBox')
        self.volumeFromCB.addItems(volumeUnits)
        self.volumeToCB = self.ui.findChild(QComboBox, 'volumeToComboBox')
        self.volumeToCB.addItems(volumeUnits)
        self.volumeFromLE = self.ui.findChild(QLineEdit, 'volumeFromLineEdit')
        self.volumeToLE = self.ui.findChild(QLineEdit, 'volumeToLineEdit')

        self.connect(self.volumeFromCB, SIGNAL('currentIndexChanged(int)'), self.convertVolume)
        self.connect(self.volumeToCB, SIGNAL('currentIndexChanged(int)'), self.convertVolume)
        self.connect(self.volumeFromLE, SIGNAL('editingFinished()'), self.convertVolume)

        # 4:  Weight Conversion
        # weightUnit = [tonUS, ozm, g, lbm, gr, quintal, tonUK, kg, t, 1000lbm]
        self.weightDict = {'kg': 1, 'pound': 2.20462, 'gr': 1000, 'ounce': 35.274}
        weightUnits = list(self.weightDict.keys())

        self.weightFromCB = self.ui.findChild(QComboBox, 'weightFromComboBox')
        self.weightFromCB.addItems(weightUnits)
        self.weightToCB = self.ui.findChild(QComboBox, 'weightToComboBox')
        self.weightToCB.addItems(weightUnits)
        self.weightFromLE = self.ui.findChild(QLineEdit, 'weightFromLineEdit')
        self.weightToLE = self.ui.findChild(QLineEdit, 'weightToLineEdit')

        self.connect(self.weightFromCB, SIGNAL('currentIndexChanged(int)'), self.convertWeight)
        self.connect(self.weightToCB, SIGNAL('currentIndexChanged(int)'), self.convertWeight)
        self.connect(self.weightFromLE, SIGNAL('editingFinished()'), self.convertWeight)

        # 5: Pressure Conversion
        # pressureUnit = [Pa, psi, atm, kPa, in Hg, MPA, ft H2O, kg/cm2, mm Hg, bar]

        # 6: Temperature Conversion
        self.tempDict = {'Celsius': [1, 0], 'Farenheit': [9 / 5, 32], 'Kelvin': [1, 273.15], 'Rankine': [9 / 5, 490]}
        tempUnits = list(self.tempDict.keys())

        self.tempFromCB = self.ui.findChild(QComboBox, 'tempFromComboBox')
        self.tempFromCB.addItems(tempUnits)
        self.tempToCB = self.ui.findChild(QComboBox, 'tempToComboBox')
        self.tempToCB.addItems(tempUnits)
        self.tempFromLE = self.ui.findChild(QLineEdit, 'tempFromLineEdit')
        self.tempToLE = self.ui.findChild(QLineEdit, 'tempToLineEdit')

        self.connect(self.tempFromCB, SIGNAL('currentIndexChanged(int)'), self.convertTemp)
        self.connect(self.tempToCB, SIGNAL('currentIndexChanged(int)'), self.convertTemp)
        self.connect(self.tempFromLE, SIGNAL('editingFinished()'), self.convertTemp)

        # 7: Density Conversion
        # densityUnit = [g/cm3, ppm, kg/m3, lbm/galUS, kg/L, S.G, gr/galUS, g/L, lbm/ft3]
        # 8: Velocity Conversion
        # velocityUnit = [Km/h, ft/min, ft/s, ft/h, mi/h, m/min, m/s]

        # 9: Flow Rate Conversion
        # flowRateUnit = [ft3/min, ft3/h, cm3/s, bbl/min, in3/s, L/min, m3/min, galUS/min, m3/s, bbl/h, L/h, m3/h, bbl/d, ft3/s]

        # 10: Power Conversion
        # powerUnit = [ho (metric), hp, W, kW, Btu/min, lbf.ft/s, lbf.ft/min]

        # ----------------------Duplex Pump----------------------
        self.dLinerDiameterCB = self.ui.findChild(QComboBox, 'dLinerDiameterCBox')
        dLinerDiameterUnit = ['in', 'cm']
        self.dLinerDiameterCB.addItems(dLinerDiameterUnit)
        self.dLinerDiameterLE = self.ui.findChild(QLineEdit, 'dLinerDiameterLineEdit')

        self.dRodDiameterCB = self.ui.findChild(QComboBox, 'dRodDiameterCBox')
        dRodDiameterUnit = ['in', 'cm']
        self.dRodDiameterCB.addItems(dRodDiameterUnit)
        self.dRodDiameterLE = self.ui.findChild(QLineEdit, 'dRodDiameterLineEdit')


        self.dStrokeLengthCB = self.ui.findChild(QComboBox, 'dStrokeLengthCBox')
        dStrokeLengthCBoxUnit = ['in', 'cm']
        self.dStrokeLengthCB.addItems(dStrokeLengthCBoxUnit)
        self.dStrokeLengthLE = self.ui.findChild(QLineEdit, 'dStrokeLengthLineEdit')


        self.dEfficiencyCB = self.ui.findChild(QComboBox, 'dEfficiencyCBox')
        dEfficiencyCBoxUnit = ['%']
        self.dEfficiencyCB.addItems(dEfficiencyCBoxUnit)
        self.dEfficiencyLE = self.ui.findChild(QLineEdit, 'dEfficiencyLineEdit')


        self.dPumpRateCB = self.ui.findChild(QComboBox, 'dPumpRateCBox')
        dPumpRateCBoxUnit = ['stk/min']
        self.dPumpRateCB.addItems(dPumpRateCBoxUnit)
        self.dPumpRateLE = self.ui.findChild(QLineEdit, 'dPumpRateLineEdit')

        self.duplexPumpLabel = self.ui.findChild(QLabel, 'duplexPumpOutputLabel')
        self.duplexFlowrateLabel = self.ui.findChild(QLabel, 'duplexFlowrateOutputLabel')
        self.duplexCalculateBtn = self.ui.findChild(QPushButton, 'dCalculateButton')
        self.duplexCalculateBtn.clicked.connect(self.duplexPump)

        # ----------------------------Triplex Pump-------------------------------
        self.tLinerDiameterCB = self.ui.findChild(QComboBox, 'tLinerDiameterCBox')
        tLinerDiameterUnit = ['in', 'cm']
        self.tLinerDiameterCB.addItems(tLinerDiameterUnit)
        self.tLinerDiameterLE = self.ui.findChild(QLineEdit, 'tLinerDiameterLineEdit')

        self.tStrokeLengthCB = self.ui.findChild(QComboBox, 'tStrokeLengthCBox')
        tStrokeLengthCBoxUnit = ['in', 'cm']
        self.tStrokeLengthCB.addItems(tStrokeLengthCBoxUnit)
        self.tStrokeLengthLE = self.ui.findChild(QLineEdit, 'tStrokeLengthLineEdit')


        self.tEfficiencyCB = self.ui.findChild(QComboBox, 'tEfficiencyCBox')
        tEfficiencyCBoxUnit = ['%']
        self.tEfficiencyCB.addItems(tEfficiencyCBoxUnit)
        self.tEfficiencyLE = self.ui.findChild(QLineEdit, 'tEfficiencyLineEdit')


        self.tPumpRateCB = self.ui.findChild(QComboBox, 'tPumpRateCBox')
        tPumpRateCBoxUnit = ['stk/min']
        self.tPumpRateCB.addItems(tPumpRateCBoxUnit)
        self.tPumpRateLE = self.ui.findChild(QLineEdit, 'tPumpRateLineEdit')

        self.triplexPumpLabel = self.ui.findChild(QLabel, 'triplexPumpOutputLabel')
        self.triplexFlowrateLabel = self.ui.findChild(QLabel, 'triplexFlowrateOutputLabel')
        self.triplexCalculateBtn = self.ui.findChild(QPushButton, 'tCalculateButton')
        self.triplexCalculateBtn.clicked.connect(self.triplexPump)
    
        # ----------------------Import Data Widgets---------------------------
        self.importBtn = self.ui.findChild(QPushButton, 'importDataButton')
        #self.fileName = 'new50.csv'

        self.model = QStandardItemModel(self)
        self.tableView = self.ui.findChild(QTableView, 'dataTableView')

        self.tableView.setModel(self.model)
        self.tableView.horizontalHeader().setStretchLastSection(True)
        self.importBtn.clicked.connect(self.importData)

        # ----------------------Hole Cleaning Widgets---------------------------
        self.holeIDWidget = self.ui.findChild(QLineEdit, 'holeIDLineEdit')
        self.pipeODWidget = self.ui.findChild(QLineEdit, 'pipeODLineEdit')
        self.mudWeightWidget1 = self.ui.findChild(QLineEdit, 'mudWeightLineEdit1')
        self.plasticViscosityWidget = self.ui.findChild(QLineEdit, 'plasticViscosityLineEdit')
        self.yieldPointWidget = self.ui.findChild(QLineEdit, 'yieldPointLineEdit')
        self.flowrateWidget1 = self.ui.findChild(QLineEdit, 'flowrateLineEdit1')

        self.calculateBtn = self.ui.findChild(QPushButton, 'calculateButton')
        self.calculateBtn.clicked.connect(self.holeCleaningCalculation)

        self.CCIWidget = self.ui.findChild(QLabel, 'CCIOutputLabel')
        self.holeCleaningWidget = self.ui.findChild(QLabel, 'holeCleaningOutputLabel')

        # -----------------------Optimization Widgets-----------------------------
        self.WOBWidget = self.ui.findChild(QLineEdit, 'WOBLineEdit')
        # self.WOBInput.setValidator()
        self.RPMWidget = self.ui.findChild(QLineEdit, 'RPMLineEdit')
        self.impactForceWidget = self.ui.findChild(QLineEdit, 'impactForceLineEdit')
        self.mudWeightWidget2 = self.ui.findChild(QLineEdit, 'mudWeightLineEdit2')
        self.flowrateWidget2 = self.ui.findChild(QLineEdit, 'flowrateLineEdit2')

        self.optimizeBtn = self.ui.findChild(QPushButton, 'optimizeButton')
        self.optimizeBtn.clicked.connect(self.optimization)

        self.ROPOutputWidget = self.ui.findChild(QLabel, 'ROPOutputLabel')


    # ---------------------------------Drill Cost--------------------------------------
        self.bitCostLE = self.ui.findChild(QLineEdit, 'bitCostLineEdit')

        self.costROPLE= self.ui.findChild(QLineEdit, 'ROPLineEdit')

        self.footageDrilledLE = self.ui.findChild(QLineEdit, 'footageDrilledLineEdit')

        self.roundTripLE = self.ui.findChild(QLineEdit, 'roundTripLineEdit')

        self.rigCostLE = self.ui.findChild(QLineEdit, 'rigCostLineEdit')

        self.drilledCostOutput = self.ui.findChild(QLabel, 'drilledCostOutputLabel')
        self.calculateCostBtn = self.ui.findChild(QPushButton, 'calculateCostButton')
        self.calculateCostBtn.clicked.connect(self.drilledCost)

    
    def convertLength(self):
        inputedLength = self.lengthFromLE.text()
        length = float(inputedLength)
        from_ = self.lengthFromCB.currentText()
        to = self.lengthToCB.currentText()
        try:
            convertedLength = length * (self.lengthUnit[to] / self.lengthUnit[from_])
            return self.lengthToLE.setText(f'{amount:.3f}')
        except KeyError:
            return "Not a valid unit!"



    def convertArea(self):
        inputedArea = self.areaFromLE.text()
        area = float(inputedArea)
        from_ = self.areaFromCB.currentText()
        to = self.areaToCB.currentText()
        try:
            convertedArea = area * (self.areaUnits[to] / self.areaUnits[from_])
            return self.areaToLE.setText(f'{convertedArea:.3f}')
        except KeyError:
            return "Not a valid unit!"

    
    def convertVolume(self):
        inputedVolume = self.volumeFromLE.text()
        volume = float(inputedVolume)
        from_ = self.volumeFromCB.currentText()
        to = self.volumeToCB.currentText()
        try:
            convertedVolume = volume * (self.volumeDict[to] / self.volumeDict[from_])
            return self.volumeToLE.setText(f'{convertedVolume:.3f}')
        except KeyError:
            return "Not a valid unit!"

    def convertWeight(self):
        inputedWeight = self.weightFromLE.text()
        weight = float(inputedWeight)
        from_ = self.weightFromCB.currentText()
        to = self.weightToCB.currentText()
        try:
            convertedWeight = weight * (self.weightDict[to] / self.weightDict[from_])
            return self.weightToLE.setText(f'{convertedWeight:.3f}')
        except KeyError:
            return "Not a valid unit!"


    def convertTemp(self):
        inputedTemp = self.tempFromLE.text()
        temp = float(inputedTemp)
        from_ = self.tempFromCB.currentText()
        to = self.tempToCB.currentText()
        try:
            #convertedTemp = temp * (self.tempDict[to] / self.tempDict[from_])
            convertedTemp = temp * self.tempDict[to][0] / self.tempDict[from_][0] + self.tempDict[to][1] - (self.tempDict[to][0] /self.tempDict[from_][0]) * self.tempDict[from_][1]
            return self.tempToLE.setText(f'{convertedTemp:.3f}')
        except KeyError:
            return "Not a valid unit!"

    def duplexPump(self):
        dLinerDiameterInput = self.dLinerDiameterLE.text()
        dLinerDiameter = float(dLinerDiameterInput)

        dRodDiameterInput = self.dRodDiameterLE.text()
        dRodDiameter = float(dRodDiameterInput)

        dStrokeLengthInput = self.dStrokeLengthLE.text()
        dStrokeLength = float(dStrokeLengthInput)

        dEfficiencyInput = self.dEfficiencyLE.text()
        dEfficiency = float(dEfficiencyInput)


        dPumpRateInput = self.dPumpRateLE.text()
        dPumpRate = float(dPumpRateInput)


        if dLinerDiameter and dRodDiameter and dStrokeLength and dEfficiency and dPumpRate:
            duplexPumpOutput = 0.000162 * dStrokeLength * (2 * (dLinerDiameter**2 - dRodDiameter**2) * dEfficiency)
            duplexFlowrate = 42 * duplexPumpOutput * dPumpRate
        else:
            print(f'Invalid Input')
        self.duplexPumpLabel.setText(f'{duplexPumpOutput:.5f}')
        self.duplexFlowrateLabel.setText(f'{duplexFlowrate:.5f}')
        


    def triplexPump(self):
        tLinerDiameterInput = self.tLinerDiameterLE.text()
        tLinerDiameter = float(tLinerDiameterInput)

        tStrokeLengthInput = self.tStrokeLengthLE.text()
        tStrokeLength = float(tStrokeLengthInput)

        tEfficiencyInput = self.tEfficiencyLE.text()
        tEfficiency = float(tEfficiencyInput)


        tPumpRateInput = self.tPumpRateLE.text()
        tPumpRate = float(tPumpRateInput)


        if tLinerDiameter and tStrokeLength and tEfficiency and tPumpRate:
            triplexPumpOutput = 0.000243 * tLinerDiameter**2 * tStrokeLength * tEfficiency
            triplexFlowrate = 42 * triplexPumpOutput * tPumpRate
        else:
            print(f'Invalid Input')
        self.triplexPumpLabel.setText(f'{triplexPumpOutput:.5f}')
        self.triplexFlowrateLabel.setText(f'{triplexFlowrate:.5f}')
        
    # def loadCSV(self, fileName):
    #     with open(fileName, "r") as fileInput:
    #         for row in csv.reader(fileInput):    
    #             items = [
    #                 QStandardItem(field)
    #                 for field in row
    #             ]
    #             self.model.appendRow(items)


    # def on_pushButtonLoad_clicked(self):
    #     self.loadCSV(self.fileName)


    def importData(self):
        #print('CSV has been imported')
        global file
        file = QFileDialog.getOpenFileName(
            self, "Open CSV", (QDir.homePath()), "CSV (*.csv *.xlsx)")

        if file[0]:
            with open(file[0], "r") as fileInput:
                for row in csv.reader(fileInput):    
                    items = [
                        QStandardItem(field)
                        for field in row
                    ]
                    self.model.appendRow(items)

        return file[0]

    # --------------Hole Cleaning Calculation--------------------
    def holeCleaningCalculation(self):
        holeIDInput = self.holeIDWidget.text()
        holeID = float(holeIDInput)

        pipeODInput = self.pipeODWidget.text()
        pipeOD = float(pipeODInput)

        mudWeightInput1 = self.mudWeightWidget1.text()
        mudWeight1 = float(mudWeightInput1)

        plasticViscosityInput = self.plasticViscosityWidget.text()
        plasticViscosity = float(plasticViscosityInput)


        yieldPointInput = self.yieldPointWidget.text()
        yieldPoint = float(yieldPointInput)

        flowrateInput1 = self.flowrateWidget1.text()
        flowrate1 = float(flowrateInput1)

        if holeID and pipeOD and mudWeight1 and plasticViscosityInput and yieldPoint and flowrate1:
            # Annular Capacity (gal/ft)
            annularCapacity = ((holeID**2) - (pipeOD**2)) / 24.51
            # Annular Velocity (ft/min)
            annularVelocity = flowrate1 / annularCapacity
            # Flow behaviour index - n
            n = 3.322 * np.log10(((2 * plasticViscosity) + yieldPoint)/ (plasticViscosity + yieldPoint))
            # Power Law Constant - K (unitless)
            K = (511)**(1-n) * (plasticViscosity + yieldPoint)
            # Carrying Capacity Index - CCI (unitless)
            CCI = (K * annularVelocity * mudWeight1) / 400000
        else:
            print(f'Invalid Input')

        self.CCIWidget.setText(f'{CCI:.5f}')
        if CCI >= 1.0:
            self.holeCleaningWidget.setText('Hole Cleaning is Good')
        elif CCI <= 0.5:
            self.holeCleaningWidget.setText('Hole Cleaning is Poor')


    def optimization(self):
        WOBInput = self.WOBWidget.text()
        WOB = float(WOBInput)

        RPMInput = self.RPMWidget.text()
        RPM = float(RPMInput)

        impactForceInput = self.impactForceWidget.text()
        impactForce = float(impactForceInput)

        mudWeightInput2 = self.mudWeightWidget2.text()
        mudWeight2 = float(mudWeightInput2)

        flowrateInput2 = self.flowrateWidget2.text()
        flowrate2 = float(flowrateInput2)
        data = file[0]

        df = pd.read_csv(data)
        df.sort_values(by='ROP (ft/hr)', inplace=True)

        X = df.iloc[:, 1:6]
        y = df.iloc[:, 0]
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=0)

        #X_train = np.array([22090.29, 233, 852.52, 693.6, 13.11])
        #y_train = np.array([233, 22090.29, 660.69, 7.92, 1.6])

        #X_train = X_train.reshape(-1, 1)
        #y_train = y_train.reshape(-1, 1)

        polynomial_features = PolynomialFeatures(degree=2)
        x_poly = polynomial_features.fit_transform(X_train)
        model = LinearRegression()
        model.fit(x_poly, y_train)

        if WOB and RPM and flowrate2 and impactForce and mudWeight2:
            input = polynomial_features.fit_transform(
                [[WOB, RPM, flowrate2, impactForce, mudWeight2]])
            ROP = model.predict(input)
        else:
            print(f'Invalid Input')

        self.ROPOutputWidget.setText(f'{ROP[0]:.3f}')


    def drilledCost(self):
        bitCostLEInput = self.bitCostLE.text()
        bitCost = float(bitCostLEInput)

        ROPInput = self.costROPLE.text()
        costROP = float(ROPInput)

        rigCostPerHourInput = self.rigCostLE.text()
        rigCostPerHour = float(rigCostPerHourInput)

        footageDrilledLEInput = self.footageDrilledLE.text()
        footageDrilled = float(footageDrilledLEInput)

        roundTripInput = self.roundTripLE.text()
        roundTrip = float(roundTripInput)

        if bitCost and costROP and rigCostPerHour and footageDrilled and roundTrip:
            drilledCostPerFoot = (bitCost + rigCostPerHour * ((footageDrilled / costROP) + roundTrip)) / footageDrilled
        else:
            print(f'Invalid Input')
        self.drilledCostOutput.setText(f'{drilledCostPerFoot:.5f} $/ft')
  
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

  
if __name__ == '__main__':
    QCoreApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
    #qt_app = QtWidgets.QApplication(sys.argv)
    app = QApplication(sys.argv)
    window = MainWindow()
    #window.show()
    sys.exit(app.exec_())
