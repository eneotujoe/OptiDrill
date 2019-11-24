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

# from PySide2.QtUiTools import QUiLoader

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvas,  NavigationToolbar2QT as NavigationToolbar

# pyside2-uic -x main.ui -o mainui.py
# export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"

from matplotlib.figure import Figure
from matplotlib import style
style.use('ggplot')


import mainui

class MainWindow(QMainWindow, mainui.Ui_MainWindow):

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        self.setupUi(self)
        
        self.center()
        #----------------------Home----------------------------
        self.logoLabel.setPixmap(QPixmap('logo.png').scaled(600, 400))

        self.versionLabel.setText('Version 0.1.1 \n All Right Reserved')

        # 1: Length conversion:
        #lengthUnit = ['--select unit--', 'mm', 'cm', 'm', 'km', 'in', '1/64in', 'ft', 'yard']
        self.lengthDict = {'m': 1, 'mm': 1000, 'cm': 100, 'km': 0.001}
        lengthUnits = sorted(self.lengthDict.keys())
        # lengthUnits = list(self.lengthDict.keys())

        self.lengthFromComboBox.addItems(lengthUnits)
        self.lengthToComboBox.addItems(lengthUnits)

        self.connect(self.lengthFromComboBox, SIGNAL('currentIndexChanged(int)'), self.convertLength)
        self.connect(self.lengthToComboBox, SIGNAL('currentIndexChanged(int)'), self.convertLength)
        self.connect(self.lengthFromLineEdit, SIGNAL('editingFinished()'), self.convertLength)

        # 2: Area Conversion:
        # areaUnit = [u'ft\u00b2', u'cm\u00b2', 'acre', u'm\u00b2', 'ha', u'in\u00b2', u'cm\u00b2', u'mm\u00b2']
        self.areaUnits = {'sq.cm': 10000, 'sq. meter': 1, 'sq. feet': 10.7639, 'ha': 1e-4}
        #sortedAreaUnit = sorted(self.areaUnits.keys())
        areaUnits = list(self.areaUnits.keys())

        self.areaFromComboBox.addItems(areaUnits)
        self.areaToComboBox.addItems(areaUnits)

        self.connect(self.areaFromComboBox, SIGNAL('currentIndexChanged(int)'), self.convertArea)
        self.connect(self.areaToComboBox, SIGNAL('currentIndexChanged(int)'), self.convertArea)
        self.connect(self.areaFromLineEdit, SIGNAL('editingFinished()'), self.convertArea)

        # 3: Volume Calculation
        # volumeUnit = [Bcf, Tcf, acre.ft, MMcf, ft3, yard3, galUK, cm3, L, dm3, bbl, m3, quart, in3, MCF, galUS]
        self.volumeDict = {'ml': 1, 'm3': 1, 'ft3': 35.3147, 'litre': 1000}
        #sortedVolumeUnits = sorted(self.areaUnits.keys())
        volumeUnits = list(self.volumeDict.keys())

        self.volumeFromComboBox.addItems(volumeUnits)
        self.volumeToComboBox.addItems(volumeUnits)

        self.connect(self.volumeFromComboBox, SIGNAL('currentIndexChanged(int)'), self.convertVolume)
        self.connect(self.volumeToComboBox, SIGNAL('currentIndexChanged(int)'), self.convertVolume)
        self.connect(self.volumeFromLineEdit, SIGNAL('editingFinished()'), self.convertVolume)

        # 4:  Weight Conversion
        # weightUnit = [tonUS, ozm, g, lbm, gr, quintal, tonUK, kg, t, 1000lbm]
        self.weightDict = {'kg': 1, 'pound': 2.20462, 'gr': 1000, 'ounce': 35.274}
        weightUnits = list(self.weightDict.keys())

        self.weightFromComboBox.addItems(weightUnits)
        self.weightToComboBox.addItems(weightUnits)

        self.connect(self.weightFromComboBox, SIGNAL('currentIndexChanged(int)'), self.convertWeight)
        self.connect(self.weightToComboBox, SIGNAL('currentIndexChanged(int)'), self.convertWeight)
        self.connect(self.weightFromLineEdit, SIGNAL('editingFinished()'), self.convertWeight)

        # 5: Pressure Conversion
        # pressureUnit = [Pa, psi, atm, kPa, in Hg, MPA, ft H2O, kg/cm2, mm Hg, bar]

        # 6: Temperature Conversion
        self.tempDict = {'Celsius': [1, 0], 'Farenheit': [9 / 5, 32], 'Kelvin': [1, 273.15], 'Rankine': [9 / 5, 490]}
        tempUnits = list(self.tempDict.keys())

        self.tempFromComboBox.addItems(tempUnits)
        self.tempToComboBox.addItems(tempUnits)

        self.connect(self.tempFromComboBox, SIGNAL('currentIndexChanged(int)'), self.convertTemp)
        self.connect(self.tempToComboBox, SIGNAL('currentIndexChanged(int)'), self.convertTemp)
        self.connect(self.tempFromLineEdit, SIGNAL('editingFinished()'), self.convertTemp)

        # 7: Density Conversion
        # densityUnit = [g/cm3, ppm, kg/m3, lbm/galUS, kg/L, S.G, gr/galUS, g/L, lbm/ft3]
        # 8: Velocity Conversion
        # velocityUnit = [Km/h, ft/min, ft/s, ft/h, mi/h, m/min, m/s]

        # 9: Flow Rate Conversion
        # flowRateUnit = [ft3/min, ft3/h, cm3/s, bbl/min, in3/s, L/min, m3/min, galUS/min, m3/s, bbl/h, L/h, m3/h, bbl/d, ft3/s]

        # 10: Power Conversion
        # powerUnit = [ho (metric), hp, W, kW, Btu/min, lbf.ft/s, lbf.ft/min]

        # ----------------------Duplex Pump----------------------
        dLinerDiameterUnit = ['in', 'cm']
        self.dLinerDiameterCBox.addItems(dLinerDiameterUnit)

        dRodDiameterUnit = ['in', 'cm']
        self.dRodDiameterCBox.addItems(dRodDiameterUnit)

        dStrokeLengthCBoxUnit = ['in', 'cm']
        self.dStrokeLengthCBox.addItems(dStrokeLengthCBoxUnit)

        dEfficiencyCBoxUnit = ['%']
        self.dEfficiencyCBox.addItems(dEfficiencyCBoxUnit)

        dPumpRateCBoxUnit = ['stk/min']
        self.dPumpRateCBox.addItems(dPumpRateCBoxUnit)

        self.dCalculateButton.clicked.connect(self.duplexPump)

        # ----------------------------Triplex Pump-------------------------------
        tLinerDiameterUnit = ['in', 'cm']
        self.tLinerDiameterCBox.addItems(tLinerDiameterUnit)

        tStrokeLengthCBoxUnit = ['in', 'cm']
        self.tStrokeLengthCBox.addItems(tStrokeLengthCBoxUnit)

        tEfficiencyCBoxUnit = ['%']
        self.tEfficiencyCBox.addItems(tEfficiencyCBoxUnit)

        tPumpRateCBoxUnit = ['stk/min']
        self.tPumpRateCBox.addItems(tPumpRateCBoxUnit)

        self.tCalculateButton.clicked.connect(self.triplexPump)
    
        # ------------Import Data Widgets--------------
        #self.fileName = 'new50.csv'

        self.model = QStandardItemModel(self)

        self.dataTableView.setModel(self.model)
        self.dataTableView.horizontalHeader().setStretchLastSection(True)
        self.importDataButton.clicked.connect(self.importData)

        # -------------------Hole Cleaning Widgets-------------------
        self.calculateButton.clicked.connect(self.holeCleaningCalculation)
        # --------------------Optimization Widgets------------------
        self.optimizeButton.clicked.connect(self.optimization)

        # -------------------Drill Cost-------------------
        self.calculateCostButton.clicked.connect(self.drilledCost)

    
    def convertLength(self):
        inputedLength = self.lengthFromLineEdit.text()
        length = float(inputedLength)
        from_ = self.lengthFromComboBox.currentText()
        to = self.lengthToComboBox.currentText()
        try:
            convertedLength = length * (self.lengthDict[to] / self.lengthDict[from_])
            return self.lengthToLineEdit.setText(f'{convertedLength:.3f}')
        except KeyError:
            return "Not a valid unit!"



    def convertArea(self):
        inputedArea = self.areaFromLineEdit.text()
        area = float(inputedArea)
        from_ = self.areaFromComboBox.currentText()
        to = self.areaToComboBox.currentText()
        try:
            convertedArea = area * (self.areaUnits[to] / self.areaUnits[from_])
            return self.areaToLineEdit.setText(f'{convertedArea:.3f}')
        except KeyError:
            return "Not a valid unit!"

    
    def convertVolume(self):
        inputedVolume = self.volumeFromLineEdit.text()
        volume = float(inputedVolume)
        from_ = self.volumeFromComboBox.currentText()
        to = self.volumeToComboBox.currentText()
        try:
            convertedVolume = volume * (self.volumeDict[to] / self.volumeDict[from_])
            return self.volumeToLineEdit.setText(f'{convertedVolume:.3f}')
        except KeyError:
            return "Not a valid unit!"

    def convertWeight(self):
        inputedWeight = self.weightFromLineEdit.text()
        weight = float(inputedWeight)
        from_ = self.weightFromComboBox.currentText()
        to = self.weightToComboBox.currentText()
        try:
            convertedWeight = weight * (self.weightDict[to] / self.weightDict[from_])
            return self.weightToLineEdit.setText(f'{convertedWeight:.3f}')
        except KeyError:
            return "Not a valid unit!"


    def convertTemp(self):
        inputedTemp = self.tempFromLineEdit.text()
        temp = float(inputedTemp)
        from_ = self.tempFromComboBox.currentText()
        to = self.tempToComboBox.currentText()
        try:
            #convertedTemp = temp * (self.tempDict[to] / self.tempDict[from_])
            convertedTemp = temp * self.tempDict[to][0] / self.tempDict[from_][0] + self.tempDict[to][1] - (self.tempDict[to][0] /self.tempDict[from_][0]) * self.tempDict[from_][1]
            return self.tempToLineEdit.setText(f'{convertedTemp:.3f}')
        except KeyError:
            return "Not a valid unit!"

    def duplexPump(self):
        dLinerDiameterInput = self.dLinerDiameterLineEdit.text()
        dLinerDiameter = float(dLinerDiameterInput)

        dRodDiameterInput = self.dRodDiameterLineEdit.text()
        dRodDiameter = float(dRodDiameterInput)

        dStrokeLengthInput = self.dStrokeLengthLineEdit.text()
        dStrokeLength = float(dStrokeLengthInput)

        dEfficiencyInput = self.dEfficiencyLineEdit.text()
        dEfficiency = float(dEfficiencyInput)


        dPumpRateInput = self.dPumpRateLineEdit.text()
        dPumpRate = float(dPumpRateInput)


        if dLinerDiameter and dRodDiameter and dStrokeLength and dEfficiency and dPumpRate:
            duplexPumpOutput = 0.000162 * dStrokeLength * (2 * (dLinerDiameter**2 - dRodDiameter**2) * dEfficiency)
            duplexFlowrate = 42 * duplexPumpOutput * dPumpRate
        else:
            print(f'Invalid Input')
        self.duplexPumpOutputLabel.setText(f'{duplexPumpOutput:.5f}')
        self.duplexFlowrateOutputLabel.setText(f'{duplexFlowrate:.5f}')
        


    def triplexPump(self):
        tLinerDiameterInput = self.tLinerDiameterLineEdit.text()
        tLinerDiameter = float(tLinerDiameterInput)

        tStrokeLengthInput = self.tStrokeLengthLineEdit.text()
        tStrokeLength = float(tStrokeLengthInput)

        tEfficiencyInput = self.tEfficiencyLineEdit.text()
        tEfficiency = float(tEfficiencyInput)


        tPumpRateInput = self.tPumpRateLineEdit.text()
        tPumpRate = float(tPumpRateInput)


        if tLinerDiameter and tStrokeLength and tEfficiency and tPumpRate:
            triplexPumpOutput = 0.000243 * tLinerDiameter**2 * tStrokeLength * tEfficiency
            triplexFlowrate = 42 * triplexPumpOutput * tPumpRate
        else:
            print(f'Invalid Input')
        self.triplexPumpOutputLabel.setText(f'{triplexPumpOutput:.5f}')
        self.triplexFlowrateOutputLabel.setText(f'{triplexFlowrate:.5f}')
        
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
        holeIDInput = self.holeIDLineEdit.text()
        holeID = float(holeIDInput)

        pipeODInput = self.pipeODLineEdit.text()
        pipeOD = float(pipeODInput)

        mudWeightInput1 = self.mudWeightLineEdit1.text()
        mudWeight1 = float(mudWeightInput1)

        plasticViscosityInput = self.plasticViscosityLineEdit.text()
        plasticViscosity = float(plasticViscosityInput)


        yieldPointInput = self.yieldPointLineEdit.text()
        yieldPoint = float(yieldPointInput)

        flowrateInput1 = self.flowrateLineEdit1.text()
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

        self.CCIOutputLabel.setText(f'{CCI:.5f}')
        if CCI >= 1.0:
            self.holeCleaningOutputLabel.setText('Hole Cleaning is Good')
        elif CCI <= 0.5:
            self.holeCleaningOutputLabel.setText('Hole Cleaning is Poor')


    def optimization(self):
        WOBInput = self.WOBLineEdit.text()
        WOB = float(WOBInput)

        RPMInput = self.RPMLineEdit.text()
        RPM = float(RPMInput)

        impactForceInput = self.impactForceLineEdit.text()
        impactForce = float(impactForceInput)

        mudWeightInput2 = self.mudWeightLineEdit2.text()
        mudWeight2 = float(mudWeightInput2)

        flowrateInput2 = self.flowrateLineEdit2.text()
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

        self.ROPOutputLabel.setText(f'{ROP[0]:.3f}')


    def drilledCost(self):
        bitCostLEInput = self.bitCostLineEdit.text()
        bitCost = float(bitCostLEInput)

        ROPInput = self.ROPLineEdit.text()
        costROP = float(ROPInput)

        rigCostPerHourInput = self.rigCostLineEdit.text()
        rigCostPerHour = float(rigCostPerHourInput)

        footageDrilledLEInput = self.footageDrilledLineEdit.text()
        footageDrilled = float(footageDrilledLEInput)

        roundTripInput = self.roundTripLineEdit.text()
        roundTrip = float(roundTripInput)

        if bitCost and costROP and rigCostPerHour and footageDrilled and roundTrip:
            drilledCostPerFoot = (bitCost + rigCostPerHour * ((footageDrilled / costROP) + roundTrip)) / footageDrilled
        else:
            print(f'Invalid Input')
        self.drilledCostOutputLabel.setText(f'{drilledCostPerFoot:.5f} $/ft')
  
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
    window.show()
    sys.exit(app.exec_())
