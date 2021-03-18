"""
Module to assign peaks
Responds to current.peaks

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2021"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2021-03-18 13:10:44 +0000 (Thu, March 18, 2021) $"
__version__ = "$Revision: 3.0.3 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import typing
import numpy as np
from functools import partial
from collections import OrderedDict
from PyQt5 import QtGui, QtWidgets, QtCore
from ccpn.core.NmrAtom import NmrAtom
from ccpn.core.NmrResidue import NmrResidue, _getNmrResidue
from ccpn.core.Peak import Peak
from ccpn.core.lib import CcpnSorting
from ccpn.core.lib.AssignmentLib import nmrAtomsForPeaks, peaksAreOnLine, sameAxisCodes
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.ButtonList import ButtonList, Button
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.Frame import Frame, ScrollableFrame
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.HLine import HLine
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.GuiTable import GuiTable
from ccpn.ui.gui.widgets.Column import ColumnClass
from ccpn.ui.gui.widgets.MessageDialog import showYesNoWarning
from ccpn.ui.gui.widgets.Splitter import Splitter
from ccpn.ui.gui.widgets.Icon import Icon
from ccpn.ui.gui.guiSettings import getColours, DIVIDER, LABEL_WARNINGFOREGROUND
from ccpn.util.Logging import getLogger
from ccpn.util.Common import greekKey, _truncateText, getIsotopeListFromCode
from ccpn.ui.gui.widgets.MessageDialog import showWarning, showYesNo
from ccpnmodel.ccpncore.lib.Constants import defaultNmrChainCode
from ccpn.core.lib.Notifiers import Notifier
from ccpn.ui.gui.widgets.ScrollArea import ScrollArea
from ccpn.ui.gui.widgets.Widget import Widget
from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget
from ccpn.ui.gui.widgets.Font import getFontHeight, TABLEFONT
from ccpn.core.lib.ContextManagers import undoBlock
from ccpn.ui.gui.guiSettings import BORDERNOFOCUS_COLOUR


logger = getLogger()

allowedResidueTypes = [('', '', ''),
                       ('Alanine', 'ALA', 'A'),
                       ('Arginine', 'ARG', 'R'),
                       ('Asparagine', 'ASN', 'N'),
                       ('Aspartic acid', 'ASP', 'D'),
                       ('ASP/ASN ambiguous', 'ASX', 'B'),
                       ('Cysteine', 'CYS', 'C'),
                       ('Glutamine', 'GLN', 'Q'),
                       ('Glutamic acid', 'GLU', 'E'),
                       ('GLU/GLN ambiguous', 'GLX', 'Z'),
                       ('Glycine', 'GLY', 'G'),
                       ('Histidine', 'HIS', 'H'),
                       ('Isoleucine', 'ILE', 'I'),
                       ('Leucine', 'LEU', 'L'),
                       ('Lysine', 'LYS', 'K'),
                       ('Methionine', 'MET', 'M'),
                       ('Phenylalanine', 'PHE', 'F'),
                       ('Proline', 'PRO', 'P'),
                       ('Serine', 'SER', 'S'),
                       ('Threonine', 'THR', 'T'),
                       ('Tryptophan', 'TRP', 'W'),
                       ('Tyrosine', 'TYR', 'Y'),
                       ('Unknown', 'UNK', ''),
                       ('Valine', 'VAL', 'V')]

MSG = '<Not-defined. Select any to start>'

ROWDEFAULT = 300
ROWSIZES = {7 : 3000,
            4 : 2500,
            0 : 2000,
            -1: 1200,
            -2: ROWDEFAULT,
            }


class PeakAssigner(CcpnModule):
    """Module for assignment of nmrAtoms to the different axes of a peak.
    Module responds to current.peak
    """

    # override in specific module implementations
    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    settingsPosition = 'top'
    className = 'PeakAssigner'


    class _emptyObject():
        def __init__(self):
            pass


    def __init__(self, mainWindow, name="Peak Assigner"):

        CcpnModule.__init__(self, mainWindow=mainWindow, name=name)

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current

        # settings
        row = 0
        self.doubleToleranceCheckbox = CheckBox(self.settingsWidget, checked=False,
                                                callback=self._updateInterface,
                                                grid=(row, 1))
        doubleToleranceCheckboxLabel = Label(self.settingsWidget, text="Double Tolerances ", grid=(row, 0))

        row += 1
        self.intraCheckbox = CheckBox(self.settingsWidget, checked=False,
                                      callback=self._updateInterface,
                                      grid=(row, 1))
        intraCheckboxLabel = Label(self.settingsWidget, text="Only Intra-residual ", grid=(row, 0))

        row += 1
        self.multiCheckbox = CheckBox(self.settingsWidget, checked=True,
                                      callback=self._updateInterface,
                                      grid=(row, 1))
        multiCheckboxLabel = Label(self.settingsWidget, text="Allow Multiple Peaks ", grid=(row, 0))

        row += 1
        self.allChainCheckBoxLabel = CheckBox(self.settingsWidget, checked=False,
                                              callback=self._updateInterface,
                                              grid=(row, 1))
        allChainCheckBoxLabel = Label(self.settingsWidget, "Peak Selection from Table", grid=(row, 0))

        self._height = getFontHeight()
        self._tableHeight = getFontHeight(name=TABLEFONT)
        self.settingsWidget.setContentsMargins(5, 5, 5, 5)
        self.settingsWidget.getLayout().setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        self.settingsWidget.setScrollBarPolicies(scrollBarPolicies=('asNeeded', 'never'))

        # setup a scroll area
        self.axisFrameWidget = ScrollableFrame(parent=self.mainWidget, showBorder=False, setLayout=True,
                                               acceptDrops=True, grid=(0, 0), gridSpan=(1, 1), spacing=(5, 5))
        self._axisFrameScrollArea = self.axisFrameWidget._scrollArea

        row = 0
        self.spacer = Spacer(self.axisFrameWidget, 5, 5,
                             QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed,
                             grid=(row, 0), gridSpan=(1, 1))

        row += 1
        # put a label and axisFrame into the axisFrameWidget, this will be wrapped in scroll bars
        self.peakLabel = Label(parent=self.axisFrameWidget, setLayout=True, spacing=(0, 0),
                               text='Current Peak: ' + MSG, bold=True,
                               grid=(row, 0),  #margins=(1, 3, 1, 3),
                               hAlign='left', vAlign='t',
                               hPolicy='ignored', vPolicy='fixed'
                               )

        row += 1
        self.axisFrame = Splitter(self, grid=(row, 0), horizontal=False)
        self.axisFrameWidget.getLayout().addWidget(self.axisFrame, row, 0, 1, 1)  # MUST be added like this

        self.axisTables = []
        self.axisDivergeLabels = []
        self.NDims = 0
        self.currentAtoms = None

        # respond to peaks
        self._registerNotifiers()

        self.closeModule = self._closeModule

        self._updateInterface()

    def _registerNotifiers(self):
        # without a tableSelection specified in the table callback, this nmrAtom callback is needed
        # to update the table
        self.setNotifier(self.current, [Notifier.CURRENT],
                         targetName=Peak._pluralLinkName,
                         callback=self._updateInterface,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.DELETE, Notifier.CREATE],
                         targetName=Peak.__name__,
                         callback=self._updateInterface,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.CHANGE, Notifier.RENAME, Notifier.CREATE],
                         targetName=NmrAtom.__name__,
                         callback=self._updateNmrAtom,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.CHANGE],
                         targetName=Peak.__name__,
                         callback=self._updateNmrResidue,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.DELETE, Notifier.CREATE],
                         targetName=NmrResidue.__name__,
                         callback=self._updateNmrResidue,
                         onceOnly=True)

    def _updateNmrAtom(self, data):
        self._updateInterface(action=data[Notifier.TRIGGER])

    def _updateNmrResidue(self, data):
        self._updateInterface(action=data[Notifier.TRIGGER])

    def _updateInterface(self, peaks: typing.List[Peak] = None,
                         enableDeleteButton=False, enableDeassignButton=False, enableAssignButton=False,
                         action=None):
        """Updates the whole module, including recalculation
           of which nmrAtoms fit to the peaks.
        """
        # self._emptyAllTablesAndLists()
        if not self.current.peaks or not self._peaksAreCompatible():
            self.axisFrame.hide()
            self.peakLabel.setText('Current Peak: ' + MSG)
        else:

            Ndimensions = len(self.current.peak.position)
            # _sizes = [1000] * Ndimensions

            if Ndimensions > self.NDims:  # len(self.axisTables):
                for addNew in range(len(self.axisTables), Ndimensions):
                    # add a new axis item to the end of the list
                    _frame = Frame(self, setLayout=True, hAlign='l', vAlign='t')
                    _newAxis = AxisAssignmentObject(self, index=addNew,
                                                    parent=_frame,
                                                    mainWindow=self.mainWindow,
                                                    grid=(0, 0), gridSpan=(1, 1))
                    self.axisTables.append(_newAxis)

                    # make a small label that appears when there is nothing to display
                    self.tempFrame = Frame(_frame, setLayout=True, grid=(1, 0))
                    self.tempDivider = None  #HLine(self.tempFrame, grid=(0, 0), gridSpan=(1, 3), colour=getColours()[DIVIDER], height=15)
                    self.tempLabel = Label(self.tempFrame, text='', grid=(1, 0), hPolicy='ignored', textColour=getColours()[LABEL_WARNINGFOREGROUND], )
                    self.tempLabel.setFixedHeight(self._height * 3)

                    self.axisDivergeLabels.append([self.tempFrame, self.tempDivider, self.tempLabel])

                    self.axisFrame.addWidget(_frame)

                for showNew in range(self.NDims, Ndimensions):
                    self.axisTables[showNew].setVisible(True)
                    self.axisDivergeLabels[showNew][0].hide()
                    self.axisFrame.widget(showNew).setVisible(True)

            elif Ndimensions < len(self.axisTables):
                for delOld in range(Ndimensions, len(self.axisTables)):
                    # self.axisTables[delOld].setVisible(False)
                    # self.axisDivergeLabels[delOld][0].hide()
                    # _sizes.append(1)
                    self.axisFrame.widget(delOld).setVisible(False)

            self.NDims = Ndimensions

            # and enable the frame
            self.axisFrame.show()

            peaksIds = ' , '.join([str(pp.id) for pp in self.current.peaks])
            if len(self.current.peaks) < 2:
                self.peakLabel.setText('Current Peak: %s' % self.current.peak.id)
            else:
                self.peakLabel.setText('Current Peaks: %s' % _truncateText(peaksIds, maxWords=6))
                self.peakLabel.setToolTip(peaksIds)

            _sizes = self._updateNewTable(enableDeleteButton=enableDeleteButton,
                                          enableDeassignButton=enableDeassignButton,
                                          enableAssignButton=enableAssignButton,
                                          action=action)

            if all(val == ROWDEFAULT for val in _sizes):
                self.axisFrame.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
            else:
                self.axisFrame.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Expanding)

            self.axisFrame.setSizes(_sizes)

    def _updateNewTable(self, enableDeleteButton=False,
                        enableDeassignButton=False,
                        enableAssignButton=False,
                        action=None):
        """
        update Assigned and alternatives tables showing which nmrAtoms
        are assigned to which peak dimensions. If multiple
        peaks are selected, only the assignment that they
        have in common are shown. Maybe this should be all
        assignments. You can see that at the peak annotation
        though.
        """
        peaks = self.current.peaks
        doubleTolerance = self.doubleToleranceCheckbox.isChecked()
        intraResidual = self.intraCheckbox.isChecked()
        validNmrAtoms = [nmrAtom for nmrAtom in self.project.nmrAtoms if not (nmrAtom.nmrResidue.isDeleted or nmrAtom.nmrResidue._flaggedForDelete)]
        nmrAtomsForTables = nmrAtomsForPeaks(peaks, validNmrAtoms,
                                             doubleTolerance=doubleTolerance,
                                             intraResidual=intraResidual)

        Ndimensions = len(nmrAtomsForTables)
        self.currentList = []

        self._tables = [self._emptyObject()] * Ndimensions

        _sizes = []
        for dim, nmrAtoms in zip(range(Ndimensions),
                                 nmrAtomsForTables):
            self.axisTables[dim].show()
            self.axisDivergeLabels[dim][0].hide()

            ll = [set(peak.dimensionNmrAtoms[dim]) for peak in self.current.peaks]
            self.nmrAtoms = list(sorted(set.intersection(*ll)))  # was intersection
            self.nmrAtoms = [nmrAtom for nmrAtom in self.nmrAtoms if not (nmrAtom.nmrResidue.isDeleted or nmrAtom.nmrResidue._flaggedForDelete)]

            self.currentList.append([str(a.pid) for a in self.nmrAtoms])  # ejb - keep another list
            self.axisTables[dim].setAssignedTable(self.nmrAtoms)
            _rows = self.axisTables[dim].tables[0].rowCount()
            # _rows = self.axisTables[dim].tables[0].sizeHint().height()

            nmrAtomsForTables[dim] = [nmr for nmr in nmrAtomsForTables[dim] if nmr not in self.nmrAtoms]

            if peaksAreOnLine(peaks, dim):
                self.axisTables[dim].setAlternativesTable(nmrAtomsForTables[dim])
                _rows = max(_rows, self.axisTables[dim].tables[1].rowCount())
                for k, val in ROWSIZES.items():
                    if _rows > k:
                        _sizes.append(val)
                        break
                else:
                    _sizes.append(ROWDEFAULT)
            else:
                self.axisTables[dim].setAlternativesTable(None)
                _sizes.append(ROWDEFAULT)

                # hide as this is not a valid table
                if not self.nmrAtoms:
                    self.axisTables[dim].setVisible(False)
                    self.axisDivergeLabels[dim][0].show()

            positions = [peak.position[dim] for peak in self.current.peaks]
            avgPos = round(sum(positions) / len(positions), 3)
            axisCode = self.current.peak.peakList.spectrum.axisCodes[dim]
            text = '%s: %.3f' % (axisCode, avgPos)
            self.axisTables[dim].axisLabel.setText(text)
            self.axisDivergeLabels[dim][2].setText(axisCode + ': peaks diverge')

            # check whether the buttons can be enabled/disabled
            currentNmrAtomSelected = (self.axisTables[dim].chainPulldown.currentText(),
                                      self.axisTables[dim].seqCodePulldown.currentText(),
                                      self.axisTables[dim].resTypePulldown.currentText(),
                                      self.axisTables[dim].atomTypePulldown.currentText())

            enable = False
            for nmrAtom in self.nmrAtoms:
                nmrChain = str(nmrAtom.nmrResidue.nmrChain.id)
                sequenceCode = str(nmrAtom.nmrResidue.sequenceCode)
                residueType = str(nmrAtom.nmrResidue.residueType)
                atomType = str(nmrAtom.name)

                item = (nmrChain, sequenceCode, residueType, atomType)
                enable = enable or (False not in self.axisTables[dim]._atomCompare(item, currentNmrAtomSelected))

            self.axisTables[dim].buttonList.setButtonEnabled('Deassign', enable)

            enable = False
            for nmrAtom in nmrAtomsForTables[dim]:
                nmrChain = str(nmrAtom.nmrResidue.nmrChain.id)
                sequenceCode = str(nmrAtom.nmrResidue.sequenceCode)
                residueType = str(nmrAtom.nmrResidue.residueType)
                atomType = str(nmrAtom.name)

                item = (nmrChain, sequenceCode, residueType, atomType)
                enable = enable or (False not in self.axisTables[dim]._atomCompare(item, currentNmrAtomSelected))

            self.axisTables[dim].buttonList.setButtonEnabled('Assign', enable)

        return _sizes

    def _getDeltaShift(self, nmrAtom: NmrAtom, dim: int) -> typing.Union[float, str]:
        """
        Calculation of delta shift to add to the table.
        """
        if (not self.current.peaks) or nmrAtom is NOL:
            return ''

        deltas = []
        for peak in self.current.peaks:
            shiftList = peak.peakList.spectrum.chemicalShiftList
            if shiftList:
                shift = shiftList.getChemicalShift(nmrAtom.id)
                if shift:
                    position = peak.position[dim]
                    deltas.append(abs(shift.value - position))
        # average = sum(deltas)/len(deltas) #Bug: ZERO DIVISION!

        if len(deltas) > 0:
            return float(np.mean(deltas))  #'%6.3f' % np.mean(deltas) - handled by table
        else:
            return ''

    def _getShift(self, nmrAtom: NmrAtom) -> typing.Union[float, str]:
        """
        Calculation of chemical shift value to add to the table.
        """
        if (not self.current.peaks) or nmrAtom is NOL:
            return ''

        for peak in self.current.peaks:
            shiftList = peak.peakList.spectrum.chemicalShiftList
            if shiftList:
                shift = shiftList.getChemicalShift(nmrAtom.id)
                if shift:
                    return shift.value  # '%8.3f' % shift.value

    def _peaksAreCompatible(self) -> bool:
        """
        If multiple peaks are selected, a check is performed
        to determine whether assignment of corresponding
        dimensions of a peak allowed.
        """
        if len(self.current.peaks) == 1:
            return True
        if not self.multiCheckbox.isChecked():
            self.project._logger.warning("Multiple peaks selected, not allowed.")
            return False
        dimensionalities = set([len(peak.position) for peak in self.current.peaks])
        if len(dimensionalities) > 1:
            self.project._logger.warning('Not all peaks have the same number of dimensions.')
            return False
        for dim in range(len(self.current.peak.position)):
            if not sameAxisCodes(self.current.peaks, dim):
                self.project._logger.warning('''The combination of axiscodes is different for multiple
                 selected peaks.''')
                return False
        return True

    def _emptyAllTablesAndLists(self):
        """
        Quick erase of all present information in ListWidgets and ObjectTables.
        """
        self.peakLabel.setText('Current Peak: ' + MSG)
        for label in self.labels:
            label.setText('')
        for objectTable in self.objectTables:
            objectTable.setObjects([])
        for listWidget in self.listWidgets:
            listWidget.clear()

    def _updatePulldownLists(self, dim: int, row: int = None, col: int = None, obj: object = None):
        objectTable = self.objectTables[dim]
        nmrAtom = objectTable.getCurrentObject()
        self._updateAssignmentWidget(dim, nmrAtom)

    def _closeModule(self):
        """
        CCPN-INTERNAL: used to close the module
        """
        # self._unRegisterNotifiers()
        for axisTable in self.axisTables:
            axisTable._close()
        self.axisTables = None
        super()._closeModule()

    def close(self):
        """
        Close the table from the commandline
        """
        self._closeModule()


class NotOnLine(object):
    """
    Small 'fake' object to get a message the user in the assignment
    Table that a specific dimension can not be assigned in one go
    since the frequencies of the peaks in this dimension are not on
    one line (i.e. the C frequencies of the CA and CB in a strip for
    instance).
    """

    def __init__(self):
        self.pid = 'Multiple selected peaks not on line.'
        self.id = 'Multiple selected peaks not on line.'


NOL = NotOnLine()


class AxisAssignmentObject(Frame):
    """
    Create a new frame for displaying information in 1 axis of peakassigner
    """

    def __init__(self, parentModule, index=None, parent=None, mainWindow=None, grid=None, gridSpan=None):
        super(AxisAssignmentObject, self).__init__(parent=parent,
                                                   setLayout=True,
                                                   spacing=(5, 0), grid=grid, gridSpan=gridSpan
                                                   )

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current
        self.currentAtoms = None
        self._clickedNmrAtom = None

        self.splitter = Splitter(self)
        self.getLayout().addWidget(self.splitter, 0, 0)
        self._assignmentsFrame = Frame(self.splitter, setLayout=True)
        self._alternativesFrame = Frame(self.splitter, setLayout=True)
        self.splitter.setSizes([1000, 1000])
        self.splitter.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        self._assignmentsFrame.getLayout().setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        self._alternativesFrame.getLayout().setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)

        # self.divider = HLine(self, grid=(0, 0), colour=getColours()[DIVIDER], height=15)
        self._alternativesLabel = Label(self._alternativesFrame, 'Alternatives', hAlign='l', grid=(0, 0))

        row = 0
        self.axisLabel = Label(self._assignmentsFrame, 'Axis', hAlign='l', grid=(row, 0), bold=True)

        # add two tables - left is current assignments, right is alternatives
        row += 1
        self.tables = [GuiTable(parent=self._assignmentsFrame,
                                mainWindow=mainWindow,
                                dataFrameObject=None,
                                setLayout=True,
                                autoResize=True, multiSelect=False,
                                actionCallback=partial(self._assignDeassignNmrAtom, 0),
                                selectionCallback=partial(self._updatePulldownLists, 0),
                                grid=(1, 0), gridSpan=(1, 1),
                                stretchLastSection=True,
                                enableSearch=False,
                                acceptDrops=True),

                       GuiTable(parent=self._alternativesFrame,
                                mainWindow=mainWindow,
                                dataFrameObject=None,
                                setLayout=True,
                                autoResize=True, multiSelect=False,
                                actionCallback=partial(self._assignDeassignNmrAtom, 1),
                                selectionCallback=partial(self._updatePulldownLists, 1),
                                grid=(1, 0), gridSpan=(1, 1),
                                stretchLastSection=True,
                                enableSearch=False,
                                acceptDrops=True)
                       ]

        # set up notifiers to changes to peaks, nmrAtoms and assignments
        self.tables[0].setTableNotifiers(tableClass=Peak,
                                         rowClass=NmrAtom,
                                         cellClassNames=None,
                                         tableName='assignedPeaks', rowName='nmrAtom',
                                         changeFunc=parentModule._updateInterface,
                                         className='peakLists',
                                         updateFunc=parentModule._updateInterface,
                                         tableSelection=None,
                                         pullDownWidget=None,
                                         callBackClass=NmrAtom,
                                         moduleParent=self)
        self.tables[1].setTableNotifiers(tableClass=Peak,
                                         rowClass=NmrAtom,
                                         cellClassNames=None,
                                         tableName='assignedPeaks', rowName='nmrAtom',
                                         changeFunc=parentModule._updateInterface,
                                         className='peakLists',
                                         updateFunc=parentModule._updateInterface,
                                         tableSelection=None,
                                         pullDownWidget=None,
                                         callBackClass=NmrAtom,
                                         moduleParent=self.tables)  # just to give a unique id

        self._assignmentsFrame.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self._alternativesFrame.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.tables[0].setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.tables[1].setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        self._bottomFrame = Frame(self, setLayout=True, showBorder=False, grid=(1, 0), gridSpan=(1, 1), margins=(3, 3, 3, 3), vPolicy='minimum', hPolicy='minimum', )

        # add pulldowns for editing new assignment
        bRow = 0
        # self._clickedFrame = Frame(self._bottomFrame, setLayout=True, grid=(bRow, 0), gridSpan=(1, 2),
        #                            vPolicy='fixed', hPolicy='minimum', hAlign='l')
        # self._clickedLabel = Label(self._clickedFrame, text='Current NmrAtom: ', grid=(0, 0))
        # self._clickedClear = Button(self._clickedFrame, grid=(0, 1),
        #                             callback=self._clearClicked, icon=Icon('icons/reset'),
        #                             vPolicy='ignored', hPolicy='Fixed', hAlign='l', vAlign='t')
        # self._clickedClear.setFlat(True)
        # self._clickedClear.setVisible(False)
        #
        # bRow += 1
        self.pulldownFrame = Frame(parent=self._bottomFrame, setLayout=True,
                                   showBorder=False, fShape='noFrame',
                                   vAlign='top', hAlign='l',
                                   vPolicy='fixed', hPolicy='minimum',
                                   # grid=(row, 0), gridSpan=(1, 1)
                                   grid=(bRow, 0), gridSpan=(1, 2)
                                   )

        self.chainPulldown = self._createChainPulldown(parent=self.pulldownFrame,
                                                       grid=(1, 0), gridSpan=(1, 1),
                                                       tipText='Chain code')
        self.seqCodePulldown = self._createPulldown(parent=self.pulldownFrame,
                                                    grid=(1, 2), gridSpan=(1, 1),
                                                    tipText='Sequence code')
        self.resTypePulldown = self._createPulldown(parent=self.pulldownFrame,
                                                    grid=(1, 4), gridSpan=(1, 1),
                                                    tipText='Residue type')
        self.atomTypePulldown = self._createPulldown(parent=self.pulldownFrame,
                                                     grid=(1, 6), gridSpan=(1, 1),
                                                     tipText='Atom type')

        _width = parentModule._height * 6
        self.chainPulldown.setMinimumWidth(_width)
        self.seqCodePulldown.setMinimumWidth(_width)
        self.resTypePulldown.setMinimumWidth(_width)
        self.atomTypePulldown.setMinimumWidth(_width)

        # add a buttonlist
        bRow += 1
        self.buttonList = ButtonList(parent=self._bottomFrame, texts=['New', 'Delete', 'Deassign', 'Assign', 'Rename'],
                                     callbacks=[partial(self._createNewNmrAtom, index),
                                                partial(self._deleteNmrAtom, index),
                                                partial(self._deassignNmrAtom, index),
                                                partial(self._assignNmrAtom, index),
                                                partial(self._reassignNmrAtom, index),
                                                ],
                                     grid=(bRow, 0), gridSpan=(1, 1),
                                     vAlign='c', hAlign='l')

        self.buttonList.setFixedHeight(parentModule._height * 1.5)
        self.buttonList.setButtonEnabled('Delete', False)
        self.buttonList.setButtonEnabled('Deassign', False)
        self.buttonList.setButtonEnabled('Assign', False)
        self.buttonList.setButtonEnabled('Rename', True)

        # self.createNew = CheckBoxCompoundWidget(
        #         self._bottomFrame,
        #         grid=(bRow, 1), vAlign='c', hAlign='left',
        #         #minimumWidths=(colwidth, 0),
        #         # fixedWidths=(None, None),
        #         orientation='right',
        #         labelText='Create New',
        #         checked=True
        #         )
        Spacer(self._bottomFrame, 5, 5, QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed,
               grid=(bRow, 2), gridSpan=(1, 1))

        # initialise axis information
        self.index = index
        self._parent = parentModule
        self.dataFrameAssigned = None
        self.dataFrameAlternatives = None
        self.lastTableSelected = None
        self.lastNmrAtomSelected = None

        # set column definitions and hidden columns for each table
        self.columnDefs = ColumnClass([('NmrAtom', lambda nmrAtom: str(nmrAtom.id), 'NmrAtom identifier', None, None),
                                       ('Pid', lambda nmrAtom: str(nmrAtom.pid), 'Pid of the nmrAtom', None, None),
                                       ('_object', lambda nmrAtom: nmrAtom, 'Object', None, None),
                                       ('Shift', lambda nmrAtom: parentModule._getShift(nmrAtom), 'Chemical shift',
                                        None, '%8.3f'),
                                       ('Delta', lambda nmrAtom: parentModule._getDeltaShift(nmrAtom, index),
                                        'Delta shift', None, '%6.3f')])
        self._hiddenColumns = [['Pid', 'Shift'], ['Pid', 'Shift']]

        self.tables[0]._hiddenColumns = ['Pid', 'Shift']
        self.tables[1]._hiddenColumns = ['Pid', 'Shift']

        self._setDefaultPulldowns()
        # self._minWidth = (self.buttonList.sizeHint() + self.createNew.sizeHint()).width()
        self._minWidth = (self.buttonList.sizeHint()).width()

    def sizeHint(self) -> QtCore.QSize:
        _size = super().sizeHint()
        _width = max(self._minWidth, self._parent.width() - 30)
        t0 = self._parent._axisFrameScrollArea.verticalScrollBar()
        if t0.isVisible():
            _width -= t0.width()
        return QtCore.QSize(_width, _size.height())

    def _close(self):
        self.tables[0]._close()
        self.tables[1]._close()
        self.tables = None

    # def _clearClicked(self, val):
    #     self._clickedNmrAtom = None
    #     self._clickedLabel.setText('Current NmrAtom: <None>')
    #     self._clickedClear.setVisible(False)

    def _assignDeassignNmrAtom(self, tableNum: int, data):
        """
        Assign/Deassign the nmrAtom that is double clicked to the
        the corresponding dimension of the selected
        peaks.
        """
        if tableNum == 0:
            # deassign from left to right
            self._deassignNmrAtom(self.index)
        elif tableNum == 1:
            # assign from right to left
            self._assignNmrAtom(self.index, action=True)

    def _updatePulldownLists(self, tableNum, data):
        self.lastTableSelected = tableNum
        obj = data[Notifier.OBJECT]
        if obj:
            self._clickedNmrAtom = obj[0]
            # self._clickedLabel.setText('Current NmrAtom: {}'.format(obj[0].pid))
            # self._clickedClear.setVisible(True)
            if tableNum == 0:
                self._updateAssignmentWidget(tableNum, obj[0])
                # self.tables[1].clearSelection()
                self.buttonList.setButtonEnabled('Delete', True)
                self.buttonList.setButtonEnabled('Deassign', True)
                self.buttonList.setButtonEnabled('Assign', False)
            elif tableNum == 1:
                self._updateAssignmentWidget(tableNum, obj[0])
                # self.tables[0].clearSelection()
                self.buttonList.setButtonEnabled('Delete', True)
                self.buttonList.setButtonEnabled('Deassign', False)
                self.buttonList.setButtonEnabled('Assign', True)

    def _createChainPulldown(self, parent=None, grid=(0, 0), gridSpan=(1, 1), tipText='') -> PulldownList:
        """Creates a PulldownList with callback, editable.
        """
        pulldownList = PulldownList(parent=parent, grid=grid, backgroundText=tipText, editable=True, gridSpan=gridSpan,
                                    tipText=tipText)
        pulldownList.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        # pulldownList.setEditable(True)
        pulldownList.lineEdit().textChanged.connect(partial(self._chainEdited, pulldownList))
        return pulldownList

    def _chainEdited(self, pulldownList):
        text = pulldownList.currentText()
        chains = [chain.id for chain in self.project.nmrChains]

        if text in chains:
            index = chains.index(text)
            thisChain = self.project.nmrChains[index]
            self._setSequenceCodes(thisChain)
            self._setResidueTypes(thisChain)
            self._setAtomNames()

        self._pulldownEdited(None)

    def _createPulldown(self, parent=None, grid=(0, 0), gridSpan=(1, 1), tipText='') -> PulldownList:
        """Creates a PulldownList with callback, editable.
        """
        pulldownList = PulldownList(parent=parent, grid=grid, backgroundText=tipText, editable=True, gridSpan=gridSpan,
                                    tipText=tipText)
        pulldownList.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        # pulldownList.setEditable(True)
        pulldownList.lineEdit().textChanged.connect(partial(self._pulldownEdited, pulldownList))
        return pulldownList

    def _createNewNmrAtom(self, dim):
        isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[dim]
        nmrAtom = self.project.fetchNmrChain(shortName=defaultNmrChainCode
                                             ).newNmrResidue().newNmrAtom(isotopeCode=isotopeCode)

        with undoBlock():
            try:

                for peak in self.current.peaks:
                    if nmrAtom not in peak.dimensionNmrAtoms[dim]:
                        # newAssignments = peak.dimensionNmrAtoms[dim] + [nmrAtom]

                        newAssignments = list(peak.dimensionNmrAtoms[dim]) + [nmrAtom]  # ejb - changed to list
                        axisCode = peak.peakList.spectrum.axisCodes[dim]
                        peak.assignDimension(axisCode, newAssignments)

                self._parent._updateInterface(enableDeleteButton=True,
                                              enableDeassignButton=True,
                                              enableAssignButton=False)

                # highlight on the table and populate the pulldowns
                self.tables[0].selectObjects([nmrAtom], setUpdatesEnabled=False)
                self.tables[1].clearSelection()
                self._updateAssignmentWidget(0, nmrAtom)

                self.lastTableSelected = 0
                self.buttonList.setButtonEnabled('Delete', True)
                self.buttonList.setButtonEnabled('Deassign', True)
                self.buttonList.setButtonEnabled('Assign', False)

            except Exception as es:
                showWarning(str(self.windowTitle()), str(es))

    def _reassignNmrAtom(self, dim: int):
        """
        Assigns dimensionNmrAtoms to peak dimension when called using Assign Button in assignment widget.
        :param dim - axis dimension of the atom:
        :param action - True if callback is action from the table:
        """
        try:
            nmrChainName = self.chainPulldown.currentText()
            seqCode = self.seqCodePulldown.currentText()
            newResType = self.resTypePulldown.currentText()
            nmrAtomName = self.atomTypePulldown.currentText()

            # get options from the pulldowns
            currentNmrAtomSelected = (nmrChainName,
                                      seqCode,
                                      newResType,
                                      nmrAtomName)
            atomCompare = self._atomCompare(self.lastNmrAtomSelected, currentNmrAtomSelected)
            nmrAtom = None

            if not self._clickedNmrAtom:
                showWarning("Rename NmrAtom", "Please select an NmrAtom from the tables")
                return

            # wrap all actions in a single undo block
            with undoBlock():

                _chainPid = 'NC:{}'.format(nmrChainName)
                # if create and not action:
                #     # get the current chain (but may create a new one)
                #     _nmrChain = self.project.fetchNmrChain(nmrChainName)
                # else:
                #     # find the existing nmrChain
                #     _nmrChain = self.project.getByPid(_chainPid)
                #     if not _nmrChain:
                #         # raise error to notify popup
                #         raise ValueError("NmrChain doesn't exists")

                _nmrChain = self.project.fetchNmrChain(nmrChainName)
                nmrResidue = _getNmrResidue(_nmrChain, seqCode, )
                nmrAtom = nmrResidue.getNmrAtom(nmrAtomName) if nmrResidue else None

                # edit existing
                if nmrResidue and self._clickedNmrAtom.nmrResidue != nmrResidue:
                    # existing different nmrResidue
                    nmrAtom = nmrResidue.getNmrAtom(nmrAtomName)
                    if nmrAtom:
                        yesNo = showYesNo('Merge NmrAtom', "Do you want to merge\n\n"
                                                           "{}   into   {}".format(self._clickedNmrAtom.id,
                                                                                   nmrAtom.id))
                        if yesNo:
                            # merge into the new nmrAtom
                            nmrAtom.mergeNmrAtoms(self._clickedNmrAtom)

                    else:
                        # assign to a new nmrAtom
                        self._clickedNmrAtom.assignTo(chainCode=nmrChainName,
                                                      sequenceCode=seqCode,
                                                      residueType=newResType,
                                                      name=nmrAtomName,
                                                      mergeToExisting=False)

                elif nmrResidue and self._clickedNmrAtom.nmrResidue == nmrResidue:
                    # rename the same nmrAtom
                    if newResType != nmrResidue.residueType:
                        nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)

                    if nmrAtomName != self._clickedNmrAtom.name:
                        nmrAtom = nmrResidue.getNmrAtom(nmrAtomName)
                        if nmrAtom:
                            raise ValueError('NmrAtom already exists {}'.format(nmrAtom))
                        self._clickedNmrAtom.rename(nmrAtomName)

                else:
                    # nmrResidue doesn't exists
                    self._clickedNmrAtom.assignTo(chainCode=nmrChainName,
                                                  sequenceCode=seqCode,
                                                  residueType=newResType,
                                                  name=nmrAtomName,
                                                  mergeToExisting=False)

            self._parent._updateInterface()

            # update the module
            self.update()

        except Exception as es:
            showWarning('Rename NmrAtom', str(es))

    def _assignNmrAtom(self, dim: int, action: bool = False, create: bool = True):
        """
        Assigns dimensionNmrAtoms to peak dimension when called using Assign Button in assignment widget.
        :param dim - axis dimension of the atom:
        :param action - True if callback is action from the table:
        """
        # FIXME Potential Bug: no error checks for dim. It can give easily an IndexError

        # return if no peaks selected
        if not self.current.peaks:
            return

        try:
            nmrChainName = self.chainPulldown.currentText()
            seqCode = self.seqCodePulldown.currentText()
            newResType = self.resTypePulldown.currentText()
            nmrAtomName = self.atomTypePulldown.currentText()

            # get options from the pulldowns
            currentNmrAtomSelected = (nmrChainName,
                                      seqCode,
                                      newResType,
                                      nmrAtomName)
            atomCompare = self._atomCompare(self.lastNmrAtomSelected, currentNmrAtomSelected)

            # create = self.createNew.isChecked()

            nmrAtom = None

            # wrap all actions in a single undo block
            with undoBlock():

                _chainPid = 'NC:{}'.format(nmrChainName)
                if create and not action:
                    # get the current chain (but may create a new one)
                    _nmrChain = self.project.fetchNmrChain(nmrChainName)
                else:
                    # find the existing nmrChain
                    _nmrChain = self.project.getByPid(_chainPid)
                    if not _nmrChain:
                        # raise error to notify popup
                        raise ValueError("NmrChain doesn't exists")

                nmrResidue = _getNmrResidue(_nmrChain, seqCode, )
                nmrAtom = nmrResidue.getNmrAtom(nmrAtomName) if nmrResidue else None

                if not action:
                    if create:
                        if nmrResidue:
                            if self._clickedNmrAtom and self._clickedNmrAtom.nmrResidue == nmrResidue and nmrResidue.residueType != newResType:
                                if len(nmrResidue.nmrAtoms) > 1:
                                    yes = showYesNoWarning('Assigning nmrAtoms',
                                                           'This will change all nmrAtoms to the residueType {}, continue?'.format(newResType))
                                    if yes:
                                        nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)
                                else:
                                    nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)

                        else:
                            # can do a residueType rename
                            nmrResidue = _nmrChain.fetchNmrResidue(seqCode, newResType)

                        nmrAtom = nmrResidue.fetchNmrAtom(nmrAtomName)

                    else:
                        pass

                        # if not self._clickedNmrAtom:
                        #     showWarning("Rename NmrAtom", "Please select an NmrAtom from the tables")
                        #     return
                        #
                        # # edit existing
                        # if nmrResidue and self._clickedNmrAtom.nmrResidue != nmrResidue:
                        #     # existing different nmrResidue
                        #     nmrAtom = nmrResidue.getNmrAtom(nmrAtomName)
                        #     if nmrAtom:
                        #         yesNo = showYesNo('Merge NmrAtom', "Do you want to merge\n\n"
                        #                                             "{}   into   {}".format(self._clickedNmrAtom.id,
                        #                                                                     nmrAtom.id))
                        #         if yesNo:
                        #             # merge into the new nmrAtom
                        #             nmrAtom.mergeNmrAtoms(self._clickedNmrAtom)
                        #
                        #     else:
                        #         # assign to a new nmrAtom
                        #         self._clickedNmrAtom.assignTo(chainCode=nmrChainName,
                        #                                       sequenceCode=seqCode,
                        #                                       residueType=newResType,
                        #                                       name=nmrAtomName,
                        #                                       mergeToExisting=False)
                        #
                        # elif nmrResidue and self._clickedNmrAtom.nmrResidue == nmrResidue:
                        #     # rename the same nmrAtom
                        #     if newResType != nmrResidue.residueType:
                        #         nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)
                        #
                        #     if nmrAtomName != self._clickedNmrAtom.name:
                        #         nmrAtom = nmrResidue.getNmrAtom(nmrAtomName)
                        #         if nmrAtom:
                        #             raise ValueError('NmrAtom already exists {}'.format(nmrAtom))
                        #         self._clickedNmrAtom.rename(nmrAtomName)
                        #
                        # else:
                        #     # nmrResidue doesn't exists
                        #     self._clickedNmrAtom.assignTo(chainCode=nmrChainName,
                        #                                   sequenceCode=seqCode,
                        #                                   residueType=newResType,
                        #                                   name=nmrAtomName,
                        #                                   mergeToExisting=False)

                try:

                    for peak in self.current.peaks:

                        dimNmrAtoms = list(peak.dimensionNmrAtoms[dim])

                        currentObject = nmrAtom
                        if nmrAtom not in dimNmrAtoms:
                            dimNmrAtoms.append(nmrAtom)

                            toAssign = dimNmrAtoms.index(currentObject)

                            dimNmrAtoms[toAssign] = nmrAtom
                            allAtoms = list(peak.dimensionNmrAtoms)
                            allAtoms[dim] = dimNmrAtoms
                            peak.dimensionNmrAtoms = allAtoms

                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

            self._parent._updateInterface()

            self.tables[0].selectObjects([nmrAtom], setUpdatesEnabled=False)

            if nmrAtom:
                self._updateAssignmentWidget(0, nmrAtom)

                self.lastTableSelected = 0
                self.buttonList.setButtonEnabled('Delete', True)
                self.buttonList.setButtonEnabled('Deassign', True)
                self.buttonList.setButtonEnabled('Assign', False)

            else:
                self._updateAssignmentWidget(0, None)

                self.lastTableSelected = 0
                self.buttonList.setButtonEnabled('Delete', False)
                self.buttonList.setButtonEnabled('Deassign', False)
                self.buttonList.setButtonEnabled('Assign', False)

            # update the module
            self.update()

        except Exception as es:
            showWarning('Assign NmrAtom', str(es))

    def _deassignNmrAtom(self, dim: int):
        """
        remove nmrAtom from peak assignment
        """

        # return if no peaks selected
        if not self.current.peaks:
            return

        try:
            currentObject = self.tables[0].getSelectedObjects()

            if currentObject:
                try:
                    with undoBlock():
                        for peak in self.current.peaks:
                            peakDimNmrAtoms = peak.dimensionNmrAtoms
                            dimNmrAtoms = list(peakDimNmrAtoms[dim])  # ejb - changed to list
                            dimNmrAtoms.remove(currentObject[0])

                            allAtoms = list(peakDimNmrAtoms)
                            allAtoms[dim] = dimNmrAtoms
                            peak.dimensionNmrAtoms = allAtoms

                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

                self._parent._updateInterface()
                self.tables[1].selectObjects([currentObject[0]], setUpdatesEnabled=False)
                nextAtom = self.tables[1].getSelectedObjects()
                if nextAtom:
                    self._updateAssignmentWidget(1, currentObject[0])

                    self.lastTableSelected = 1
                    self.buttonList.setButtonEnabled('Delete', True)
                    self.buttonList.setButtonEnabled('Deassign', False)
                    self.buttonList.setButtonEnabled('Assign', True)

                else:
                    self._updateAssignmentWidget(1, None)

                    self.lastTableSelected = 1
                    self.buttonList.setButtonEnabled('Delete', False)
                    self.buttonList.setButtonEnabled('Deassign', False)
                    self.buttonList.setButtonEnabled('Assign', False)

        except Exception as es:
            showWarning('Deassign NmrAtom', str(es))

    def setAssignedTable(self, atomList: list):

        self.tables[0].populateTable(rowObjects=atomList,
                                     columnDefs=self.columnDefs
                                     )
        self.tables[0].sortByColumn(4, QtCore.Qt.AscendingOrder)
        objs = self.tables[0].getFirstObject()
        if objs:
            objPid = objs.get('Pid')
            if objPid:
                nmrAtom = self.project.getByPid(objPid)
                if nmrAtom:
                    self._updatePulldownLists(0, {Notifier.OBJECT: [nmrAtom]})

    def setAlternativesTable(self, atomList: list):

        self.tables[1].populateTable(rowObjects=atomList,
                                     columnDefs=self.columnDefs
                                     )
        self.tables[1].sortByColumn(4, QtCore.Qt.AscendingOrder)

    def _updateAssignmentWidget(self, tableNum: int, item: object):
        """
        Update all information in assignment widget when NmrAtom is selected in list widget of that
        assignment widget.
        """
        nmrAtom = item

        if nmrAtom:
            nmrChain = nmrAtom.nmrResidue.nmrChain
            sequenceCode = nmrAtom.nmrResidue.sequenceCode
            residueType = nmrAtom.nmrResidue.residueType

            if not self._parent.allChainCheckBoxLabel.isChecked():
                self._setChains(nmrChain)
                # self.chainPulldown.setIndex(self.chainPulldown.texts.index(nmrChain.id) if nmrChain.id in self.chainPulldown.texts else 0)

                self._setSequenceCodes(nmrChain)
                self.seqCodePulldown.setIndex(self.seqCodePulldown.texts.index(sequenceCode) if sequenceCode in self.seqCodePulldown.texts else 0)

                self._setResidueTypes(nmrChain)
                self.resTypePulldown.setIndex(self.resTypePulldown.texts.index(residueType) if residueType in self.resTypePulldown.texts else 0)

                self._setAtomNames(nmrAtom)
                self.atomTypePulldown.setIndex(self.atomTypePulldown.texts.index(nmrAtom.name) if nmrAtom.name in self.atomTypePulldown.texts else 0)
            else:

                # only allow selection of peaks from the table
                # atoms = self.objectTables[dim].getObjects()
                atoms = self.tables[tableNum]._dataFrameObject.objects
                if atoms:
                    options = [[''], [''], [''], ['']]  #'[None] * 4  # 4 empty lists
                    for atom in atoms:
                        thisOpt = atom.id.split('.')

                        for optionNum in range(0, len(thisOpt)):
                            if options[optionNum]:
                                if thisOpt[optionNum] not in options[optionNum]:
                                    options[optionNum].append(thisOpt[optionNum])
                            else:
                                options[optionNum] = [thisOpt[optionNum]]

                    self.chainPulldown.setData(options[0])
                    self.seqCodePulldown.setData(options[1])
                    self.resTypePulldown.setData(options[2])
                    self.atomTypePulldown.setData(options[3])
                else:
                    self._setDefaultPulldowns()

                self.chainPulldown.setIndex(self.chainPulldown.texts.index(nmrChain.id) if nmrChain.id in self.chainPulldown.texts else 0)
                self.seqCodePulldown.setIndex(self.seqCodePulldown.texts.index(sequenceCode) if sequenceCode in self.seqCodePulldown.texts else 0)
                self.resTypePulldown.setIndex(self.resTypePulldown.texts.index(residueType) if residueType in self.resTypePulldown.texts else 0)
                self.atomTypePulldown.setIndex(self.atomTypePulldown.texts.index(nmrAtom.name) if nmrAtom.name in self.atomTypePulldown.texts else 0)

            self.lastNmrAtomSelected = (self.chainPulldown.currentText(),
                                        self.seqCodePulldown.currentText(),
                                        self.resTypePulldown.currentText(),
                                        self.atomTypePulldown.currentText())
        else:
            self._setDefaultPulldowns()
            self.lastNmrAtomSelected = None

    def _setDefaultPulldowns(self):
        """Clear the contents of the pullDowns
        """
        self.chainPulldown.clear()
        self.seqCodePulldown.clear()
        self.resTypePulldown.clear()
        self.atomTypePulldown.clear()

        self._setChains()
        self._setResidueTypes()
        self._setAtomNames()

    def _setChains(self, nmrChain=None):
        """Populate the chain pulldown from the project
        """
        thisChain = self.chainPulldown.currentText()
        chains = ['']
        chains.extend([chain.id for chain in self.project.nmrChains])
        if nmrChain:
            thisChain = nmrChain.id

        self.chainPulldown.setData(chains)
        self.chainPulldown.setIndex(self.chainPulldown.texts.index(thisChain) if thisChain in self.chainPulldown.texts else 0)

    def _setSequenceCodes(self, nmrChain=None):
        """Populate the sequenceCode pulldown from the nmrChain or project
        """
        thisSeq = self.seqCodePulldown.currentText()
        sequenceCodes = ['']
        if nmrChain:
            sequenceCodes.extend([nmrResidue.sequenceCode for nmrResidue in nmrChain.nmrResidues])
        else:
            sequenceCodes.extend([nmrResidue.sequenceCode for nmrResidue in self.project.nmrResidues])

        self.seqCodePulldown.setData(sorted(sequenceCodes, key=CcpnSorting.stringSortKey))
        self.seqCodePulldown.setIndex(self.seqCodePulldown.texts.index(thisSeq) if thisSeq in self.seqCodePulldown.texts else 0)

    def _setResidueTypes(self, nmrChain=None):
        """Populate the residueTypes pulldown from the nmrChain or project
        """
        thisRes = self.resTypePulldown.currentText()
        residueTypes = ['']
        if nmrChain:
            residueTypes.extend([nmrResidue.residueType for nmrResidue in nmrChain.nmrResidues])
        else:
            residueTypes.extend([nmrResidue.residueType for nmrResidue in self.project.nmrResidues])

        residueTypes.extend([nmrResidue[1] for nmrResidue in allowedResidueTypes])  # self.project.nmrResidues]
        residueTypes = list(set(OrderedDict.fromkeys(residueTypes)))

        self.resTypePulldown.setData(sorted(residueTypes, key=CcpnSorting.stringSortKey))
        self.resTypePulldown.setIndex(self.resTypePulldown.texts.index(thisRes) if thisRes in self.resTypePulldown.texts else 0)

    def _setAtomNames(self, nmrAtom=None):
        """Populate the atomNames pulldown from the project
        """
        thisAtom = self.atomTypePulldown.currentText()
        atomNames = ['']
        if self.current.peak:
            isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[self.index]
            atomNames += getIsotopeListFromCode(isotopeCode)

        if nmrAtom:
            atomNames.insert(0, nmrAtom.name)
            thisAtom = nmrAtom.name  # set only if nmrAtom defined

        if self.lastNmrAtomSelected:
            atomNames.extend([self.lastNmrAtomSelected[3]])

        self.atomTypePulldown.setData(sorted(list(set(atomNames)), key=greekKey))
        self.atomTypePulldown.setIndex(self.atomTypePulldown.texts.index(thisAtom) if thisAtom in self.atomTypePulldown.texts else 0)

    def _deleteNmrAtom(self, dim: int):
        """
        delete selected nmrAtom from project
        """
        if self.lastTableSelected is not None:

            # # deassign if assigned
            # self._deassignNmrAtom(dim)

            # remove from the table
            deleted = self.tables[self.lastTableSelected].deleteObjFromTable()
            if deleted:
                nextAtoms = self.tables[self.lastTableSelected].getSelectedObjects()

                # reset buttons
                if not nextAtoms:

                    self.buttonList.setButtonEnabled('Delete', False)
                    self.buttonList.setButtonEnabled('Deassign', False)
                    self.buttonList.setButtonEnabled('Assign', False)

                    self._updateAssignmentWidget(self.lastTableSelected, None)
                else:
                    self._updateAssignmentWidget(self.lastTableSelected, nextAtoms[0])

    def _pulldownEdited(self, pulldown: object):
        """
        Enable the assignment button if the text has changed in the pulldown
        """
        currentNmrAtomSelected = (self.chainPulldown.currentText(),
                                  self.seqCodePulldown.currentText(),
                                  self.resTypePulldown.currentText(),
                                  self.atomTypePulldown.currentText())
        enable = False in self._atomCompare(self.lastNmrAtomSelected, currentNmrAtomSelected)
        self.buttonList.setButtonEnabled('Assign', enable)

    def _atomCompare(self, atom1: tuple, atom2: tuple):
        """
        check whether the selection has changed from being clicked
        """
        if atom1 and atom2:
            return [True if a == b else False for a, b in zip(atom1, atom2)]
        else:
            return [False]
