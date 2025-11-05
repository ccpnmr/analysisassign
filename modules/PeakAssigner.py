"""
Module to assign peaks
Responds to current.peaks

"""
from __future__ import annotations


#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (https://www.ccpn.ac.uk) 2014 - 2025"
__credits__ = ("Ed Brooksbank, Morgan Hayward, Victoria A Higman, Luca Mureddu, Eliza Płoskoń",
               "Timothy J Ragan, Brian O Smith, Daniel Thompson",
               "Gary S Thompson & Geerten W Vuister")
__licence__ = ("CCPN licence. See https://ccpn.ac.uk/software/licensing/")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, https://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: djt540 $"
__dateModified__ = "$dateModified: 2025-11-05 11:31:58 +0000 (Wed, November 05, 2025) $"
__version__ = "$Revision: 3.3.3 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

__all__ = ["PeakAssigner"]

import numpy as np
import pandas as pd
from functools import partial
from collections import OrderedDict, Counter
from PyQt5 import QtGui, QtCore, QtWidgets
from time import time_ns

from ccpn.core.NmrAtom import NmrAtom, UnknownIsotopeCode
from ccpn.core.NmrResidue import NmrResidue, _getNmrResidue, MoveToEnd
from ccpn.core.Peak import Peak
from ccpn.core.lib import CcpnSorting
from ccpn.core.lib.AssignmentLib import nmrAtomsForPeaks, peaksAreOnLine, PROTEIN_NEF_ATOM_NAMES, NEF_ATOM_NAMES
from ccpn.core.lib.ContextManagers import undoBlock, undoBlockWithoutSideBar
from ccpn.core.lib.Notifiers import Notifier, _removeDuplicatedNotifiers
from ccpn.core.lib.DataFrameObject import DataFrameObject
from ccpn.core.lib.WeakRefLib import WeakRefDescriptor
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.ButtonList import ButtonList
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.Frame import Frame, ScrollableFrame
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.HLine import LabeledHLine
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.Splitter import Splitter, SplitterGroup
from ccpn.ui.gui.widgets.table._ProjectTable import _ProjectTableABC
from ccpn.ui.gui.widgets.Column import ColumnClass, Column
from ccpn.ui.gui.widgets.SpeechBalloon import SpeechBalloon
from ccpn.ui.gui.widgets.MessageDialog import showWarning, showYesNo
from ccpn.ui.gui.widgets.Font import getFontHeight, TABLEFONT, setWidgetFont
from ccpn.ui.gui.widgets.DropBase import DropBase
from ccpn.ui.gui.lib.GuiNotifier import GuiNotifier
from ccpn.ui.gui.guiSettings import getColours, DIVIDER, LABEL_WARNINGFOREGROUND
from ccpn.util.Logging import getLogger
from ccpn.util.Common import greekKey, _truncateText, getIsotopeListFromCode, makeIterableList
from ccpn.util.UpdateScheduler import UpdateScheduler
from ccpn.util.UpdateQueue import UpdateQueue
from ccpn.util.OrderedSet import OrderedSet
from ccpnmodel.ccpncore.lib.Constants import defaultNmrChainCode


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

MSG = 'Not-defined >Select any to start<'

ROWDEFAULT = 300
ROWSIZES = {7 : 3000,
            4 : 2500,
            0 : 2000,
            -1: 1200,
            -2: ROWDEFAULT,
            }

_showBorders = False  # for debugging of layout's
_margins = (2, 2, 2, 2)
ASSIGNEDROWS = 3
ALTERNATIVEROWS = 5
MINTABLEWIDTH = 150
DEFAULT_COLOR = QtGui.QColor('black')

PULLDOWNPREFIX = '--'
OtherByIC = f'{PULLDOWNPREFIX} Amino/Isotope specific {PULLDOWNPREFIX}'
OtherByResType = f'{PULLDOWNPREFIX} In this nmrResidue {PULLDOWNPREFIX}'
OtherNames = f'{PULLDOWNPREFIX} All other atom-types {PULLDOWNPREFIX}'


#=========================================================================================
# PeakAssigner
#=========================================================================================

class PeakAssigner(CcpnModule):
    """Module for assignment of nmrAtoms to the different axes of a peak.
    Module responds to current.peak
    """
    className = 'PeakAssigner'
    # override in specific module implementations
    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    settingsPosition = 'left'
    activePulldownClass = None

    # set the queue handling parameters
    _maximumQueueLength = 10
    _logQueue = False

    def __init__(self, mainWindow, name="Peak Assigner"):
        """
        Initialise the Module widgets
        """
        super().__init__(mainWindow=mainWindow, name=name)

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current

        self.Ndims = 1
        self.maxDims = 8
        self.dimensionTabs = []
        self.currentAtoms = None
        self._visibleDims = 0
        self._chemShifts = {}

        # add widgets to the module
        self._setWidgets()

        # populate the tables
        self._updateInterface(self.current.peaks)

        # notifier queue handling
        self._queuePending = UpdateQueue()
        self._queueActive = None
        self._lock = QtCore.QMutex()
        self._scheduler = UpdateScheduler(self.project, self._queueProcess, name='PeakAssigner',
                                          log=False, completeCallback=self.update)
        # set notifiers to respond to peaks
        self._registerNotifiers()

    def _setWidgets(self):
        """Add the widgets to the module
        """
        self.blockSignals(True)

        if self.includeSettingsWidget:
            # settings
            row = 0
            self.doubleToleranceCheckbox = CheckBox(self.settingsWidget, checked=False,
                                                    callback=self._updateInterface,
                                                    grid=(row, 1))
            Label(self.settingsWidget, text="Double Tolerances ", grid=(row, 0))

            row += 1
            self.intraCheckbox = CheckBox(self.settingsWidget, checked=False,
                                          callback=self._updateInterface,
                                          grid=(row, 1))
            Label(self.settingsWidget, text="Only Intra-residual ", grid=(row, 0))

            row += 1
            self.multiCheckbox = CheckBox(self.settingsWidget, checked=True,
                                          callback=self._updateInterface,
                                          grid=(row, 1))
            Label(self.settingsWidget, text="Allow Multiple Peaks ", grid=(row, 0))

            row += 1
            self.allChainCheckBoxLabel = CheckBox(self.settingsWidget, checked=False,
                                                  callback=self._updateInterface,
                                                  grid=(row, 1))
            Label(self.settingsWidget, "Peak Selection from Table", grid=(row, 0))

            self._height = getFontHeight()
            self._tableHeight = getFontHeight(name=TABLEFONT)
            self.settingsWidget.setContentsMargins(5, 5, 5, 5)
            self.settingsWidget.getLayout().setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
            self.settingsWidget.setScrollBarPolicies(scrollBarPolicies=('asNeeded', 'never'))

        row = 0
        # add a label for the selected peaks
        self.peakLabel = Label(parent=self.mainWidget, setLayout=True, spacing=(0, 0),
                               text='Current Peak: ' + MSG, bold=True,
                               grid=(row, 0), margins=_margins,
                               hAlign='left',
                               hPolicy='ignored',
                               vPolicy='fixed'
                               )
        row += 1
        # set up a frame for the dimension frames - scrollable frame not resizing correctly
        self.axisFrameWidget = ScrollableFrame(parent=self.mainWidget, showBorder=False, setLayout=True,
                                               acceptDrops=True, grid=(row, 0),
                                               scrollBarPolicies=('asNeeded', 'never')
                                               )
        # Set the constraint so that the frame can shrink to a minimum before
        # the scrollbars become active
        self.axisFrameWidget.layout().setSizeConstraint(QtWidgets.QLayout.SetMinimumSize)

        splitterList = []
        row += 1
        colIndex = 0
        for dimIndex in range(self.maxDims):
            # 4 axes per row
            if dimIndex == 4:
                row += 1
                colIndex = 0
            dimTab = AxisAssignmentObject(parent=self.axisFrameWidget, grid=(row, colIndex),
                                          parentModule=self, dimIndex=dimIndex,
                                          mainWindow=self.mainWindow,
                                          )
            splitterList.append(dimTab.split)
            self.dimensionTabs.append(dimTab)
            colIndex += 1

        self.splitterGroup = SplitterGroup(splitterList)

        self.blockSignals(False)

    def _registerNotifiers(self):
        # without a tableSelection specified in the table callback, this nmrAtom callback is needed
        # to update the table
        self.setNotifier(self.current, [Notifier.CURRENT],
                         targetName=Peak._pluralLinkName,
                         callback=self._updateCurrent,  # self._updateInterface,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.DELETE, Notifier.CREATE],
                         targetName=Peak.__name__,
                         callback=self._updatePeak,  # self._updateInterface,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.CHANGE, Notifier.RENAME, Notifier.CREATE, Notifier.DELETE],
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

    #-----------------------------------------------------------------------------------------
    # Notifier queue handling
    #-----------------------------------------------------------------------------------------

    def queueFull(self):
        """Method that is called when the queue is deemed to be too big.
        Apply overall operation instead of all individual notifiers.
        """
        if self._logQueue:
            # log the queue-time if required
            getLogger().debug(f'queueFull  {self.__class__.__name__}')
        self._updateInterface()

    def _queueProcess(self):
        """Process current items in the queue
        """
        with QtCore.QMutexLocker(self._lock):
            # protect the queue switching
            self._queueActive = self._queuePending
            self._queuePending = UpdateQueue()

        startTime = 0.0
        useQueueFull = (self._maximumQueueLength not in [0, None] and len(self._queueActive) > self._maximumQueueLength)
        if self._logQueue:
            # log the queue-time if required
            startTime = time_ns()
            getLogger().debug(
                    f'_queueProcess  {self.__class__.__name__}  len: {len(self._queueActive)}  useQueueFull: {useQueueFull}')

        if useQueueFull:
            # rebuild from scratch if the queue is too big
            try:
                self._queueActive = None
                self.queueFull()
            except Exception as es:
                getLogger().debug(f'Error in {self.__class__.__name__} update queueFull: {es}')

        else:
            executeQueue = _removeDuplicatedNotifiers(self._queueActive)
            if self._logQueue:
                getLogger().debug(f'execute-queue {self.__class__.__name__} len: {len(executeQueue)}')

            for itm in executeQueue:
                # process item if different from previous
                if self.application and self.application._disableQueueException:
                    func, data = itm
                    func(data)

                else:
                    try:
                        func, data = itm
                        func(data)
                    except Exception as es:
                        getLogger().debug(f'Error in {self.__class__.__name__} update - {es}')

        if self._logQueue:
            getLogger().debug(f'elapsed time {self.__class__.__name__} - {(time_ns() - startTime) / 1e9}')

    def _queueAppend(self, itm):
        """Append a new item to the queue
        """
        self._queuePending.put(itm)
        if not self._scheduler.isActive and not self._scheduler.isBusy:
            self._scheduler.start()

        elif self._scheduler.isBusy:
            # caught during the queue processing event, need to restart
            self._scheduler.signalRestart()

    #-----------------------------------------------------------------------------------------
    # Notifier queue handling
    #-----------------------------------------------------------------------------------------

    def _updateCurrent(self, data):
        # not a very efficient way of doing this
        # self._updateInterface(data, action=data[Notifier.TRIGGER])
        self._queueAppend([self._updateInterface, data])

    def _updatePeak(self, data):
        # not a very efficient way of doing this
        # self._updateInterface(data, action=data[Notifier.TRIGGER])
        self._queueAppend([self._updateInterface, data])

    def _updateNmrAtom(self, data):
        # not a very efficient way of doing this
        # self._updateInterface(data, action=data[Notifier.TRIGGER])
        self._queueAppend([self._updateInterface, data])

    def _updateNmrResidue(self, data):
        # not a very efficient way of doing this
        # self._updateInterface(data, action=data[Notifier.TRIGGER])
        self._queueAppend([self._updateInterface, data])

    def _updateInterface(self, data=None, action=None):
        """Updates the whole module, including recalculation
           of which nmrAtoms fit to the peaks.
        """
        peaks = self.current.peaks
        self._cachedShifts = {}
        self._cachedTableShifts = {}
        self._cachedTableDeltas = {}

        if not peaks or not self._peaksAreCompatible(peaks):
            self.axisFrameWidget.hide()
            self.peakLabel.setText('Current Peak: ' + MSG)
        else:

            if len(self.current.peaks) < 2:
                self.peakLabel.setText(f'Current Peak: {self.current.peak.id if self.current.peak else ""}')

            else:
                # update the peaksLabel
                peaksIds = ' , '.join([str(pp.id) for pp in self.current.peaks])
                self.peakLabel.setText(f'Current Peaks: {_truncateText(peaksIds, maxWords=6)}')
                self.peakLabel.setToolTip(peaksIds)

            self.Ndims = self.current.peak.spectrum.dimensionCount

            # show/hide the dimension tabs
            for dimIndex, dimTab in enumerate(self.dimensionTabs):
                if dimIndex < self.Ndims:
                    dimTab.show()
                else:
                    dimTab.hide()

            for dimIndex, dimTab in enumerate(self.dimensionTabs[:self.Ndims]):
                # show/hide if peaks are aligned
                aligned = peaksAreOnLine(peaks=peaks, dimIndex=dimIndex)
                if not aligned:
                    txt = '%s: peaks\nnot aligned' % (peaks[0].axisCodes[dimIndex],)
                    dimTab.setNotAlignedText(txt)
                else:
                    ppmValues = np.array([pk.ppmPositions[dimIndex] for pk in peaks])
                    txt = '%s: %.3f' % (peaks[0].axisCodes[dimIndex], ppmValues.mean())
                    dimTab.setHLineText(txt)
                dimTab.showNotAligned(not aligned)

            self._updateTables(peaks=peaks)
            self.axisFrameWidget.show()

    def _updateTables(self, peaks):
        """
        update Assigned and alternatives tables showing which nmrAtoms
        are assigned to which peak dimensions. If multiple peaks are selected,
        only the assignment that they have in common are shown. Maybe this should be all
        assignments. You can see that at the peak annotation though.
        """

        doubleTolerance = self.doubleToleranceCheckbox.isChecked()
        intraResidual = self.intraCheckbox.isChecked()

        validNmrAtoms = [nmrAtom for nmrAtom in self.project.nmrAtoms if
                         not nmrAtom.nmrResidue.isDeleted or not nmrAtom.isDeleted]
        nmrAtomsForTables = nmrAtomsForPeaks(peaks, validNmrAtoms,
                                             doubleTolerance=doubleTolerance,
                                             intraResidual=intraResidual)

        Ndimensions = self.Ndims
        self.currentList = []

        _sizes = []
        for dim, nmrAtoms in zip(range(Ndimensions), nmrAtomsForTables):
            ll = [set(peak.dimensionNmrAtoms[dim]) for peak in peaks]
            self.nmrAtoms = list(sorted(set.intersection(*ll)))
            self.nmrAtoms = [nmrAtom for nmrAtom in self.nmrAtoms if
                             not nmrAtom.nmrResidue.isDeleted or not nmrAtom.isDeleted]

            self.currentList.append([str(a.pid) for a in self.nmrAtoms])  # ejb - keep another list
            self.dimensionTabs[dim].setAssignedTable(self.nmrAtoms)

            nmrAtomsForTables[dim] = [nmr for nmr in nmrAtomsForTables[dim] if nmr not in self.nmrAtoms]
            self.dimensionTabs[dim].setAlternativesTable(nmrAtomsForTables[dim])

        return _sizes

    def _getCachedShift(self, shiftList, nmrAtom):
        """Get the chemicalShift or the cached if exists
        """
        if (shiftList, nmrAtom) in self._cachedShifts:
            return self._cachedShifts.get((shiftList, nmrAtom))

        sh = shiftList.getChemicalShift(nmrAtom)
        self._cachedShifts[(shiftList, nmrAtom)] = sh

        return sh

    def _getDeltaShift(self, nmrAtom: NmrAtom, dim: int) -> float | str:
        """
        Calculation of delta shift to add to the table.
        """
        if (not self.current.peaks) or nmrAtom is NOL:
            return ''

        if nmrAtom in self._cachedTableDeltas:
            return self._cachedTableDeltas[nmrAtom]

        deltas = []
        for peak in self.current.peaks:
            if (shiftList := peak.peakList.spectrum.chemicalShiftList):
                if (shift := self._getCachedShift(shiftList, nmrAtom)):  # shiftList.getChemicalShift(nmrAtom)
                    _value = shift.value
                    if _value is not None:
                        position = peak.position[dim]
                        deltas.append(abs(shift.value - position))
        # average = sum(deltas)/len(deltas) #Bug: ZERO DIVISION!

        if len(deltas) > 0:
            _val = float(np.mean(deltas))  #'%6.3f' % np.mean(deltas) - handled by table
        else:
            _val = ''

        self._cachedTableDeltas[nmrAtom] = _val
        return _val

    def _getShift(self, nmrAtom: NmrAtom) -> float | str | None:
        """
        Calculation of chemical shift value to add to the table.
        """
        if (not self.current.peaks) or nmrAtom is NOL:
            return ''

        if nmrAtom in self._cachedTableShifts:
            return self._cachedTableShifts[nmrAtom]

        for peak in self.current.peaks:
            if (shiftList := peak.peakList.spectrum.chemicalShiftList):
                if (shift := self._getCachedShift(shiftList, nmrAtom)):  # shiftList.getChemicalShift(nmrAtom)
                    _val = shift.value  # '%8.3f' % shift.value
                    self._cachedTableShifts[nmrAtom] = _val
                    return _val
                    # return shift.value  # '%8.3f' % shift.value

    def _peaksAreCompatible(self, peaks) -> bool:
        """
        If multiple peaks are selected, a check is performed
        to determine whether assignment of corresponding
        dimensions of a peak allowed.
        """
        if len(peaks) == 1:
            return True

        if not self.multiCheckbox.isChecked():
            getLogger().warning("Selecting multiple peaks is currently not active (change in settings pane)")
            return False

        dimensionalities = set(peak.spectrum.dimensionCount for peak in peaks)
        if len(dimensionalities) > 1:
            getLogger().warning('Not all selected peaks have the same number of dimensions')
            return False

        for dimIndex in range(peaks[0].spectrum.dimensionCount):
            isotopeCodes = set(peak.spectrum.isotopeCodes[dimIndex] for peak in peaks)
            if len(isotopeCodes) > 1:
                getLogger().warning('Selected peaks have different isotopeCodes along dimension %d' % (dimIndex + 1))
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


#=========================================================================================
# NotOnLine
#=========================================================================================

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
_EDIT_OPTION = 'Edit nmrAtom'
_NEW_OPTION = 'New nmrAtom'

#=========================================================================================
# AssignmentTable
#=========================================================================================

# The two AssignmentTables for each dimension
_ASSIGNED_TABLE = 0
_ALTERNATIVES_TABLE = 1


class AssignmentTable(_ProjectTableABC):
    """Subclassed for some added functionality"""

    # define the notifiers that are required for the specific table-type
    tableClass = None
    rowClass = None
    cellClass = None
    tableName = 'assignedPeaks'
    rowName = None
    cellClassNames = None
    selectCurrent = True
    callBackClass = NmrAtom
    search = False

    _enableSelectionCallback = True
    _enableActionCallback = True

    # set the queue handling parameters
    _maximumQueueLength = 10  # shouldn't be responding to any notifiers

    defaultHidden = ['Pid']
    _internalColumns = ['_object']

    _dim = None
    _enableSearch = False
    _owner: WeakRefDescriptor[AxisAssignmentObject] = WeakRefDescriptor()
    atomList: list[NmrAtom] | None = []

    defaultSortColumn = 'Delta'
    defaultSortOrder = QtCore.Qt.AscendingOrder

    def __init__(self, parent, dim=0, dimIndex=None, owner=None, *args, **kwds):
        """Initialise the table and store as top- (dim=0) or-bottom (dim=1) table
        :param dimIndex: the dimension index for the assignments of a peak.
        """
        self._dim = dim
        if dimIndex is None or dimIndex < 0:
            raise ValueError(f'Initialising AssignmentTable: invalid {dimIndex = }')
        self._dimIndex = dimIndex
        self._owner = owner
        
        super().__init__(parent, *args, **kwds)

    #-----------------------------------------------------------------------------------------
    # Build the dataFrame for the table
    #-----------------------------------------------------------------------------------------

    def buildTableDataFrame(self):
        """Return a Pandas dataFrame from an internal list of objects
        """
        """Return a Pandas dataFrame from an internal list of objects.
        The columns are based on the 'func' functions in the columnDefinitions.
        :return pandas dataFrame
        """
        allItems = []
        objects = []

        if self.atomList:
            self._columnDefs = self._getTableColumns(self.atomList)

            for col, obj in enumerate(self.atomList):
                listItem = OrderedDict()
                for header in self._columnDefs.columns:
                    try:
                        listItem[header.headerText] = header.getValue(obj)
                    except Exception as es:
                        getLogger().debug2(f'Error creating table information {es}')
                        listItem[header.headerText] = None

                allItems.append(listItem)
                objects.append(obj)

            df = pd.DataFrame(allItems, columns=self._columnDefs.headings)

        else:
            self._columnDefs = self._getTableColumns()
            df = pd.DataFrame(columns=self._columnDefs.headings)

        # use the object as the index, object always exists even if isDeleted
        df.set_index(df[self.OBJECTCOLUMN], inplace=True, )

        _dfObject = DataFrameObject(dataFrame=df,
                                    columnDefs=self._columnDefs or [],
                                    table=self)

        return _dfObject

    #-----------------------------------------------------------------------------------------
    # Table functions
    #-----------------------------------------------------------------------------------------

    def _getTableColumns(self, nmrAtoms=None):
        """Add default columns plus the ones according to peakList.spectrum dimension
        format of column = ( Header Name, value, tipText, editOption)
        editOption allows the user to modify the value content by doubleclick
        """

        # set column definitions and hidden columns for each table
        self._columnDefs = ColumnClass([])
        self._columnDefs._columns = [
            Column('NmrAtom', lambda nmrAtom: str(nmrAtom.id), tipText='NmrAtom identifier'),
            Column('Pid', lambda nmrAtom: str(nmrAtom.pid), tipText='Pid of the nmrAtom'),
            Column('_object', lambda nmrAtom: nmrAtom, tipText='Object'),
            Column('Delta',
                   lambda nmrAtom: self.moduleParent._getDeltaShift(nmrAtom, self._dimIndex),
                   tipText='Delta-shift', format='%0.3f'),
            Column('Shift', lambda nmrAtom: self.moduleParent._getShift(nmrAtom), tipText='Chemical-shift',
                   format='%8.3f'),
            ]
        return self._columnDefs

    def clearSelection(self):
        """Clear the current selection in the table
        """
        # core-object selection is not required here
        self.selectionModel().clearSelection()

    def addTableMenuOptions(self, menu):
        """Add options to the right-mouse menu
        """
        super().addTableMenuOptions(menu)

        if self._dim == 0:
            self._peakMenuAction = menu.addAction(f'Deassign from Peak', self._peakActionCallback)
        else:
            self._peakMenuAction = menu.addAction(f'Assign to Peak', self._peakActionCallback)

        self._editMenuAction = menu.addAction(f'{_EDIT_OPTION}...', self._editNmrAtom)
        self._newMenuAction = menu.addAction(_NEW_OPTION, self._newNmrAtom)

        if (_actions := menu.actions()):
            _topMenuItem = _actions[0]
            _topSeparator = menu.insertSeparator(_topMenuItem)
            # move new actions to the top of the list
            menu.insertAction(_topSeparator, self._newMenuAction)
            menu.insertAction(_topSeparator, self._editMenuAction)
            menu.insertAction(self._newMenuAction, self._peakMenuAction)

    def setTableMenuOptions(self, menu):
        """Update options in the right-mouse menu
        """
        super().setTableMenuOptions(menu)

        selection = self.getSelectedObjects()
        data = self.getRightMouseItem()
        if data is not None and not data.empty and selection:
            # add more information to the edit nmrAtom option in the menu
            currentNmrAtom = selection[0]
            self._editMenuAction.setText(f'{_EDIT_OPTION}{" " + currentNmrAtom.id if currentNmrAtom else "..."}')
            self._editMenuAction.setEnabled(True if currentNmrAtom else False)
            self._peakMenuAction.setEnabled(True if currentNmrAtom else False)

        else:
            # disabled but visible lets user know that menu items exist
            self._editMenuAction.setText(f'{_EDIT_OPTION}...')
            self._editMenuAction.setEnabled(False)
            self._peakMenuAction.setEnabled(False)

        # hide the previous edit balloon (looks a little cleaner)
        self._owner.setEditPopupVisible(False)

    def _editNmrAtom(self):
        """Edit the nmrAtom from the parent widget
        """
        selection = self.getSelectedObjects()
        data = self.getRightMouseItem()
        if data is not None and not data.empty and selection:
            self._owner.lastTableSelected = self._dim
            # call the edit popup balloon
            self._owner._reassignNmrAtomPopup(mode=1)

    def _newNmrAtom(self):
        """Create new nmrAtom from the parent widget
        """
        self._owner.lastTableSelected = self._dim
        # call the new popup balloon
        self._owner._newNmrAtomPopup(mode=1)

    def _peakActionCallback(self):
        """Assign/deassign the peak
        """
        self._owner.lastTableSelected = self._dim
        if self._dim == 0:
            # deAssign from top to bottom
            self._owner._deassignNmrAtom(self._owner.dimIndex)
        elif self._dim == 1:
            # assign bottom - up
            self._owner._assignNmrAtom(self._owner.dimIndex, action=True)

    #-----------------------------------------------------------------------------------------
    # Selection/action callbacks
    #-----------------------------------------------------------------------------------------

    def actionCallback(self, selection, lastItem):
        """Notifier DoubleClick action on item in table. Mark a chemicalShift based on all attached nmrAtoms
        """
        try:
            objs = list(lastItem[self._OBJECT])

        except Exception as es:
            getLogger().debug2(f'{self.__class__.__name__}.actionCallback: No selection\n{es}')

        else:
            # nmrAtom = objs[0] if isinstance(objs, (list, tuple)) else objs

            if self._dim == _ASSIGNED_TABLE:
                # deAssign from top to bottom
                self._owner._deassignNmrAtom(self._owner.dimIndex)
            elif self._dim == _ALTERNATIVES_TABLE:
                # assign bottom - up
                self._owner._assignNmrAtom(self._owner.dimIndex, action=True)

    def selectionCallback(self, selected, deselected, selection, lastItem):
        """Notifier Callback for selecting rows in the table
        """
        try:
            objs = list(selection[self._OBJECT])

        except Exception as es:
            getLogger().debug2(f'{self.__class__.__name__}.selectionCallback: No selection\n{es}')

        else:
            # enable the edit-button
            self._owner._clickedTableCallback(self._dim, {Notifier.OBJECT: objs})

    def _selectCurrentCallBack(self, data):
        """Callback from a current changed notifier to highlight the current objects
        :param data
        """
        pass


#=========================================================================================
# EditNmrAtomBalloon
#=========================================================================================

class EditNmrAtomBalloon(SpeechBalloon):
    """Balloon to hold the pulldown lists for editing the nmrAtom
    """

    def __init__(self, mainWindow=None, project=None, *args, **kwds):
        super().__init__(*args, **kwds)
        self._mainWindow = mainWindow
        self._project = project
        self.setWindowModality(QtCore.Qt.NonModal)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setWindowFlags(int(self.windowFlags()) | QtCore.Qt.Popup)
        setWidgetFont(self)
        self.setStyleSheet('QToolTip {{ background-color: {TOOLTIP_BACKGROUND}; '
                           'color: {TOOLTIP_FOREGROUND}; '
                           'font-size: {_size}pt ; }}'.format(_size=self.font().pointSize(), **getColours()))


#=========================================================================================
# AxisAssignmentObject
#=========================================================================================

class AxisAssignmentObject(Frame):
    """
    Create a new frame for displaying information in 1 axis of peakassigner
    """

    # soft-links to external classes
    mainWindow = WeakRefDescriptor()
    application = WeakRefDescriptor()
    project = WeakRefDescriptor()
    current = WeakRefDescriptor()
    _parent = WeakRefDescriptor()

    def __init__(self, parent, parentModule, dimIndex, mainWindow, grid=None, **kwds):

        super().__init__(parent=parent,
                         setLayout=True, showBorder=_showBorders,
                         grid=grid, **kwds
                         )

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        if self.mainWindow:
            self.application = mainWindow.application
            self.project = mainWindow.application.project
            self.current = mainWindow.application.current
        self.currentAtoms = None
        self._clickedNmrAtom = None
        self._blockEscapeFlag = None

        # initialise axis information
        self.dimIndex = dimIndex
        self._parent = parentModule
        self.dataFrameAssigned = None
        self.dataFrameAlternatives = None
        self.lastTableSelected = None
        self.lastNmrAtomSelected = None
        self.tables: list[AssignmentTable | None] = [None, None]  # The two tables (assignment and alternatives)

        height = 20
        # self._minWidth = 150
        _minTabWidth = 100
        _tabHeight = 100
        _pullDownWidth = 65

        # assignments
        asRow = 0
        self._assignmentsFrame = Frame(self, setLayout=True, showBorder=_showBorders,
                                       grid=(asRow, 0), margins=_margins, acceptDrops=True)
        self._parent.setGuiNotifier(self._assignmentsFrame, [GuiNotifier.DROPEVENT], [DropBase.PIDS],
                                    callback=self._handleDropsFromSideBar)

        self.split = Splitter(parent=self._assignmentsFrame, horizontal=False)

        self.getLayout().addWidget(self.split)

        self.topSplit = Frame(None, setLayout=True, showBorder=True)
        self.bottomSplit = Frame(None, setLayout=True, showBorder=True)

        self.split.addWidget(self.topSplit)
        self.split.addWidget(self.bottomSplit)

        row = 0
        self.hLine = LabeledHLine(self.topSplit, text='axis', grid=(row, 0), height=16,
                                  colour=getColours()[DIVIDER])

        row += 1
        tt = self.tables[_ASSIGNED_TABLE] = AssignmentTable(parent=self.topSplit,
                                                            mainWindow=mainWindow,
                                                            grid=(row, 0), gridSpan=(1, 1),
                                                            # tipText='Click to select; double-click to de-assign'
                                                            showVerticalHeader=False,
                                                            multiSelect=False,
                                                            dim=_ASSIGNED_TABLE,
                                                            dimIndex=dimIndex,
                                                            owner=self
                                                            )
        tt.moduleParent = self._parent
        # Slight priority to the upper table
        # self._assignmentsFrame.layout().setRowStretch(row, 5)

        row += 1
        _buttons = ButtonList(self.topSplit, texts=['Edit', 'New'],
                              tipTexts=['Rename selected nmrAtom', 'Create new nmrAtom'],
                              callbacks=[self._reassignNmrAtomPopup,
                                         self._newNmrAtomPopup],
                              grid=(row, 0),
                              hAlign='l'
                              )
        row += 1
        self._alternativesLabel = Label(self.bottomSplit, 'Alternatives', hAlign='l', grid=(row, 0))
        self._alternativesLabel.setMinimumHeight(height)
        row += 1
        tt = self.tables[_ALTERNATIVES_TABLE] = AssignmentTable(parent=self.bottomSplit,
                                                                mainWindow=mainWindow,
                                                                grid=(row, 0), gridSpan=(1, 1),
                                                                # tipText='Click to select; double-click to assign'
                                                                showVerticalHeader=False,
                                                                multiSelect=False,
                                                                dim=_ALTERNATIVES_TABLE,
                                                                dimIndex=dimIndex,
                                                                owner=self
                                                                )
        tt.moduleParent = self._parent
        # self._assignmentsFrame.layout().setRowStretch(row, 4)

        self.editButton = _buttons.getButton('Edit')
        self.newNmrAtomButton = _buttons.getButton('New')

        #-----------------------------------------------------------------------------------------==
        # Not-aligned frame
        self.notAlignedFrame = Frame(self, setLayout=True, showBorder=_showBorders, grid=(asRow, 0),
                                     margins=_margins, )
        self.notAlignedLabel = Label(parent=self.notAlignedFrame, text='peaks\nnot aligned', grid=(0, 0),
                                     hAlign='centre',
                                     textColour=getColours()[LABEL_WARNINGFOREGROUND])

        self._assignmentWidget = self._nmrAtomWidget(parent=self.topSplit, minWidth=_pullDownWidth,
                                                     setLayout=True, showBorder=_showBorders, grid=(0, 0))
        self.editPopup = EditNmrAtomBalloon(mainWindow=self.mainWindow, project=self.project, on_top=True)
        self.editPopup.setCentralWidget(self._assignmentWidget)
        self.editPopup.hide()

        self.editButton.enableWidget(False)

    def _nmrAtomWidget(self, parent, minWidth, **kwds):
        """Make Frame with the nmrAtom Pulldown widgets
        :return Frame instance
        """
        _frame = Frame(parent=parent, **kwds)
        self.chainPulldown = self._createChainPulldown(parent=_frame,
                                                       grid=(0, 0), gridSpan=(1, 1),
                                                       tipText='Chain code', minWidth=minWidth)
        self.seqCodePulldown = self._createPulldown(parent=_frame,
                                                    grid=(0, 1), gridSpan=(1, 1),
                                                    tipText='Sequence code', minWidth=minWidth)
        self.resTypePulldown = self._createPulldown(parent=_frame,
                                                    grid=(0, 2), gridSpan=(1, 1),
                                                    tipText='Residue type', minWidth=minWidth)
        self.atomTypePulldown = self._createPulldown(parent=_frame,
                                                     grid=(0, 3), gridSpan=(1, 1),
                                                     tipText='Atom type', minWidth=minWidth)
        _innerFrame = Frame(parent=_frame, setLayout=True, grid=(1, 0), gridSpan=(1, 4))
        self._acceptMode = 0
        self._acceptFuncs = [self._reassignAccept, self._assignNewAccept]
        self._acceptButton = Button(parent=_innerFrame, text='Accept', grid=(1, 0), hAlign='r',
                                    callback=self._acceptNmrAtomCallback)
        # activate return/enter on the button when focussed
        self._acceptButton.setAutoDefault(True)
        self.chainPulldown.activated.connect(self._userSelectChainFromPulldown)
        self.seqCodePulldown.activated.connect(self._userSelectSeqCodeFromPulldown)
        self.resTypePulldown.activated.connect(self._userSelectResTypeFromPulldown)
        return _frame

    def _userSelectChainFromPulldown(self, *args):
        """Check the chain/sequenceCode and update the residueType/atomNames if set
        """
        # set residues
        nmrChain = self.project.getNmrChain(self.chainPulldown.currentText())
        nmrResidues = [f'{nmrResidue.sequenceCode}' for nmrResidue in nmrChain.nmrResidues]
        self.seqCodePulldown.setData(texts=nmrResidues)
        # just clarify that they are from different pulldowns
        self._userSelectSeqCodeFromPulldown(*args)

    def _userSelectSeqCodeFromPulldown(self, *args):
        """Check the chain/sequenceCode and update the residueType/atomNames if set
        """
        nmrChain = self.chainPulldown.currentText()
        seqCode = self.seqCodePulldown.currentText()

        _nmrChain = self.project.getNmrChain(nmrChain)
        if nmrResidue := _getNmrResidue(_nmrChain, seqCode):
            # set the residueType pulldown
            resType = nmrResidue.residueType
            _ind = self.resTypePulldown.texts.index(resType) if resType in self.resTypePulldown.texts else 0
            self.resTypePulldown.setIndex(_ind)
            self.resTypePulldown.repaint()
            # set the atom-names
            self._setAtomNames(nmrResidue=nmrResidue)
            # colour as required
            self._setPulldownColours(nmrResidue)

        else:
            self._setAtomNames()
            self._resetPulldownColours()

    def _userSelectResTypeFromPulldown(self, *args):
        """Check the residueType and update the atomname list as necessary
        """
        nmrChain = self.chainPulldown.currentText()
        seqCode = self.seqCodePulldown.currentText()
        _nmrChain = self.project.getNmrChain(nmrChain)
        nmrResidue = _getNmrResidue(_nmrChain, seqCode)

        self._setAtomNames(nmrResidue=nmrResidue or None)
        if nmrResidue:
            self._setPulldownColours(nmrResidue)
        else:
            self._resetPulldownColours()

    def _acceptNmrAtomCallback(self, pulldown=None, *args):
        """Perform different acceptFunc depending on the mode
        """
        if pulldown:
            QtCore.QTimer.singleShot(0, pulldown.setFocus)
        else:
            _func = self._acceptFuncs[self._acceptMode]
            _func()

    def setEditPopupVisible(self, visible):
        """Hide the edit popup balloon
        """
        self.editPopup.setVisible(visible)

    def showNotAligned(self, flag):
        """Show/hide of notAligned and assignmentFrame"""
        self.notAlignedFrame.setVisible(flag)
        self.topSplit.setVisible(not flag)
        self.bottomSplit.setVisible(not flag)
        self._assignmentsFrame.setVisible(not flag)

    def setHLineText(self, text):
        """Set the text of the top horizontal Line"""
        self.hLine.setText(text)

    def setNotAlignedText(self, text):
        """Set the text of the notAligned widget"""
        self.notAlignedLabel.setText(text)

    def _handleDropsFromSideBar(self, dataDict):
        """
        Handle drops from SideBar. If NmrAtoms, then assign to the selected peaks.
        """
        objs = self.project.getObjectsByPids(dataDict.get(DropBase.PIDS))
        nmrAtoms = [x for x in objs if isinstance(x, NmrAtom)]

        if self.current.peak:
            failedNmrAtoms = []
            isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[self.dimIndex]
            for nmrAtom in nmrAtoms:
                if isotopeCode == nmrAtom.isotopeCode or nmrAtom.isotopeCode in [UnknownIsotopeCode, None]:
                    self._assignNmrAtom(self.dimIndex, nmrAtoms=[nmrAtom])
                else:
                    failedNmrAtoms.append(nmrAtom)
            if failedNmrAtoms:
                showWarning('Incompatible IsotopeCode Error',
                            f'Cannot assign NmrAtoms: {nmrAtoms} to peaks with IsotopeCode {isotopeCode} ')

    def _handleDroppedItems(self, droppingToTableNum: int, dataDict, ):
        """
        Notifier callback activated upon a DropEvent of an object.
        Note, the source of the drag can be from anywhere, therefore here is limited only if the source is
        within the module and right tables pairs. The correct instance of the dropped object is checked afterwards.
        """
        assignmentTableNum = 0
        alternativeTableNum = 1
        sourceTable = dataDict.get('source')
        nmrAtoms = self.project.getObjectsByPids(dataDict.get(DropBase.PIDS))

        ## Action 0, Assignment: dropping to Assignment (Table-0) from Alternative (Table-1)
        if droppingToTableNum == assignmentTableNum and sourceTable == self.tables[alternativeTableNum]:
            self._assignNmrAtom(self.dimIndex, nmrAtoms=nmrAtoms)
            return
        ## Action 1, DeAssign from top to bottom: dropping to Alternative (Table-1) from Assignment (Table-0)
        if droppingToTableNum == alternativeTableNum and sourceTable == self.tables[assignmentTableNum]:
            self._deassignNmrAtom(self.dimIndex, nmrAtoms=nmrAtoms)
            return

    def _assignDeassignNmrAtom(self, tableNum: int, data):
        """
        Assign/Deassign the nmrAtom that is double-clicked to
        the corresponding dimension of the selected
        peaks.
        """
        if tableNum == 0:
            # deAssign from top to bottom
            self._deassignNmrAtom(self.dimIndex)
        elif tableNum == 1:
            # assign bottom - up
            self._assignNmrAtom(self.dimIndex, action=True)

    def _clickedTableCallback(self, tableNum, data):
        if obj := data[Notifier.OBJECT]:
            self.lastTableSelected = tableNum
            self._clickedNmrAtom = obj[0]

            if tableNum == 0:
                # this will clear the other table and fire its selection
                # which will first disable the edit-button below
                self.tables[1].clearSelection()
            elif tableNum == 1:
                self.tables[0].clearSelection()

            # re-enable the button
            self.editButton.enableWidget(True)

        else:
            self.editButton.enableWidget(False)

    def _clearTableCallback(self, tableNum, data):
        self._clickedNmrAtom = None
        self.editButton.enableWidget(False)
        self.tables[0].clearSelection()
        self.tables[1].clearSelection()

    def _createChainPulldown(self, parent=None, grid=(0, 0), gridSpan=(1, 1), tipText='', minWidth=50) -> PulldownList:
        """Creates a PulldownList with callback, editable.
        """
        pulldownList = PulldownList(parent=parent, grid=grid, backgroundText=tipText, editable=True, gridSpan=gridSpan,
                                    tipText=tipText)
        # pulldownList.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        # pulldownList.setEditable(True)
        # pulldownList.lineEdit().returnPressed.connect(partial(self._acceptNmrAtomCallback, pulldown=pulldownList))
        pulldownList.lineEdit().textChanged.connect(partial(self._chainEdited, pulldownList))
        pulldownList.setMinimumWidth(minWidth)
        self._updateCompleter(pulldownList)
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

    def _createPulldown(self, parent=None, grid=(0, 0), gridSpan=(1, 1), tipText='', minWidth=50,
                        popupMode=QtWidgets.QCompleter.PopupCompletion) -> PulldownList:
        """Creates a PulldownList with callback, editable.
        """
        pulldownList = PulldownList(parent=parent, grid=grid, backgroundText=tipText, editable=True, gridSpan=gridSpan,
                                    tipText=tipText)
        pulldownList.setMinimumWidth(minWidth)
        self._updateCompleter(pulldownList, popupMode=popupMode)
        return pulldownList

    def _updateCompleter(self, pulldownList, popupMode=QtWidgets.QCompleter.PopupCompletion):
        pulldownList.lineEdit().returnPressed.connect(partial(self._acceptNmrAtomCallback, pulldown=pulldownList))
        completer = pulldownList.completer()
        completer.setCompletionMode(popupMode)
        completer.setMaxVisibleItems(16)
        popup = completer.popup()
        popup.setObjectName('_COMPLETER')
        popup.installEventFilter(self)
        setWidgetFont(popup)
        pulldownList.view().setObjectName('_PULLDOWNVIEW')
        pulldownList.view().installEventFilter(self)
        pulldownList.setObjectName('_PULLDOWN')
        pulldownList.installEventFilter(self)

    def eventFilter(self, source: 'QObject', event: 'QEvent') -> bool:
        if source.objectName() in {'_PULLDOWN'}:
            if event.type() == QtCore.QEvent.KeyPress and event.key() in {QtCore.Qt.Key_Escape}:
                QtCore.QTimer.singleShot(0, source.setFocus)
                if source.view().isVisible():
                    source.hidePopup()
                    return True
                if self._blockEscapeFlag:  # == ('setfocus', source):
                    self._blockEscapeFlag = None
                    return True
            elif event.type() == QtCore.QEvent.KeyPress and event.key() in {QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter}:
                if self._blockEscapeFlag:
                    self._blockEscapeFlag = None
                else:
                    self._acceptNmrAtomCallback()
        elif source.objectName() in {'_PULLDOWNVIEW'}:
            if event.type() in {QtCore.QEvent.Hide}:
                self._blockEscapeFlag = True  # ('setfocus', source._pulldown)
                QtCore.QTimer.singleShot(0, source._pulldown.lineEdit().setFocus)
            elif event.type() == QtCore.QEvent.KeyPress and event.key() in {QtCore.Qt.Key_Escape}:
                return True
        elif source.objectName() in {'_COMPLETER'}:
            if event.type() in {QtCore.QEvent.Hide}:
                self._blockEscapeFlag = None
            elif event.type() == QtCore.QEvent.KeyPress and event.key() in {QtCore.Qt.Key_Escape}:
                source.hide()
                return True
            elif event.type() == QtCore.QEvent.KeyPress and event.key() in {QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter}:
                self._blockEscapeFlag = None

        return super().eventFilter(source, event)

    def _newNmrAtomPopup(self, mode=0):
        """Callback for the newNmrAtom button
        """
        from ccpn.core.NmrAtom import NmrAtom
        from ccpn.ui.gui.widgets.BalloonMetrics import Side

        # get the next available Ids
        _id = NmrAtom._nextId()
        nmrChain, sequenceCode, residueType, atomName = _id.split('.')

        # add the names to the pulldowns
        self._updateAssignmentWidget(0, None)
        self.chainPulldown.select('@-')
        self.seqCodePulldown.addItem(sequenceCode)
        self.seqCodePulldown.select(sequenceCode)

        if atomName not in self.atomTypePulldown.texts:
            self.atomTypePulldown.insertText(0, atomName)
            self.atomTypePulldown.select(atomName)

        # set the accept mode
        self._acceptMode = 1

        # show the popup
        pos = QtGui.QCursor().pos()
        if mode == 0:
            # if called from the button then set the pointer size - otherwise hide it
            global_rect = QtCore.QRect(self.newNmrAtomButton.mapToGlobal(QtCore.QPoint(0, 0)),
                                       self.newNmrAtomButton.geometry().size())
            self.editPopup.pointerHeight = 10
        else:
            global_rect = pos
            self.editPopup.pointerHeight = 0

        mouse_screen = next((screen for screen in QtGui.QGuiApplication.screens() if screen.geometry().contains(pos)),
                            None)
        for combo in {self.chainPulldown, self.seqCodePulldown, self.resTypePulldown, self.atomTypePulldown}:
            # need to change from editable->non-editable->editable to force stylesheets to update correctly
            combo.setEditable(False)
            combo.setEditable(True)
        self.editPopup.showAt(global_rect, preferred_side=Side.TOP,
                              side_priority=(Side.TOP, Side.BOTTOM, Side.RIGHT, Side.LEFT),
                              target_screen=mouse_screen)
        QtCore.QTimer.singleShot(0, self.chainPulldown.setFocus)

    def _assignNewAccept(self):
        _undoIndex = self.project._undo.nextIndex
        try:
            dim = self.dimIndex

            with undoBlock():
                # get the isotope code for the current dimension
                isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[dim]

                nmrChainName = self.chainPulldown.currentText()
                seqCode = self.seqCodePulldown.currentText()
                resType = self.resTypePulldown.currentText()
                nmrAtomName = self.atomTypePulldown.currentText()

                # search for an existing nmrAtom or create a new one
                nmrChain = self.project.fetchNmrChain(shortName=nmrChainName or defaultNmrChainCode)
                if not (nmrResidue := _getNmrResidue(nmrChain, seqCode)):
                    nmrResidue = nmrChain.fetchNmrResidue(sequenceCode=seqCode, residueType=resType)

                    mainRess = nmrChain.mainNmrResidues
                    if nmrChain.isConnected and len(mainRess) > 1:
                        assignedPks = makeIterableList(self.current.peak.assignedNmrAtoms)
                        assignedRess = {nmrAt.nmrResidue for nmrAt in assignedPks}

                        foundInds = Counter(mainRess.index(nmr) for nmr in assignedRess if nmr in mainRess)
                        if 0 in foundInds and showYesNo('Assign New NmrAtom',
                                                        'A new nmrResidue has been created for the new nmrAtom.\n'
                                                        'By default this will be placed at the end of the connected nmrChain.\n\n'
                                                        'Do you want to move the new nmrResidue to the front of the nmrChain instead?'):
                            nmrResidue.moveToEnd(MoveToEnd.HEAD)

                elif nmrResidue.residueType != resType:
                    # if existing then check the residueType matches
                    raise ValueError(f'residueType does not match existing nmrResidue {nmrResidue.id}')

                self._acceptNmrAtom = nmrResidue.fetchNmrAtom(name=nmrAtomName, isotopeCode=isotopeCode)
                nmrAtom = self._acceptNmrAtom

                for peak in self.current.peaks:
                    if nmrAtom not in peak.dimensionNmrAtoms[dim]:
                        newAssignments = list(peak.dimensionNmrAtoms[dim]) + [nmrAtom]
                        axisCode = peak.spectrum.axisCodes[dim]
                        peak.assignDimension(axisCode, newAssignments)

                # highlight on the table and populate the pulldowns
                self.tables[0].highlightObjects([nmrAtom])  #, setUpdatesEnabled=False)
                self.tables[1].clearSelection()

                # No need for update, as this will be done by the callback on the newNmrAtom/newNmrResidue
                # self._parent._updateInterface()
                # self._updateAssignmentWidget(0, nmrAtom)

                self.lastTableSelected = 0
                self.lastNmrAtomSelected = nmrAtom

                self._clickedNmrAtom = nmrAtom

        except Exception as es:
            while (ni := self.project._undo.nextIndex) > _undoIndex:
                self.project._undo.undo()
                if self.project._undo.nextIndex == ni:
                    getLogger().debug(f'*** {self.__class__.__name__}:_assignNewAccept - error processing undo stack')
                    break
            showWarning(str(self.windowTitle()), str(es))
        else:
            self._reassignNmrAtom()
        finally:
            self.editPopup.setVisible(False)

    def _showNmrAtomPopup(self, nmrAtom, tableNum, mode=0):
        """Call the popup with the supplied nmrAtom
        """
        from ccpn.ui.gui.widgets.BalloonMetrics import Side

        # if nmrAtom:
        self._updateAssignmentWidget(tableNum, nmrAtom)

        pos = QtGui.QCursor().pos()
        if mode == 0:
            # if called from the button then set the pointer size - otherwise hide it
            global_rect = QtCore.QRect(self.editButton.mapToGlobal(QtCore.QPoint(0, 0)),
                                       self.editButton.geometry().size())
            self.editPopup.pointerHeight = 10
        else:
            global_rect = pos
            self.editPopup.pointerHeight = 0

        mouse_screen = next((screen for screen in QtGui.QGuiApplication.screens() if screen.geometry().contains(pos)),
                            None)
        for combo in {self.chainPulldown, self.seqCodePulldown, self.resTypePulldown, self.atomTypePulldown}:
            # need to change from editable->non-editable->editable to force stylesheets to update correctly
            combo.setEditable(False)
            combo.setEditable(True)
        self.editPopup.showAt(global_rect, preferred_side=Side.TOP,
                              side_priority=(Side.TOP, Side.BOTTOM, Side.RIGHT, Side.LEFT),
                              target_screen=mouse_screen)

        # give the popup time to appear
        QtCore.QTimer.singleShot(0, self.chainPulldown.setFocus)

    def _reassignNmrAtomPopup(self, mode=0):
        """Show the edit popup
        """
        _tableNum = self.lastTableSelected
        if (nextAtom := self.tables[_tableNum].getSelectedObjects()):
            self._clickedNmrAtom = nextAtom[0]

        self._acceptMode = 0
        self._showNmrAtomPopup(nextAtom[0] if nextAtom else None, _tableNum, mode)

    def _reassignAccept(self):
        """Handle the accept button in the balloon popup
        """
        self._reassignNmrAtom()
        # self.editPopup.setVisible(False)

    def _reassignNmrAtom(self):
        """
        Assigns dimensionNmrAtoms to peak dimension when called using Assign Button in assignment widget.
        """
        _undoIndex = self.project._undo.nextIndex
        try:
            nmrChainName = self.chainPulldown.currentText()
            seqCode = self.seqCodePulldown.currentText()
            newResType = self.resTypePulldown.currentText()
            nmrAtomName = self.atomTypePulldown.currentText()

            if not self._clickedNmrAtom:
                showWarning("Rename NmrAtom", "Please select an NmrAtom from the tables")
                return

            # wrap all actions in a single undo block
            with undoBlock():
                _chainPid = f'NC:{nmrChainName}'
                _nmrChain = self.project.fetchNmrChain(nmrChainName)
                nmrResidue = _getNmrResidue(_nmrChain, seqCode, )
                # nmrAtom = nmrResidue.getNmrAtom(nmrAtomName) if nmrResidue else None

                # edit existing
                if nmrResidue and self._clickedNmrAtom.nmrResidue != nmrResidue:
                    if nmrAtom := nmrResidue.getNmrAtom(nmrAtomName):
                        if showYesNo('Merge NmrAtom',
                                     f'Do you want to merge\n'
                                     f'{self._clickedNmrAtom.id}   into   {nmrAtom.id}',
                                     dontShowEnabled=True,
                                     defaultResponse=True,
                                     popupId=f'{self.__class__.__name__}Merge'):
                            # merge into the existing nmrAtom
                            nmrAtom.mergeNmrAtoms(self._clickedNmrAtom)
                    else:
                        # assign to a new nmrAtom
                        self._clickedNmrAtom.assignTo(chainCode=nmrChainName,
                                                      sequenceCode=seqCode,
                                                      residueType=newResType,
                                                      name=nmrAtomName,
                                                      mergeToExisting=False)

                elif nmrResidue:
                    # rename the same nmrAtom
                    if newResType != nmrResidue.residueType:
                        nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)
                    if nmrAtomName != self._clickedNmrAtom.name:
                        if nmrAtom := nmrResidue.getNmrAtom(nmrAtomName):
                            if showYesNo('NmrAtom already exists',
                                         f'Do you want to merge\n'
                                         f'{self._clickedNmrAtom.id}   into   {nmrAtom.id}',
                                         dontShowEnabled=True,
                                         defaultResponse=True,
                                         popupId=f'{self.__class__.__name__}MergeExist'):
                                # merge into the existing nmrAtom
                                nmrAtom.mergeNmrAtoms(self._clickedNmrAtom)
                        else:
                            # rename the nmrAtom
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
            while (ni := self.project._undo.nextIndex) > _undoIndex:
                self.project._undo.undo()
                if self.project._undo.nextIndex == ni:
                    getLogger().debug(f'*** {self.__class__.__name__}:_reassignNmrAtom - error processing undo stack')
                    break
            showWarning('Rename NmrAtom', str(es))
        finally:
            # okay, close the popup
            self.editPopup.setVisible(False)

    def _assignNmrAtom(self, dim: int, action: bool = False, create: bool = True, nmrAtoms=None):
        """
        Assigns dimensionNmrAtoms to peak dimension when called using Assign Button in assignment widget.
        :param dim: axis dimension of the atom
        :type dim: int
        :param action: True if callback is action from the table
        :type action: bool
        """
        # FIXME Potential Bug: no error checks for dim. It can give easily an IndexError

        # return if no peaks selected
        if not self.current.peaks:
            return

        try:
            selectedObjects = nmrAtoms or self.tables[1].getSelectedObjects()
            if not (selectedObjects and selectedObjects[0]):
                return
            nmrAtom = selectedObjects[0]
            if not isinstance(nmrAtom, NmrAtom) or nmrAtom.isDeleted:
                return

            # nmrAtom = None

            # wrap all actions in a single undo block
            with undoBlockWithoutSideBar():

                # NOTE:ED need to keep for the minute
                # _chainPid = 'NC:{}'.format(nmrChainName)
                # if create and not action:
                #     # get the current chain (but may create a new one)
                #     _nmrChain = self.project.fetchNmrChain(nmrChainName)
                # else:
                #     # find the existing nmrChain
                #     _nmrChain = self.project.getByPid(_chainPid)
                #     if not _nmrChain:
                #         # raise error to notify popup
                #         raise ValueError("NmrChain doesn't exists")
                #
                # nmrResidue = _getNmrResidue(_nmrChain, seqCode, )
                # nmrAtom = nmrResidue.getNmrAtom(nmrAtomName) if nmrResidue else None
                #
                # if not action:
                #     if create:
                #         if nmrResidue:
                #             if self._clickedNmrAtom and self._clickedNmrAtom.nmrResidue == nmrResidue and nmrResidue.residueType != newResType:
                #                 if len(nmrResidue.nmrAtoms) > 1:
                #                     yes = showYesNoWarning('Assigning nmrAtoms',
                #                                            'This will change all nmrAtoms to the residueType {}, continue?'.format(newResType))
                #                     if yes:
                #                         nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)
                #                 else:
                #                     nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)
                #
                #         else:
                #             # can do a residueType rename
                #             nmrResidue = _nmrChain.fetchNmrResidue(seqCode, newResType)
                #
                #         nmrAtom = nmrResidue.fetchNmrAtom(nmrAtomName)
                #
                #     else:
                #         pass

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

                    ## Set the isotopeCode here if was not defined yet
                    if not nmrAtom.isotopeCode:
                        isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[dim]
                        nmrAtom._setIsotopeCode(isotopeCode)


                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

            self._parent._updateInterface()

            self.tables[0].highlightObjects([nmrAtom])  #, setUpdatesEnabled=False)
            self.lastTableSelected = 0

            # if nmrAtom:
            #     # self._updateAssignmentWidget(0, nmrAtom)
            #     self.lastTableSelected = 0
            #
            # else:
            #     # self._updateAssignmentWidget(0, None)
            #     self.lastTableSelected = 0

            # update the module
            self.update()

        except Exception as es:
            showWarning('Assign NmrAtom', str(es))

    def _deassignNmrAtom(self, dim: int, nmrAtoms=None):
        """
        remove nmrAtom from peak assignment
        """
        # return if no peaks selected
        if not self.current.peaks:
            return

        try:
            currentObjects = nmrAtoms or self.tables[0].getSelectedObjects()
            if not currentObjects:
                return
            nmrAtom = currentObjects[0]
            if not isinstance(nmrAtom, NmrAtom) or nmrAtom.isDeleted:
                return

            try:
                with undoBlockWithoutSideBar():
                    for peak in self.current.peaks:
                        peakDimNmrAtoms = peak.dimensionNmrAtoms
                        dimNmrAtoms = list(peakDimNmrAtoms[dim])  # ejb - changed to list
                        if nmrAtom in dimNmrAtoms:
                            dimNmrAtoms.remove(nmrAtom)

                        allAtoms = list(peakDimNmrAtoms)
                        allAtoms[dim] = dimNmrAtoms
                        peak.dimensionNmrAtoms = allAtoms

            except Exception as es:
                showWarning(str(self.windowTitle()), str(es))

            self._parent._updateInterface()
            self.tables[1].highlightObjects([nmrAtom])  #, setUpdatesEnabled=False)
            nextAtom = self.tables[1].getSelectedObjects()
            self.lastTableSelected = 1

            # if nextAtom:
            #     # self._updateAssignmentWidget(1, currentObject[0])
            #
            #     self.lastTableSelected = 1
            #     # self.buttonList.setButtonEnabled('Delete', True)
            #     # self.buttonList.setButtonEnabled('Deassign', False)
            #     # self.buttonList.setButtonEnabled('Assign', True)
            #
            # else:
            #     # self._updateAssignmentWidget(1, None)
            #
            #     self.lastTableSelected = 1
            #     # self.buttonList.setButtonEnabled('Delete', False)
            #     # self.buttonList.setButtonEnabled('Deassign', False)
            #     # self.buttonList.setButtonEnabled('Assign', True) #False)

        except Exception as es:
            showWarning('Deassign NmrAtom', str(es))

    def setAssignedTable(self, atomList: list):

        self.tables[0].atomList = atomList
        self.tables[0].populateTable()

    def setAlternativesTable(self, atomList: list):

        self.tables[1].atomList = atomList
        self.tables[1].populateTable()

    def _updateAssignmentWidget(self, tableNum: int, item: NmrAtom | None):
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
                self.chainPulldown.selectValue(nmrChain.id)

                self._setSequenceCodes(nmrChain)
                idx = self.seqCodePulldown.selectValue(sequenceCode)

                self._setResidueTypes(nmrChain)
                idx = self.resTypePulldown.selectValue(residueType)

                self._setAtomNames(nmrAtom)
                idx = self.atomTypePulldown.selectValue(nmrAtom.name)

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

                self.chainPulldown.setIndex(
                        self.chainPulldown.texts.index(nmrChain.id) if nmrChain.id in self.chainPulldown.texts else 0)
                self.seqCodePulldown.setIndex(
                        self.seqCodePulldown.texts.index(
                                sequenceCode) if sequenceCode in self.seqCodePulldown.texts else 0)
                self.resTypePulldown.setIndex(
                        self.resTypePulldown.texts.index(
                                residueType) if residueType in self.resTypePulldown.texts else 0)
                self.atomTypePulldown.setIndex(self.atomTypePulldown.texts.index(
                        nmrAtom.name) if nmrAtom.name in self.atomTypePulldown.texts else 0)
                self._setPulldownColours(nmrAtom.nmrResidue)

            self.lastNmrAtomSelected = nmrAtom
        else:
            self._setDefaultPulldowns()
            self.lastNmrAtomSelected = None

    def _resetPulldownColours(self):
        # reset all items to the palette foreground colour
        for combo in (self.resTypePulldown, self.atomTypePulldown):
            model = combo.model()
            for ii in range(len(combo.texts)):
                itm = model.item(ii)
                if PULLDOWNPREFIX not in itm.text():
                    # clear the colour and revert to palette.Text
                    itm.setData(None, QtCore.Qt.ForegroundRole)
            combo.repaint()

    def _setPulldownColours(self, nmrResidue: NmrResidue):
        if not nmrResidue:
            return

        residueType = nmrResidue.residueType
        blueCol = QtGui.QColor('dodgerblue')
        greenCol = QtGui.QColor('seagreen')
        combo = self.resTypePulldown
        model = combo.model()
        _inds = [ii for ii, val in enumerate(self.resTypePulldown.texts) if val and val == residueType]
        for ind in range(len(combo.texts)):
            itm = model.item(ind)
            if PULLDOWNPREFIX not in itm.text() and ind in _inds:
                itm.setData(blueCol, QtCore.Qt.ForegroundRole)
        combo.repaint()

        combo = self.atomTypePulldown
        model = combo.model()
        _inds = {ii for ii, val in enumerate(self.atomTypePulldown.texts) if
                 val in [nmrAt.name for nmrAt in nmrResidue.nmrAtoms]}
        for ind in range(len(combo.texts)):
            itm = model.item(ind)
            if PULLDOWNPREFIX not in itm.text() and ind in _inds:
                itm.setData(greenCol, QtCore.Qt.ForegroundRole)

        if self._clickedNmrAtom:
            _inds = {ii for ii, val in enumerate(self.atomTypePulldown.texts) if
                     val and val == self._clickedNmrAtom.name}
            for ind in _inds:
                itm = model.item(ind)
                if PULLDOWNPREFIX not in itm.text():
                    itm.setData(blueCol, QtCore.Qt.ForegroundRole)
        combo.repaint()

    def _setDefaultPulldowns(self):
        """Clear the contents of the pullDowns
        """
        self.chainPulldown.clear()
        self.seqCodePulldown.clear()
        self.resTypePulldown.clear()
        self.atomTypePulldown.clear()

        self._setChains()
        self._setSequenceCodes()
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
        self.chainPulldown.setIndex(
                self.chainPulldown.texts.index(thisChain) if thisChain in self.chainPulldown.texts else 0)

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
        self.seqCodePulldown.setIndex(
                self.seqCodePulldown.texts.index(thisSeq) if thisSeq in self.seqCodePulldown.texts else 0)

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
        self.resTypePulldown.setIndex(self.resTypePulldown.texts.index(thisRes)
                                      if thisRes in self.resTypePulldown.texts else 0)
        self.resTypePulldown.repaint()

    def _setAtomNames(self, nmrAtom=None, nmrResidue=None):
        """Populate the atomNames pulldown from the project
        """
        thisAtom = self.atomTypePulldown.currentText()
        thisNmrResidueType = self.resTypePulldown.currentText()
        # get atom-names from the isotope-list for the selected dimension (NEF_ATOM_NAMES)
        try:
            isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[self.dimIndex]
        except Exception:
            # should always exist, but will cause following to get all atom-names
            isotopeCode = '?'
        isotopeCodeAtoms = OrderedSet(sorted(getIsotopeListFromCode(isotopeCode), key=greekKey))

        nmrAtomName = None
        if nmrAtom:
            nmrAtomName = nmrAtom.name  # set only if nmrAtom defined, find parent nmrResidue
            nmrResidue = nmrAtom.nmrResidue
        # start with the last typed in value
        _atomNameOptions = [thisAtom, ] if (thisAtom and not thisAtom.startswith(PULLDOWNPREFIX)) else []

        if thisNmrResidueType:
            # get the defined nmrAtoms on the selected nmrResidue
            thisNmrResAtoms = OrderedSet(sorted([atm.name for atm in nmrResidue.nmrAtoms] if nmrResidue else [],
                                                key=greekKey))
            # get the list of specific codes based on residueType
            if atomsByResType := PROTEIN_NEF_ATOM_NAMES.get(thisNmrResidueType, []):
                isotopeCodeAtoms = OrderedSet(sorted(atomsByResType, key=greekKey))
            else:
                isotopeCodeAtoms = OrderedSet(sorted(getIsotopeListFromCode(None),
                                                     key=greekKey))

            if isotopeCode in NEF_ATOM_NAMES:
                # isotope-code is valid from the spectrum dimension
                atomsByIsotopeCode = OrderedSet(sorted(getIsotopeListFromCode(isotopeCode or nmrAtom.isotopeCode),
                                                       key=greekKey))
                atomOfSameIsotopeCode = isotopeCodeAtoms & atomsByIsotopeCode
                atomNotOfSameIsotopeCode = isotopeCodeAtoms - atomsByIsotopeCode
                if atomOfSameIsotopeCode:
                    _atomNameOptions += ([OtherByIC] + list(atomOfSameIsotopeCode))
                # if thisNmrResAtoms:
                #     _atomNameOptions += ([OtherByResType] + list(thisNmrResAtoms - atomOfSameIsotopeCode))
                # if atomNotOfSameIsotopeCode:
                #     _atomNameOptions += ([OtherNames] + list(atomNotOfSameIsotopeCode - thisNmrResAtoms))

            elif thisNmrResAtoms:
                _atomNameOptions += ([OtherByResType] +
                                     list(thisNmrResAtoms) +
                                     [OtherNames] +
                                     list(isotopeCodeAtoms - thisNmrResAtoms))
            else:
                if thisNmrResAtoms:
                    _atomNameOptions += ([OtherByResType] + list(thisNmrResAtoms))
                if isotopeCodeAtoms:
                    _atomNameOptions += ([OtherByIC] + list(isotopeCodeAtoms - thisNmrResAtoms))

        elif isotopeCodeAtoms:
            _atomNameOptions += ([OtherByIC] + list(isotopeCodeAtoms))

        if self.lastNmrAtomSelected:
            # add the last typed in value
            val = self.lastNmrAtomSelected.pid.fields[3]
            # if val not in self.atomTypePulldown.texts:
            if val not in _atomNameOptions:
                _atomNameOptions.insert(0, val)

        self.atomTypePulldown.setData(_atomNameOptions)
        self.atomTypePulldown.disableLabelsOnPullDown([OtherNames, OtherByIC, OtherByResType])

        if nmrResidue:
            self._setPulldownColours(nmrResidue)
        if nmrAtomName:
            self.atomTypePulldown.setIndex(self.atomTypePulldown.texts.index(thisAtom)
                                           if thisAtom in self.atomTypePulldown.texts else 0)
        self.atomTypePulldown.repaint()
        self.resTypePulldown.repaint()

        self.atomTypePulldown.disableLabelsOnPullDown([OtherNames, OtherByIC, OtherByResType])

    @staticmethod
    def _atomCompare(atom1: tuple, atom2: tuple):
        """
        check whether the selection has changed from being clicked
        """
        if atom1 and atom2:
            return [True if a == b else False for a, b in zip(atom1, atom2)]
        else:
            return [False]
