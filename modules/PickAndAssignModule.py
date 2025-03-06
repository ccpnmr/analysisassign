"""
Do a restricted peak pick along the 'y-axis' of (a set of) spectra.
Use settings to define the spectral displays, the active spectra and the tolerances for peak picking

This module closely works with the Atom Selector module

First version by SS
Refactored by GWV to be responsible to active state on start; proper callbacks 
and to include "Restricted pick and assign" button.
Refactored by GST to allow multiple NMRResidues to be picked and provide
progress dialogs during longer operations

TODO:
Assign selected should be deactivated if there are no peaks selected
All buttons should be deactivated if there are no NMRResidues selected
Button deassign selected
Button restricted Assign
Button delete peaks
Add tool tips for buttons
Meta A on table should select all and should be in the table menu
Enable up and down arrows
Command up and down should be up and down on table, up and down should be on last table with focus if still focused...
deleting a selection on a residue table doesn't update the pick and assign residue table and also doesn't trigger
enabling or disabling buttons...
can current nmrResidues be strings?
editing NMRResidue table not reflected in current table
deleting NMRResidue in NMRResidue table pops up a dialog!
deleting delete key doesn't delete current selected = NMRResidues

"""
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
__modifiedBy__ = "$modifiedBy: Daniel Thompson $"
__dateModified__ = "$dateModified: 2025-03-06 14:30:42 +0000 (Thu, March 06, 2025) $"
__version__ = "$Revision: 3.3.1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Geerten Vuister $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

from functools import partial
from typing import Iterator, Iterable

from OpenGL.logs import getLog
from PyQt5.QtWidgets import QStackedWidget
from PyQt5 import QtWidgets, QtCore
# from icecream import ic

from ccpn.core.Peak import Peak
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.lib.AssignmentLib import copyAssignmentsFromReference, propagateAssignments, copyAssignments
from ccpn.ui.gui.lib import PeakListLib
from ccpn.ui.gui.lib import StripLib
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTableModule, _NewNmrResidueTableWidget, NmrResidueTableFrame
from ccpn.ui.gui.modules.PeakTable import _NewPeakTableWidget, _PeakTableFrame
from ccpn.ui.gui.widgets.Base import Base
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.Font import getFontHeight
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.PulldownListsForObjects import PeakPulldown, NmrChainPulldown
from ccpn.ui.gui.widgets.SettingsWidgets import PickAndAssignSettings
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.Widget import Widget
from ccpn.util.OrderedSet import OrderedSet
from ccpn.util.Logging import getLogger
from ccpnmodel.ccpncore.lib.Io.PyMMLibPDB import KEYWDS


logger = getLogger()

ALL = '<Use all>'
SelectToAdd = '> select-to-add <'


# class PickAndAssignModuleO(NmrResidueTableModule):
#     """
#     Do a restricted peak pick along the 'y-axis' of (a set of) spectra.
#     Use settings to define the spectral displays, the active spectra and the tolerances for peak picking
#
#     # GST is this true any more?
#     This module closely works with the Atom Selector module
#     """
#     className = 'PickAndAssignModule'
#
#     includeSettingsWidget = True
#     maxSettingsState = 2
#     settingsPosition = 'left'
#     settingsMinimumSizes = (500, 200)
#
#     includePeakLists = False
#     includeNmrChains = False
#     includeSpectrumTable = True
#
#     includeDisplaySettings = True
#
#     def __init__(self, mainWindow, name='Pick and Assign'):
#
#         super().__init__(mainWindow=mainWindow, name=name, selectFirstItem=True)  # ejb ='Pick And Assign')
#
#         # Derive application, project, and current from mainWindow
#         self.mainWindow = mainWindow
#         self.application = mainWindow.application
#         self.project = mainWindow.application.project
#         self.current = mainWindow.application.current
#
#         # Main widget
#         restrictedPickAndAssignWithAssignFalse = partial(self.restrictedPickAndAssign, assign=False)
#         self.restrictedPickButton = Button(text='Restricted\nPick', callback=restrictedPickAndAssignWithAssignFalse)
#         self.tableFrame.addWidgetToPos(self.restrictedPickButton, row=0, col=2)
#
#         self.assignSelectedButton = Button(text='Assign\nSelected', callback=self.assignSelected)
#         self.tableFrame.addWidgetToPos(self.assignSelectedButton, row=0, col=3)
#
#         restrictedPickAndAssignWithAssignTrue = partial(self.restrictedPickAndAssign, assign=True)
#         self.restrictedPickAndAssignButton = Button(text='Restricted\nPick and Assign',
#                                                     callback=restrictedPickAndAssignWithAssignTrue)
#         self.tableFrame.addWidgetToPos(self.restrictedPickAndAssignButton, row=0, col=4)
#
#         self.restrictedPickButton.setEnabled(True)
#         self.assignSelectedButton.setEnabled(True)
#         self.restrictedPickAndAssignButton.setEnabled(True)
#
#         # change default-settings inherited from NmrResidueTableModule
#         self.nmrResidueTableSettings.sequentialStripsWidget.checkBox.setChecked(False)
#
#         if self.nmrResidueTableSettings.displaysWidget:
#             self.nmrResidueTableSettings.displaysWidget.addPulldownItem(0)  # select the <all> option
#
#         self.nmrResidueTableSettings.setLabelText('Navigate to\nDisplay(s)')
#
#         # need to feedback to current.nmrResidueTable
#         self._registerNotifiers()
#
#         # # these need to change whenever different spectrumDisplays are selected
#         if self.nmrResidueTableSettings.axisCodeOptions:
#             self.nmrResidueTableSettings.axisCodeOptions.selectAll()
#
#             # just clear the 'C' axes - this is the usual configuration
#             for ii, box in enumerate(self.nmrResidueTableSettings.axisCodeOptions.checkBoxes):
#                 if box.text().upper().startswith('C'):
#                     self.nmrResidueTableSettings.axisCodeOptions.clearIndex(ii)
#
#         # fix the second column to stop extra widgets flickering
#         alignWidgets(self.nmrResidueTableSettings, columnScale=1.2)


class PickAndAssignModule(CcpnModule):
    className = 'PickAndAssignModule'

    includeSettingsWidget = True
    maxSettingsState = 2
    settingsPosition = 'left'
    settingsMinimumSizes = (500, 200)

    includePeakLists = False
    includeNmrChains = False
    includeSpectrumTable = True

    includeDisplaySettings = True
    pickAndAssignSettings = True

    def __init__(self, mainWindow, name='Pick and Assign'):
        super().__init__(mainWindow=mainWindow, name=name)

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current

        self._settings = PickAndAssignSettings(parent=self.settingsWidget, mainWindow=mainWindow)
        self.stackedTableWidget = StackedTableFrameWidget(parent=self.mainWidget,
                                                          grid=(0, 0), moduleParent=self)

        self.tables = []
        self._setupTables()
        self._setupWidgets()

        # need to feedback to current.nmrResidueTable
        self._registerNotifiers()

    @property
    def currentTable(self):
        return self.stackedTableWidget.currentTable

    def _setupWidgets(self):
        for table in self.tables:
            # TODO Re-add button functionality
            # Main widget
            restrictedPickAndAssignWithAssignFalse = partial(self.restrictedPickAndAssign, assign=False)
            restrictedPickAndAssignWithAssignTrue = partial(self.restrictedPickAndAssign, assign=True)

            self.restrictedPickButton = Button(text='Restricted\nPick',
                                               callback=restrictedPickAndAssignWithAssignFalse)
            table.addWidgetToPos(self.restrictedPickButton, row=0, col=2)

            self.assignSelectedButton = Button(text='Assign\nSelected',
                                               callback=self.assignSelected)
            table.addWidgetToPos(self.assignSelectedButton, row=0, col=3)

            self.restrictedPickAndAssignButton = Button(text='Restricted\nPick and Assign',
                                                        callback=restrictedPickAndAssignWithAssignTrue)
            table.addWidgetToPos(self.restrictedPickAndAssignButton, row=0, col=4)

            self.restrictedPickButton.setEnabled(True)
            self.assignSelectedButton.setEnabled(True)
            self.restrictedPickAndAssignButton.setEnabled(True)

    def _setupTables(self):
        self.nmrChainTable = NmrResidueTableFrame(parent=self.stackedTableWidget, mainWindow=self.mainWindow,
                                                  moduleParent=self, grid=(0, 0))
        self.peakTable = _PeakTableFrame(parent=self.stackedTableWidget, mainWindow=self.mainWindow,
                                         moduleParent=self, grid=(0, 0))

        self.nmrChainTable.nmrResidueTableSettings = self._settings.nmrResidueTableSettings
        self.nmrResidueTableSettings = self.nmrChainTable.nmrResidueTableSettings
        self.peakTable._settings = self._settings.peakTableSettings

        # set existing widgets to false.
        self.peakTable.posUnitPulldownLabel.setEnabled(False)
        self.peakTable.posUnitPulldownLabel.setVisible(False)
        self.peakTable.posUnitPulldown.setEnabled(False)
        self.peakTable.posUnitPulldown.setVisible(False)

        self.tables = [self.nmrChainTable, self.peakTable]
        self.stackedTableWidget.addTablesToFrame(self.tables)

    def _registerNotifiers(self):
        """
        set up the notifiers
        """
        self.setNotifier(self.current,
                         [Notifier.CURRENT],
                         targetName=NmrResidue._pluralLinkName,
                         callback=self._selectionCallback)

    def _selectionCallback(self, data):
        """
        enable/disable the pick buttons
        """
        selected = data[Notifier.OBJECT].nmrResidue

        if selected:
            self.restrictedPickButton.setEnabled(True)
            self.assignSelectedButton.setEnabled(True)
            self.restrictedPickAndAssignButton.setEnabled(True)
        else:
            self.restrictedPickButton.setEnabled(False)
            self.assignSelectedButton.setEnabled(False)
            self.restrictedPickAndAssignButton.setEnabled(False)

    def _getDisplay(self):
        """Get the current selected spectrum-display from the pulldown
        """
        if (pulldown := self.nmrResidueTableSettings.spectrumDisplayPulldown) and \
                (texts := self.nmrResidueTableSettings.spectrumDisplayPulldown.getTexts()):
            if ALL in texts:
                gids = [self.application.getByGid(gid) for gid in pulldown.pulldownList.texts
                        if gid not in [ALL, SelectToAdd]]
            else:
                gids = [self.application.getByGid(gid) for gid in texts if gid not in [ALL, SelectToAdd]]
            return gids

    def _getMsgIfSetupInvalid(self):
        msg = None
        if self.currentTable is self.nmrChainTable:
            if not self.current.nmrResidues:
                # check that is defined and of the correct type
                msg = 'no NmrResidues selected, please pick one or more NmrResidues'
        elif self.currentTable is self.peakTable:
            if not self.current.peaks:
                msg = 'no Peaks selected, please pick one or more Peaks'

        if not msg and not self._getDisplay():
            # check the selected display
            msg = 'Undefined display;\nselect display in gearbox settings before proceeding'

        return msg

    def _getSelected(self):
        if self.currentTable is self.nmrChainTable:
            return list(self.current.nmrResidues)
        elif self.currentTable is self.peakTable:
            return list(self.currentTable.current.peaks)

    # def _getNmrResidues(self) -> list[NmrResidue] | list[Peak]:
    #     """ get the current selected NmrResidues
    #     """
    #
    #     nmrResidues = list(self.current.nmrResidues)
    #
    #     # GST: not sure if this needed - can current.nmrResidue[s] be a string?
    #     for i, nmrResidue in enumerate(nmrResidues):
    #         nmrResidues[i] = (self.project.getByPid(nmrResidue)) if isinstance(nmrResidues, str) else nmrResidue
    #
    #     return nmrResidues

    @staticmethod
    def _getValidPeakListViews(displays):
        """Get the list of valid peakListViews
        """
        validPeakListViews = {}
        # loop through all the selected displays/spectrumViews/peakListViews that are visible
        for dp in displays:

            # ignore undefined displays
            if not dp:
                continue

            if dp.strips:
                for sv in dp.strips[0].spectrumViews:
                    for plv in sv.peakListViews:
                        if plv.isDisplayed and sv.isDisplayed:
                            if plv.peakList not in validPeakListViews:
                                validPeakListViews[plv.peakList] = (sv.spectrum, plv)
                            else:
                                # skip for now, only one valid peakListView needed per peakList
                                pass
        return validPeakListViews

    # def assignSelected(self):
    #     """Assign current.peaks on the bases of nmrAtoms of current.nmrResidues
    #     """
    #     nmrResidues = self._getNmrResidues()
    #
    #     peaks = self.current.peaks
    #     if len(peaks) == 0:
    #         showWarning('Pick and Assign', 'No peaks currently selected')
    #         return
    #
    #     with undoBlockWithoutSideBar():
    #         self._assignPeaks(peaks, nmrResidues)

    def assignSelected(self):
        """Assign current.peaks on the bases of nmrAtoms of current.nmrResidues
        TODO: improve docstring
        """
        peaks = self.current.peaks

        if len(peaks) == 0:
            showWarning('Pick and Assign', 'No peaks currently selected')
            return

        with undoBlockWithoutSideBar():
            if self.currentTable is self.peakTable:
                self._assignSelectedPeaks(peaks)
            elif self.currentTable is self.nmrChainTable:
                nmrResidues = self._getSelected()
                self._assignSelectedResidues(peaks, nmrResidues)

    def _assignSelectedPeaks(self, peaks):
        copyAssignments(peaks)

    # convert to be an iterator...
    def _assignSelectedResidues(self, peaks, nmrResidues):
        displays = self._getDisplay()
        for display in displays:
            currentAxisCodeIndexes = self.nmrResidueTableSettings.axisCodeOptionsDict.get(f'{display}')

            for nmrResidue in nmrResidues:
                shiftDict = {}
                for atom in nmrResidue.nmrAtoms:
                    shiftDict[atom.isotopeCode] = []

                for specInd in self.nmrResidueTableSettings.spectrumIndex:
                    for peak in peaks:
                        if (spectrum := peak.peakList.spectrum) not in specInd:
                            continue

                        shiftList = peak.peakList.spectrum.chemicalShiftList
                        for nmrAtom in nmrResidue.nmrAtoms:
                            if nmrAtom.isotopeCode in shiftDict.keys():
                                cShift = shiftList.getChemicalShift(nmrAtom)
                                if cShift:
                                    shiftDict[nmrAtom.isotopeCode].append((nmrAtom, cShift.value))

                        for ii, isotopeCode in enumerate(spectrum.isotopeCodes):
                            if ii in specInd[spectrum]:
                                _restrictedIdx = specInd[spectrum].index(ii)
                                if (_restrictedIdx not in currentAxisCodeIndexes):
                                    continue
                            pValue = peak.position[ii]
                            if isotopeCode in shiftDict.keys():
                                shiftList = set()
                                for shift in shiftDict[isotopeCode]:
                                    sValue = shift[1]
                                    if abs(sValue - pValue) <= spectrum.assignmentTolerances[ii]:
                                        shiftList.add(shift[0])
                                if shiftList:
                                    peak.assignDimension(spectrum.axisCodes[ii], list(shiftList))

    @staticmethod
    def _getActionMsg(assign):
        return 'Restricted Pick and Assign' if assign else 'Restricted Pick'

    def restrictedPickAndAssign(self, assign=True):
        """
        Takes the selected NmrResidues from current NmrResidues feeds them into restricted pick lib functions
        and picks peaks for all spectrum displays specified in the settings tab. Pick uses X and Z axes for each
        spectrumView as centre points with tolerances and the y as the long axis to pick the whole region.
        """

        # if invalidMsg := self._getMsgIfSetupInvalid():
        #     showWarning(self._getActionMsg(assign), invalidMsg)
        # else:
        #     nmrResidues = self._getNmrResidues()
        #     self._doPickAndAssignOnSelectedNmrResidues(nmrResidues, assign)

        if invalidMsg := self._getMsgIfSetupInvalid():
            showWarning(self._getActionMsg(assign), invalidMsg)
            return

        if self.currentTable is self.peakTable:
            peaks = self._getSelected()
            self._doPickAndAssignOnSelectedObjs(peaks, assign)
        elif self.currentTable is self.nmrChainTable:
            nmrResidues = self._getSelected()
            self._doPickAndAssignOnSelectedObjs(nmrResidues, assign)

    def _doPickAndAssignOnSelectedNmrResidues(self, nmrResidues, assign):
        from ccpn.core.lib.ContextManagers import progressHandler

        undoStack = self.application._getUndo()
        originalUndoState = undoStack.undoList
        # ic('orig', originalUndoState)

        with undoBlockWithoutSideBar():
            msg = "Picking and Assigning Peaks..." if assign else "Picking peaks..."
            stopButtonText = 'Stop Pick and Assign' if assign else "Stop Picking"

            numResidues = len(nmrResidues)
            curPeaks = set()
            self.current.peaks = []  # option to do this?
            with progressHandler(text=msg, cancelButtonText=stopButtonText,
                                 maximum=numResidues) as progress:
                for i, nmrResidue, errorMsg, peaks in self._restrictedPeakPickIterator(nmrResidues):
                    progress.checkCancel()
                    if errorMsg:
                        showWarning(self._getActionMsg(assign), errorMsg)
                        progress.cancel()
                    progress.setValue(i)
                    if peaks and assign:
                        self._assignSelectedResidues(peaks, [nmrResidue, ])
                    curPeaks |= OrderedSet(peaks)

            self.current.peaks = list(OrderedSet(self.current.peaks) | curPeaks)
            if progress.cancelled:
                while undoStack.undoList != originalUndoState and undoStack.nextIndex > 0:
                    undoStack.undo()

    def _doPickAndAssignOnSelectedObjs(self, objs, assign):
        from ccpn.core.lib.ContextManagers import progressHandler

        undoStack = self.application._getUndo()
        originalUndoState = undoStack.undoList

        with undoBlockWithoutSideBar():
            msg = "Picking and Assigning Peaks..." if assign else "Picking peaks..."
            stopButtonText = 'Stop Pick and Assign' if assign else "Stop Picking"

            curPeaks = set()
            self.current.peaks = []  # option to do this?
            with progressHandler(text=msg, cancelButtonText=stopButtonText,
                                 maximum=len(objs)) as progress:
                for i, obj, errorMsg, peaks in self._restrictedPeakPickIterator(objs):
                    progress.checkCancel()
                    if errorMsg:
                        showWarning(self._getActionMsg(assign), errorMsg)
                        progress.cancel()
                    progress.setValue(i)
                    if peaks and assign:
                        # assign based on object type
                        if isinstance(obj, Peak):
                            self._assignSelectedPeaks([obj])
                        if isinstance(obj, NmrResidue):
                            self._assignSelectedResidues(peaks, [obj, ])
                    curPeaks |= OrderedSet(peaks)

            self.current.peaks = list(OrderedSet(self.current.peaks) | curPeaks)
            if progress.cancelled:
                while undoStack.undoList != originalUndoState and undoStack.nextIndex > 0:
                    undoStack.undo()

    # def _restrictedPeakPickIteratorOld(self, nmrResidues: Iterable[NmrResidue]) \
    #         -> Iterator[tuple[int | None, NmrResidue, str | None, list[Peak] | None]]:
    #
    #     displays = self._getDisplay()
    #     for display in displays:
    #         validPeakListViews = self._getValidPeakListViews([display])
    #         currentAxisCodeIndexes = self.nmrResidueTableSettings.axisCodeOptionsDict.get(f'{display}')
    #         for specInd in self.nmrResidueTableSettings.spectrumIndex:
    #             try:
    #
    #                 specAxisCodes = [[spectrum.axisCodes[specInd[spectrum].index(ii)]
    #                                   for ii in currentAxisCodeIndexes
    #                                   if ii in specInd[spectrum]]
    #                                  for spectrum, peakListView in validPeakListViews.values()
    #                                  if spectrum in specInd]
    #             except Exception:
    #                 # TODO: this should be a DataClass or named tuple for clarity
    #                 continue
    #                 # return None, None, badAxisCodeMsg, None
    #
    #             for i, nmrResidue in enumerate(nmrResidues):
    #
    #                 peaks = []
    #                 try:
    #                     for (spectrum, peakListView), axisCodes in zip(validPeakListViews.values(), specAxisCodes):
    #
    #                         # axis-codes should be valid at this point
    #                         peakList, pks = PeakListLib.restrictedPick(peakListView=peakListView,
    #                                                                    axisCodes=axisCodes, nmrResidue=nmrResidue)
    #                         if pks:
    #                             peaks += list(pks)
    #
    #                 except Exception:
    #                     continue
    #                     # return None, nmrResidue, badAxisCodeMsg, None
    #
    #                 yield i, nmrResidue, None, list(peaks)

    def _restrictedPeakPickIterator(self, iterObjs: Iterable[NmrResidue] | Iterable[Peak]):
        displays = self._getDisplay()
        for display in displays:
            validPeakListViews = self._getValidPeakListViews([display])
            currentAxisCodeIndexes = self.nmrResidueTableSettings.axisCodeOptionsDict.get(f'{display}')
            for specInd in self.nmrResidueTableSettings.spectrumIndex:
                try:
                    specAxisCodes = [[spectrum.axisCodes[specInd[spectrum].index(ii)]
                                      for ii in currentAxisCodeIndexes
                                      if ii in specInd[spectrum]]
                                     for spectrum, peakListView in validPeakListViews.values()
                                     if spectrum in specInd]
                except Exception:
                    continue

                for i, iterObj in enumerate(iterObjs):
                    pks = peaks = []
                    try:
                        for (spectrum, peakListView), axisCodes in zip(validPeakListViews.values(), specAxisCodes):
                            if isinstance(iterObj, Peak):
                                peakList, pks = PeakListLib.restrictedPick(peakListView=peakListView,
                                                                           axisCodes=axisCodes, peak=iterObj)
                            elif isinstance(iterObj, NmrResidue):
                                peakList, pks = PeakListLib.restrictedPick(peakListView=peakListView,
                                                                           axisCodes=axisCodes, nmrResidue=iterObj)
                            if pks:
                                peaks += list(pks)
                    except Exception as e:
                        getLogger().warning(f'{e.__traceback__}')
                    yield i, iterObj, None, list(peaks)


class StackedWidget(QStackedWidget, Base):
    def __init__(self, parent=None, **kwds):
        super().__init__(parent)
        Base._init(self, **kwds)


class StackedTableFrameWidget(Frame):
    """A frame that contains multiple stacked tables and a pulldown to control them."""

    def __init__(self, parent=None, mainWindow=None, moduleParent=None, **kwds):
        super().__init__(parent, setLayout=True, **kwds)

        self.mainWindow = mainWindow
        if mainWindow:
            self.application = mainWindow.application
            self.project = mainWindow.application.project
            self.current = mainWindow.application.current
        else:
            self.application = self.project = self.current = None

        self.moduleParent = moduleParent
        self.tableNameDict = dict()

        self.tablesWidget = StackedWidget(parent=self, grid=(0, 0), gridSpan=(2, 1), )
        self.currentTablePulldown = PulldownList(parent=self,
                                                 grid=(0, 0), hAlign='right', vAlign='t',
                                                 callback=self._switchTableCallback,
                                                 sizeAdjustPolicy=QtWidgets.QComboBox.AdjustToContents,
                                                 minimumWidths=(0, 100))

    @property
    def currentTable(self):
        return self.tablesWidget.currentWidget()

    @currentTable.setter
    def currentTable(self, table):
        if isinstance(table, str):
            try:
                table = self.tableNameDict.get(table)
            except KeyError:
                getLogger().error(f'{self.__class__} _switchTableCallback KeyError, table not found in nameDict')
                return

        self.tablesWidget.setCurrentWidget(table)

    def addTablesToFrame(self, tableFrames: list() = None):
        if tableFrames is None:
            getLogger().warning('No table frames given to initialise')
            return

        for tableFrame in tableFrames:
            self.tablesWidget.addWidget(tableFrame)
            self.addToControlPulldown(tableFrame)

        if self.tablesWidget.currentWidget() is None:
            self.tablesWidget.setCurrentIndex(0)

    def addToControlPulldown(self, table):
        tableName = table.guiTable.attributeName
        self.tableNameDict.update({tableName: table})
        self.currentTablePulldown.addItem(tableName)

    def removeFromControlPulldown(self, table):
        tableName = table.guiTable.attributeName
        self.tableNameDict.pop({tableName: table})
        self.currentTablePulldown.removeItem(tableName)

    def _switchTableCallback(self, value: None = None):
        self.currentTable = self.currentTablePulldown.getText()


# class PickAndAssignModuleNEW(CcpnModule):
#     className = 'PickAndAssignModule'
#
#     includeSettingsWidget = True
#     maxSettingsState = 2
#     settingsPosition = 'left'
#     settingsMinimumSizes = (500, 200)
#
#     includePeakLists = False
#     includeNmrChains = False
#     includeSpectrumTable = True
#
#     includeDisplaySettings = True
#     pickAndAssignSettings = True
#
#     def __init__(self, mainWindow, name='Pick and Assign'):
#         super().__init__(mainWindow=mainWindow, name=name)
#
#         # Derive application, project, and current from mainWindow
#         self.mainWindow = mainWindow
#         self.application = mainWindow.application
#         self.project = mainWindow.application.project
#         self.current = mainWindow.application.current
#
#         self._settings = PickAndAssignSettings(parent=self.settingsWidget)
#         self.stackedTableWidget = StackedTableFrameWidget(parent=self.mainWidget,
#                                                           grid=(0, 0), moduleParent=self)
#
#         self.tables = []
#         self._setupTables()
#         self._setupWidgets()
#
#     def _setupWidgets(self):
#         for table in self.tables:
#             # TODO Re-add button functionality
#             # Main widget
#             self.restrictedPickButton = Button(text='Restricted\nPick', callback=None)
#             table.addWidgetToPos(self.restrictedPickButton, row=0, col=2)
#
#             self.assignSelectedButton = Button(text='Assign\nSelected', callback=None)
#             table.addWidgetToPos(self.assignSelectedButton, row=0, col=3)
#
#             self.restrictedPickAndAssignButton = Button(text='Restricted\nPick and Assign', callback=None)
#             table.addWidgetToPos(self.restrictedPickAndAssignButton, row=0, col=4)
#
#             self.restrictedPickButton.setEnabled(True)
#             self.assignSelectedButton.setEnabled(True)
#             self.restrictedPickAndAssignButton.setEnabled(True)
#
#     def _setupTables(self):
#         self.nmrChainTable = NmrResidueTableFrame(parent=self.stackedTableWidget, mainWindow=self.mainWindow,
#                                                   moduleParent=self, grid=(0, 0))
#         self.peakTable = _PeakTableFrame(parent=self.stackedTableWidget, mainWindow=self.mainWindow,
#                                          moduleParent=self, grid=(0, 0))
#
#         self.nmrChainTable.nmrResidueTableSettings = self._settings.nmrResidueTableSettings
#         self.peakTable._settings = self._settings.peakTableSettings
#
#         self.tables = [self.nmrChainTable, self.peakTable]
#         self.stackedTableWidget.addTablesToFrame(self.tables)
#
#
#     def restrictedPickAndAssign(self, assign=True):
#         """
#         Takes the selected NmrResidues from current NmrResidues feeds them into restricted pick lib functions
#         and picks peaks for all spectrum displays specified in the settings tab. Pick uses X and Z axes for each
#         spectrumView as centre points with tolerances and the y as the long axis to pick the whole region.
#         """
#
#         if invalidMsg := self._getMsgIfSetupInvalid():
#             showWarning(self._getActionMsg(assign), invalidMsg)
#         else:
#             nmrResidues = self._getNmrResidues()
#             self._doPickAndAssignOnSelectedNmrResidues(nmrResidues, assign)
#
#     def goToPositionInModules(self, nmrResidue=None, row=None, col=None):
#         """Go to the positions defined my NmrAtoms of nmrResidue in the active displays"""
#
#         nmrResidue = self.project.getByPid(nmrResidue) if isinstance(nmrResidue, str) else nmrResidue
#
#         activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()
#
#         with undoBlockWithoutSideBar():
#
#             if nmrResidue is not None:
#                 mainWindow = self.application.ui.mainWindow
#                 mainWindow.clearMarks()
#                 for display in activeDisplays:
#                     strip = display.strips[0]
#                     n = len(strip.axisCodes)
#                     if n == 2:
#                         widths = ['default', 'default']
#                     else:
#                         widths = ['default', 'full'] + (n - 2) * ['']
#
#                     StripLib.navigateToNmrAtomsInStrip(strip=strip,
#                                                        nmrAtoms=nmrResidue.nmrAtoms,
#                                                        widths=strip._getCurrentZoomRatio(strip.viewRange()),
#                                                        markPositions=(n == 2))
#                 self.current.nmrResidue = nmrResidue


