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
__dateModified__ = "$dateModified: 2025-03-14 12:34:36 +0000 (Fri, March 14, 2025) $"
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
from statistics import mean
from collections import defaultdict
# from icecream import ic

from ccpn.core.Peak import Peak
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.PeakList import PeakList
from ccpn.core.lib.AssignmentLib import copyAssignmentsFromReference, propagateAssignments, copyAssignments
from ccpn.ui.gui.lib import PeakListLib
from ccpn.ui.gui.lib import StripLib
from ccpn.ui.gui.lib.StripLib import navigateToNmrAtomsInStrip
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
from ccpn.ui.gui.widgets.Tabs import Tabs
from ccpn.ui.gui.widgets.Widget import Widget
from ccpn.util.OrderedSet import OrderedSet
from ccpn.util.Logging import getLogger
from ccpnmodel.ccpncore.lib.Io.PyMMLibPDB import KEYWDS


logger = getLogger()

ALL = '<Use all>'
SelectToAdd = '> select-to-add <'

# ------------------------------------ BbAssign settings ------------------------------------ #

assignIsotope = '13C'  # Isotope of dimension to be assigned
rootIsotope = '1H'  # Isotope of one of the root dimensions which has already been assigned

iSpectra = ['H[N[CA]]', 'H[N[ca[CO]]]', 'H[N[{CA|ca[Cali]}]]']

exptTypeFilter = ['H[N[CA]]', 'H[N[co[CA]]]', 'H[N[ca[CO]]]', 'H[N[CO]]',
                  'H[N[{CA|ca[Cali]}]]', 'h{CA|Cca}coNH', 'H[N[co[{CA|ca[C]}]]]']

# ------------------------------------------------------------------------------------------- #

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

ZEROMARGINS = (0, 0, 0, 0) # l, t, r, b


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
        self.nmrResidueTableSettings = self._settings.nmrResidueTableSettings
        self.tabWidget = Tabs(parent=self.mainWidget, grid=(0, 0), gridSpan=(1, 3))
        self.tabWidget.setContentsMargins(*ZEROMARGINS)

        self.tables = []
        self._setupTables()
        self._setupWidgets()

        # need to feedback to current.nmrResidueTable
        self._registerNotifiers()
        self.tabWidget.setTabClickCallback(self.tabCallback)

    @property
    def currentTable(self):
        """Returns the current table widget."""
        return self.tabWidget.currentWidget()

    def tabCallback(self, data):
        """Callback for changing tabs

        Enables/disables different widgets based on the selected table.
        """
        # this seems inverted because the callback is before the tab change.
        if self.currentTable is self.peakTable:
            self._settings.nmrResidueTableSettings.sequentialStripsWidget.setEnabled(True)
            self._settings.nmrResidueTableSettings.linkToPulldownClass.setEnabled(True)
        else:
            self._settings.nmrResidueTableSettings.sequentialStripsWidget.setEnabled(False)
            self._settings.nmrResidueTableSettings.linkToPulldownClass.setEnabled(False)

    def _setupWidgets(self):
        """Sets up the table widgets

        Handles the additions of buttons and callbacks to the tables
        """
        for table in self.tables:
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

            # ensure all buttons are enabled
            self.restrictedPickButton.setEnabled(True)
            self.assignSelectedButton.setEnabled(True)
            self.restrictedPickAndAssignButton.setEnabled(True)

    def _setupTables(self):
        """Creates the table frames and adds them to the tabs widget

        This also ensures the settings are set correctly for each table.
        """
        self.nmrChainTable = NmrResidueTableFrame(parent=self.mainWidget, mainWindow=self.mainWindow,
                                                  moduleParent=self, grid=(0, 0))
        self.peakTable = _PeakTableFrame(parent=self.mainWidget, mainWindow=self.mainWindow,
                                         moduleParent=self, grid=(0, 0))

        self.tabWidget.addTab(self.nmrChainTable, 'NmrResidue Table')
        self.tabWidget.addTab(self.peakTable, 'Peak Table')

        self.nmrChainTable.nmrResidueTableSettings = self.nmrResidueTableSettings
        self.peakTable._settings = self._settings.peakTableSettings
        self.peakTable._tableWidget.setActionCallback(self.peakTableActionCallback)

        # set existing widgets to false.
        self.peakTable.posUnitPulldownLabel.setEnabled(False)
        self.peakTable.posUnitPulldownLabel.setVisible(False)
        self.peakTable.posUnitPulldown.setEnabled(False)
        self.peakTable.posUnitPulldown.setVisible(False)

        self.tables = [self.nmrChainTable, self.peakTable]

    def peakTableActionCallback(self, selection, lastItem):
        """Use nmrResidueTableSettings to navigate to nmrAtoms in all displays based off peak.
        """
        from ccpn.ui.gui.lib.StripLib import navigateToPositionInStrip, _getCurrentZoomRatio

        try:
            if not (objs := list(lastItem[self.peakTable._tableWidget._OBJECT])):
                return
        except Exception as es:
            getLogger().debug2(f'{self.__class__.__name__}.actionCallback: No selection\n{es}')
            return

        peak = objs[0] if isinstance(objs, (tuple, list)) else objs

        markPositionsBool = self.nmrResidueTableSettings.markPositionsWidget.checkBox.isChecked()

        if self.nmrResidueTableSettings.displaysWidget:
            displays = self.nmrResidueTableSettings.displaysWidget.getDisplays()
        elif self.current.strip:
            displays = [self.current.strip.spectrumDisplay]

        if not displays and self.nmrResidueTableSettings.displaysWidget:
            logger.warning('Undefined display module(s); select in settings first')
            showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
            return

        with undoBlockWithoutSideBar():
            if self.nmrResidueTableSettings.autoClearMarksWidget.checkBox.isChecked():
                self.application.ui.mainWindow.clearMarks()

            for display in displays:
                for strip in display.strips:
                    if ((optDict := self.nmrResidueTableSettings.axisCodeOptionsDict) and
                            (options := optDict.get(f'{display}')) and
                            display.axes):
                        axisMask = [True if num in options else None for num, axis in enumerate(display.axes)]
                    else:
                        axisMask = None

                    flattenedAssignedNmrAtoms = [atom for axis in peak.assignedNmrAtoms
                                                 for atom in axis if atom is not None]

                    navigateToNmrAtomsInStrip(strip,
                                              flattenedAssignedNmrAtoms,
                                              widths=[],
                                              markPositions=markPositionsBool,
                                              axisMask=axisMask
                                              )

    @property
    def automaticBbNmrAtomAssignment(self):
        return self.nmrResidueTableSettings.automaticBbNmrAtomAssignment.isChecked()

    @property
    def glyHasCaSign(self):
        return self.nmrResidueTableSettings.glyHasCaSign.isChecked()

    @property
    def casPosCbsNeg(self):
        radInd = self.nmrResidueTableSettings.casPosCbsNeg.getIndex()
        if radInd == 0:
            return True
        elif radInd == 1:
            return False

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
        if self.nmrResidueTableSettings.spectrumDisplayPulldown and \
                (texts := self.nmrResidueTableSettings.spectrumDisplayPulldown.getTexts()):
            if ALL in texts:
                gids = self.project.spectrumDisplays
            else:
                gids = [self.application.getByGid(gid) for gid in texts if gid not in [ALL, SelectToAdd]]
            return gids

    def _getMsgIfSetupInvalid(self):
        """Returns an error message based on table and project current"""
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

    def _getSelected(self) -> list:
        """Returns current peaks/nmrResidues depending on the current tab"""
        if self.currentTable is self.nmrChainTable:
            return list(self.current.nmrResidues)
        elif self.currentTable is self.peakTable:
            return list(self.current.peaks)

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

    def assignSelected(self):
        """Assign the currently selected peaks/nmrResidues

        For NmrChainTable: current.peaks on the bases of nmrAtoms of current.nmrResidues
        For PeakTable: copy assignments across peaks
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

            if self.automaticBbNmrAtomAssignment:
                self.bbAssignCarbonNmrAtoms(currentPeaks=peaks)

    def _assignSelectedPeaks(self, peaks=None):
        """Unifies assignments across all selected peaks

        :param peaks: Peaks to unify assignments across
        """
        if peaks is None:
            getLogger().warning('No peaks given to assign')
            return

        for peak in self.current.peaks:
            copyAssignmentsFromReference(peaks, peak)

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
        Takes the selected NmrResidues/Peaks from current NmrResidue/Peak feeds them into restricted pick lib functions
        and picks peaks for all spectrum displays specified in the settings tab. Pick uses X and Z axes for each
        spectrumView as centre points with tolerances and the y as the long axis to pick the whole region.
        """
        if invalidMsg := self._getMsgIfSetupInvalid():
            showWarning(self._getActionMsg(assign), invalidMsg)
            return

        if self.currentTable is self.peakTable:
            peaks = self._getSelected()
            self._doPickAndAssignOnSelectedObjs(peaks, assign)
        elif self.currentTable is self.nmrChainTable:
            nmrResidues = self._getSelected()
            self._doPickAndAssignOnSelectedObjs(nmrResidues, assign)

    def _doPickAndAssignOnSelectedObjs(self, objs: list[NmrResidue,] | list[Peak], assign: bool):
        """Picks and assigns based on given objects.

        :param objs:
        :param assign: Whether to assign the picked peaks.
        :return:
        """
        from ccpn.core.lib.ContextManagers import progressHandler

        undoStack = self.application._getUndo()
        originalUndoState = undoStack.undoList

        with undoBlockWithoutSideBar():
            msg = "Picking and Assigning Peaks..." if assign else "Picking peaks..."
            stopButtonText = 'Stop Pick and Assign' if assign else "Stop Picking"

            curPeaks = set()
            if self.currentTable is self.nmrChainTable:
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
                            self._assignSelectedPeaks(peaks)
                        if isinstance(obj, NmrResidue):
                            self._assignSelectedResidues(peaks, [obj, ])

                        if self.automaticBbNmrAtomAssignment:
                            self.bbAssignCarbonNmrAtoms(currentPeaks=peaks)

                    curPeaks |= OrderedSet(peaks)

            self.current.peaks = list(OrderedSet(self.current.peaks) | curPeaks)
            if progress.cancelled:
                while undoStack.undoList != originalUndoState and undoStack.nextIndex > 0:
                    undoStack.undo()

    def _restrictedPeakPickIterator(self, iterObjs: Iterable[NmrResidue] | Iterable[Peak]):
        """For each display do restricted picks on iterObjs using PeakListLib.restrictedPick
        """
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

    def bbAssignCarbonNmrAtoms(self, currentPeaks: list[Peak] | None = None):
        if len(currentPeaks) == 0:
            showWarning('No Peaks selected', 'Please make sure you have selected some peaks with '
                                             'assigned root (NH) resonances.')
            return

        with undoBlockWithoutSideBar():
            pkDict = defaultdict(list)
            GlyCheck = False
            glyCheckDict = {'CA-1': {'shifts': [], 'peaks': []},
                            'CB-1': {'shifts': [], 'peaks': []},
                            'CA0' : {'shifts': [], 'peaks': []},
                            'CB0' : {'shifts': [], 'peaks': []}}
            GSTCheck = False
            gstCheckDict = {'CA-1': {'shifts': [], 'peaks': []},
                            'CB-1': {'shifts': [], 'peaks': []}}

            for peak in currentPeaks:
                rootDim = \
                    [ind for ind, value in enumerate(peak.peakList.spectrum.isotopeCodes) if value == rootIsotope][0]

                # Check peak root dim is assigned
                if peak.assignmentsByDimensions[rootDim]:
                    peakNmrRes = peak.assignmentsByDimensions[rootDim][0].nmrResidue
                    peakNmrChain = peak.assignmentsByDimensions[rootDim][0].nmrResidue.nmrChain
                else:
                    showWarning('Missing Root Assignment', 'Please make sure all your peaks have '
                                                           'their root NH NmrAtoms assigned')
                    return
                # Put peaks into pkDict
                peakExptType = peak.peakList.spectrum.experimentType

                if peakExptType not in exptTypeFilter:
                    getLogger().warning('Spectrum Experiment Type not valid for automatic C/CA/CB NmrAtom '
                                        'assignment...skipping')
                    continue

                if peakExptType is None:
                    showWarning('Missing Experiment Type', 'Please make sure all '
                                                           'your spectra have an Experiment Type '
                                                           'associated with them (use shortcut ET '
                                                           'to set these)')
                    return

                if peakExptType not in pkDict:
                    pkDict[peakExptType] = [peak]
                else:
                    pkDict[peakExptType].append(peak)

            for expt in pkDict:
                if expt in iSpectra:
                    allPeaks = sorted(pkDict[expt], key=lambda x: x.height if x.height else 0, reverse=True)
                    highestPeak = allPeaks[0]
                    lowestPeak = allPeaks[-1]
                for peak in pkDict[expt]:
                    assignDim = getAssignDim(peak)
                    peakNmrRes = peak.assignmentsByDimensions[rootDim][0].nmrResidue
                    peakSeqCode = peakNmrRes.sequenceCode
                    peakShift = peak.ppmPositions[assignDim]

                    if expt == 'H[N[CA]]':
                        if peak == highestPeak:
                            assignNmrAtom(peakSeqCode, atomName='CA', offset=0, pkNmrChain=peakNmrChain, pk=peak,
                                          pkNmrRes=peakNmrRes)
                        else:
                            assignNmrAtom(peakSeqCode, atomName='CA', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                          pkNmrRes=peakNmrRes)
                    elif expt == 'H[N[co[CA]]]':
                        assignNmrAtom(peakSeqCode, atomName='CA', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                      pkNmrRes=peakNmrRes)
                    elif expt == 'H[N[ca[CO]]]':
                        if peak == highestPeak:
                            assignNmrAtom(peakSeqCode, atomName='C', offset=0, pkNmrChain=peakNmrChain, pk=peak,
                                          pkNmrRes=peakNmrRes)
                        else:
                            assignNmrAtom(peakSeqCode, atomName='C', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                          pkNmrRes=peakNmrRes)
                    elif expt == 'H[N[CO]]':
                        assignNmrAtom(peakSeqCode, atomName='C', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                      pkNmrRes=peakNmrRes)
                    elif expt == 'H[N[{CA|ca[Cali]}]]':
                        if peak == highestPeak:
                            if self.casPosCbsNeg:
                                assignNmrAtom(peakSeqCode, atomName='CA', offset=0, pkNmrChain=peakNmrChain, pk=peak,
                                              pkNmrRes=peakNmrRes)
                                storeDataForGlyCheck(peakShift, peak, atomType='CA0', glyCheckDict=glyCheckDict)
                            else:
                                assignNmrAtom(peakSeqCode, atomName='CB', offset=0, pkNmrChain=peakNmrChain, pk=peak,
                                              pkNmrRes=peakNmrRes)
                                storeDataForGlyCheck(peakShift, peak, atomType='CB0', glyCheckDict=glyCheckDict)
                        elif peak == lowestPeak:
                            if self.casPosCbsNeg:
                                assignNmrAtom(peakSeqCode, atomName='CB', offset=0, pkNmrChain=peakNmrChain, pk=peak,
                                              pkNmrRes=peakNmrRes)
                                storeDataForGlyCheck(peakShift, peak, atomType='CB0', glyCheckDict=glyCheckDict)
                            else:
                                assignNmrAtom(peakSeqCode, atomName='CA', offset=0, pkNmrChain=peakNmrChain, pk=peak,
                                              pkNmrRes=peakNmrRes)
                                storeDataForGlyCheck(peakShift, peak, atomType='CA0', glyCheckDict=glyCheckDict)
                        elif peak.height > 0:
                            if self.casPosCbsNeg:
                                assignNmrAtom(peakSeqCode, atomName='CA', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                              pkNmrRes=peakNmrRes)
                                storeDataForGlyCheck(peakShift, peak, atomType='CA-1', glyCheckDict=glyCheckDict)
                            else:
                                assignNmrAtom(peakSeqCode, atomName='CB', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                              pkNmrRes=peakNmrRes)
                                storeDataForGlyCheck(peakShift, peak, atomType='CB-1', glyCheckDict=glyCheckDict)
                        elif peak.height < 0:
                            if self.casPosCbsNeg:
                                assignNmrAtom(peakSeqCode, atomName='CB', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                              pkNmrRes=peakNmrRes)
                                storeDataForGlyCheck(peakShift, peak, atomType='CB-1', glyCheckDict=glyCheckDict)
                            else:
                                assignNmrAtom(peakSeqCode, atomName='CA', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                              pkNmrRes=peakNmrRes)
                                storeDataForGlyCheck(peakShift, peak, atomType='CA-1', glyCheckDict=glyCheckDict)
                        GlyCheck = True
                    elif expt == 'H[N[co[{CA|ca[C]}]]]' or expt == 'h{CA|Cca}coNH':
                        if peakShift >= 47.0:
                            assignNmrAtom(peakSeqCode, atomName='CA', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                          pkNmrRes=peakNmrRes)
                            storeDataForGSTCheck(peakShift, peak, atomType='CA-1', checkDict=gstCheckDict)
                        elif peakShift < 47.0:
                            assignNmrAtom(peakSeqCode, atomName='CB', offset=-1, pkNmrChain=peakNmrChain, pk=peak,
                                          pkNmrRes=peakNmrRes)
                            storeDataForGSTCheck(peakShift, peak, atomType='CB-1', checkDict=gstCheckDict)
                        GSTCheck = True

            if GlyCheck:
                checkForGly(glyCheckDict, self.glyHasCaSign)
            if GSTCheck:
                checkForGST(gstCheckDict)


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


def getAssignDim(peak):
    return [ind for ind, value in enumerate(peak.peakList.spectrum.isotopeCodes) if value == assignIsotope][0]


def getAssignAxisCode(peak):
    assignDim = getAssignDim(peak)
    assignAxCde = peak.peakList.spectrum.axisCodes[assignDim]
    return assignAxCde


def assignNmrAtom(seqCode, atomName, offset, pkNmrChain, pk, pkNmrRes):
    if offset == -1:
        newsc = ''.join((seqCode, '-1'))
        newnr = pkNmrChain.fetchNmrResidue(sequenceCode=newsc, residueType=None)
        newna = newnr.fetchNmrAtom(name=atomName, isotopeCode=assignIsotope)
        assignAxCde = getAssignAxisCode(pk)
        pk.assignDimension(axisCode=assignAxCde, value=newna)
    elif offset == 0:
        newna = pkNmrRes.fetchNmrAtom(name=atomName, isotopeCode=assignIsotope)
        assignAxCde = getAssignAxisCode(pk)
        pk.assignDimension(axisCode=assignAxCde, value=newna)


def storeDataForGlyCheck(peakShift, peak, atomType, glyCheckDict):
    glyCheckDict[atomType]['shifts'].append(peakShift)
    glyCheckDict[atomType]['peaks'].append(peak)


def checkForGly(glyDict, hasCaSign):
    cas_1 = glyDict['CA-1']['shifts']
    cbs_1 = glyDict['CB-1']['shifts']
    cas0 = glyDict['CA0']['shifts']
    cbs0 = glyDict['CB0']['shifts']
    if hasCaSign:
        if len(cbs_1) == 0 and len(cas0) >= 1:
            if 48.5 > mean(cas0) > 40.0:
                # this is an i Glycine
                for pk in glyDict['CB0']['peaks']:
                    assignDim = getAssignDim(pk)
                    nr = pk.assignmentsByDimensions[assignDim][0].nmrResidue.getOffsetNmrResidue(-1)
                    na = nr.fetchNmrAtom(name='CB', isotopeCode=assignIsotope)
                    assignAxCde = getAssignAxisCode(pk)
                    pk.assignDimension(axisCode=assignAxCde, value=na)
    else:
        if len(cas_1) == 0 and len(cas0) >= 1:
            if 48.5 > mean(cbs_1) > 40.0:
                # this is an i-1 Glycine
                for pk in glyDict['CB-1']['peaks']:
                    assignDim = getAssignDim(pk)
                    nr = pk.assignmentsByDimensions[assignDim][0].nmrResidue.mainNmrResidue
                    na = nr.fetchNmrAtom(name='CA', isotopeCode=assignIsotope)
                    assignAxCde = getAssignAxisCode(pk)
                    pk.assignDimension(axisCode=assignAxCde, value=na)
        elif len(cas_1) == 0 and len(cas0) >= 1:
            if 48.5 > mean(cbs0) > 40.0:
                # this is an i Glycine
                for pk in glyDict['CB0']['peaks']:
                    assignDim = getAssignDim(pk)
                    nr = pk.assignmentsByDimensions[assignDim][0].nmrResidue
                    na = nr.fetchNmrAtom(name='CA', isotopeCode=assignIsotope)
                    assignAxCde = getAssignAxisCode(pk)
                    pk.assignDimension(axisCode=assignAxCde, value=na)
                for pk in glyDict['CA0']['peaks']:
                    assignDim = getAssignDim(pk)
                    nr = pk.assignmentsByDimensions[assignDim][0].nmrResidue.getOffsetNmrResidue(-1)
                    na = nr.fetchNmrAtom(name='CA', isotopeCode=assignIsotope)
                    assignAxCde = getAssignAxisCode(pk)
                    pk.assignDimension(axisCode=assignAxCde, value=na)


def storeDataForGSTCheck(peakShift, peak, atomType, checkDict):
    checkDict[atomType]['shifts'].append(peakShift)
    checkDict[atomType]['peaks'].append(peak)


def checkForGST(gstDict):
    cas = gstDict['CA-1']['shifts']
    cbs = gstDict['CB-1']['shifts']
    peaks = gstDict['CA-1']['peaks'] + gstDict['CB-1']['peaks']
    if len(cas) == 0 and len(cbs) != 0:
        # this is a Glycine
        for pk in peaks:
            assignDim = getAssignDim(pk)
            nr = pk.assignmentsByDimensions[assignDim][0].nmrResidue
            if 48.5 > pk.ppmPositions[assignDim] > 40.0:
                na = nr.fetchNmrAtom(name='CA', isotopeCode=assignIsotope)
                assignAxCde = getAssignAxisCode(pk)
                pk.assignDimension(axisCode=assignAxCde, value=na)
            na = nr.fetchNmrAtom(name='CB', isotopeCode=assignIsotope)
            if not na.assignedPeaks:
                na.delete()
    elif len(cbs) == 0 and len(cas) >= 2:
        # this is a Serine or Threonine
        for pk in peaks:
            assignDim = getAssignDim(pk)
            if pk.ppmPositions[assignDim] > mean(cas):
                nr = pk.assignmentsByDimensions[assignDim][0].nmrResidue
                na = nr.fetchNmrAtom(name='CB', isotopeCode=assignIsotope)
                assignAxCde = getAssignAxisCode(pk)
                pk.assignDimension(axisCode=assignAxCde, value=na)
