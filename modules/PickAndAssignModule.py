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
__dateModified__ = "$dateModified: 2025-04-03 16:12:04 +0100 (Thu, April 03, 2025) $"
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
from typing import Iterable

from statistics import mean, StatisticsError
from collections import defaultdict

from OpenGL.logs import getLog

from ccpn.core.Peak import Peak
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.PeakList import PeakList
from ccpn.core.lib.AssignmentLib import propagateAssignmentsFromReference
from ccpn.ui.gui.lib import PeakListLib

from ccpn.ui.gui.lib.StripLib import navigateToPositionInStrip
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTableFrame
from ccpn.ui.gui.modules.PeakTable import _PeakTableFrame
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar, logCommandManager, progressHandler, \
    notificationEchoBlocking
from ccpn.ui.gui.widgets.SettingsWidgets import PickAndAssignSettings
from ccpn.ui.gui.widgets.Tabs import Tabs
from ccpn.util.OrderedSet import OrderedSet
from ccpn.util.Logging import getLogger


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
        """Initialise both tables and settings to default states.
        """
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
        self._tableButtons = {}
        self._setupTables()
        self._setupWidgets()

        # need to feedback to current.nmrResidueTable
        self._registerNotifiers()
        self.tabWidget.setTabClickCallback(self.tabCallback)

    @property
    def currentTable(self):
        """Returns the current table widget."""
        return self.tabWidget.currentWidget()

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
        restrictedPickAndAssignWithAssignFalse = partial(self.restrictedPickAndAssign, assign=False)
        restrictedPickAndAssignWithAssignTrue = partial(self.restrictedPickAndAssign, assign=True)
        for table in self.tables:

            restrictedPickButton = Button(text='Restricted\nPick',
                                          callback=restrictedPickAndAssignWithAssignFalse)
            assignSelectedButton = Button(text='Assign\nSelected',
                                          callback=self.assignSelected)
            restrictedPickAndAssignButton = Button(text='Restricted\nPick and Assign',
                                                        callback=restrictedPickAndAssignWithAssignTrue)

            self._tableButtons[table] = [restrictedPickButton, assignSelectedButton, restrictedPickAndAssignButton]

            table.addWidgetToPos(restrictedPickButton, row=0, col=2)
            table.addWidgetToPos(assignSelectedButton, row=0, col=3)
            table.addWidgetToPos(restrictedPickAndAssignButton, row=0, col=4)

            # ensure all buttons are enabled
            restrictedPickButton.setEnabled(True)
            assignSelectedButton.setEnabled(True)
            restrictedPickAndAssignButton.setEnabled(True)

    def _setupTables(self):
        """Creates the table frames and adds them to the tabs widget

        This also ensures the settings are set correctly for each table.
        """
        self.nmrChainTable = NmrResidueTableFrame(parent=self.mainWidget, mainWindow=self.mainWindow,
                                                  moduleParent=self, grid=(0, 0), selectFirstItem=True)
        self.peakTable = _PeakTableFrame(parent=self.mainWidget, mainWindow=self.mainWindow,
                                         moduleParent=self, grid=(0, 0), selectFirstItem=True)

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
        """Navigate to and mark peaks based on the settings widget.
        """

        try:
            if not (objs := list(lastItem[self.peakTable._tableWidget._OBJECT])):
                return
        except Exception as es:
            getLogger().debug2(f'{self.__class__.__name__}.actionCallback: No selection\n{es}')
            return

        peak = objs[0] if isinstance(objs, (tuple, list)) else objs

        markPositionsBool = self.nmrResidueTableSettings.markPositionsWidget.checkBox.isChecked()
        clearMarksBool = self.nmrResidueTableSettings.autoClearMarksWidget.checkBox.isChecked()

        if self.nmrResidueTableSettings.displaysWidget:
            displays = self.nmrResidueTableSettings.displaysWidget.getDisplays()
        elif self.current.strip:
            displays = [self.current.strip.spectrumDisplay]

        if not displays and self.nmrResidueTableSettings.displaysWidget:
            logger.warning('Undefined display module(s); select in settings first')
            showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
            return

        with undoBlockWithoutSideBar():
            if clearMarksBool:
                self.application.ui.mainWindow.clearMarks()

            for display in displays:
                for strip in display.strips:
                    navigateToPositionInStrip(strip=strip,
                                              positions=peak.position,
                                              axisCodes=peak.axisCodes,
                                              markPositions=markPositionsBool)

    def _registerNotifiers(self):
        """
        set up the notifiers
        """
        self.setNotifier(self.current,
                         [Notifier.CURRENT],
                         targetName=NmrResidue._pluralLinkName,
                         callback=self._selectionCallback)

        self.setNotifier(self.current,
                         [Notifier.CURRENT],
                         targetName=Peak._pluralLinkName,
                         callback=self._selectionCallback)

    def _selectionCallback(self, data):
        """enable/disable the pick buttons
        """
        if self.currentTable is (table := self.peakTable):
            selected = data[Notifier.OBJECT].peak
        elif self.currentTable is (table := self.nmrChainTable):
            selected = data[Notifier.OBJECT].nmrResidue
        else:
            return

        if selected:
            self._tableButtons[table][0].setEnabled(True)
            self._tableButtons[table][1].setEnabled(True)
            self._tableButtons[table][2].setEnabled(True)
        else:
            self._tableButtons[table][0].setEnabled(False)
            self._tableButtons[table][1].setEnabled(False)
            self._tableButtons[table][2].setEnabled(False)

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
        with logCommandManager(f'{self.__class__.__name__}', funcName='assignSelected'):
            with undoBlockWithoutSideBar():
                if self.nmrResidueTableSettings.nmrChainPeakListRadioButton.getIndex() == 1:
                    peakLists = self._settings.peakListPulldownTexts
                    assignees = [peak for peakList in peakLists for peak in peakList.peaks]
                    references = [peak for peak in self.current.peaks if peak in self.peakTable.table.peaks]
                    if assignees:
                        for reference in references:
                            self._assignSelectedPeaks(assignees, reference)
                elif self.currentTable is self.peakTable:
                    # split out based on assigning table
                    assignees = [peak for peak in peaks if peak not in self.peakTable.table.peaks]
                    references = [peak for peak in self.current.peaks if peak in self.peakTable.table.peaks]
                    if assignees:
                        for reference in references:
                            self._assignSelectedPeaks(assignees, reference)
                    else:
                        # if all peaks are from the same table.
                        self._assignSelectedPeaks(peaks[:-1], peaks[-1])
                elif self.currentTable is self.nmrChainTable:
                    nmrResidues = self._getSelected()
                    self._assignSelectedResidues(peaks, nmrResidues)

                if self.automaticBbNmrAtomAssignment:
                    if self.checkDisplayForExptType():
                        self.bbAssignCarbonNmrAtoms(currentPeaks=peaks)

    @staticmethod
    def _assignSelectedPeaks(peaks: list[Peak] = None, refPeak: Peak = None):
        """Assign peaks based on

        :param peaks: peaks gain assignments
        :param refPeak: reference peak for the assignment
        """
        if refPeak is None:
            getLogger().warning('No reference peak given')
            return

        if peaks is None:
            getLogger().warning('No peaks given to assign')
            return

        propagateAssignmentsFromReference(peaks=peaks, referencePeak=refPeak,
                                          tolerancesByAxisCode={})

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
        with logCommandManager(f'{self.__class__.__name__}', funcName='restrictedPickAndAssign', assign=assign):
            if self.nmrResidueTableSettings.nmrChainPeakListRadioButton.getIndex() == 1:
                peakLists = self._settings.peakListPulldownTexts
                self.pickFromRootAssignOnPeaks(peakLists=peakLists, assign=assign)
            elif self.currentTable is self.peakTable:
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
        # originalUndoState = undoStack.undoList

        if self.automaticBbNmrAtomAssignment:
            exptTypeValid = self.checkDisplayForExptType()

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
                            self._assignSelectedPeaks(peaks, obj)
                        if isinstance(obj, NmrResidue):
                            self._assignSelectedResidues(peaks, [obj, ])

                        if self.automaticBbNmrAtomAssignment and exptTypeValid:
                            self.bbAssignCarbonNmrAtoms(currentPeaks=peaks)

                    curPeaks |= OrderedSet(peaks)

            self.current.peaks = list(OrderedSet(self.current.peaks) | curPeaks)
            if progress.cancelled:
                # while undoStack.undoList != originalUndoState and undoStack.nextIndex > 0:
                undoStack.undo()
                undoStack.clearRedoItems()
                # self.application.undo()

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

    def pickFromRootAssignOnPeaks(self, peakLists: list[PeakList] = None, assign: bool = False):

        msg = "Picking and Assigning Peaks..." if assign else "Picking peaks..."
        stopButtonText = 'Stop Pick and Assign' if assign else "Stop Picking"

        undoStack = self.application._getUndo()
        with notificationEchoBlocking():
            with progressHandler(text=msg, cancelButtonText=stopButtonText,
                                 maximum=len(peakLists)) as progress:

                with undoBlockWithoutSideBar():
                    for i, peakList in enumerate(peakLists):
                        progress.checkCancel()
                        progress.setValue(i)

                        for peak in self.current.peaks:
                            _positionCodeDict = dict(zip(peak.axisCodes, peak.position))
                            peaks = peakList.restrictedPick(positionCodeDict=_positionCodeDict, doPos=True, doNeg=False)

                            if assign and peaks:
                                self._assignSelectedPeaks(peaks, peak)

            if progress.cancelled:
                undoStack.undo()
                undoStack.clearRedoItems()


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
                    if self.currentTable is self.nmrChainTable:
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

    def checkDisplayForExptType(self, enableWarning : bool = True) -> bool:
        displays = self.nmrResidueTableSettings.displaysWidget.getDisplays()

        validExptFromDisplay = [specView.spectrum.experimentType for display in displays
                                for specView in display.spectrumViews
                                if specView.spectrum.experimentType in exptTypeFilter]

        if enableWarning and not validExptFromDisplay:
            showWarning('Automatic C/CA/CB NmrAtom Assigment',
                        'Spectrum Experiment Types not valid for automatic C/CA/CB NmrAtom assignment..\n'
                        'Skipping automatic C/CA/CB assignment...')

        return bool(validExptFromDisplay)


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
    try:
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
            if len(cas_1) == 0 and len(cbs_1) >= 1:
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
    except StatisticsError as e:
        getLogger().warning(e)


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
            try:
                if pk.ppmPositions[assignDim] > mean(cas):
                    nr = pk.assignmentsByDimensions[assignDim][0].nmrResidue
                    na = nr.fetchNmrAtom(name='CB', isotopeCode=assignIsotope)
                    assignAxCde = getAssignAxisCode(pk)
                    pk.assignDimension(axisCode=assignAxCde, value=na)
            except StatisticsError as e:
                getLogger().warning(e)
