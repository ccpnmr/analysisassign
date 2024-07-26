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
__copyright__ = "Copyright (C) CCPN project (https://www.ccpn.ac.uk) 2014 - 2024"
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
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2024-07-09 11:52:17 +0100 (Tue, July 09, 2024) $"
__version__ = "$Revision: 3.2.5 $"
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
# from icecream import ic

from ccpn.core import Peak
from ccpn.core.NmrResidue import NmrResidue
from ccpn.ui.gui.lib import PeakListLib
from ccpn.ui.gui.lib import StripLib
from ccpn.ui.gui.lib.alignWidgets import alignWidgets
from ccpn.ui.gui.lib.StripLib import getZoomRatio
from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTableModule
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.core.lib.Notifiers import Notifier, CurrentNotifier
from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar
from ccpn.util.OrderedSet import OrderedSet
from ccpn.util.Logging import getLogger


logger = getLogger()


class PickAndAssignModule(NmrResidueTableModule):
    """
    Do a restricted peak pick along the 'y-axis' of (a set of) spectra.
    Use settings to define the spectral displays, the active spectra and the tolerances for peak picking

    # GST is this true anymore?
    This module closely works with the Atom Selector module
    """
    className = 'PickAndAssignModule'

    includeSettingsWidget = True
    maxSettingsState = 2
    settingsPosition = 'top'
    settingsMinimumSizes = (500, 200)

    includePeakLists = False
    includeNmrChains = False
    includeSpectrumTable = True

    includeDisplaySettings = True

    def __init__(self, mainWindow, name='Pick and Assign'):

        super().__init__(mainWindow=mainWindow, name=name, selectFirstItem=True)  # ejb ='Pick And Assign')

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current

        # Main widget
        restrictedPickAndAssignWithAssignFalse = partial(self.restrictedPickAndAssign, assign=False)
        self.restrictedPickButton = Button(text='Restricted\nPick', callback=restrictedPickAndAssignWithAssignFalse)
        self.tableFrame.addWidgetToPos(self.restrictedPickButton, row=0, col=2)

        self.assignSelectedButton = Button(text='Assign\nSelected', callback=self.assignSelected)
        self.tableFrame.addWidgetToPos(self.assignSelectedButton, row=0, col=3)

        restrictedPickAndAssignWithAssignTrue = partial(self.restrictedPickAndAssign, assign=True)
        self.restrictedPickAndAssignButton = Button(text='Restricted\nPick and Assign',
                                                    callback=restrictedPickAndAssignWithAssignTrue)
        self.tableFrame.addWidgetToPos(self.restrictedPickAndAssignButton, row=0, col=4)

        self.restrictedPickButton.setEnabled(True)
        self.assignSelectedButton.setEnabled(True)
        self.restrictedPickAndAssignButton.setEnabled(True)

        # change default-settings inherited from NmrResidueTableModule
        self.nmrResidueTableSettings.sequentialStripsWidget.checkBox.setChecked(False)

        if self.nmrResidueTableSettings.displaysWidget:
            self.nmrResidueTableSettings.displaysWidget.addPulldownItem(0)  # select the <all> option

        self.nmrResidueTableSettings.setLabelText('Navigate to\nDisplay(s)')

        # need to feedback to current.nmrResidueTable
        self._selectOnTableCurrentNmrResiduesNotifier = None
        self._registerNotifiers()

        # # these need to change whenever different spectrumDisplays are selected
        if self.nmrResidueTableSettings.axisCodeOptions:
            self.nmrResidueTableSettings.axisCodeOptions.selectAll()

            # just clear the 'C' axes - this is the usual configuration
            for ii, box in enumerate(self.nmrResidueTableSettings.axisCodeOptions.checkBoxes):
                if box.text().upper().startswith('C'):
                    self.nmrResidueTableSettings.axisCodeOptions.clearIndex(ii)

        # fix the second column to stop extra widgets flickering
        alignWidgets(self.nmrResidueTableSettings, columnScale=1.2)

    def _registerNotifiers(self):
        """
        set up the notifiers
        """
        self._selectOnTableCurrentNmrResiduesNotifier = CurrentNotifier(targetName=NmrResidue._pluralLinkName,
                                                                        callback=self._selectionCallback)

    def _unRegisterNotifiers(self):
        """
        clean up the notifiers
        """
        if self._selectOnTableCurrentNmrResiduesNotifier is not None:
            self._selectOnTableCurrentNmrResiduesNotifier.unRegisterNotifier()

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

    def _closeModule(self):
        """
        Unregister notifiers and close module.
        """
        self._unRegisterNotifiers()
        super()._closeModule()

    def _getDisplay(self):
        """Get the current selected spectrum-display from the pulldown
        """
        if self.nmrResidueTableSettings.spectrumDisplayPulldown and \
                (gid := self.nmrResidueTableSettings.spectrumDisplayPulldown.getText()):
            return self.application.getByGid(gid)

    def _getMsgIfSetupInvalid(self):
        msg = None
        if not self.current.nmrResidues:
            # check that is defined and of the correct type
            msg = 'no NmrResidues selected, please pick one or more NmrResidues'

        if not msg and not self._getDisplay():
            # check the selected display
            msg = 'Undefined display;\nselect display in gearbox settings before proceeding'

        if not msg and not self.nmrResidueTableSettings.axisCodeOptions:
            # check that the settings have been populated correctly
            msg = 'Undefined display;\nselect display in gearbox settings before proceeding'

        return msg

    def _getNmrResidues(self) -> list[NmrResidue]:
        """ get the current selected NmrResidues
        """

        nmrResidues = list(self.current.nmrResidues)

        # GST: not sure if this needed - can current.nmrResidue[s] be a string?
        for i, nmrResidue in enumerate(nmrResidues):
            nmrResidues[i] = (self.project.getByPid(nmrResidue)) if isinstance(nmrResidues, str) else nmrResidue

        return nmrResidues

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
                                # validPeakListViews[plv.peakList] += (plv,)
                                pass
        return validPeakListViews

    def assignSelected(self):
        """Assign current.peaks on the bases of nmrAtoms of current.nmrResidues
        """

        nmrResidues = self._getNmrResidues()

        peaks = self.current.peaks
        if len(peaks) == 0:
            return

        with undoBlockWithoutSideBar():
            self._assignPeaks(peaks, nmrResidues)

    # convert to be an iterator...
    def _assignPeaks(self, peaks, nmrResidues):

        currentAxisCodeIndexes = self.nmrResidueTableSettings.axisCodeOptions.getSelectedIndexes()

        for nmrResidue in nmrResidues:
            shiftDict = {}
            for atom in nmrResidue.nmrAtoms:
                shiftDict[atom.isotopeCode] = []

            for peak in peaks:
                if (spectrum := peak.peakList.spectrum) not in self.nmrResidueTableSettings.spectrumIndex:
                    continue

                shiftList = peak.peakList.spectrum.chemicalShiftList
                for nmrAtom in nmrResidue.nmrAtoms:
                    if nmrAtom.isotopeCode in shiftDict.keys():
                        cShift = shiftList.getChemicalShift(nmrAtom)
                        if cShift:
                            shiftDict[nmrAtom.isotopeCode].append((nmrAtom, cShift.value))

                for ii, isotopeCode in enumerate(spectrum.isotopeCodes):
                    if ii in self.nmrResidueTableSettings.spectrumIndex[spectrum]:
                        _restrictedIdx = self.nmrResidueTableSettings.spectrumIndex[spectrum].index(ii)
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

        if invalidMsg := self._getMsgIfSetupInvalid():
            showWarning(self._getActionMsg(assign), invalidMsg)
        else:
            nmrResidues = self._getNmrResidues()
            self._doPickAndAssignOnSelectedNmrResidues(nmrResidues, assign)

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
                        self._assignPeaks(peaks, [nmrResidue, ])
                    curPeaks |= OrderedSet(peaks)

            self.current.peaks = list(OrderedSet(self.current.peaks) | curPeaks)
            if progress.cancelled:
                while undoStack.undoList != originalUndoState and undoStack.nextIndex > 0:
                    undoStack.undo()

    def _restrictedPeakPickIterator(self, nmrResidues: Iterable[NmrResidue]) \
            -> Iterator[tuple[int | None, NmrResidue, str | None, list[Peak] | None]]:

        currentAxisCodeIndexes = self.nmrResidueTableSettings.axisCodeOptions.getSelectedIndexes()

        displays = [self._getDisplay()]
        validPeakListViews = self._getValidPeakListViews(displays)

        badAxisCodeMsg = """\
            Cannot pick some or peaks all peaks; check selected spectrumDisplay
            possibly missing axis-codes or one of the selected nmrResidues has no matching axis-codes
        """

        try:

            specAxisCodes = [[spectrum.axisCodes[self.nmrResidueTableSettings.spectrumIndex[spectrum].index(ii)]
                              for ii in currentAxisCodeIndexes
                              if ii in self.nmrResidueTableSettings.spectrumIndex[spectrum]]
                             for spectrum, peakListView in validPeakListViews.values()]
        except Exception:
            # TODO: this should be a DataClass or named tuple for clarity,,,
            return None, None, badAxisCodeMsg, None

        for i, nmrResidue in enumerate(nmrResidues):

            peaks = []
            try:
                for (spectrum, peakListView), axisCodes in zip(validPeakListViews.values(), specAxisCodes):

                    # axis-codes should be valid at this point
                    peakList, pks = PeakListLib.restrictedPick(peakListView=peakListView,
                                                               axisCodes=axisCodes, nmrResidue=nmrResidue)
                    if pks:
                        peaks += list(pks)

            except Exception:
                return None, nmrResidue, badAxisCodeMsg, None

            yield i, nmrResidue, None, list(peaks)

    def goToPositionInModules(self, nmrResidue=None, row=None, col=None):
        """Go to the positions defined my NmrAtoms of nmrResidue in the active displays"""

        nmrResidue = self.project.getByPid(nmrResidue) if isinstance(nmrResidue, str) else nmrResidue

        activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()

        with undoBlockWithoutSideBar():

            if nmrResidue is not None:
                mainWindow = self.application.ui.mainWindow
                mainWindow.clearMarks()
                for display in activeDisplays:
                    strip = display.strips[0]
                    n = len(strip.axisCodes)
                    if n == 2:
                        widths = ['default', 'default']
                    else:
                        widths = ['default', 'full'] + (n - 2) * ['']

                    StripLib.navigateToNmrAtomsInStrip(strip=strip,
                                                       nmrAtoms=nmrResidue.nmrAtoms,
                                                       widths=getZoomRatio(strip.viewRange()),
                                                       markPositions=(n == 2))
                self.current.nmrResidue = nmrResidue
