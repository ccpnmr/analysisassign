"""
Do a restricted peak pick along the 'y-axis' of (a set of) spectra.
Use settings to define the spectral displays, the active spectra and the tolerances for peak picking

This module closely works with the Atom Selector module

First version by SS
Refactored by GWV to be responsible to active state on start; proper callbacks 
and to include "Restricted pick and assign" button.

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (https://www.ccpn.ac.uk) 2014 - 2024"
__credits__ = ("Ed Brooksbank, Joanna Fox, Morgan Hayward, Victoria A Higman, Luca Mureddu",
               "Eliza Płoskoń, Timothy J Ragan, Brian O Smith, Gary S Thompson & Geerten W Vuister")
__licence__ = ("CCPN licence. See https://ccpn.ac.uk/software/licensing/")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, https://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Geerten Vuister $"
__dateModified__ = "$dateModified: 2024-03-20 13:39:43 +0000 (Wed, March 20, 2024) $"
__version__ = "$Revision: 3.2.2 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Geerten Vuister $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

from ccpn.ui.gui.lib import PeakListLib
from ccpn.ui.gui.lib import StripLib
from ccpn.ui.gui.lib.alignWidgets import alignWidgets
from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTableModule
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.core.lib.Notifiers import Notifier, CurrentNotifier
from ccpn.core.NmrResidue import NmrResidue
from ccpn.util.Logging import getLogger
from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar


logger = getLogger()


class PickAndAssignModule(NmrResidueTableModule):
    """
    Do a restricted peak pick along the 'y-axis' of (a set of) spectra.
    Use settings to define the spectral displays, the active spectra and the tolerances for peak picking
  
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
        self.restrictedPickButton = Button(text='Restricted\nPick', callback=self.restrictedPick, )
        self.tableFrame.addWidgetToPos(self.restrictedPickButton, row=0, col=2)

        self.assignSelectedButton = Button(text='Assign\nSelected', callback=self.assignSelected)
        self.tableFrame.addWidgetToPos(self.assignSelectedButton, row=0, col=3)

        self.restrictedPickAndAssignButton = Button(text='Restricted\nPick and Assign', callback=self.restrictedPickAndAssign)
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
        # self._setAxisCodes()
        if self.nmrResidueTableSettings.axisCodeOptions:
            self.nmrResidueTableSettings.axisCodeOptions.selectAll()

            # just clear the 'C' axes - this is the usual configuration
            for ii, box in enumerate(self.nmrResidueTableSettings.axisCodeOptions.checkBoxes):
                if box.text().upper().startswith('C'):
                    self.nmrResidueTableSettings.axisCodeOptions.clearIndex(ii)

        # fix the second column to stop extra widgets flickering
        alignWidgets(self.nmrResidueTableSettings, columnScale=1.2)
        # if self.nmrResidueTableSettings.displaysWidget:
        #     alignWidgets(self.nmrResidueTableSettings.displaysWidget)

    def _registerNotifiers(self):
        """
        set up the notifiers
        """
        self._selectOnTableCurrentNmrResiduesNotifier = CurrentNotifier(
                                                                 targetName=NmrResidue._pluralLinkName,
                                                                 callback=self._selectionCallback)

    def _unRegisterNotifiers(self):
        """
        clean up the notifiers
        """
        if self._selectOnTableCurrentNmrResiduesNotifier is not None:
            self._selectOnTableCurrentNmrResiduesNotifier.unRegister()

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

    def _verify(self, msgHeader, nmrResidue):
        """Verify that the settings are valid
        """
        nmrResidue = self.project.getByPid(nmrResidue) if isinstance(nmrResidue, str) else nmrResidue
        if not nmrResidue:
            # use current if not set
            nmrResidue = self.application.current.nmrResidue
        if nmrResidue is None or not isinstance(nmrResidue, NmrResidue):
            # check that is defined and of the correct type
            showWarning(msgHeader, 'Undefined nmrResidue; select one first before proceeding')
            return

        if not self._getDisplay():
            # check the selected display
            showWarning(msgHeader, 'Undefined display;\nselect display in gearbox settings before proceeding')
            return

        if not self.nmrResidueTableSettings.axisCodeOptions:
            # check that the settings have been populated correctly
            showWarning(msgHeader, 'Undefined display;\nselect display in gearbox settings before proceeding')
            return

        return nmrResidue

    @staticmethod
    def _getValidPeakListViews(displays):
        """Get tehist of valid peakListViews
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

    def assignSelected(self, nmrResidue=None, msgHeader=None):
        """Assign current.peaks on the bases of nmrAtoms of current.nmrResidue
        """
        msgHeader = msgHeader or 'Assign Selected'
        if not (nmrResidue := self._verify(msgHeader, nmrResidue)):
            return

        if len(self.application.current.peaks) == 0:
            showWarning(msgHeader, 'Undefined peak(s); select one or more before proceeding')
            return

        with undoBlockWithoutSideBar():
            lastNmrResidue = nmrResidue  # self.application.current.nmrResidue
            currentAxisCodeIndexes = self.nmrResidueTableSettings.axisCodeOptions.getSelectedIndexes()

            shiftDict = {}
            for atom in self.application.current.nmrResidue.nmrAtoms:
                shiftDict[atom.isotopeCode] = []

            for peak in self.application.current.peaks:
                shiftList = peak.peakList.spectrum.chemicalShiftList
                spectrum = peak.peakList.spectrum

                for nmrAtom in self.application.current.nmrResidue.nmrAtoms:
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
                            # pValue = peak.position[ii]
                            if abs(sValue - pValue) <= spectrum.assignmentTolerances[ii]:
                                # peak.assignDimension(spectrum.axisCodes[ii], [shift[0]])
                                shiftList.add(shift[0])

                        if shiftList:
                            peak.assignDimension(spectrum.axisCodes[ii], list(shiftList))

            # self.application.current.peaks = []
            # update the NmrResidue table
            self.tableWidget._table = self.application.current.nmrResidue.nmrChain
            self.tableWidget._update()

            # reset to the last selected nmrResidue - stops other tables messing up
            self.application.current.nmrResidue = lastNmrResidue

    #TODO:GEERTEN: compact the two routines
    def restrictedPick(self, nmrResidue=None, msgHeader=None):
        """
        Routine refactored in revision 9381.
     
        Takes an NmrResidue feeds it into restricted pick lib functions and picks peaks for all
        spectrum displays specified in the settings tab. Pick uses X and Z axes for each spectrumView as
        centre points with tolerances and the y as the long axis to pick the whole region.
        """
        msgHeader = msgHeader or 'Restricted Pick'
        if not (nmrResidue := self._verify(msgHeader, nmrResidue)):
            return

        currentAxisCodeIndexes = self.nmrResidueTableSettings.axisCodeOptions.getSelectedIndexes()

        with undoBlockWithoutSideBar():
            peaks = []

            # displays = self._getDisplays()
            # gid = self.nmrResidueTableSettings.spectrumDisplayPulldown.getText()
            displays = [self._getDisplay()]  # [self.application.getByGid(gid)]

            validPeakListViews = self._getValidPeakListViews(displays)
            try:
                specAxisCodes = [[spectrum.axisCodes[self.nmrResidueTableSettings.spectrumIndex[spectrum].index(ii)]
                                  for ii in currentAxisCodeIndexes if ii in self.nmrResidueTableSettings.spectrumIndex[spectrum]]
                                 for spectrum, peakListView in validPeakListViews.values()]
            except Exception:
                showWarning(msgHeader, 'Cannot pick peaks; check selected spectrumDisplay,\n'
                                       'possibly missing axis-codes or selected nmrResidue has no matching axis-codes')

            else:
                try:
                    for (spectrum, peakListView), axisCodes in zip(validPeakListViews.values(), specAxisCodes):
                        # axisCodes = [spectrum.axisCodes[self.nmrResidueTableSettings.spectrumIndex[spectrum].index(ii)]
                        #              for ii in currentAxisCodeIndexes if ii in self.nmrResidueTableSettings.spectrumIndex[spectrum]]

                        # axis-codes should be valid this time
                        peakList, pks = PeakListLib.restrictedPick(peakListView=peakListView,
                                                                   axisCodes=axisCodes, nmrResidue=nmrResidue)
                        if pks:
                            peaks += list(pks)
                except Exception:
                    showWarning(msgHeader, 'Cannot pick peaks; check selected spectrumDisplay,\n'
                                           'possibly missing axis-codes or selected nmrResidue has no matching axis-codes')

                # for module in self.application.project.spectrumDisplays:
                #     if len(module.axisCodes) >= 2:
                #         for spectrumView in module.strips[0].spectrumViews:
                #
                #             visiblePeakListViews = [peakListView for peakListView in spectrumView.peakListViews
                #                                     if peakListView.isVisible()]
                #
                #             # if len(visiblePeakListViews) == 0:
                #             #     continue
                #             # else:
                #             #     peakList, pks = PeakList.restrictedPick(peakListView=visiblePeakListViews[0],
                #             #                                             axisCodes=module.axisCodes[0::2], nmrResidue=nmrResidue)
                #             #     peaks = peaks + pks
                #
                #             # if len(visiblePeakListViews) == 0:
                #             #     spectrum = spectrumView.spectrum
                #             #
                #             #     axisCodes = [axis for ]
                else:
                    # set the current peaks - may need intermediate list here
                    self.application.current.peaks = peaks

                    # update the NmrResidue table
                    self.tableWidget._table = nmrResidue.nmrChain
                    self.tableWidget._update()

                    return True  # pick was successful

    # from ccpn.util.decorators import profile
    # @profile
    def restrictedPickAndAssign(self, nmrResidue=None):
        """
        Functionality for beta2 to include the Assign part
         
        Takes an NmrResidue feeds it into restricted pick lib functions and picks peaks for all
        spectrum displays specified in the settings tab. Pick uses X and Z axes for each spectrumView as
        centre points with tolerances and the y as the long axis to pick the whole region.
        
        Calls assignSelected to assign
        """
        msgHeader = 'Restricted Pick and Assign'
        if not (nmrResidue := self._verify(msgHeader, nmrResidue)):
            return

        with undoBlockWithoutSideBar():

            if self.restrictedPick(nmrResidue, msgHeader=msgHeader) and self.application.current.peaks:
                # if peaks have been selected then assign them
                self.assignSelected(msgHeader=msgHeader)

                # notifier for other modules
                nmrResidue._finaliseAction('change')

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
                    #Strip.navigateToNmrAtomsInStrip(strip=strip, nmrAtoms=nmrResidue.nmrAtoms, widths=['default', 'full', ''])
                    if n == 2:
                        widths = ['default', 'default']
                    else:
                        widths = ['default', 'full'] + (n - 2) * ['']

                    StripLib.navigateToNmrAtomsInStrip(strip=strip,
                                                       nmrAtoms=nmrResidue.nmrAtoms,
                                                       widths=strip._getCurrentZoomRatio(strip.viewRange()),
                                                       markPositions=(n == 2))
                self.application.current.nmrResidue = nmrResidue
