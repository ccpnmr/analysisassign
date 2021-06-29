"""Module Documentation here
"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2021"
__credits__ = ("Ed Brooksbank, Joanna Fox, Victoria A Higman, Luca Mureddu, Eliza Płoskoń",
               "Timothy J Ragan, Brian O Smith, Gary S Thompson & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2021-06-29 14:27:29 +0100 (Tue, June 29, 2021) $"
__version__ = "$Revision: 3.0.4 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Geerten Vuister $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import ccpn.util.isotopes as Isotopes
from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule
from ccpn.ui.gui.lib.SpectrumDisplay import makeStripPlot, makeStripPlotFromSingles
from ccpn.core.lib.Notifiers import Notifier

from ccpn.ui.gui.lib.Strip import matchAxesAndNmrAtoms


class SideChainAssignmentModule(PickAndAssignModule):

    className = 'SideChainAssignmentModule'

    def __init__(self, mainWindow, name='Sidechain Assignment'):

        super().__init__(mainWindow=mainWindow, name=name)

        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current

        self._notifier = self.setNotifier(self.project,
                                          [Notifier.RENAME, Notifier.CREATE, Notifier.CHANGE, Notifier.DELETE],
                                          targetName='NmrAtom', callback=self._updateModules,
                                          )

        # self.refreshButton.show()
        # self.refreshButton.setCallback(self._startAssignment)
        # self.spectrumSelectionWidget.refreshBox.setCallback(self._mediateRefresh)
        # self.nmrResidueTable.nmrResidueTable.setTableCallback(self._startAssignment)
        # self.mode = 'pairs'

    def _mediateRefresh(self):
        """
        Activate/de-activate notifiers depending on the state of the auto-refresh checkbox
        in the spectrum selection widget.
        """
        checked = self.spectrumSelectionWidget.refreshBox.isChecked()
        self._notifier.setBlanking(not checked)

    def _updateModules(self, nmrAtom):
        """
        Convenience function called by notifiers to refresh strip plots when an NmrAtom is created, deleted,
        modified or rename. Calls _startAssignment as to carry out changes.
        """
        if not nmrAtom.nmrResidue is self.current.nmrResidue:
            return
        else:
            self._startAssignment()

    # def _closeModule(self):
    #   super()._closeModule()

    def _startAssignment(self):
        self.mainWindow.clearMarks()
        if self.mode == 'singles':
            self._startAssignmentFromSingles()
        elif self.mode == 'pairs':
            self._startAssignmentFromPairs()

    def _startAssignmentFromPairs(self):
        from ccpn.core.lib.AssignmentLib import getBoundNmrAtomPairs
        activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()

        for display in activeDisplays:
            axisCodes = display.strips[0].axisCodes
            nmrAtomPairs = getBoundNmrAtomPairs(self.current.nmrResidue.nmrAtoms, axisCodes[-1][0])
            displayIsotopeCodes = [Isotopes.name2IsotopeCode(code) for code in axisCodes]
            pairsToRemove = []
            for nmrAtomPair in nmrAtomPairs:
                pairIsotopeCodes = [nap.isotopeCode for nap in nmrAtomPair]
                nmrAtoms = set()
                if (displayIsotopeCodes[1] in pairIsotopeCodes
                        and displayIsotopeCodes[0] not in pairIsotopeCodes):
                    pairsToRemove.append(nmrAtomPair)
                    nmrAtoms.add(nmrAtomPair[0])
                    nmrAtoms.add(nmrAtomPair[1])
                if not all(x.isotopeCode in displayIsotopeCodes for x in nmrAtomPair):
                    pairsToRemove.append(nmrAtomPair)
                    nmrAtoms.add(nmrAtomPair[0])
                    nmrAtoms.add(nmrAtomPair[1])
                elif nmrAtomPair[0].isotopeCode == nmrAtomPair[1].isotopeCode and not \
                        any(displayIsotopeCodes.count(x) > 1 for x in displayIsotopeCodes):
                    pairsToRemove.append(nmrAtomPair)
                    nmrAtoms.add(nmrAtomPair[0])
                    nmrAtoms.add(nmrAtomPair[1])
                if len(displayIsotopeCodes) > 2:
                    if nmrAtomPair[0].isotopeCode == nmrAtomPair[1].isotopeCode and displayIsotopeCodes[0] != \
                            displayIsotopeCodes[2]:
                        if displayIsotopeCodes.count(nmrAtomPair[0].isotopeCode) != 2:
                            nmrAtoms.add(nmrAtomPair[0])
                            nmrAtoms.add(nmrAtomPair[1])
                            pairsToRemove.append(nmrAtomPair)
            for pair in pairsToRemove:
                if pair in nmrAtomPairs:
                    nmrAtomPairs.remove(pair)
            if len(nmrAtomPairs) > 1:
                sortedNmrAtomPairs = self.sortNmrAtomPairs(nmrAtomPairs)

            else:
                sortedNmrAtomPairs = nmrAtomPairs

            if len(display.strips[0].axisCodes) > 2:
                makeStripPlot(display, sortedNmrAtomPairs, autoWidth=False)
            nmrAtoms = [x for y in nmrAtomPairs for x in y]
            axisCodePositionDict = matchAxesAndNmrAtoms(display.strips[0], nmrAtoms)
            self.mainWindow.markPositions(self.project, list(axisCodePositionDict.keys()), list(axisCodePositionDict.values()))

    def sortNmrAtomPairs(self, nmrAtomPairs):
        """
        Sorts pairs of NmrAtoms into 'greek' order. Used in _startAssignmentFromPairs to pass correctly
        ordered lists to makeStripPlot so strips are in the correct order.
        """
        order = ['C', 'CA', 'CB', 'CG', 'CG1', 'CG2', 'CD1', 'CD2', 'CE', 'CZ', 'N', 'ND', 'NE', 'NZ', 'NH',
                 'H', 'HA', 'HB', 'HG', 'HD', 'HE', 'HZ', 'HH']
        ordering = []
        for p in nmrAtomPairs:
            if p[0].name[:len(p[0].name)] in order:
                ordering.append((order.index(p[0].name[:len(p[0].name)]), p))

        if len(nmrAtomPairs) > 1:
            sortedNmrAtomPairs = [x[1] for x in sorted(ordering, key=lambda x: x[0])]
        else:
            sortedNmrAtomPairs = nmrAtomPairs
        return sortedNmrAtomPairs

    def _startAssignmentFromSingles(self):
        activeDisplays = self.spectrumSelectionWidget.getActiveDisplays()

        for display in activeDisplays:
            axisCodes = display.strips[0].axisCodes
            nmrAtoms = set()
            displayIsotopeCodes = [Isotopes.name2IsotopeCode(code) for code in axisCodes]

            for nmrAtom in self.current.nmrResidue.nmrAtoms:
                if nmrAtom.isotopeCode in displayIsotopeCodes and nmrAtom.isotopeCode == displayIsotopeCodes[2]:
                    nmrAtoms.add(nmrAtom)

            makeStripPlotFromSingles(display, list(nmrAtoms))
            axisCodePositionDict = matchAxesAndNmrAtoms(display.strips[0], list(nmrAtoms))
            self.mainWindow.markPositions(self.project, list(axisCodePositionDict.keys()), list(axisCodePositionDict.values()))

            display.setColumnStretches(stretchValue=True)

