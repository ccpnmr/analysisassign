"""
AnalysisAssign Program
"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2022"
__credits__ = ("Ed Brooksbank, Joanna Fox, Victoria A Higman, Luca Mureddu, Eliza Płoskoń",
               "Timothy J Ragan, Brian O Smith, Gary S Thompson & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Luca Mureddu $"
__dateModified__ = "$dateModified: 2022-01-28 14:50:06 +0000 (Fri, January 28, 2022) $"
__version__ = "$Revision: 3.0.4 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

from ccpn.framework.Framework import Framework
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets import MessageDialog
from ccpn.util.Logging import getLogger
from ccpn.util.decorators import logCommand


class Assign(Framework):
    """Root class for Assign application"""

    def __init__(self, applicationName, applicationVersion, commandLineArguments):
        Framework.__init__(self, applicationName, applicationVersion, commandLineArguments)
        # self.components.add('Assignment')

    def _setupMenus(self):
        super()._setupMenus()
        menuSpec = ('Assign', [("Set up NmrResidues", self.showSetupNmrResiduesPopup, [('shortcut', 'sn')]),
                               ("Pick and Assign", self.showPickAndAssignModule, [('shortcut', 'pa')]),
                               (),
                               ("Backbone Assignment", self.showBackboneAssignmentModule, [('shortcut', 'bb')]),
                               ("Sidechain Assignment", self.showSidechainAssignmentModule, [('shortcut', 'sc'), ('enabled', False)]),
                               (),
                               ("Peak Assigner", self.showPeakAssigner, [('shortcut', 'ap')]),
                               ("NmrAtom Assigner", self.showAtomSelector, [('shortcut', 'an')]),
                               ("Assignment Inspector", self.showAssignmentInspectorModule, [('shortcut', 'ai')]),
                               # ("Residue Information", self.showResidueInformation, [('shortcut', 'ri')]),
                               ])
        self._addApplicationMenuSpec(menuSpec)

        viewMenuItems = [("Sequence Graph", self.showSequenceGraph, [('shortcut', 'sg')]),
                         # ("NmrAtom Assigner", self.showAtomSelector, [('shortcut', 'as')]),
                         ()
                         ]
        self._addApplicationMenuItems('View', viewMenuItems, position=9)

    # overrides superclass
    def _closeExtraWindows(self):

        # remove links to modules when closing them
        for attr in ('sequenceGraph', 'backboneModule', 'sidechainAssignmentModule'):
            if hasattr(self, attr):
                delattr(self, attr)

        Framework._closeExtraWindows(self)

    def showSetupNmrResiduesPopup(self):
        if not self.project.peakLists:
            getLogger().warning('No peaklists in project. Cannot assign peaklists.')
            MessageDialog.showWarning('No peaklists in project.', 'Cannot assign peaklists.')
        else:
            from ccpn.ui.gui.popups.SetupNmrResiduesPopup import SetupNmrResiduesPopup

            popup = SetupNmrResiduesPopup(parent=self.ui.mainWindow, mainWindow=self.ui.mainWindow)
            popup.exec_()

    @logCommand('application.')
    def showPickAndAssignModule(self, position: str = 'bottom', relativeTo: CcpnModule = None):
        """Display the Pick and Assign module.
        """
        from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea
        pickAndAssignModule = PickAndAssignModule(mainWindow=mainWindow)
        mainWindow.moduleArea.addModule(pickAndAssignModule, position=position, relativeTo=relativeTo)
        return pickAndAssignModule

    @logCommand('application.')
    def showBackboneAssignmentModule(self, position: str = 'bottom', relativeTo: CcpnModule = None):
        """Display the Backbone Assignment module.
        """
        from ccpn.AnalysisAssign.modules.BackboneAssignmentModule import BackboneAssignmentModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea
        backboneModule = BackboneAssignmentModule(mainWindow=mainWindow)
        mainWindow.moduleArea.addModule(backboneModule, position=position, relativeTo=relativeTo)
        return backboneModule

    @logCommand('application.')
    def showSidechainAssignmentModule(self, position: str = 'bottom', relativeTo: CcpnModule = None):
        """Display the SideChain module.
        """
        MessageDialog.showWarning('Not implemented',
                                  'Sidechain Assignment Module\n'
                                  'is not implemented yet')

    @logCommand('application.')
    def showPeakAssigner(self, position='bottom', relativeTo=None):
        """Display the Peak Assigner module.
        """
        from ccpn.AnalysisAssign.modules.PeakAssigner import PeakAssigner

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea
        assignmentModule = PeakAssigner(mainWindow=mainWindow)
        mainWindow.moduleArea.addModule(assignmentModule, position=position, relativeTo=relativeTo)
        return assignmentModule

    @logCommand('application.')
    def showAssignmentInspectorModule(self, nmrAtom=None, position: str = 'bottom', relativeTo: CcpnModule = None):
        """Display the Assignment Inspector module.
        """
        from ccpn.AnalysisAssign.modules.AssignmentInspectorModule import AssignmentInspectorModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea
        assignmentInspectorModule = AssignmentInspectorModule(mainWindow=mainWindow, selectFirstItem=True)
        mainWindow.moduleArea.addModule(assignmentInspectorModule, position=position, relativeTo=relativeTo)
        return assignmentInspectorModule

    @logCommand('application.')
    def showSequenceGraph(self, position: str = 'bottom', relativeTo: CcpnModule = None, nmrChain=None):
        """Displays Sequence Graph at the bottom of the screen, relative to another module if nextTo is specified.
        """
        from ccpn.AnalysisAssign.modules.SequenceGraph import SequenceGraphModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea
        sequenceGraphModule = SequenceGraphModule(mainWindow=mainWindow, nmrChain=nmrChain)
        mainWindow.moduleArea.addModule(sequenceGraphModule, position=position, relativeTo=relativeTo)
        return sequenceGraphModule

    @logCommand('application.')
    def showAtomSelector(self, position: str = 'bottom', relativeTo: CcpnModule = None, nmrAtom=None):
        """Displays Atom Selector module.
        """
        from ccpn.AnalysisAssign.modules.NmrAtomAssigner import NmrAtomAssignerModule

        mainWindow = self.ui.mainWindow

        if not relativeTo:
            relativeTo = mainWindow.moduleArea
        nmrAtomAssigner = NmrAtomAssignerModule(mainWindow=mainWindow, nmrAtom=nmrAtom)
        mainWindow.moduleArea.addModule(nmrAtomAssigner, position=position, relativeTo=relativeTo)
        return nmrAtomAssigner

    @logCommand('application.')
    def showPipeline(self, position='bottom', relativeTo=None):
        """Display the Screening pipeLine Module
        """
        from ccpn.pipes import loadedPipes
        from ccpn.ui.gui.modules.PipelineModule import GuiPipeline
        guiPipeline = GuiPipeline(mainWindow=self.ui.mainWindow, pipes=loadedPipes, templates=None)
        self.ui.mainWindow.moduleArea.addModule(guiPipeline, position=position)
        return guiPipeline
