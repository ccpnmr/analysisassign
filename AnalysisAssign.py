"""
AnalysisAssign Program
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
__dateModified__ = "$dateModified: 2024-06-26 14:52:13 +0100 (Wed, June 26, 2024) $"
__version__ = "$Revision: 3.2.4 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

from ccpn.framework.Version import applicationVersion
from ccpn.framework.Framework import Framework
from ccpn.framework.Application import ANALYSIS_ASSIGN

from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets import MessageDialog
from ccpn.util.Logging import getLogger
from ccpn.util.decorators import logCommand


class Assign(Framework):
    """Root class for AnalysisAssign application
    """
    applicationName = ANALYSIS_ASSIGN
    applicationVersion = applicationVersion

    def _getUI(self):
        """Get the user interface
        :return a Ui instance
        """
        if self.args.interface == 'Gui':
            from ccpn.AnalysisAssign.gui.Gui import AnalysisAssignGui
            ui = AnalysisAssignGui(application=self)

        else:
            from ccpn.ui.Ui import NoUi
            ui = NoUi(application=self)

        return ui

    # GWV 9/2/24: Now in Gui.py
    # def _setupMenus(self):
    #     """Augment the menu's
    #     :return the MenuDefs (i.e. a list) instance
    #     """
    #     from ccpn.ui.gui.Menus import getMenuDefs, VIEW_MENU
    #
    #     menuDefs = getMenuDefs()
    #
    #     assignDef = ('Assign', [("Set up NmrResidues", self.showSetupNmrResiduesPopup, [('shortcut', 'sn')]),
    #                            ("Pick and Assign", self.showPickAndAssignModule, [('shortcut', 'pa')]),
    #                            (),
    #                            ("Backbone Assignment", self.showBackboneAssignmentModule, [('shortcut', 'bb')]),
    #                            # ("Sidechain Assignment", self.showSidechainAssignmentModule, [('shortcut', 'sc'), ('enabled', False)]),
    #                            (),
    #                            ("Peak Assigner", self.showPeakAssigner, [('shortcut', 'ap')]),
    #                            ("NmrAtom Assigner", self.showAtomSelector, [('shortcut', 'an')]),
    #                            ("Assignment Inspector", self.showAssignmentInspectorModule, [('shortcut', 'ai')]),
    #                            # ("Residue Information", self.showResidueInformation, [('shortcut', 'ri')]),
    #                            ])
    #     menuDefs._addMenuDef(assignDef, position=4)
    #
    #     viewMenuItems = [("Sequence Graph", self.showSequenceGraph, [('shortcut', 'sg')]),
    #                     ]
    #     menuDefs._addMenuItems(VIEW_MENU, viewMenuItems, position=11)
    #     return menuDefs

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
        mainWindow._addModule(pickAndAssignModule, position=position, relativeTo=relativeTo)
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
        mainWindow._addModule(backboneModule, position=position, relativeTo=relativeTo)
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
        mainWindow._addModule(assignmentModule, position=position, relativeTo=relativeTo)
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
        mainWindow._addModule(assignmentInspectorModule, position=position, relativeTo=relativeTo)
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
        mainWindow._addModule(sequenceGraphModule, position=position, relativeTo=relativeTo)
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
        mainWindow._addModule(nmrAtomAssigner, position=position, relativeTo=relativeTo)
        return nmrAtomAssigner

    @logCommand('application.')
    def showPipeline(self, position='bottom', relativeTo=None):
        """Display the Screening pipeLine Module
        """
        from ccpn.pipes import loadedPipes
        from ccpn.ui.gui.modules.PipelineModule import GuiPipeline
        guiPipeline = GuiPipeline(mainWindow=self.ui.mainWindow, pipes=loadedPipes, templates=None)
        self.ui.mainWindow._addModule(guiPipeline, position=position)
        return guiPipeline

    def propagateAssignments(self):
        """Propagate assignments across selected peaks.
        """
        self.ui.mainWindow.propagateAssignments()

    def copyAssignments(self):
        """Copy assignments across selected peaks.
        """
        self.ui.mainWindow.copyAssignments()
