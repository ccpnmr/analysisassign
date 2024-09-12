"""
AnalysisAssign Gui
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
__modifiedBy__ = "$modifiedBy: Geerten Vuister $"
__dateModified__ = "$dateModified: 2024-09-12 08:50:01 +0100 (Thu, September 12, 2024) $"
__version__ = "$Revision: 3.2.5 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: gvuister $"
__date__ = "$Date: 2024-02-09 10:28:40 +0000 (Fri, Feb 09, 2024) $"
#=========================================================================================
# Start of code
#=========================================================================================

from ccpn.ui.gui.Gui import Gui
from ccpn.ui.gui.menus.MenuDefs import \
    VIEW_MENU, MACRO_MENU, VIEW_CHEMICAL_SHIFT_MAPPING, \
    _projectHasSpectra, _projectHasPeaks
from ccpn.ui.gui.menus._MenuItems import Menu, Action, Separator

from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets import MessageDialog

from ccpn.util.Logging import getLogger
from ccpn.util.decorators import logCommand


class AnalysisAssignGui(Gui):
    """Extended Gui interface for AnalysisAssign
    """

    def _getMenuDefs(self):
        """:return the MenuDefs instance; modified with AnalysisAssign menu additions
        """
        menuDefs = super()._getMenuDefs()

        app = self.application
        _assignMenu = \
    Menu("Assign",

        Action("Set up NmrResidues", self.showSetupNmrResiduesPopup, shortcut = 'sn', checkEnabled=_projectHasSpectra),
        Action("Pick and Assign", self.showPickAndAssignModule, shortcut = 'pa', checkEnabled=_projectHasSpectra),

        Separator(),
        Action("Backbone Assignment", self.showBackboneAssignmentModule, shortcut = 'bb', checkEnabled=_projectHasSpectra),
        # Action("Sidechain Assignment", self.showSidechainAssignmentModule, shortcut = 'sc', enabled = False),

        Separator(),
        Action("Peak Assigner", self.showPeakAssigner, shortcut = 'ap', checkEnabled=_projectHasPeaks),
        Action("NmrAtom Assigner", self.showAtomSelector, shortcut = 'an', checkEnabled=_projectHasPeaks),
        Action("Assignment Inspector", self.showAssignmentInspectorModule, shortcut = 'ai', checkEnabled=_projectHasPeaks),
        # Action("Residue Information", self.showResidueInformation, shortcut = 'ri'),

    )  # end menu Assign

        # put it before the MACRO_MENU
        menuDefs.insertBefore([MACRO_MENU], menuDef=_assignMenu)

        # Add sequence graph to VIEW menu, before CHEMICAL_SHIFT_MAPPING
        _seqGraphMenu = Action("Sequence Graph", self.showSequenceGraph, shortcut = 'sg')
        menuDefs.insertBefore([VIEW_MENU, VIEW_CHEMICAL_SHIFT_MAPPING], menuDef=_seqGraphMenu)

        return menuDefs

    def showSetupNmrResiduesPopup(self):
        if not self.project.peakLists:
            getLogger().warning('No peaklists in project. Cannot assign peaklists.')
            MessageDialog.showWarning('No peaklists in project.', 'Cannot assign peaklists.')
        else:
            from ccpn.ui.gui.popups.SetupNmrResiduesPopup import SetupNmrResiduesPopup
            popup = SetupNmrResiduesPopup(parent=self.mainWindow, mainWindow=self.mainWindow)
            popup.exec_()

    @logCommand('ui.')
    def showPickAndAssignModule(self, position: str = 'bottom', relativeTo: CcpnModule = None):
        """Display the Pick and Assign module.
        """
        from ccpn.AnalysisAssign.modules.PickAndAssignModule import PickAndAssignModule

        mainWindow = self.mainWindow
        pickAndAssignModule = PickAndAssignModule(mainWindow=mainWindow)
        mainWindow._addModule(pickAndAssignModule, position=position, relativeTo=relativeTo)
        return pickAndAssignModule

    @logCommand('ui.')
    def showBackboneAssignmentModule(self, position: str = 'bottom', relativeTo: CcpnModule = None):
        """Display the Backbone Assignment module.
        """
        from ccpn.AnalysisAssign.modules.BackboneAssignmentModule import BackboneAssignmentModule

        mainWindow = self.mainWindow
        backboneModule = BackboneAssignmentModule(mainWindow=mainWindow)
        mainWindow._addModule(backboneModule, position=position, relativeTo=relativeTo)
        return backboneModule

    @logCommand('ui.')
    def showSidechainAssignmentModule(self, position: str = 'bottom', relativeTo: CcpnModule = None):
        """Display the SideChain module.
        """
        MessageDialog.showWarning('Not implemented',
                                  'Sidechain Assignment Module\n'
                                  'is not implemented yet')

    @logCommand('ui.')
    def showPeakAssigner(self, position='bottom', relativeTo=None):
        """Display the Peak Assigner module.
        """
        from ccpn.AnalysisAssign.modules.PeakAssigner import PeakAssigner

        mainWindow = self.mainWindow
        assignmentModule = PeakAssigner(mainWindow=mainWindow)
        mainWindow._addModule(assignmentModule, position=position, relativeTo=relativeTo)
        return assignmentModule

    @logCommand('ui.')
    def showAssignmentInspectorModule(self, nmrAtom=None, position: str = 'bottom', relativeTo: CcpnModule = None):
        """Display the Assignment Inspector module.
        """
        from ccpn.AnalysisAssign.modules.AssignmentInspectorModule import AssignmentInspectorModule

        mainWindow = self.mainWindow
        assignmentInspectorModule = AssignmentInspectorModule(mainWindow=mainWindow, selectFirstItem=True)
        mainWindow._addModule(assignmentInspectorModule, position=position, relativeTo=relativeTo)
        return assignmentInspectorModule

    @logCommand('ui.')
    def showSequenceGraph(self, position: str = 'bottom', relativeTo: CcpnModule = None, nmrChain=None):
        """Displays Sequence Graph at the bottom of the screen, relative to another module if nextTo is specified.
        """
        from ccpn.AnalysisAssign.modules.SequenceGraph import SequenceGraphModule

        mainWindow = self.mainWindow
        sequenceGraphModule = SequenceGraphModule(mainWindow=mainWindow, nmrChain=nmrChain)
        mainWindow._addModule(sequenceGraphModule, position=position, relativeTo=relativeTo)
        return sequenceGraphModule

    @logCommand('ui.')
    def showAtomSelector(self, position: str = 'bottom', relativeTo: CcpnModule = None, nmrAtom=None):
        """Displays Atom Selector module.
        """
        from ccpn.AnalysisAssign.modules.NmrAtomAssigner import NmrAtomAssignerModule

        mainWindow = self.mainWindow
        nmrAtomAssigner = NmrAtomAssignerModule(mainWindow=mainWindow, nmrAtom=nmrAtom)
        mainWindow._addModule(nmrAtomAssigner, position=position, relativeTo=relativeTo)
        return nmrAtomAssigner

    @logCommand('ui.')
    def showPipeline(self, position='bottom', relativeTo=None):
        """Display the Screening pipeLine Module
        """
        from ccpn.pipes import loadedPipes
        from ccpn.ui.gui.modules.PipelineModule import GuiPipeline

        mainWindow = self.mainWindow
        guiPipeline = GuiPipeline(mainWindow=mainWindow, pipes=loadedPipes, templates=None)
        mainWindow._addModule(guiPipeline, position=position, relativeTo=relativeTo)
        return guiPipeline
