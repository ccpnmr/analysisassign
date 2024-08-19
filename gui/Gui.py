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
__dateModified__ = "$dateModified: 2024-08-19 15:20:00 +0100 (Mon, August 19, 2024) $"
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
from ccpn.ui.gui.menus.MenuDefs import VIEW_MENU, MACRO_MENU, VIEW_CHEMICAL_SHIFT_MAPPING, \
    Menu, Action, Separator, _projectHasSpectra, _projectHasPeaks


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
        Action("Set up NmrResidues", app.showSetupNmrResiduesPopup, shortcut = 'sn', checkEnabled=_projectHasSpectra),
        Action("Pick and Assign", app.showPickAndAssignModule, shortcut = 'pa', checkEnabled=_projectHasSpectra),

        Separator(),
        Action("Backbone Assignment", app.showBackboneAssignmentModule, shortcut = 'bb', checkEnabled=_projectHasSpectra),
        # Action("Sidechain Assignment", self.showSidechainAssignmentModule, shortcut = 'sc', enabled = False),

        Separator(),
        Action("Peak Assigner", app.showPeakAssigner, shortcut = 'ap', checkEnabled=_projectHasPeaks),
        Action("NmrAtom Assigner", app.showAtomSelector, shortcut = 'an', checkEnabled=_projectHasPeaks),
        Action("Assignment Inspector", app.showAssignmentInspectorModule, shortcut = 'ai', checkEnabled=_projectHasPeaks),
        # Action("Residue Information", app.showResidueInformation, shortcut = 'ri'),
)  # end _assignMenu

        # put it before the MACRO_MENU
        menuDefs.insertBefore([MACRO_MENU], menuDef=_assignMenu)

        # Add sequence graph to VIEW menu, before CHEMICAL_SHIFT_MAPPING
        _seqGraphMenu = Action("Sequence Graph", app.showSequenceGraph, shortcut = 'sg')
        menuDefs.insertBefore([VIEW_MENU, VIEW_CHEMICAL_SHIFT_MAPPING], menuDef=_seqGraphMenu)

        return menuDefs


