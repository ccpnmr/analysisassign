"""
AnalysisAssign Gui
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
__dateModified__ = "$dateModified: 2024-02-09 12:14:30 +0000 (Fri, February 09, 2024) $"
__version__ = "$Revision: 3.2.2 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: gvuister $"
__date__ = "$Date: 2024-02-09 10:28:40 +0000 (Fri, Feb 09, 2024) $"
#=========================================================================================
# Start of code
#=========================================================================================

from ccpn.ui.gui.Gui import Gui
from ccpn.ui.gui.Menus import MACRO_MENU
from ccpn.ui.gui.widgets import MessageDialog
from ccpn.util.Logging import getLogger
from ccpn.util.decorators import logCommand


class AnalysisAssignGui(Gui):
    """Extended Gui interface for AnalysisAssign
    """

    def _getMenuDefs(self):
        """:return the MenuDefs instance; modified with AnalysisAssign menu additions
        """
        from ccpn.ui.gui.Menus import VIEW_MENU
        menuDefs = super()._getMenuDefs()

        app = self.application
        menuDef = ('Assign',  [("Set up NmrResidues", app.showSetupNmrResiduesPopup, [('shortcut', 'sn')]),
                               ("Pick and Assign", app.showPickAndAssignModule, [('shortcut', 'pa')]),
                               (),
                               ("Backbone Assignment", app.showBackboneAssignmentModule, [('shortcut', 'bb')]),
                               # ("Sidechain Assignment", self.showSidechainAssignmentModule, [('shortcut', 'sc'), ('enabled', False)]),
                               (),
                               ("Peak Assigner", app.showPeakAssigner, [('shortcut', 'ap')]),
                               ("NmrAtom Assigner", app.showAtomSelector, [('shortcut', 'an')]),
                               ("Assignment Inspector", app.showAssignmentInspectorModule, [('shortcut', 'ai')]),
                               # ("Residue Information", app.showResidueInformation, [('shortcut', 'ri')]),
                               ])
        # put it before the MACRO_MENU
        _position = menuDefs._getMenuIndex(MACRO_MENU)
        menuDefs._addMenuDef(menuDef, position=_position)

        viewMenuItems = [("Sequence Graph", app.showSequenceGraph, [('shortcut', 'sg')]),
                        ]
        menuDefs._addMenuItems(VIEW_MENU, viewMenuItems, position=11)

        return menuDefs


