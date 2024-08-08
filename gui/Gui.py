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
__dateModified__ = "$dateModified: 2024-08-08 15:43:19 +0100 (Thu, August 08, 2024) $"
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
from ccpn.ui.gui.Menus import VIEW_MENU, MACRO_MENU, VIEW_CHEMICAL_SHIFT_MAPPING, \
    Separator, options, _projectHasSpectra, _projectHasPeaks
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
        menuDef = ('Assign',  [("Set up NmrResidues", app.showSetupNmrResiduesPopup, options(shortcut = 'sn'), _projectHasSpectra),
                               ("Pick and Assign", app.showPickAndAssignModule, options(shortcut = 'pa'), _projectHasSpectra),

                               Separator(),
                               ("Backbone Assignment", app.showBackboneAssignmentModule, options(shortcut = 'bb'), _projectHasSpectra),
                               # ("Sidechain Assignment", self.showSidechainAssignmentModule, options(shortcut = 'sc', enabled = False)),

                               Separator(),
                               ("Peak Assigner", app.showPeakAssigner, options(shortcut = 'ap'), _projectHasPeaks),
                               ("NmrAtom Assigner", app.showAtomSelector, options(shortcut = 'an'), _projectHasPeaks),
                               ("Assignment Inspector", app.showAssignmentInspectorModule, options(shortcut = 'ai'), _projectHasPeaks),
                               # ("Residue Information", app.showResidueInformation, options(shortcut = 'ri')),
                               ]
                   )

        # put it before the MACRO_MENU
        menuDefs.insertBefore([MACRO_MENU], menuDef=menuDef)

        # Add sequence graph to VIEW menu
        _seqGraphMenu = ("Sequence Graph", app.showSequenceGraph, options(shortcut = 'sg'))
        menuDefs.insertBefore([VIEW_MENU, VIEW_CHEMICAL_SHIFT_MAPPING], menuDef=_seqGraphMenu)

        return menuDefs


