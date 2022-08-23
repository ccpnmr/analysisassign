"""

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (https://www.ccpn.ac.uk) 2014 - 2022"
__credits__ = ("Ed Brooksbank, Joanna Fox, Victoria A Higman, Luca Mureddu, Eliza Płoskoń",
               "Timothy J Ragan, Brian O Smith, Gary S Thompson & Geerten W Vuister")
__licence__ = ("CCPN licence. See https://ccpn.ac.uk/software/licensing/")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Luca Mureddu $"
__dateModified__ = "$dateModified: 2022-06-29 11:57:45 +0100 (Wed, June 29, 2022) $"
__version__ = "$Revision: 3.1.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Harrison Wang $"
__date__ = "$Date: 2022-05-20 12:59:02 +0100 (Fri, May 20, 2022) $"
#=========================================================================================
# Start of code
#=========================================================================================

from ccpn.AnalysisAssign.modules.backboneExtensions.BackboneAssignmentExtensionABC import BackboneAssignmentExtensionFrame

from ccpn.ui.gui.widgets.RadioButtons import RadioButtons
from ccpn.ui.gui.widgets.Label import Label

class ExampleBackboneExtensionFrame(BackboneAssignmentExtensionFrame):
    """
    Base class for Gui Extension panels.
    """

    NAME = 'Simple Example'

    def __init__(self, guiModule, *args, **Framekwargs):
        BackboneAssignmentExtensionFrame.__init__(self, guiModule, **Framekwargs)


    def registerNotifiers(self):
        pass

    def initWidgets(self):

        row = 0
        l = Label(self, text='This is an example', grid=(row,0))
        r = RadioButtons(self, texts=['Test1', 'Test2'], grid=(row,1))


    def onInstall(self):
        pass

    def updatePanel(self, *args, **kwargs):
        pass

    def close(self):
        """ de-register anything left or close table etc"""
        pass

## Register the Extension in the BackboneAssignmentModule
from ccpn.AnalysisAssign.modules.BackboneAssignmentModule import BackboneAssignmentModule
BackboneAssignmentModule.registerExtension(BackboneAssignmentModule, ExampleBackboneExtensionFrame)