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
__modifiedBy__ = "$modifiedBy: Geerten Vuister $"
__dateModified__ = "$dateModified: 2022-02-01 15:30:12 +0000 (Tue, February 01, 2022) $"
__version__ = "$Revision: 3.0.4 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: TJ $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import os
import sys
from PyQt5 import QtGui, QtWidgets

from ccpn.framework.Application import defineProgramArguments
from ccpn.AnalysisAssign.AnalysisAssign import Assign as Application


def main():

    # argument parser
    parser = defineProgramArguments()

    # add any additional commandline argument here
    commandLineArguments = parser.parse_args()

    # viewportFormat = QtGui.QSurfaceFormat()
    # viewportFormat.setSwapInterval(0)  #disable VSync - this works here!
    # QtGui.QSurfaceFormat().setDefaultFormat(viewportFormat)
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--enable-logging --log-level=3"

    application = Application(commandLineArguments)

    application._startApplication()
    QtWidgets.QApplication.quit()

    if sys.platform.startswith('win'):
        os._exit(0)


if __name__ == '__main__':
    main()
