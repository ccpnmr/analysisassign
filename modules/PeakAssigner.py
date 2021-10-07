"""
Module to assign peaks
Responds to current.peaks

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
__dateModified__ = "$dateModified: 2021-10-05 11:34:11 +0100 (Tue, October 05, 2021) $"
__version__ = "$Revision: 3.0.4 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import typing
import numpy as np
from functools import partial
from collections import OrderedDict
from PyQt5 import QtGui, QtWidgets, QtCore
from ccpn.core.NmrAtom import NmrAtom
from ccpn.core.NmrResidue import NmrResidue, _getNmrResidue
from ccpn.core.Peak import Peak
from ccpn.core.lib import CcpnSorting
from ccpn.core.lib.AssignmentLib import nmrAtomsForPeaks, peaksAreOnLine, sameAxisCodes
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.ButtonList import ButtonList, Button
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.Frame import Frame, ScrollableFrame
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.HLine import HLine, LabeledHLine
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.GuiTable import GuiTable
from ccpn.ui.gui.widgets.Column import ColumnClass
from ccpn.ui.gui.widgets.SpeechBalloon import SpeechBalloon
from ccpn.ui.gui.guiSettings import getColours, DIVIDER, LABEL_WARNINGFOREGROUND
from ccpn.util.Logging import getLogger
from ccpn.util.Common import greekKey, _truncateText, getIsotopeListFromCode
from ccpn.ui.gui.widgets.MessageDialog import showWarning, showYesNo
from ccpnmodel.ccpncore.lib.Constants import defaultNmrChainCode
from ccpn.core.lib.Notifiers import Notifier
from ccpn.ui.gui.widgets.Font import getFontHeight, TABLEFONT
from ccpn.core.lib.ContextManagers import undoBlock, undoBlockWithoutSideBar
from ccpn.ui.gui.widgets.DropBase import DropBase
from ccpn.ui.gui.lib.GuiNotifier import GuiNotifier


allowedResidueTypes = [('', '', ''),
                       ('Alanine', 'ALA', 'A'),
                       ('Arginine', 'ARG', 'R'),
                       ('Asparagine', 'ASN', 'N'),
                       ('Aspartic acid', 'ASP', 'D'),
                       ('ASP/ASN ambiguous', 'ASX', 'B'),
                       ('Cysteine', 'CYS', 'C'),
                       ('Glutamine', 'GLN', 'Q'),
                       ('Glutamic acid', 'GLU', 'E'),
                       ('GLU/GLN ambiguous', 'GLX', 'Z'),
                       ('Glycine', 'GLY', 'G'),
                       ('Histidine', 'HIS', 'H'),
                       ('Isoleucine', 'ILE', 'I'),
                       ('Leucine', 'LEU', 'L'),
                       ('Lysine', 'LYS', 'K'),
                       ('Methionine', 'MET', 'M'),
                       ('Phenylalanine', 'PHE', 'F'),
                       ('Proline', 'PRO', 'P'),
                       ('Serine', 'SER', 'S'),
                       ('Threonine', 'THR', 'T'),
                       ('Tryptophan', 'TRP', 'W'),
                       ('Tyrosine', 'TYR', 'Y'),
                       ('Unknown', 'UNK', ''),
                       ('Valine', 'VAL', 'V')]

MSG = 'Not-defined >Select any to start<'

ROWDEFAULT = 300
ROWSIZES = {7 : 3000,
            4 : 2500,
            0 : 2000,
            -1: 1200,
            -2: ROWDEFAULT,
            }

_showBorders = False  # for debugging of layout's
_margins = (2, 2, 2, 2)
ASSIGNEDROWS = 3
ALTERNATIVEROWS = 5
MINTABLEWIDTH = 150


class PeakAssigner(CcpnModule):
    """Module for assignment of nmrAtoms to the different axes of a peak.
    Module responds to current.peak
    """

    # override in specific module implementations
    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    settingsPosition = 'left'

    className = 'PeakAssigner'

    activePulldownClass = None


    class _emptyObject():
        def __init__(self):
            pass


    def __init__(self, mainWindow, name="Peak Assigner"):
        """
        Initialise the Module widgets
        """
        super().__init__(mainWindow=mainWindow, name=name)

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current

        self.Ndims = 1
        self.maxDims = 8
        self.dimensionTabs = []
        self.currentAtoms = None
        self._visibleDims = 0

        # add widgets to the module
        self._setWidgets()

        # set notifiers to respond to peaks
        self._registerNotifiers()

        # install event filter to track changes in width
        self.mainWidget.installEventFilter(self)

        # populate the tables
        self._updateInterface(self.current.peaks)

        self.installMaximiseEventHandler(self._maximise, self._closeModule)

    def eventFilter(self, target, event):
        """Event filter to handle a mainWidget resizing
        """
        # tables are acting very strange in this module - caused by setStretchLastColumn
        if event.type() == QtCore.QEvent.Resize:
            self._resize(event.size().width())
        return super().eventFilter(target, event)

    def _maximise(self):
        """
        Maximise the attached table
        """
        self.chemicalShiftTable._maximise()

    def _setWidgets(self):
        """Add the widgets to the module
        """
        self.blockSignals(True)

        if self.includeSettingsWidget:
            # settings
            row = 0
            self.doubleToleranceCheckbox = CheckBox(self.settingsWidget, checked=False,
                                                    callback=self._updateInterface,
                                                    grid=(row, 1))
            Label(self.settingsWidget, text="Double Tolerances ", grid=(row, 0))

            row += 1
            self.intraCheckbox = CheckBox(self.settingsWidget, checked=False,
                                          callback=self._updateInterface,
                                          grid=(row, 1))
            Label(self.settingsWidget, text="Only Intra-residual ", grid=(row, 0))

            row += 1
            self.multiCheckbox = CheckBox(self.settingsWidget, checked=True,
                                          callback=self._updateInterface,
                                          grid=(row, 1))
            Label(self.settingsWidget, text="Allow Multiple Peaks ", grid=(row, 0))

            row += 1
            self.allChainCheckBoxLabel = CheckBox(self.settingsWidget, checked=False,
                                                  callback=self._updateInterface,
                                                  grid=(row, 1))
            Label(self.settingsWidget, "Peak Selection from Table", grid=(row, 0))

            self._height = getFontHeight()
            self._tableHeight = getFontHeight(name=TABLEFONT)
            self.settingsWidget.setContentsMargins(5, 5, 5, 5)
            self.settingsWidget.getLayout().setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
            self.settingsWidget.setScrollBarPolicies(scrollBarPolicies=('asNeeded', 'never'))

        row = 0
        # add a label for the selected peaks
        self.peakLabel = Label(parent=self.mainWidget, setLayout=True, spacing=(0, 0),
                               text='Current Peak: ' + MSG, bold=True,
                               grid=(row, 0), margins=_margins,
                               hAlign='left', vAlign='t',
                               hPolicy='ignored', vPolicy='fixed'
                               )

        row += 1
        # setup a frame for the dimension frames - scrollable frame not resizing correctly
        self.axisFrameWidget = ScrollableFrame(parent=self.mainWidget, showBorder=False, setLayout=True,
                                               acceptDrops=True, grid=(row, 0),
                                               )
        # self._axisFrameScrollArea = self.axisFrameWidget._scrollArea
        # row = -1

        # self.Ndims = 0
        # self.maxDims = 4
        # self.dimensionTabs = []
        # self.currentAtoms = None
        # self._visibleDims = 0

        row += 1
        # _frame = Frame(self.axisFrameWidget, grid=(row,0), setLayout=True, showBorder=_showBorders,
        #                # hAlign='left',
        #                # hPolicy='minimumExpanding'
        #                )
        # self.axisFrameWidget.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.MinimumExpanding)
        # self.axisFrameWidget.getLayout().setSizeConstraint(QtWidgets.QLayout.SetFixedSize)

        colIndex = 0
        for dimIndex in range(self.maxDims):
            # 4 axes per row
            if dimIndex == 4:
                row += 1
                colIndex = 0
            dimTab = AxisAssignmentObject(parent=self.axisFrameWidget, grid=(row, colIndex),
                                          parentModule=self, dimIndex=dimIndex,
                                          mainWindow=self.mainWindow,
                                          )
            self.dimensionTabs.append(dimTab)
            colIndex += 1

        self.mainWidget.getLayout().setAlignment(QtCore.Qt.AlignTop)
        self.axisFrameWidget.setVisible(False)

        # # testing
        # self.dimensionTabs[1].setNotAligned(True)
        # self.dimensionTabs[0].setHLineText('H: 10.7')
        # self.dimensionTabs[2].hide()

        self.blockSignals(False)

    def _resize(self, width):
        if self.Ndims:
            try:
                wid = self.axisFrameWidget.scrollArea.verticalScrollBar()
                visible = wid.isVisible()
                offset = wid.width() if visible else 0
            except:
                offset = 0
            w = (width - 6 - offset) / min(self.Ndims, 4)
            for tab in self.dimensionTabs:
                tab.setFixedWidth(max(w, MINTABLEWIDTH))

    def _registerNotifiers(self):
        # without a tableSelection specified in the table callback, this nmrAtom callback is needed
        # to update the table
        self.setNotifier(self.current, [Notifier.CURRENT],
                         targetName=Peak._pluralLinkName,
                         callback=self._updateInterface,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.DELETE, Notifier.CREATE],
                         targetName=Peak.__name__,
                         callback=self._updateInterface,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.CHANGE, Notifier.RENAME, Notifier.CREATE],
                         targetName=NmrAtom.__name__,
                         callback=self._updateNmrAtom,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.CHANGE],
                         targetName=Peak.__name__,
                         callback=self._updateNmrResidue,
                         onceOnly=True)
        self.setNotifier(self.project, [Notifier.DELETE, Notifier.CREATE],
                         targetName=NmrResidue.__name__,
                         callback=self._updateNmrResidue,
                         onceOnly=True)

    def _updateNmrAtom(self, data):
        # not a very efficient way of doing this
        self._updateInterface(data, action=data[Notifier.TRIGGER])

    def _updateNmrResidue(self, data):
        # not a very efficient way of doing this
        self._updateInterface(data, action=data[Notifier.TRIGGER])

    def _updateInterface(self, data=None, action=None):
        """Updates the whole module, including recalculation
           of which nmrAtoms fit to the peaks.
        """
        peaks = self.current.peaks

        if not peaks or not self._peaksAreCompatible(peaks):
            self.axisFrameWidget.hide()
            self.peakLabel.setText('Current Peak: ' + MSG)
        else:

            # update the peaksLabel
            peaksIds = ' , '.join([str(pp.id) for pp in self.current.peaks])
            if len(self.current.peaks) < 2:
                self.peakLabel.setText('Current Peak: %s' % self.current.peak.id)
            else:
                self.peakLabel.setText('Current Peaks: %s' % _truncateText(peaksIds, maxWords=6))
                self.peakLabel.setToolTip(peaksIds)

            self.Ndims = self.current.peak.spectrum.dimensionCount

            # show/hide the dimension tabs
            for dimIndex, dimTab in enumerate(self.dimensionTabs):
                if dimIndex < self.Ndims:
                    dimTab.show()
                else:
                    dimTab.hide()

            for dimIndex, dimTab in enumerate(self.dimensionTabs[:self.Ndims]):
                # show/hide if peaks are aligned
                aligned = peaksAreOnLine(peaks=peaks, dimIndex=dimIndex)
                if not aligned:
                    txt = '%s: peaks\nnot aligned' % (peaks[0].axisCodes[dimIndex],)
                    dimTab.setNotAlignedText(txt)
                else:
                    ppmValues = np.array([pk.ppmPositions[dimIndex] for pk in peaks])
                    txt = '%s: %.3f' % (peaks[0].axisCodes[dimIndex], ppmValues.mean())
                    dimTab.setHLineText(txt)
                dimTab.showNotAligned(not aligned)

            self._updateTables(peaks=peaks)
            self._resize(self.mainWidget.width())
            self.axisFrameWidget.show()

    def _updateTables(self, peaks):
        """
        update Assigned and alternatives tables showing which nmrAtoms
        are assigned to which peak dimensions. If multiple peaks are selected,
        only the assignment that they have in common are shown. Maybe this should be all
        assignments. You can see that at the peak annotation though.
        """

        doubleTolerance = self.doubleToleranceCheckbox.isChecked()
        intraResidual = self.intraCheckbox.isChecked()

        validNmrAtoms = [nmrAtom for nmrAtom in self.project.nmrAtoms if not (nmrAtom.nmrResidue.isDeleted or nmrAtom.nmrResidue._flaggedForDelete)]
        nmrAtomsForTables = nmrAtomsForPeaks(peaks, validNmrAtoms,
                                             doubleTolerance=doubleTolerance,
                                             intraResidual=intraResidual)

        Ndimensions = self.Ndims
        self.currentList = []

        self._tables = [self._emptyObject()] * Ndimensions

        _sizes = []
        for dim, nmrAtoms in zip(range(Ndimensions), nmrAtomsForTables):
            # self.axisTables[dim].show()
            # self.axisDivergeLabels[dim][0].hide()

            ll = [set(peak.dimensionNmrAtoms[dim]) for peak in peaks]
            self.nmrAtoms = list(sorted(set.intersection(*ll)))  # was intersection
            self.nmrAtoms = [nmrAtom for nmrAtom in self.nmrAtoms if not (nmrAtom.nmrResidue.isDeleted or nmrAtom.nmrResidue._flaggedForDelete)]

            self.currentList.append([str(a.pid) for a in self.nmrAtoms])  # ejb - keep another list
            self.dimensionTabs[dim].setAssignedTable(self.nmrAtoms)
            # _rows = self.axisTables[dim].tables[0].rowCount()
            # _rows = self.axisTables[dim].tables[0].sizeHint().height()

            nmrAtomsForTables[dim] = [nmr for nmr in nmrAtomsForTables[dim] if nmr not in self.nmrAtoms]
            self.dimensionTabs[dim].setAlternativesTable(nmrAtomsForTables[dim])
            # if peaksAreOnLine(peaks, dim):
            #     self.axisTables[dim].setAlternativesTable(nmrAtomsForTables[dim])
            #     # _rows = max(_rows, self.axisTables[dim].tables[1].rowCount())
            #     # for k, val in ROWSIZES.items():
            #     #     if _rows > k:
            #     #         _sizes.append(val)
            #     #         break
            #     # else:
            #     #     _sizes.append(ROWDEFAULT)
            # else:
            #     self.axisTables[dim].setAlternativesTable(None)
            #     # _sizes.append(ROWDEFAULT)
            #     #

        return _sizes

    def _getDeltaShift(self, nmrAtom: NmrAtom, dim: int) -> typing.Union[float, str]:
        """
        Calculation of delta shift to add to the table.
        """
        if (not self.current.peaks) or nmrAtom is NOL:
            return ''

        deltas = []
        for peak in self.current.peaks:
            shiftList = peak.peakList.spectrum.chemicalShiftList
            if shiftList:
                shift = shiftList.getChemicalShift(nmrAtom)
                if shift:
                    _value = shift.value
                    if _value is not None:
                        position = peak.position[dim]
                        deltas.append(abs(shift.value - position))
        # average = sum(deltas)/len(deltas) #Bug: ZERO DIVISION!

        if len(deltas) > 0:
            return float(np.mean(deltas))  #'%6.3f' % np.mean(deltas) - handled by table
        else:
            return ''

    def _getShift(self, nmrAtom: NmrAtom) -> typing.Union[float, str]:
        """
        Calculation of chemical shift value to add to the table.
        """
        if (not self.current.peaks) or nmrAtom is NOL:
            return ''

        for peak in self.current.peaks:
            shiftList = peak.peakList.spectrum.chemicalShiftList
            if shiftList:
                shift = shiftList.getChemicalShift(nmrAtom)
                if shift:
                    return shift.value  # '%8.3f' % shift.value

    def _peaksAreCompatible(self, peaks) -> bool:
        """
        If multiple peaks are selected, a check is performed
        to determine whether assignment of corresponding
        dimensions of a peak allowed.
        """
        if len(peaks) == 1:
            return True

        if not self.multiCheckbox.isChecked():
            getLogger().warning("Selecting multiple peaks is currently not active (change in settings pane)")
            return False

        dimensionalities = set(peak.spectrum.dimensionCount for peak in peaks)
        if len(dimensionalities) > 1:
            getLogger().warning('Not all selected peaks have the same number of dimensions')
            return False

        for dimIndex in range(peaks[0].spectrum.dimensionCount):
            isotopeCodes = set(peak.spectrum.isotopeCodes[dimIndex] for peak in peaks)
            if len(isotopeCodes) > 1:
                getLogger().warning('Selected peaks have different isotopeCodes along dimension %d' % dimIndex + 1)
                return False

        return True

    def _emptyAllTablesAndLists(self):
        """
        Quick erase of all present information in ListWidgets and ObjectTables.
        """
        self.peakLabel.setText('Current Peak: ' + MSG)
        for label in self.labels:
            label.setText('')
        for objectTable in self.objectTables:
            objectTable.setObjects([])
        for listWidget in self.listWidgets:
            listWidget.clear()

    def _updatePulldownLists(self, dim: int, row: int = None, col: int = None, obj: object = None):
        objectTable = self.objectTables[dim]
        nmrAtom = objectTable.getCurrentObject()
        self._updateAssignmentWidget(dim, nmrAtom)

    def _closeModule(self):
        """
        CCPN-INTERNAL: used to close the module
        """
        for dimTab in self.dimensionTabs:
            dimTab._close()
        self.dimensionTabs = []
        super()._closeModule()


class NotOnLine(object):
    """
    Small 'fake' object to get a message the user in the assignment
    Table that a specific dimension can not be assigned in one go
    since the frequencies of the peaks in this dimension are not on
    one line (i.e. the C frequencies of the CA and CB in a strip for
    instance).
    """

    def __init__(self):
        self.pid = 'Multiple selected peaks not on line.'
        self.id = 'Multiple selected peaks not on line.'


NOL = NotOnLine()


class AssignmentTable(GuiTable):
    """Subclassed for some added functionality"""

    def __init__(self, parent=None,
                 mainWindow=None,
                 dataFrameObject=None,  # collate into a single object that can be changed quickly
                 actionCallback=None,
                 selectionCallback=None,
                 checkBoxCallback=None,
                 clearSelectionCallback=None,
                 _pulldownKwds=None, enableMouseMoveEvent=True,
                 multiSelect=False, selectRows=True, numberRows=False, autoResize=False,
                 enableExport=True, enableDelete=True, enableSearch=True, allowRowDragAndDrop=True,
                 hideIndex=True, stretchLastSection=True,
                 **kwds):

        super().__init__(parent=parent,
                         mainWindow=mainWindow,
                         dataFrameObject=dataFrameObject,
                         actionCallback=actionCallback,
                         selectionCallback=selectionCallback,
                         checkBoxCallback=checkBoxCallback,
                         _pulldownKwds=_pulldownKwds, enableMouseMoveEvent=enableMouseMoveEvent,
                         allowRowDragAndDrop=allowRowDragAndDrop,
                         multiSelect=multiSelect, selectRows=selectRows,
                         numberRows=numberRows, autoResize=autoResize,
                         enableExport=enableExport, enableDelete=enableDelete, enableSearch=enableSearch,
                         hideIndex=hideIndex, stretchLastSection=stretchLastSection,
                         **kwds)

        self._clearSelectionCallbackFunction = clearSelectionCallback

    def _clearSelectionCallback(self):
        super(AssignmentTable, self)._clearSelectionCallback()
        if self._clearSelectionCallbackFunction:
            data = {}
            self._clearSelectionCallbackFunction(data)

    def _setContextMenu(self, enableExport=True, enableDelete=True):
        """Subclass guiTable to add new item to top of context menu
        """
        super()._setContextMenu(enableExport=enableExport, enableDelete=enableDelete)
        _actions = self.tableMenu.actions()
        if _actions:
            _topMenuItem = _actions[0]
            _topSeparator = self.tableMenu.insertSeparator(_topMenuItem)
            self._editMenuAction = self.tableMenu.addAction('Edit nmrAtom ...', self._editNmrAtom)
            # move new actions to the top of the list
            self.tableMenu.insertAction(_topSeparator, self._editMenuAction)

    def _raiseTableContextMenu(self, pos):
        """Create a new menu and popup at cursor position
        Update text for edit nmrAtom
        """
        selection = self.getSelectedObjects()
        data = self.getRightMouseItem()
        if data and selection:
            # add more information to the edit nmrAtom option in the menu
            currentNmrAtom = selection[0]
            self._editMenuAction.setText('Edit NmrAtom {}'.format(currentNmrAtom.id if currentNmrAtom else '...'))
            self._editMenuAction.setEnabled(True if currentNmrAtom else False)

        else:
            # disabled but visible lets user know that menu items exist
            self._editMenuAction.setText('Edit NmrAtom ...')
            self._editMenuAction.setEnabled(False)

        # hide the previous edit balloon (looks a little cleaner)
        self._owner.setEditPopupVisible(False)
        super()._raiseTableContextMenu(pos)

    def _editNmrAtom(self):
        """Edit the nmrAtom from the parent widget
        """
        # NOTE:ED - should the edit functionality be here or in the AxisAssignmentObject?
        selection = self.getSelectedObjects()
        data = self.getRightMouseItem()
        if data and selection:
            # call the edit popup balloon
            self._owner._reassignNmrAtomPopup()


class EditNmrAtomBalloon(SpeechBalloon):
    """Balloon to hold the pulldown lists for editing the nmrAtom
    """

    def __init__(self, mainWindow=None, *args, **kwds):
        super().__init__(*args, **kwds)

        # simplest way to make the popup function as modal and disappear as required
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.Popup)

    #     # attach handle to focus and maximise/minimise signals from app and parent window
    #     QtWidgets.qApp.focusChanged.connect(self._focusChanged)
    #     if mainWindow:
    #         mainWindow.WindowMaximiseMinimise.connect(self._stateChanged)
    #
    # @QtCore.pyqtSlot('QWidget*', 'QWidget*')
    # def _focusChanged(self, oldWidget, newWidget):
    #     """Hide the balloon if the focus has changed to a widget outside the parent widget
    #     """
    #     # Check whether the focus is still inside the popup
    #     if self.isVisible() and newWidget not in self.findChildren(QtWidgets.QWidget):
    #         self.setVisible(False)
    #         getLogger().debug2(f'Closing _focusChanged')
    #
    # @QtCore.pyqtSlot(bool)
    # def _stateChanged(self, state):
    #     """Hide the balloon if the parent window is minimised
    #     """
    #     # Check whether widget is visible
    #     if state is False and self.isVisible() is True:
    #         self.setVisible(False)
    #         getLogger().debug2(f'Closing _stateChanged')

    # NOTE:ED - test method, attach eventFilter to all children - above method much more flexible
    # def _attachFocusOut(self, mainWindow=None):
    #     """Attach eventFilter to child widgets to capture focusOut event
    #     Used to close parent if lost focus to outside widget
    #     """
    #     children = self.findChildren(QtWidgets.QWidget)
    #     for ch in children:
    #         ch.installEventFilter(self)
    #
    # def eventFilter(self, target, event):
    #     """Event filter to handle a child of a widget losing focus
    #     Calls parent check to see if another child still has focus
    #     """
    #     if event.type() == QtCore.QEvent.FocusOut:
    #         # use singleShot to ensure that all widgets are up-to-date before checking focus from parent
    #         QtCore.QTimer.singleShot(0, self._checkFocus)
    #     return super().eventFilter(target, event)
    #
    # def _checkFocus(self):
    #     """Check whether any of the children still have focus
    #     """
    #     children = self.findChildren(QtWidgets.QWidget)
    #     for ch in children:
    #         if ch.hasFocus():
    #             # editPopup is okay, stay visible
    #             break
    #     else:
    #         # no children have focus, so need to hide - caused by any click outside the editBalloon
    #         self.setVisible(False)
    #         getLogger().debug2(f'Closing nmrAtom editor')


class AxisAssignmentObject(Frame):
    """
    Create a new frame for displaying information in 1 axis of peakassigner
    """

    def __init__(self, parent, parentModule, dimIndex, mainWindow, grid=None, **kwds):

        # settings = dict(hPolicy = 'minimum', hAlign='left', vPolicy = 'expanding', vAlign='top')
        settings = dict(vAlign='top', )

        super(AxisAssignmentObject, self).__init__(parent=parent,
                                                   setLayout=True, showBorder=_showBorders,
                                                   grid=grid, **settings, **kwds
                                                   )

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current
        self.currentAtoms = None
        self._clickedNmrAtom = None

        # initialise axis information
        self.dimIndex = dimIndex
        self._parent = parentModule
        self.dataFrameAssigned = None
        self.dataFrameAlternatives = None
        self.lastTableSelected = None
        self.lastNmrAtomSelected = None
        self.tables = [None, None]  # The two tables (assignment and alternatives)

        height = 20
        # self._minWidth = 150
        _minTabWidth = 100
        _tabHeight = 100
        _pullDownWidth = 65

        aRow = -1  # Toplevel row in the widget
        #=========================================
        # divider line
        #=========================================
        # aRow += 1
        # self.hLine = LabeledHLine(self, text='axis', grid=(aRow,0), height=10, colour=getColours()[DIVIDER])

        #=========================================
        # assignments
        #=========================================
        aRow += 1
        self._assignmentsFrame = Frame(self, setLayout=True, showBorder=_showBorders,
                                       grid=(aRow, 0), margins=_margins, acceptDrops=True, **settings)
        # self._assignmentsFrame.getLayout().setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        # row = -1
        self._parent.setGuiNotifier(self._assignmentsFrame, [GuiNotifier.DROPEVENT], [DropBase.PIDS],
                                    callback=self._handleDropsFromSideBar)

        row = 0
        self.hLine = LabeledHLine(self._assignmentsFrame, text='axis', grid=(row, 0), height=16, colour=getColours()[DIVIDER])

        # row += 1
        # self.axisLabel = Label(self._assignmentsFrame, 'Axis', hAlign='l', grid=(row, 0))
        # self.axisLabel.setMinimumHeight(height)
        row += 1
        self.tables[0] = AssignmentTable(parent=self._assignmentsFrame,
                                         mainWindow=mainWindow,
                                         dataFrameObject=None,
                                         setLayout=True,
                                         autoResize=False, multiSelect=False,
                                         actionCallback=partial(self._assignDeassignNmrAtom, 0),
                                         selectionCallback=partial(self._clickedTableCallback, 0),
                                         clearSelectionCallback=partial(self._clearTableCallback, 0),
                                         grid=(row, 0), gridSpan=(1, 1),
                                         # **settings,
                                         stretchLastSection=True,
                                         enableSearch=False,
                                         acceptDrops=True,
                                         enableExport=False,
                                         tipText='Click to select; double-click to de-assign')
        self.tables[0]._owner = self
        self.tables[0].setFixedHeight((ASSIGNEDROWS + 1) * getFontHeight() * 1.5)
        self._parent.setGuiNotifier(self.tables[0], [GuiNotifier.DROPEVENT], [DropBase.PIDS],
                                    callback=partial(self._handleDroppedItems, 0))
        self._parent.setGuiNotifier(self.tables[0], [GuiNotifier.DRAGMOVEEVENT], [DropBase.PIDS],
                                    callback=partial(self._handleDragMoveEvent, 0))

        row += 1
        self._alternativesLabel = Label(self._assignmentsFrame, 'Alternatives', hAlign='l', grid=(row, 0))
        self._alternativesLabel.setMinimumHeight(height)
        row += 1
        self.tables[1] = AssignmentTable(parent=self._assignmentsFrame,
                                         mainWindow=mainWindow,
                                         dataFrameObject=None,
                                         setLayout=True,
                                         autoResize=False, multiSelect=False,
                                         actionCallback=partial(self._assignDeassignNmrAtom, 1),
                                         selectionCallback=partial(self._clickedTableCallback, 1),
                                         clearSelectionCallback=partial(self._clearTableCallback, 1),
                                         grid=(row, 0), gridSpan=(1, 1),
                                         # **settings,
                                         stretchLastSection=True,
                                         enableSearch=False,
                                         enableExport=False,
                                         acceptDrops=True,
                                         tipText='Click to select; double-click to assign')
        self._parent.setGuiNotifier(self.tables[1], [GuiNotifier.DROPEVENT],
                                    [DropBase.PIDS], callback=partial(self._handleDroppedItems, 1))
        self._parent.setGuiNotifier(self.tables[1], [GuiNotifier.DRAGMOVEEVENT],
                                    [DropBase.PIDS], callback=partial(self._handleDragMoveEvent, 1))
        self.tables[1]._owner = self
        self.tables[1].setFixedHeight((ALTERNATIVEROWS + 1) * getFontHeight() * 1.5)

        # self.tables[1].setMinimumWidth(_minTabWidth)
        # self.tables[1].setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding)
        # row += 1

        # self._assignmentWidget = self._nmrAtomWidget(parent=self._assignmentsFrame, minWidth=_pullDownWidth,
        #                                              setLayout=True, showBorder=_showBorders, grid=(0, 0))

        row += 1
        _buttons = ButtonList(self._assignmentsFrame, texts=['Edit', 'New'],
                              tipTexts=['Rename selected nmrAtom', 'Create new nmrAtom'],
                              callbacks=[self._reassignNmrAtomPopup,
                                         partial(self._createNewNmrAtom, self.dimIndex)],
                              grid=(row, 0), hAlign='l'
                              )
        self.editButton = _buttons.getButton('Edit')
        self.newNmrAtomButton = _buttons.getButton('New')
        # _frame = Frame(parent=self._assignmentsFrame, grid=(row,0), setLayout=True, showBorder=_showBorders, ) #**settings)
        # self.editButton = Button(  parent=_frame, text='Edit',
        #                            callback=partial(self._reassignNmrAtom, self.dimIndex),
        #                            grid=(0,0), hAlign='centre',
        #                            tipText='Rename selected nmrAtom')
        #
        # self.newNmrAtomButton = Button(parent=_frame, text='New',
        #                                callback=partial(self._createNewNmrAtom, self.dimIndex),
        #                                grid=(0,1), hAlign='centre',
        #                                tipText='Create new nmrAtom')

        # row += 1
        # self._assignmentsFrame.addSpacer(5, 5, grid=(row,0), expandX=True, expandY=True)

        #===========================================
        # Not-aligned frame
        #===========================================
        # aRow += 1
        self.notAlignedFrame = Frame(self, setLayout=True, showBorder=_showBorders, grid=(aRow, 0), margins=_margins, )  #**settings)
        self.notAlignedLabel = Label(parent=self.notAlignedFrame, text='peaks\nnot aligned', grid=(0, 0),
                                     hPolicy='minimal', hAlign='centre',
                                     textColour=getColours()[LABEL_WARNINGFOREGROUND])

        #===========================================
        # set up notifiers to changes to peaks, nmrAtoms and assignments
        #===========================================
        self.tables[0].setTableNotifiers(tableClass=Peak,
                                         rowClass=NmrAtom,
                                         cellClassNames=None,
                                         tableName='assignedPeaks', rowName='nmrAtom',
                                         changeFunc=parentModule._updateInterface,
                                         className='peakLists',
                                         updateFunc=parentModule._updateInterface,
                                         tableSelection=None,
                                         pullDownWidget=None,
                                         callBackClass=NmrAtom,
                                         moduleParent=self)
        self.tables[1].setTableNotifiers(tableClass=Peak,
                                         rowClass=NmrAtom,
                                         cellClassNames=None,
                                         tableName='assignedPeaks', rowName='nmrAtom',
                                         changeFunc=parentModule._updateInterface,
                                         className='peakLists',
                                         updateFunc=parentModule._updateInterface,
                                         tableSelection=None,
                                         pullDownWidget=None,
                                         callBackClass=NmrAtom,
                                         moduleParent=self.tables)  # just to give a unique id

        # set column definitions and hidden columns for each table
        self.columnDefs = ColumnClass([('NmrAtom', lambda nmrAtom: str(nmrAtom.id), 'NmrAtom identifier', None, None),
                                       ('Pid', lambda nmrAtom: str(nmrAtom.pid), 'Pid of the nmrAtom', None, None),
                                       ('_object', lambda nmrAtom: nmrAtom, 'Object', None, None),
                                       ('Shift', lambda nmrAtom: parentModule._getShift(nmrAtom), 'Chemical shift',
                                        None, '%8.3f'),
                                       ('Delta', lambda nmrAtom: parentModule._getDeltaShift(nmrAtom, self.dimIndex),
                                        'Delta shift', None, '%6.3f')])
        self._hiddenColumns = [['Pid', 'Shift'], ['Pid', 'Shift']]

        self.tables[0]._hiddenColumns = ['Pid', 'Shift']
        self.tables[1]._hiddenColumns = ['Pid', 'Shift']

        self._assignmentWidget = self._nmrAtomWidget(parent=self._assignmentsFrame, minWidth=_pullDownWidth,
                                                     setLayout=True, showBorder=_showBorders, grid=(0, 0))
        self.editPopup = EditNmrAtomBalloon(mainWindow=self.mainWindow, on_top=True)
        self.editPopup.setCentralWidget(self._assignmentWidget)
        self.editPopup.hide()

        self.editButton.enableWidget(False)

    def _nmrAtomWidget(self, parent, minWidth, **kwds):
        """Make Frame with the nmrAtom Pulldown widgets
        :return Frame instance
        """
        _frame = Frame(parent=parent, **kwds)
        self.chainPulldown = self._createChainPulldown(parent=_frame,
                                                       grid=(0, 0), gridSpan=(1, 1),
                                                       tipText='Chain code')
        self.seqCodePulldown = self._createPulldown(parent=_frame,
                                                    grid=(0, 1), gridSpan=(1, 1),
                                                    tipText='Sequence code')
        self.resTypePulldown = self._createPulldown(parent=_frame,
                                                    grid=(0, 2), gridSpan=(1, 1),
                                                    tipText='Residue type')
        self.atomTypePulldown = self._createPulldown(parent=_frame,
                                                     grid=(0, 3), gridSpan=(1, 1),
                                                     tipText='Atom type')
        _innerFrame = Frame(parent=_frame, setLayout=True, grid=(1, 0), gridSpan=(1, 4))

        _accept = partial(self._reassignAccept, self.dimIndex)
        self._acceptButton = Button(parent=_innerFrame, text='Accept', grid=(1, 0), hAlign='r',
                                    callback=_accept)

        # set the response to pressing enter/return in the popup
        self.chainPulldown.lineEdit().returnPressed.connect(_accept)
        self.seqCodePulldown.lineEdit().returnPressed.connect(_accept)
        self.resTypePulldown.lineEdit().returnPressed.connect(_accept)
        self.atomTypePulldown.lineEdit().returnPressed.connect(_accept)
        # activate return/enter on the button when focussed
        self._acceptButton.setAutoDefault(True)

        for w in [self.chainPulldown, self.seqCodePulldown, self.resTypePulldown, self.atomTypePulldown]:
            w.setMinimumWidth(minWidth)

        return _frame

    def setEditPopupVisible(self, visible):
        """Hide the edit popup balloon
        """
        self.editPopup.setVisible(visible)

    def showNotAligned(self, flag):
        """Show/hide of notAligned and assignmentFrame"""
        self.notAlignedFrame.setVisible(flag)
        self._assignmentsFrame.setVisible(not flag)

    def setHLineText(self, text):
        """Set the text of the top horizontal Line"""
        self.hLine.setText(text)

    def setNotAlignedText(self, text):
        """Set the text of the notAligned widget"""
        self.notAlignedLabel.setText(text)

    def _close(self):
        self.tables[0]._close()
        self.tables[1]._close()
        self.tables = None

    def _clearTableOveray(self):
        for table in self.tables:
           table.setStyleSheet(table._defaultStyleSheet)

    def _handleDropsFromSideBar(self, dataDict):
        """
        Handle drops from SideBar. If NmrAtoms, then assign to the selected peaks.
        """
        objs = self.project.getObjectsByPids(dataDict.get(DropBase.PIDS))
        nmrAtoms = [x for x in objs if isinstance(x, NmrAtom)]

        if self.current.peak:
            failedNmrAtoms = []
            isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[self.dimIndex]
            for nmrAtom in nmrAtoms:
                if isotopeCode == nmrAtom.isotopeCode or nmrAtom.isotopeCode is None:
                    self._assignNmrAtom(self.dimIndex, nmrAtoms=[nmrAtom])
                else:
                    failedNmrAtoms.append(nmrAtom)
            if len(failedNmrAtoms)>0:
                showWarning('Incompatible IsotopeCode Error',
                            f'Cannot assign NmrAtoms: {nmrAtoms} to peaks with IsotopeCode {isotopeCode} ')

    def _handleDragMoveEvent(self, enteringToTableNum: int, dataDict):
        """
        Notifier callback activated upon a DragEnterEvent of an object.
        Add a border overlay if the dragEnterEvent is in the permitted table
        """
        source = dataDict.get('source')
        if source == self.tables[0] and enteringToTableNum == 1:
            self.tables[1]._setDraggingStyleSheet()
        elif source == self.tables[1] and enteringToTableNum == 0:
            self.tables[0]._setDraggingStyleSheet()
        else:
           self._clearTableOveray()


    def _handleDroppedItems(self, droppingToTableNum: int, dataDict, ):
        """
        Notifier callback activated upon a DropEvent of an object.
        Note, the source of the drag can be from anywhere, therefore here is limited only if the source is
        within the module and right tables pairs. The correct instance of the dropped object is checked afterwards.
        """
        assignmentTableNum = 0
        alternativeTableNum = 1
        sourceTable = dataDict.get('source')
        nmrAtoms = self.project.getObjectsByPids(dataDict.get(DropBase.PIDS))

        ## Action 0, Assignment: dropping to Assignment (Table-0) from Alternative (Table-1)
        if droppingToTableNum == assignmentTableNum:
            if sourceTable == self.tables[alternativeTableNum]: # needs this constraint to avoid cross-table drag&drop
                self._assignNmrAtom(self.dimIndex, nmrAtoms=nmrAtoms)
                return
        ## Action 1, DeAssign from top to bottom: dropping to Alternative (Table-1) from Assignment (Table-0)
        if droppingToTableNum == alternativeTableNum:
            if sourceTable == self.tables[assignmentTableNum]: # needs this constraint to avoid cross-table drag&drop
                self._deassignNmrAtom(self.dimIndex, nmrAtoms=nmrAtoms)
                return

    # def _clearClicked(self, val):
    #     self._clickedNmrAtom = None
    #     self._clickedLabel.setText('Current NmrAtom: <None>')
    #     self._clickedClear.setVisible(False)

    def _assignDeassignNmrAtom(self, tableNum: int, data):
        """
        Assign/Deassign the nmrAtom that is double clicked to the
        the corresponding dimension of the selected
        peaks.
        """
        if tableNum == 0:
            # deAssign from top to bottom
            self._deassignNmrAtom(self.dimIndex)
        elif tableNum == 1:
            # assign bottom - up
            self._assignNmrAtom(self.dimIndex, action=True)

    def _clickedTableCallback(self, tableNum, data):
        self.lastTableSelected = tableNum
        obj = data[Notifier.OBJECT]
        if obj:
            self._clickedNmrAtom = obj[0]
            self.editButton.enableWidget(True)

            if tableNum == 0:
                self.tables[1].clearSelection()
            elif tableNum == 1:
                self.tables[0].clearSelection()

    def _clearTableCallback(self, tableNum, data):
        self._clickedNmrAtom = None
        self.editButton.enableWidget(False)

    #     def _assignNmrAtomCallback(self, data):
    #         obj = data[Notifier.OBJECT]
    #         if obj:
    #             nmrAtom = obj[0]
    #
    #
    #
    #     # def _updatePulldownLists(self, tableNum, data):
    #     #     self.lastTableSelected = tableNum
    #     #     obj = data[Notifier.OBJECT]
    #     #     if obj:
    #     #         self._clickedNmrAtom = obj[0]
    #     #         # self._clickedLabel.setText('Current NmrAtom: {}'.format(obj[0].pid))
    #     #         # self._clickedClear.setVisible(True)
    #     #     if tableNum == 0:
    #     #         # self._updateAssignmentWidget(tableNum, obj[0])
    #     #         self.tables[1].clearSelection()
    #     #
    #     #     elif tableNum == 1:
    #     #         # self._updateAssignmentWidget(tableNum, obj[0])
    #     #         self.tables[0].clearSelection()

    def _createChainPulldown(self, parent=None, grid=(0, 0), gridSpan=(1, 1), tipText='') -> PulldownList:
        """Creates a PulldownList with callback, editable.
        """
        pulldownList = PulldownList(parent=parent, grid=grid, backgroundText=tipText, editable=True, gridSpan=gridSpan,
                                    tipText=tipText)
        # pulldownList.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        # pulldownList.setEditable(True)
        pulldownList.lineEdit().textChanged.connect(partial(self._chainEdited, pulldownList))
        return pulldownList

    def _chainEdited(self, pulldownList):
        text = pulldownList.currentText()
        chains = [chain.id for chain in self.project.nmrChains]

        if text in chains:
            index = chains.index(text)
            thisChain = self.project.nmrChains[index]
            self._setSequenceCodes(thisChain)
            self._setResidueTypes(thisChain)
            self._setAtomNames()

        self._pulldownEdited(None)

    def _createPulldown(self, parent=None, grid=(0, 0), gridSpan=(1, 1), tipText='') -> PulldownList:
        """Creates a PulldownList with callback, editable.
        """
        pulldownList = PulldownList(parent=parent, grid=grid, backgroundText=tipText, editable=True, gridSpan=gridSpan,
                                    tipText=tipText)
        # pulldownList.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        # pulldownList.setEditable(True)
        pulldownList.lineEdit().textChanged.connect(partial(self._pulldownEdited, pulldownList))
        return pulldownList

    def _createNewNmrAtom(self, dim):
        """Callback for the newNmrAtom button
        """

        isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[dim]

        with undoBlockWithoutSideBar():
            nmrChain = self.project.fetchNmrChain(shortName=defaultNmrChainCode)
            nmrResidue = nmrChain.newNmrResidue()
            nmrAtom = nmrResidue.newNmrAtom(isotopeCode=isotopeCode)

            try:

                for peak in self.current.peaks:
                    if nmrAtom not in peak.dimensionNmrAtoms[dim]:
                        # newAssignments = peak.dimensionNmrAtoms[dim] + [nmrAtom]

                        newAssignments = list(peak.dimensionNmrAtoms[dim]) + [nmrAtom]  # ejb - changed to list
                        axisCode = peak.spectrum.axisCodes[dim]
                        peak.assignDimension(axisCode, newAssignments)

                # highlight on the table and populate the pulldowns
                self.tables[0].selectObjects([nmrAtom], setUpdatesEnabled=False)
                self.tables[1].clearSelection()

                # No need for update, as this will be done by the callback on the newNmrAtom/newNmrResidue
                # self._parent._updateInterface()
                # self._updateAssignmentWidget(0, nmrAtom)

                self.lastTableSelected = 0
                self.lastNmrAtomSelected = nmrAtom

                self._clickedNmrAtom = nmrAtom
                self._showNmrAtomPopup(nmrAtom, self.tables[0])

            except Exception as es:
                showWarning(str(self.windowTitle()), str(es))

    def _showNmrAtomPopup(self, nmrAtom, table):
        """Call the popup with the supplied nmrAtom
        """
        from ccpn.ui.gui.widgets.BalloonMetrics import Side

        if nmrAtom:
            self._updateAssignmentWidget(table, nmrAtom)

        pos = QtGui.QCursor().pos()
        self.editPopup.showAt(pos, preferred_side=Side.TOP, side_priority=(Side.TOP, Side.BOTTOM, Side.RIGHT, Side.LEFT))
        self.chainPulldown.setFocus()

    def _reassignNmrAtomPopup(self):
        """Show the edit popup
        """
        from ccpn.ui.gui.widgets.BalloonMetrics import Side

        _table = self.lastTableSelected
        nextAtom = self.tables[_table].getSelectedObjects()

        self._showNmrAtomPopup(nextAtom[0] if nextAtom else None, _table)

    def _reassignAccept(self, dim: int):
        """Handle the accept button in the balloon popup
        """
        self._reassignNmrAtom(self.dimIndex)
        self.editPopup.setVisible(False)

    def _reassignNmrAtom(self, dim: int):
        """
        Assigns dimensionNmrAtoms to peak dimension when called using Assign Button in assignment widget.
        :param dim - axis dimension of the atom:
        :param action - True if callback is action from the table:
        """
        try:
            nmrChainName = self.chainPulldown.currentText()
            seqCode = self.seqCodePulldown.currentText()
            newResType = self.resTypePulldown.currentText()
            nmrAtomName = self.atomTypePulldown.currentText()

            # # get options from the pulldowns
            # currentNmrAtomSelected = (nmrChainName,
            #                           seqCode,
            #                           newResType,
            #                           nmrAtomName)
            # atomCompare = self._atomCompare(self.lastNmrAtomSelected, currentNmrAtomSelected)
            nmrAtom = None

            if not self._clickedNmrAtom:
                showWarning("Rename NmrAtom", "Please select an NmrAtom from the tables")
                return

            # wrap all actions in a single undo block
            with undoBlockWithoutSideBar():

                _chainPid = 'NC:{}'.format(nmrChainName)
                # if create and not action:
                #     # get the current chain (but may create a new one)
                #     _nmrChain = self.project.fetchNmrChain(nmrChainName)
                # else:
                #     # find the existing nmrChain
                #     _nmrChain = self.project.getByPid(_chainPid)
                #     if not _nmrChain:
                #         # raise error to notify popup
                #         raise ValueError("NmrChain doesn't exists")

                _nmrChain = self.project.fetchNmrChain(nmrChainName)
                nmrResidue = _getNmrResidue(_nmrChain, seqCode, )
                nmrAtom = nmrResidue.getNmrAtom(nmrAtomName) if nmrResidue else None

                # edit existing
                if nmrResidue and self._clickedNmrAtom.nmrResidue != nmrResidue:
                    # existing different nmrResidue
                    nmrAtom = nmrResidue.getNmrAtom(nmrAtomName)
                    if nmrAtom:
                        yesNo = showYesNo('Merge NmrAtom', "Do you want to merge\n\n"
                                                           "{}   into   {}".format(self._clickedNmrAtom.id,
                                                                                   nmrAtom.id))
                        if yesNo:
                            # merge into the new nmrAtom
                            nmrAtom.mergeNmrAtoms(self._clickedNmrAtom)

                    else:
                        # assign to a new nmrAtom
                        self._clickedNmrAtom.assignTo(chainCode=nmrChainName,
                                                      sequenceCode=seqCode,
                                                      residueType=newResType,
                                                      name=nmrAtomName,
                                                      mergeToExisting=False)

                elif nmrResidue and self._clickedNmrAtom.nmrResidue == nmrResidue:
                    # rename the same nmrAtom
                    if newResType != nmrResidue.residueType:
                        nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)

                    if nmrAtomName != self._clickedNmrAtom.name:
                        nmrAtom = nmrResidue.getNmrAtom(nmrAtomName)
                        if nmrAtom:
                            raise ValueError('NmrAtom already exists {}'.format(nmrAtom))
                        self._clickedNmrAtom.rename(nmrAtomName)

                else:
                    # nmrResidue doesn't exists
                    self._clickedNmrAtom.assignTo(chainCode=nmrChainName,
                                                  sequenceCode=seqCode,
                                                  residueType=newResType,
                                                  name=nmrAtomName,
                                                  mergeToExisting=False)

            self._parent._updateInterface()

            # update the module
            self.update()

        except Exception as es:
            showWarning('Rename NmrAtom', str(es))

    def _assignNmrAtom(self, dim: int, action: bool = False, create: bool = True, nmrAtoms = None):
        """
        Assigns dimensionNmrAtoms to peak dimension when called using Assign Button in assignment widget.
        :param dim - axis dimension of the atom:
        :param action - True if callback is action from the table:
        """
        # FIXME Potential Bug: no error checks for dim. It can give easily an IndexError

        # return if no peaks selected
        if not self.current.peaks:
            return

        try:
            # nmrChainName = self.chainPulldown.currentText()
            # seqCode = self.seqCodePulldown.currentText()
            # newResType = self.resTypePulldown.currentText()
            # nmrAtomName = self.atomTypePulldown.currentText()

            # # get options from the pulldowns
            # currentNmrAtomSelected = (nmrChainName,
            #                           seqCode,
            #                           newResType,
            #                           nmrAtomName)
            # atomCompare = self._atomCompare(self.lastNmrAtomSelected, currentNmrAtomSelected)

            # create = self.createNew.isChecked()

            selectedObjects = nmrAtoms or self.tables[1].getSelectedObjects()
            if not (selectedObjects and selectedObjects[0]):
                return
            nmrAtom = selectedObjects[0]
            if not isinstance(nmrAtom, NmrAtom):
                return

            # nmrAtom = None

            # wrap all actions in a single undo block
            with undoBlockWithoutSideBar():

                # NOTE:ED need to keep for the minute
                # _chainPid = 'NC:{}'.format(nmrChainName)
                # if create and not action:
                #     # get the current chain (but may create a new one)
                #     _nmrChain = self.project.fetchNmrChain(nmrChainName)
                # else:
                #     # find the existing nmrChain
                #     _nmrChain = self.project.getByPid(_chainPid)
                #     if not _nmrChain:
                #         # raise error to notify popup
                #         raise ValueError("NmrChain doesn't exists")
                #
                # nmrResidue = _getNmrResidue(_nmrChain, seqCode, )
                # nmrAtom = nmrResidue.getNmrAtom(nmrAtomName) if nmrResidue else None
                #
                # if not action:
                #     if create:
                #         if nmrResidue:
                #             if self._clickedNmrAtom and self._clickedNmrAtom.nmrResidue == nmrResidue and nmrResidue.residueType != newResType:
                #                 if len(nmrResidue.nmrAtoms) > 1:
                #                     yes = showYesNoWarning('Assigning nmrAtoms',
                #                                            'This will change all nmrAtoms to the residueType {}, continue?'.format(newResType))
                #                     if yes:
                #                         nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)
                #                 else:
                #                     nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)
                #
                #         else:
                #             # can do a residueType rename
                #             nmrResidue = _nmrChain.fetchNmrResidue(seqCode, newResType)
                #
                #         nmrAtom = nmrResidue.fetchNmrAtom(nmrAtomName)
                #
                #     else:
                #         pass

                # if not self._clickedNmrAtom:
                #     showWarning("Rename NmrAtom", "Please select an NmrAtom from the tables")
                #     return
                #
                # # edit existing
                # if nmrResidue and self._clickedNmrAtom.nmrResidue != nmrResidue:
                #     # existing different nmrResidue
                #     nmrAtom = nmrResidue.getNmrAtom(nmrAtomName)
                #     if nmrAtom:
                #         yesNo = showYesNo('Merge NmrAtom', "Do you want to merge\n\n"
                #                                             "{}   into   {}".format(self._clickedNmrAtom.id,
                #                                                                     nmrAtom.id))
                #         if yesNo:
                #             # merge into the new nmrAtom
                #             nmrAtom.mergeNmrAtoms(self._clickedNmrAtom)
                #
                #     else:
                #         # assign to a new nmrAtom
                #         self._clickedNmrAtom.assignTo(chainCode=nmrChainName,
                #                                       sequenceCode=seqCode,
                #                                       residueType=newResType,
                #                                       name=nmrAtomName,
                #                                       mergeToExisting=False)
                #
                # elif nmrResidue and self._clickedNmrAtom.nmrResidue == nmrResidue:
                #     # rename the same nmrAtom
                #     if newResType != nmrResidue.residueType:
                #         nmrResidue.moveToNmrChain(_chainPid, seqCode, newResType)
                #
                #     if nmrAtomName != self._clickedNmrAtom.name:
                #         nmrAtom = nmrResidue.getNmrAtom(nmrAtomName)
                #         if nmrAtom:
                #             raise ValueError('NmrAtom already exists {}'.format(nmrAtom))
                #         self._clickedNmrAtom.rename(nmrAtomName)
                #
                # else:
                #     # nmrResidue doesn't exists
                #     self._clickedNmrAtom.assignTo(chainCode=nmrChainName,
                #                                   sequenceCode=seqCode,
                #                                   residueType=newResType,
                #                                   name=nmrAtomName,
                #                                   mergeToExisting=False)

                try:

                    for peak in self.current.peaks:

                        dimNmrAtoms = list(peak.dimensionNmrAtoms[dim])

                        currentObject = nmrAtom
                        if nmrAtom not in dimNmrAtoms:
                            dimNmrAtoms.append(nmrAtom)

                            toAssign = dimNmrAtoms.index(currentObject)

                            dimNmrAtoms[toAssign] = nmrAtom
                            allAtoms = list(peak.dimensionNmrAtoms)
                            allAtoms[dim] = dimNmrAtoms
                            peak.dimensionNmrAtoms = allAtoms

                    ## Set the isotopeCode here if was not defined yet
                    if not nmrAtom.isotopeCode:
                        isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[dim]
                        nmrAtom._setIsotopeCode(isotopeCode)


                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

            self._parent._updateInterface()

            self.tables[0].selectObjects([nmrAtom], setUpdatesEnabled=False)

            if nmrAtom:
                # self._updateAssignmentWidget(0, nmrAtom)
                self.lastTableSelected = 0

            else:
                # self._updateAssignmentWidget(0, None)
                self.lastTableSelected = 0

            # update the module
            self.update()

        except Exception as es:
            showWarning('Assign NmrAtom', str(es))

    def _deassignNmrAtom(self, dim: int, nmrAtoms=None):
        """
        remove nmrAtom from peak assignment
        """

        # return if no peaks selected
        if not self.current.peaks:
            return

        try:
            currentObjects = nmrAtoms or self.tables[0].getSelectedObjects()
            if not currentObjects:
                return

            nmrAtom = currentObjects[0]
            if isinstance(nmrAtom, NmrAtom):
                try:
                    with undoBlockWithoutSideBar():
                        for peak in self.current.peaks:
                            peakDimNmrAtoms = peak.dimensionNmrAtoms
                            dimNmrAtoms = list(peakDimNmrAtoms[dim])  # ejb - changed to list
                            if nmrAtom in dimNmrAtoms:
                                dimNmrAtoms.remove(nmrAtom)

                            allAtoms = list(peakDimNmrAtoms)
                            allAtoms[dim] = dimNmrAtoms
                            peak.dimensionNmrAtoms = allAtoms

                except Exception as es:
                    showWarning(str(self.windowTitle()), str(es))

                self._parent._updateInterface()
                self.tables[1].selectObjects([nmrAtom], setUpdatesEnabled=False)
                nextAtom = self.tables[1].getSelectedObjects()
                if nextAtom:
                    # self._updateAssignmentWidget(1, currentObject[0])

                    self.lastTableSelected = 1
                    # self.buttonList.setButtonEnabled('Delete', True)
                    # self.buttonList.setButtonEnabled('Deassign', False)
                    # self.buttonList.setButtonEnabled('Assign', True)

                else:
                    # self._updateAssignmentWidget(1, None)

                    self.lastTableSelected = 1
                    # self.buttonList.setButtonEnabled('Delete', False)
                    # self.buttonList.setButtonEnabled('Deassign', False)
                    # self.buttonList.setButtonEnabled('Assign', True) #False)
        except Exception as es:
            showWarning('Deassign NmrAtom', str(es))

    def setAssignedTable(self, atomList: list):

        self.tables[0].populateTable(rowObjects=atomList,
                                     columnDefs=self.columnDefs
                                     )
        self.tables[0].sortByColumn(4, QtCore.Qt.AscendingOrder)

        # # Set the pulldowns with either thelast NmrAtom selected or the first in the list
        # if (nmrAtom := self.lastNmrAtomSelected) is None:
        #     if (objs := self.tables[0].getFirstObject()) is not None:
        #         if (objPid := objs.get('Pid')) is not None:
        #             nmrAtom = self.project.getByPid(objPid)
        #
        # if nmrAtom and not nmrAtom.isDeleted:
        #     self._updatePulldownLists(0, {Notifier.OBJECT: [nmrAtom]})

    def setAlternativesTable(self, atomList: list):

        self.tables[1].populateTable(rowObjects=atomList,
                                     columnDefs=self.columnDefs
                                     )
        self.tables[1].sortByColumn(4, QtCore.Qt.AscendingOrder)

    def _updateAssignmentWidget(self, tableNum: int, item: object):
        """
        Update all information in assignment widget when NmrAtom is selected in list widget of that
        assignment widget.
        """
        nmrAtom = item

        if nmrAtom:
            nmrChain = nmrAtom.nmrResidue.nmrChain
            sequenceCode = nmrAtom.nmrResidue.sequenceCode
            residueType = nmrAtom.nmrResidue.residueType

            if not self._parent.allChainCheckBoxLabel.isChecked():
                self._setChains(nmrChain)
                self.chainPulldown.selectValue(nmrChain.id)

                self._setSequenceCodes(nmrChain)
                idx = self.seqCodePulldown.selectValue(sequenceCode)

                self._setResidueTypes(nmrChain)
                idx = self.resTypePulldown.selectValue(residueType)

                self._setAtomNames(nmrAtom)
                idx = self.atomTypePulldown.selectValue(nmrAtom.name)

            else:

                # only allow selection of peaks from the table
                # atoms = self.objectTables[dim].getObjects()
                atoms = self.tables[tableNum]._dataFrameObject.objects
                if atoms:
                    options = [[''], [''], [''], ['']]  #'[None] * 4  # 4 empty lists
                    for atom in atoms:
                        thisOpt = atom.id.split('.')

                        for optionNum in range(0, len(thisOpt)):
                            if options[optionNum]:
                                if thisOpt[optionNum] not in options[optionNum]:
                                    options[optionNum].append(thisOpt[optionNum])
                            else:
                                options[optionNum] = [thisOpt[optionNum]]

                    self.chainPulldown.setData(options[0])
                    self.seqCodePulldown.setData(options[1])
                    self.resTypePulldown.setData(options[2])
                    self.atomTypePulldown.setData(options[3])
                else:
                    self._setDefaultPulldowns()

                self.chainPulldown.setIndex(self.chainPulldown.texts.index(nmrChain.id) if nmrChain.id in self.chainPulldown.texts else 0)
                self.seqCodePulldown.setIndex(self.seqCodePulldown.texts.index(sequenceCode) if sequenceCode in self.seqCodePulldown.texts else 0)
                self.resTypePulldown.setIndex(self.resTypePulldown.texts.index(residueType) if residueType in self.resTypePulldown.texts else 0)
                self.atomTypePulldown.setIndex(self.atomTypePulldown.texts.index(nmrAtom.name) if nmrAtom.name in self.atomTypePulldown.texts else 0)

            self.lastNmrAtomSelected = nmrAtom
        else:
            self._setDefaultPulldowns()
            self.lastNmrAtomSelected = None

    def _setDefaultPulldowns(self):
        """Clear the contents of the pullDowns
        """
        self.chainPulldown.clear()
        self.seqCodePulldown.clear()
        self.resTypePulldown.clear()
        self.atomTypePulldown.clear()

        self._setChains()
        self._setSequenceCodes()
        self._setResidueTypes()
        self._setAtomNames()

    def _setChains(self, nmrChain=None):
        """Populate the chain pulldown from the project
        """
        thisChain = self.chainPulldown.currentText()
        chains = ['']
        chains.extend([chain.id for chain in self.project.nmrChains])
        if nmrChain:
            thisChain = nmrChain.id

        self.chainPulldown.setData(chains)
        self.chainPulldown.setIndex(self.chainPulldown.texts.index(thisChain) if thisChain in self.chainPulldown.texts else 0)

    def _setSequenceCodes(self, nmrChain=None):
        """Populate the sequenceCode pulldown from the nmrChain or project
        """
        thisSeq = self.seqCodePulldown.currentText()
        sequenceCodes = ['']
        if nmrChain:
            sequenceCodes.extend([nmrResidue.sequenceCode for nmrResidue in nmrChain.nmrResidues])
        else:
            sequenceCodes.extend([nmrResidue.sequenceCode for nmrResidue in self.project.nmrResidues])

        self.seqCodePulldown.setData(sorted(sequenceCodes, key=CcpnSorting.stringSortKey))
        self.seqCodePulldown.setIndex(self.seqCodePulldown.texts.index(thisSeq) if thisSeq in self.seqCodePulldown.texts else 0)

    def _setResidueTypes(self, nmrChain=None):
        """Populate the residueTypes pulldown from the nmrChain or project
        """
        thisRes = self.resTypePulldown.currentText()
        residueTypes = ['']
        if nmrChain:
            residueTypes.extend([nmrResidue.residueType for nmrResidue in nmrChain.nmrResidues])
        else:
            residueTypes.extend([nmrResidue.residueType for nmrResidue in self.project.nmrResidues])

        residueTypes.extend([nmrResidue[1] for nmrResidue in allowedResidueTypes])  # self.project.nmrResidues]
        residueTypes = list(set(OrderedDict.fromkeys(residueTypes)))

        self.resTypePulldown.setData(sorted(residueTypes, key=CcpnSorting.stringSortKey))
        self.resTypePulldown.setIndex(self.resTypePulldown.texts.index(thisRes) if thisRes in self.resTypePulldown.texts else 0)

    def _setAtomNames(self, nmrAtom=None):
        """Populate the atomNames pulldown from the project
        """
        thisAtom = self.atomTypePulldown.currentText()
        atomNames = ['']
        if self.current.peak:
            isotopeCode = self.current.peak.peakList.spectrum.isotopeCodes[self.dimIndex]
            atomNames += getIsotopeListFromCode(isotopeCode)

        if nmrAtom:
            atomNames.insert(0, nmrAtom.name)
            thisAtom = nmrAtom.name  # set only if nmrAtom defined

        if self.lastNmrAtomSelected:
            atomNames.append(self.lastNmrAtomSelected.pid.fields[3])

        self.atomTypePulldown.setData(sorted(list(set(atomNames)), key=greekKey))
        self.atomTypePulldown.setIndex(self.atomTypePulldown.texts.index(thisAtom) if thisAtom in self.atomTypePulldown.texts else 0)

    def _deleteNmrAtom(self, dim: int):
        """
        delete selected nmrAtom from project
        """
        if self.lastTableSelected is not None:

            # # deassign if assigned
            # self._deassignNmrAtom(dim)

            # remove from the table
            deleted = self.tables[self.lastTableSelected].deleteObjFromTable()
            if deleted:
                nextAtoms = self.tables[self.lastTableSelected].getSelectedObjects()

                # reset buttons
                if not nextAtoms:

                    # self.buttonList.setButtonEnabled('Delete', False)
                    # self.buttonList.setButtonEnabled('Deassign', False)
                    # self.buttonList.setButtonEnabled('Assign', True) #False)

                    self._updateAssignmentWidget(self.lastTableSelected, None)
                else:
                    self._updateAssignmentWidget(self.lastTableSelected, nextAtoms[0])

    def _pulldownEdited(self, pulldown: object):
        """
        Enable the assignment button if the text has changed in the pulldown
        """
        pass
        # currentNmrAtomSelected = (self.chainPulldown.currentText(),
        #                           self.seqCodePulldown.currentText(),
        #                           self.resTypePulldown.currentText(),
        #                           self.atomTypePulldown.currentText())
        # enable = False in self._atomCompare(self.lastNmrAtomSelected, currentNmrAtomSelected)
        # self.buttonList.setButtonEnabled('Assign', True) #enable)

    def _atomCompare(self, atom1: tuple, atom2: tuple):
        """
        check whether the selection has changed from being clicked
        """
        if atom1 and atom2:
            return [True if a == b else False for a, b in zip(atom1, atom2)]
        else:
            return [False]
