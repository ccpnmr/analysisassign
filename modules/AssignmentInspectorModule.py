"""This file contains AssignmentInspectorModule class

modified by Geerten 1-9/12/2016:
- initialisation with 'empty' settings possible,
- now responsive to current.nmrResidues
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
                 "J.Biomol.Nmr (2016), 66, 111-124, https://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2022-10-12 15:27:02 +0100 (Wed, October 12, 2022) $"
__version__ = "$Revision: 3.1.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: geertenv $"
__date__ = "$Date: 2016-07-09 14:17:30 +0100 (Sat, 09 Jul 2016) $"
#=========================================================================================
# Start of code
#=========================================================================================

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QAbstractScrollArea
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Optional

from ccpn.core.NmrAtom import NmrAtom, NmrResidue
from ccpn.core.ChemicalShiftList import ChemicalShiftList
from ccpn.ui._implementation.Strip import Strip
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.lib.CallBack import CallBack
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget
from ccpn.ui.gui.widgets.SettingsWidgets import SpectrumDisplaySelectionWidget, IncludeCurrent
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.Splitter import Splitter
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.ui.gui.widgets.PulldownListsForObjects import ChemicalShiftListPulldown
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.modules.ChemicalShiftTable import _NewChemicalShiftTable
from ccpn.ui.gui.modules.PeakTable import _NewPeakTableWidget
from ccpn.ui.gui.lib.StripLib import navigateToNmrResidueInDisplay  #, _getCurrentZoomRatio
from ccpn.ui.gui.lib.SpectrumDisplay import navigateToNmrResidueInStrip
from ccpn.util.OrderedSet import OrderedSet
from ccpn.util.Logging import getLogger
from ccpn.util.AttrDict import AttrDict


logger = getLogger()
ALL = '<Use all>'
NMRRESIDUES = 'nmrResidues'
NMRATOMS = 'nmrAtoms'


# small object to facilitate passing data to peakTable
@dataclass
class _emptyObject:
    peaks = []
    spectrum = None


#=========================================================================================
# Assignment Inspector Module
#=========================================================================================

class AssignmentInspectorModule(CcpnModule):
    """
    This Module allows inspection of the NmrAtoms of a selected NmrResidue
    It responds to current.nmrResidues, taking the last added residue to this list
    The NmrAtom listWidget allows for selection of the nmrAtom; subsequently its assignedPeaks
    are displayed.
    """

    # override in specific module implementations
    className = 'AssignmentInspectorModule'
    attributeName = 'peaks'

    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    settingsPosition = 'top'
    LARGE_STRETCH = 100000
    SETTING_PADDING = 4

    def __init__(self, mainWindow, name='Assignment Inspector', chemicalShiftList=None, selectFirstItem=False):
        super().__init__(mainWindow=mainWindow, name=name, settingsScrollBarPolicies=('asNeeded', 'never'))  # gwv

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        self.application = mainWindow.application
        self.project = mainWindow.application.project
        self.current = mainWindow.application.current

        # settings window
        self.splitter = Splitter(self.mainWidget, horizontal=False)
        self._chemicalShiftFrame = Frame(self.splitter, setLayout=True)
        self._assignmentFrame = Frame(self.splitter, setLayout=True)
        self.mainWidget.getLayout().addWidget(self.splitter)

        self.splitter.setStretchFactor(0, 3)
        self.splitter.setStretchFactor(1, 2)
        self.splitter.setChildrenCollapsible(False)
        self._assignmentFrame.setMinimumHeight(100)

        # cannot set a notifier for displays, as these are not (yet?) implemented and the Notifier routines
        # underpinning the addNotifier call do not allow for it either
        colwidth = 140

        self.settingsWidget.layout().setColumnStretch(2, self.LARGE_STRETCH)

        self._settingsScrollArea.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)
        self._settingsScrollArea.setStyleSheet(".ScrollArea {padding: %ipx}" % self.SETTING_PADDING)
        self._settingsScrollArea.setScrollBarPolicies(('asNeeded', 'never'))

        self._splitWidget = Frame(self.settingsWidget, grid=(0, 0), setLayout=True, vPolicy='minimumExpanding')

        self._tickLisWidget = Frame(self._splitWidget, grid=(0, 1), setLayout=True, vPolicy='minimum')

        self.displaysWidget = SpectrumDisplaySelectionWidget(self._splitWidget, mainWindow=self.mainWindow,
                                                             grid=(0, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                                                             vPolicy='maximum',
                                                             orientation='left',
                                                             labelText='Display(s):',
                                                             tipText='SpectrumDisplay modules to respond to double-click',
                                                             # texts=[ALL, UseCurrent] + [display.pid for display in self.application.ui.mainWindow.spectrumDisplays],
                                                             defaults=[ALL],
                                                             standardListItems=[ALL, IncludeCurrent]
                                                             )

        self.sequentialStripsWidget = CheckBoxCompoundWidget(
                self._tickLisWidget,
                grid=(0, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                #minimumWidths=(colwidth, 0),
                fixedWidths=(colwidth, 30),
                orientation='left',
                labelText='Show sequential strips:',
                checked=False
                )

        self.markPositionsWidget = CheckBoxCompoundWidget(
                self._tickLisWidget,
                grid=(1, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                #minimumWidths=(colwidth, 0),
                fixedWidths=(colwidth, 30),
                orientation='left',
                labelText='Mark positions:',
                checked=True
                )
        self.autoClearMarksWidget = CheckBoxCompoundWidget(
                self._tickLisWidget,
                grid=(2, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                #minimumWidths=(colwidth, 0),
                fixedWidths=(colwidth, 30),
                orientation='left',
                labelText='Auto clear marks:',
                checked=True
                )
        self.showNmrAtomListWidget = CheckBoxCompoundWidget(
                self._tickLisWidget,
                grid=(3, 0), vAlign='top', stretch=(0, 0), hAlign='left',
                #minimumWidths=(colwidth, 0),
                fixedWidths=(colwidth, 30),
                orientation='left',
                labelText='Show nmrAtom list:',
                checked=True,
                callback=self._setNmrAtomListVisible,
                )

        self._tickLisWidget.layout().setRowStretch(4, self.LARGE_STRETCH)

        minHeight = self._calculateMinHeight()
        self._settingsScrollArea.setMinimumSizes((self._settingsScrollArea.minimumWidth(), minHeight))
        self.nmrAtomBlocking = True
        self._nmrResidues = []

        # main window
        # AssignedPeaksTable need to be initialised before chemicalShiftTable, as the callback of the latter requires
        # the former to be present
        self._setWidgets()

        # disable current callback - not required for assignmentInspector
        # responds to changes in current nmrAtoms and nmrResidues?
        self.chemicalShiftTable.clearCurrentCallback()

        # notifier to handle deleting items
        # self.chemicalShiftTable._tableSelectionChanged.connect(self._tableSelectionCallback)

        if chemicalShiftList is not None:
            self._selectTable(chemicalShiftList)
        elif selectFirstItem:
            self._modulePulldown.selectFirstItem()

        self._setNmrAtomListVisible(True)
        self._registerNotifiers()

    @QtCore.pyqtSlot(list)
    def _tableSelectionCallback(self, shifts):

        nmrResidues = list(set(cs.nmrAtom.nmrResidue for cs in shifts if cs.nmrAtom))
        _temp = AttrDict()
        _temp.nmrResidues = nmrResidues

        self._highlightNmrResidues({CallBack.OBJECT: _temp})

    def _setWidgets(self):
        """Set up the widgets for the new chemicalShiftTable
        """
        _topWidget = self._chemicalShiftFrame

        # main widgets in the top splitter
        row = 0
        Spacer(_topWidget, 5, 5,
               QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed,
               grid=(0, 0), gridSpan=(1, 1))
        row += 1

        self._modulePulldown = ChemicalShiftListPulldown(parent=_topWidget,
                                                         mainWindow=self.mainWindow, default=None,
                                                         grid=(row, 0), gridSpan=(1, 1), minimumWidths=(0, 100),
                                                         showSelectName=True,
                                                         sizeAdjustPolicy=QtWidgets.QComboBox.AdjustToContents,
                                                         callback=self._selectionPulldownCallback,
                                                         )
        # fixed height
        self._modulePulldown.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)

        row += 1
        self.spacer = Spacer(_topWidget, 5, 5,
                             QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed,
                             grid=(2, 1), gridSpan=(1, 1))
        _topWidget.getLayout().setColumnStretch(1, 2)

        row += 1
        self.chemicalShiftTable = _AssignmentInspectorTable(parent=_topWidget,
                                                            mainWindow=self.mainWindow,
                                                            moduleParent=self,
                                                            grid=(row, 0), gridSpan=(1, 6),
                                                            )

        # bottom frame containing atomList and peakTable - bottom splitter
        _bottomWidget = self._assignmentFrame

        # frame containing the list of nmrAtoms in the chemicalShift selection
        self.nmrAtomListFrame = Frame(parent=_bottomWidget, grid=(0, 0), gridSpan=(2, 1), setLayout=True)

        self.nmrAtomLabel = Label(self.nmrAtomListFrame, 'Filter by NmrAtom(s):', bold=True,
                                  grid=(0, 0), gridSpan=(1, 1))

        self.attachedNmrAtomsList = ListWidget(self.nmrAtomListFrame,
                                               contextMenu=False,
                                               grid=(1, 0), gridSpan=(1, 1)
                                               )
        self.attachedNmrAtomsList.itemSelectionChanged.connect(self._updatePeakTableCallback)
        self.attachedNmrAtomsList.setDragEnabled(False)
        self.attachedNmrAtomsList.setSortingEnabled(True)

        self.attachedNmrAtomsList.setFixedWidth(self.nmrAtomLabel.sizeHint().width())
        self.attachedNmrAtomsList.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        self.nmrAtomListFrame.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Ignored)

        # peaks assigned to the selection nmrAtoms - from chemicalShifts
        self.peaksLabel = Label(_bottomWidget, 'Peaks assigned to NmrAtom(s):', bold=True,
                                grid=(0, 1), gridSpan=(1, 6))

        self.assignedPeaksTable = _NewPeakTableWidget(parent=_bottomWidget,
                                                      mainWindow=self.mainWindow,
                                                      moduleParent=self,
                                                      grid=(1, 1), gridSpan=(1, 6),
                                                      actionCallback=self._peakActionCallback,
                                                      )

    def _selectTable(self, chemicalShiftList=None):
        """Manually select a ChemicalShiftList from the pullDown
        """
        if chemicalShiftList is None:
            self._modulePulldown.selectFirstItem()
        else:
            if not isinstance(chemicalShiftList, ChemicalShiftList):
                getLogger().warning('select: Object is not of type ChemicalShiftList')
                raise TypeError('select: Object is not of type ChemicalShiftList')
            else:
                self._modulePulldown.select(chemicalShiftList.pid)

    def _calculateMinHeight(self):
        contentsMargins = self._settingsScrollArea.contentsMargins()
        marginsTotalVertical = contentsMargins.top() + contentsMargins.bottom()
        minHeight = max(self._tickLisWidget.sizeHint().height(), self.displaysWidget.minimumSizeHint().height()) + (
                self.SETTING_PADDING * 2) + marginsTotalVertical
        return minHeight

    def _maximise(self):
        """
        refresh the table on a maximise event
        """
        self._refreshTable()

    def _registerNotifiers(self):
        """Set up the notifiers
        """
        self.setNotifier(self.current, [Notifier.CURRENT], targetName=NmrResidue._pluralLinkName,
                         callback=self._highlightNmrResidues)
        self.setNotifier(self.project, [Notifier.RENAME, Notifier.CREATE, Notifier.DELETE],
                         NmrAtom.__name__, self._updateNmrAtoms, onceOnly=True)

        self.nmrAtomBlocking = False

    def _closeModule(self):
        """
        CCPN-INTERNAL: used to close the module
        """
        self.assignedPeaksTable._close()
        self.chemicalShiftTable._close()
        super()._closeModule()

    def _selectionPulldownCallback(self, item):
        """Notifier Callback for selecting ChemicalShiftList from the pull down menu
        """
        self._table = self._modulePulldown.getSelectedObject()
        self.chemicalShiftTable._table = self._table

        if self._table is not None:
            self.chemicalShiftTable.populateTable(rowObjects=self._table.chemicalShifts,
                                                  selectedObjects=self.current.chemicalShifts)
        else:
            self.chemicalShiftTable.populateEmptyTable()

    def _setNmrAtomListVisible(self, visible=None):
        """change the visibility of the nmrAtom list in the peak tables widget
        """
        if visible is None:
            visible = self.showNmrAtomListWidget.get()
        else:
            if not isinstance(visible, bool):
                raise TypeError('nmrList visibility must be True/False')

            self.showNmrAtomListWidget.set(visible)

        self.nmrAtomListFrame.setVisible(visible)

    @contextmanager
    def _notifierBlanking(self):
        """block nmrAtom and nmrResidue notifiers for this module only
        """
        self.setBlankingAllNotifiers(True)
        try:
            # transfer control to the calling function
            yield

        except Exception as es:
            getLogger().warning('Error in AssignmentInspectorModule', str(es))

        finally:
            # reset notifiers
            self.setBlankingAllNotifiers(False)

    def navigateToNmrResidueCallBack(self, selection, lastItem):
        """Navigate in selected displays to nmrResidue; skip if none defined
        """
        # handle a single chemicalShift - SHOULD always contain an object
        try:
            if not (objs := list(lastItem[self._OBJECT])):
                return
        except Exception as es:
            getLogger().debug2(f'{self.__class__.__name__}.navigateToNmrResidueCallBack: No selection\n{es}')
            return

        if isinstance(objs, (tuple, list)):
            chemicalShift = objs[0]
        else:
            chemicalShift = objs

        self._navigateByChemicalShift(chemicalShift)

    def _navigateByChemicalShift(self, chemicalShift):
        if not chemicalShift.nmrAtom:
            return
        nmrResidue = chemicalShift.nmrAtom.nmrResidue

        getLogger().debug('nmrResidue=%s' % (nmrResidue.id))

        displays = self.displaysWidget.getDisplays()
        if len(displays) == 0:
            logger.warning('Undefined display module(s); select in settings first')
            showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
            return

        from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar
        from ccpn.ui.gui.lib.StripLib import navigateToPositionInStrip, _getCurrentZoomRatio

        with undoBlockWithoutSideBar():
            # optionally clear the marks
            if self.autoClearMarksWidget.checkBox.isChecked():
                self.mainWindow.clearMarks()

            # navigate the displays
            for display in displays:
                if isinstance(display, Strip):
                    strip = display
                    display = strip.spectrumDisplay
                    newWidths = []  #_getCurrentZoomRatio(display.strips[0].viewBox.viewRange())
                    navigateToNmrResidueInStrip(display, strip=strip,
                                                nmrResidue=nmrResidue,
                                                widths=newWidths,  #['full'] * len(display.strips[0].axisCodes),
                                                markPositions=self.markPositionsWidget.checkBox.isChecked()
                                                )
                elif len(display.strips) > 0:
                    newWidths = []  #_getCurrentZoomRatio(display.strips[0].viewBox.viewRange())
                    navigateToNmrResidueInDisplay(nmrResidue, display, stripIndex=0,
                                                  widths=newWidths,  #['full'] * len(display.strips[0].axisCodes),
                                                  showSequentialResidues=(len(display.axisCodes) > 2) and
                                                                         self.sequentialStripsWidget.checkBox.isChecked(),
                                                  markPositions=self.markPositionsWidget.checkBox.isChecked()
                                                  )

    def _selectionCallback(self, data):
        """
        Notifier Callback for double-clicking a row in the table
        Highlight all nmrAtoms belonging to the same nmrResidue as nmrAtom in checmicalShifts
        """

        objList = data[CallBack.OBJECT]

        if objList:
            getLogger().debug('AssignmentInspectorChemicalShift>>> selection', objList)

        self._selectByChemicalShifts(objList)

    def _updateNmrAtoms(self, data):
        """
        Notifier Callback for nmrAtom change - to update assignment table and list
        """
        if self.nmrAtomBlocking:
            return

        data[CallBack.OBJECT] = []
        self._refreshNmrAtoms(data)

    def _updateChemicalShifts(self, data):
        """
        Notifier Callback for chemcialShift change - to update assignment table and list
        """
        data[CallBack.OBJECT] = []
        self._refreshNmrAtoms(data)

    def _selectByChemicalShifts(self, chemicalShifts):
        if chemicalShifts:
            nmrResidues = tuple(set(cs.nmrAtom.nmrResidue for cs in chemicalShifts if cs.nmrAtom))
        else:
            nmrResidues = []

        if nmrResidues:
            with self._notifierBlanking():
                nmrAtoms = tuple(set(nmrAtom for nmrRes in nmrResidues for nmrAtom in nmrRes.nmrAtoms))

                self.current.nmrAtoms = nmrAtoms
                self.current.nmrResidues = nmrResidues

                self._highlightChemicalShifts(nmrResidues)
                self._updateModuleCallback({NMRRESIDUES: nmrResidues,
                                            NMRATOMS   : nmrAtoms},
                                           updateFromNmrResidues=True)

    def _highlightChemicalShifts(self, nmrResidues):
        """
        Highlight chemical shifts in the table
        """
        if self.chemicalShiftTable._table:
            getLogger().debug('_highlightChemicalShifts ', nmrResidues)

            chemicalShifts = self.chemicalShiftTable._table.chemicalShifts
            residues = set(nmrResidues)
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom and not cs.nmrAtom.isDeleted and cs.nmrAtom.nmrResidue in residues]

            self.chemicalShiftTable._highLightObjs(highlightList)

    def _highlightNmrResidues(self, data):
        """
        Notifier Callback for highlighting all NmrAtoms in the table
        """
        if self.nmrAtomBlocking:
            return

        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._table:
            getLogger().debug('_highlightNmrResidues ', objList)

            chemicalShifts = self.chemicalShiftTable._table.chemicalShifts
            nmrResidues = set(objList.nmrResidues)  #        set([atom.nmrResidue for atom in self.current.nmrAtoms if atom])
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom and not cs.nmrAtom.isDeleted and cs.nmrAtom.nmrResidue in nmrResidues]

            self.chemicalShiftTable._highLightObjs(highlightList)

            # will respond to selection of nmrAtom in sequenceGraph
            self._updateModuleCallback({NMRRESIDUES: nmrResidues},
                                       updateFromNmrResidues=True)

    def _highlightNmrAtoms(self, data):
        """
        Notifier Callback for highlighting all NmrAtoms in the table
        """
        if self.nmrAtomBlocking:
            return

        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._table:
            getLogger().debug('_highlightNmrAtoms ', objList)

            chemicalShifts = self.chemicalShiftTable._table.chemicalShifts
            nmrResidues = set([atom.nmrResidue for atom in self.current.nmrAtoms if atom])
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom and not cs.nmrAtom.isDeleted and cs.nmrAtom.nmrResidue in nmrResidues]

            self.chemicalShiftTable._highLightObjs(highlightList)

            # will respond to selection of nmrAtom in sequenceGraph
            self._updateModuleCallback({NMRRESIDUES: nmrResidues,
                                        NMRATOMS   : self.current.nmrAtoms},
                                       updateFromNmrResidues=False)

    def _refreshNmrAtoms(self, data):
        """
        Notifier Callback for refreshing all NmrAtoms in the table
        """
        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._table:
            getLogger().debug('_refreshNmrAtoms ', objList)

            self._updateModuleCallback(None, updateFromNmrResidues=False)

    def _updateModuleCallback(self, data: Optional[dict], updateFromNmrResidues=True):
        """
        Callback function: Module responsive to nmrResidues; updates the list widget with nmrAtoms and updates peakTable if
        current.nmrAtom belongs to nmrResidue
        """
        if data:

            # update list from nmrResidues, show all nmrAtoms connected
            nmrResidues = data[NMRRESIDUES] if NMRRESIDUES in data else []

            if updateFromNmrResidues:
                # generate nmrAtom list from nmrResidues
                nmrAtoms = OrderedSet()
                for nmrRes in nmrResidues:
                    for nmrAtom in nmrRes.nmrAtoms:
                        nmrAtoms.add(nmrAtom)

            else:
                # get from the data dict
                nmrAtoms = data[NMRATOMS] if NMRATOMS in data else []

            self._nmrResidues = nmrResidues
            self._nmrAtoms = nmrAtoms

            _select = self.attachedNmrAtomsList.getSelectedTexts()
            # there is currently a hidden list widget containing the nmrAtom ids
            with self.attachedNmrAtomsList.blockWidgetSignals(self.attachedNmrAtomsList):
                self.attachedNmrAtomsList.clear()

                # NOTE:ED - should we only display those nmrAtoms that are in the chemicalShift table? - similarly for peaks
                self.ids = [atm.id for atm in nmrAtoms if not atm.isDeleted]

                self.attachedNmrAtomsList.addItems(self.ids)
                self.attachedNmrAtomsList.selectItems(_select)

            # populate peak table with the correct peaks
            self._updatePeakTable(nmrAtoms, messageAll=updateFromNmrResidues)

        else:
            # data is None, so update from current settings - should only be called from _refreshNmrAtoms

            _select = self.attachedNmrAtomsList.getSelectedTexts()
            # there is currently a hidden list widget containing the nmrAtom ids
            with self.attachedNmrAtomsList.blockWidgetSignals(self.attachedNmrAtomsList):
                self.attachedNmrAtomsList.clear()
                _nmrAtoms = [atm for _nmrRes in self._nmrResidues if not _nmrRes.isDeleted for atm in _nmrRes.nmrAtoms if not atm.isDeleted]
                self.ids = [atm.id for atm in _nmrAtoms]
                self.attachedNmrAtomsList.addItems(self.ids)
                self.attachedNmrAtomsList.selectItems(_select)

            # populate peak table with the correct peaks
            self._updatePeakTable(_nmrAtoms, messageAll=updateFromNmrResidues)

    def _updatePeakTableCallback(self, data=None):
        """
        Update the peakTable using item.text (which contains a NmrAtom pid or <all>)
        """
        # something has been selected, so get all selected items
        numTexts = self.attachedNmrAtomsList.count()
        selectedTexts = self.attachedNmrAtomsList.getSelectedTexts()
        nmrAtoms = [self.project.getByPid('NA:' + _id) for _id in selectedTexts]

        # populate the table with valid nmrAtoms
        self._updatePeakTable([atm for atm in nmrAtoms if atm is not None],
                              messageAll=True if numTexts == len(nmrAtoms) else False)

    def _updatePeakTable(self, nmrAtoms, messageAll=True):
        """
        Update peak table depending on value of id;
        clears peakTable if pid is None
        """
        if not nmrAtoms:
            # get all items from the table
            nmrAtoms = [self.project.getByPid('NA:' + _id) for _id in self.attachedNmrAtomsList.getTexts()]

            # # populate the table with valid nmrAtoms
            # self._updatePeakTable([atm for atm in nmrAtoms if atm is not None], messageAll=True)

        self._peakList = _emptyObject()
        allPeaks = list(set([pk for nmrAtom in nmrAtoms if nmrAtom for pk in nmrAtom.assignedPeaks]))
        chemicalShiftList = self._modulePulldown.getSelectedObject()
        if chemicalShiftList:
            spectra = chemicalShiftList.spectra  #show peaks only for spectra currently available for the selected CSL
            peaks = [peak for peak in allPeaks if peak.peakList.spectrum in spectra]
        else:
            peaks = []

        self._peakList.peaks = peaks
        self._peakList.spectrum = None
        if peaks:
            # get a spectrum with the largest dimensionCount - just take last from sorted
            specs = sorted((pk.spectrum.dimensionCount, pk.spectrum) for pk in peaks)
            _maxCount, self._peakList.spectrum = specs[-1]

        self.assignedPeaksTable._table = self._peakList
        self.assignedPeaksTable.populateTable()

        ids = [atm.id for atm in nmrAtoms]

        if messageAll:
            self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % ALL)
        else:
            atomList = ', '.join([str(_id) for _id in ids])
            self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % atomList)  # nmrAtom.id)

    #=========================================================================================
    # Peak-table callbacks
    #=========================================================================================

    def _peakActionCallback(self, selection, lastItem):
        """If current strip contains the double-clicked peak will navigateToPositionInStrip
        """
        from ccpn.ui.gui.lib.StripLib import navigateToPositionInStrip, _getCurrentZoomRatio

        try:
            if not (objs := list(lastItem[self.assignedPeaksTable._OBJECT])):
                return
        except Exception as es:
            getLogger().debug2(f'{self.__class__.__name__}._peakActionCallback: No selection\n{es}')
            return

        if isinstance(objs, (tuple, list)):
            peak = objs[0]
        else:
            peak = objs

        dpObjs = self.displaysWidget.getDisplays()
        if dpObjs:
            # check which spectrumDisplays to navigate to
            for dp in dpObjs:
                if isinstance(dp, Strip):
                    widths = None
                    if peak.peakList.spectrum.dimensionCount <= 2:
                        widths = _getCurrentZoomRatio(dp.viewRange())

                    navigateToPositionInStrip(strip=dp,
                                              positions=peak.position,
                                              axisCodes=peak.axisCodes,
                                              widths=widths
                                              )
                elif dp.strips:
                    widths = None
                    if peak.peakList.spectrum.dimensionCount <= 2:
                        widths = _getCurrentZoomRatio(dp.strips[0].viewRange())

                    navigateToPositionInStrip(strip=dp.strips[0],
                                              positions=peak.position,
                                              axisCodes=peak.axisCodes,
                                              widths=widths
                                              )

        else:
            logger.warning('Impossible to navigate to peak position. Set a current strip first or select spectrumDisplays in gearbox settings')


#=========================================================================================
# Assignment Inspector Table - subclassed from ChemicalShiftTable
#=========================================================================================

class _AssignmentInspectorTable(_NewChemicalShiftTable):
    """ChemicalShift table in AssignmentInspector module with modified behaviour
    """

    def actionCallback(self, selection, lastItem):
        """Notifier DoubleClick action on item in table. Mark a chemicalShift based on all attached nmrAtoms
        """
        from ccpn.AnalysisAssign.modules.BackboneAssignmentModule import markNmrAtoms

        cShifts = self.getSelectedObjects()
        if len(self.mainWindow.marks):
            if self.moduleParent.autoClearMarksWidget.checkBox.isChecked():
                self.mainWindow.clearMarks()
        if cShifts:
            nmrAtoms = list(set(cs.nmrAtom for cs in cShifts if cs.nmrAtom))
            markNmrAtoms(self.mainWindow, nmrAtoms)

    def selectionCallback(self, selected, deselected, selection, lastItem):
        """Notifier Callback for selecting rows in the table
        """
        try:
            if not (objs := list(selection[self._OBJECT])):
                return

        except Exception as es:
            getLogger().debug2(f'{self.__class__.__name__}.selectionCallback: No selection\n{es}')
            return

        self.current.chemicalShifts = objs or []

        if objs:
            nmrResidues = tuple(set(cs.nmrAtom.nmrResidue for cs in objs if cs.nmrAtom))
        else:
            nmrResidues = []

        if nmrResidues:
            # set the associated nmrResidue and nmrAtoms
            nmrAtoms = tuple(set(nmrAtom for nmrRes in nmrResidues for nmrAtom in nmrRes.nmrAtoms))
            self.current.nmrAtoms = nmrAtoms
            self.current.nmrResidues = nmrResidues

        else:
            self.current.nmrAtoms = []
            self.current.nmrResidues = []

        # get all the chemicalShifts linked by nmrResidue
        allShifts = list(filter(None, set(cs for nmrAt in nmrAtoms for cs in nmrAt.chemicalShifts)))

        # highlight all the chemicalShifts linked to the nmrResidues
        self._highLightObjs(allShifts, scrollToSelection=False)
