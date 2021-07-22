"""This file contains AssignmentInspectorModule class

modified by Geerten 1-9/12/2016:
- intialisation with 'empty' settings possible,
- now responsive to current.nmrResidues
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
__dateModified__ = "$dateModified: 2021-06-04 15:23:19 +0100 (Fri, June 04, 2021) $"
__version__ = "$Revision: 3.0.4 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: geertenv $"
__date__ = "$Date: 2016-07-09 14:17:30 +0100 (Sat, 09 Jul 2016) $"
#=========================================================================================
# Start of code
#=========================================================================================

from PyQt5 import QtCore, QtWidgets
from contextlib import contextmanager
from typing import Optional
from ccpn.util.OrderedSet import OrderedSet
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.ListWidget import ListWidget
from ccpn.ui.gui.widgets.GuiTable import GuiTable
from ccpn.ui.gui.widgets.Column import ColumnClass, Column
from ccpn.util.Logging import getLogger
from ccpn.ui.gui.widgets.CompoundWidgets import CheckBoxCompoundWidget
from ccpn.ui.gui.widgets.CompoundWidgets import ListCompoundWidget
from ccpn.core.lib.peakUtils import getPeakPosition, getPeakAnnotation
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.NmrAtom import NmrAtom, NmrResidue
from ccpn.core.Peak import Peak
from ccpn.ui.gui.modules.ChemicalShiftTable import ChemicalShiftTable
from ccpn.ui.gui.widgets.Splitter import Splitter
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.core.lib.CallBack import CallBack
from ccpn.ui.gui.lib.Strip import navigateToNmrAtomsInStrip, \
    _getCurrentZoomRatio, navigateToNmrResidueInDisplay
from ccpn.core.PeakList import PeakList
from ccpn.util.Common import makeIterableList
from PyQt5.QtWidgets import QAbstractScrollArea
from ccpn.ui.gui.widgets.SettingsWidgets import SpectrumDisplaySelectionWidget


logger = getLogger()
ALL = '<all>'
NMRRESIDUES = 'nmrResidues'
NMRATOMS = 'nmrAtoms'


class _emptyObject():
    # small object to facilitate passing data to simulated event
    def __init__(self):
        pass


class AssignmentInspectorModule(CcpnModule):
    """
    This Module allows inspection of the NmrAtoms of a selected NmrResidue
    It responds to current.nmrResidues, taking the last added residue to this list
    The NmrAtom listWidget allows for selection of the nmrAtom; subsequently its assignedPeaks
    are displayed.
    """

    # overide in specific module implementations
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
                                                 texts=[ALL] + [display.pid for display in self.application.ui.mainWindow.spectrumDisplays],
                                                 defaults=[ALL]
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

        # main window
        # AssignedPeaksTable need to be initialised before chemicalShiftTable, as the callback of the latter requires
        # the former to be present
        self.assignedPeaksTable = AssignmentInspectorTable(parent=self._assignmentFrame,
                                                           mainWindow=self.mainWindow,
                                                           moduleParent=self,
                                                           setLayout=True,
                                                           selectionCallback=self._setCurrentPeak,
                                                           actionCallback=self._navigateToPeak,
                                                           grid=(0, 0),
                                                           hiddenColumns=['Pid'])
        self.chemicalShiftTable = ChemicalShiftTable(parent=self._chemicalShiftFrame,
                                                     mainWindow=self.mainWindow,
                                                     moduleParent=self.assignedPeaksTable,  # just to give a unique id
                                                     setLayout=True,
                                                     actionCallback=self.navigateToNmrResidueCallBack,
                                                     selectionCallback=self._selectionCallback,
                                                     grid=(0, 0),
                                                     hiddenColumns=['Pid', 'Shift list peaks', 'All peaks'])

        # disable current callback - not required for assignmentInspector
        self.chemicalShiftTable.clearCurrentCallback()
        # notifier to handle deleting items
        self.chemicalShiftTable._tableSelectionChanged.connect(self._tableSelectionCallback)

        # settingsWidget
        if chemicalShiftList is not None:
            self.chemicalShiftTable.selectChemicalShiftList(chemicalShiftList)
        elif selectFirstItem:
            firstItemText = self.chemicalShiftTable._chemicalShiftListPulldown.getFirstItemText()
            if firstItemText:
                self.chemicalShiftTable._chemicalShiftListPulldown.selectFirstItem()
                chemicalShiftList = self.chemicalShiftTable._chemicalShiftListPulldown.getSelectedObject()
                self.chemicalShiftTable._update(chemicalShiftList)

                dataFrameObject = self.chemicalShiftTable._dataFrameObject
                if len(dataFrameObject.objects) > 0:
                    self._selectByChemicalShifts([dataFrameObject.objects[0]])

        # install the event filter to handle maximising from floated dock
        self.installMaximiseEventHandler(self._maximise, self._closeModule)
        self._setNmrAtomListVisible(True)
        self._registerNotifiers()

    @QtCore.pyqtSlot(list)
    def _tableSelectionCallback(self, shifts):

        from ccpn.util.AttrDict import AttrDict

        nmrResidues = list(set(cs.nmrAtom.nmrResidue for cs in shifts))
        _temp = AttrDict()
        _temp.nmrResidues = nmrResidues

        self._highlightNmrResidues({CallBack.OBJECT: _temp})

    def _calculateMinHeight(self):
        contentsMargins = self._settingsScrollArea.contentsMargins()
        marginsTotalVertical = contentsMargins.top() + contentsMargins.bottom()
        minHeight = max(self._tickLisWidget.sizeHint().height(), self.displaysWidget.minimumSizeHint().height()) + (
                self.SETTING_PADDING * 2) + marginsTotalVertical
        return minHeight


    def _getDisplays(self):
        """
        Return list of displays to navigate - if needed
        """
        displays = []
        # check for valid displays
        gids = self.displaysWidget.getTexts()
        if len(gids) == 0: return displays
        if ALL in gids:
            displays = self.mainWindow.spectrumDisplays
        else:
            displays = [self.application.getByGid(gid) for gid in gids if gid != ALL]
        return displays

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
        self.setNotifier(self.current, [Notifier.CURRENT], targetName=NmrAtom._pluralLinkName,
                         callback=self._highlightNmrAtoms)
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

    def _setNmrAtomListVisible(self, visible=None):
        """change the visibility of the nmrAtom list in the peak tables widget
        """
        if visible is None:
            visible = self.showNmrAtomListWidget.get()
        else:
            if not isinstance(visible, bool):
                raise TypeError('nmrList visibility must be True/False')

            self.showNmrAtomListWidget.set(visible)

        self.assignedPeaksTable.nmrAtomListFrame.setVisible(visible)

    def _setCurrentPeak(self, data):
        """
        PeakTable select callback
        """
        peaks = data[CallBack.OBJECT]
        # multiselection allowed, set current to all selected peaks
        if peaks:
            self.application.current.peaks = peaks

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

    def navigateToNmrResidueCallBack(self, data):
        """Navigate in selected displays to nmrResidue; skip if none defined
        """
        # handle a single chemicalShift - SHOULD always contain an object
        objs = data[CallBack.OBJECT]
        if not objs:
            return
        if isinstance(objs, (tuple, list)):
            chemicalShift = objs[0]
        else:
            chemicalShift = objs

        self._navigateByChemicalShift(chemicalShift)

    def _navigateByChemicalShift(self, chemicalShift):
        nmrResidue = chemicalShift.nmrAtom.nmrResidue

        getLogger().debug('nmrResidue=%s' % (nmrResidue.id))

        displays = self._getDisplays()
        if len(displays) == 0:
            logger.warning('Undefined display module(s); select in settings first')
            showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
            return

        from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar

        with undoBlockWithoutSideBar():
            # optionally clear the marks
            if self.autoClearMarksWidget.checkBox.isChecked():
                self.mainWindow.clearMarks()

            # navigate the displays
            for display in displays:
                if len(display.strips) > 0:
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
            getLogger().debug('AssignmentInspector_ChemicalShift>>> selection', objList)

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
            nmrResidues = tuple(set(cs.nmrAtom.nmrResidue for cs in chemicalShifts))
        else:
            nmrResidues = []

        if nmrResidues:
            with self._notifierBlanking():
                nmrAtoms = tuple(set(nmrAtom for nmrRes in nmrResidues for nmrAtom in nmrRes.nmrAtoms))

                self.current.nmrAtoms = nmrAtoms
                self.current.nmrResidues = nmrResidues

                self._highlightChemicalShifts(nmrResidues)
                self.assignedPeaksTable._updateModuleCallback({NMRRESIDUES: nmrResidues,
                                                               NMRATOMS   : nmrAtoms},
                                                              updateFromNmrResidues=True)

    def _highlightChemicalShifts(self, nmrResidues):
        """
        Highlight chemical shifts in the table
        """
        if self.chemicalShiftTable._dataFrameObject:
            getLogger().debug('_highlightChemicalShifts ', nmrResidues)

            chemicalShifts = self.chemicalShiftTable._dataFrameObject._objects

            residues = set(nmrResidues)
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom and not cs.nmrAtom.isDeleted and cs.nmrAtom.nmrResidue in residues]

            self.chemicalShiftTable._highLightObjs(highlightList)

    def _navigateToPeak(self, data):
        """
        PeakTable double-click callback; navigate in to peak in current.strip
        """
        displays = self._getDisplays()
        markPositions = self.markPositionsWidget.checkBox.isChecked()

        if len(displays) == 0:
            logger.warning('Undefined display module(s); select in settings first')
            showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
            return

        # get the first object from the callback
        objs = data[CallBack.OBJECT]
        if not objs:
            return
        if isinstance(objs, (tuple, list)):
            peak = objs[0]
        else:
            peak = objs

        if peak:
            self.current.peak = peak

            from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar

            with undoBlockWithoutSideBar():
                # optionally clear the marks
                if self.autoClearMarksWidget.checkBox.isChecked():
                    self.application.ui.mainWindow.clearMarks()

                # navigate the displays
                for display in displays:
                    for strip in display.strips:

                        validPeakListViews = [pp.peakList for pp in strip.peakListViews if isinstance(pp.peakList, PeakList)]

                        if peak.peakList in validPeakListViews:
                            widths = None
                            if peak.peakList.spectrum.dimensionCount <= 2:
                                widths = _getCurrentZoomRatio(strip.viewRange())

                            navigateToNmrAtomsInStrip(strip, makeIterableList(peak.assignedNmrAtoms),
                                                      widths=widths, markPositions=markPositions,
                                                      setNmrResidueLabel=False)

    def _highlightNmrResidues(self, data):
        """
        Notifier Callback for highlighting all NmrAtoms in the table
        """
        if self.nmrAtomBlocking:
            return

        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._dataFrameObject:
            getLogger().debug('_highlightNmrResidues ', objList)

            chemicalShifts = self.chemicalShiftTable._dataFrameObject._objects
            nmrResidues = set(objList.nmrResidues)  #        set([atom.nmrResidue for atom in self.current.nmrAtoms if atom])
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom and not cs.nmrAtom.isDeleted and cs.nmrAtom.nmrResidue in nmrResidues]

            self.chemicalShiftTable._highLightObjs(highlightList)

            # will respond to selection of nmrAtom in sequenceGraph
            self.assignedPeaksTable._updateModuleCallback({NMRRESIDUES: nmrResidues},
                                                          updateFromNmrResidues=True)

    def _highlightNmrAtoms(self, data):
        """
        Notifier Callback for highlighting all NmrAtoms in the table
        """
        if self.nmrAtomBlocking:
            return

        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._dataFrameObject:
            getLogger().debug('_highlightNmrAtoms ', objList)

            chemicalShifts = self.chemicalShiftTable._dataFrameObject._objects
            nmrResidues = set([atom.nmrResidue for atom in self.current.nmrAtoms if atom])
            highlightList = [cs for cs in chemicalShifts if cs.nmrAtom and not cs.nmrAtom.isDeleted and cs.nmrAtom.nmrResidue in nmrResidues]

            self.chemicalShiftTable._highLightObjs(highlightList)

            # will respond to selection of nmrAtom in sequenceGraph
            self.assignedPeaksTable._updateModuleCallback({NMRRESIDUES: nmrResidues,
                                                           NMRATOMS   : self.current.nmrAtoms},
                                                          updateFromNmrResidues=False)

    def _refreshNmrAtoms(self, data):
        """
        Notifier Callback for refreshing all NmrAtoms in the table
        """
        objList = data[CallBack.OBJECT]

        if self.chemicalShiftTable._dataFrameObject:
            getLogger().debug('_refreshNmrAtoms ', objList)

            self.assignedPeaksTable._updateModuleCallback(None, updateFromNmrResidues=False)


class AssignmentInspectorTable(GuiTable):
    """
    Class to present a NmrResidue Table and a NmrChain pulldown list, wrapped in a widget
    """
    className = 'AssignmentInspectorTable'
    attributeName = 'chemicalShifts'

    OBJECT = 'object'
    TABLE = 'table'

    def __init__(self, parent=None, mainWindow=None, moduleParent=None, actionCallback=None, selectionCallback=None, **kwds):
        """
        Initialise the widgets for the module.
        """
        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        if mainWindow:
            self.application = mainWindow.application
            self.project = mainWindow.application.project
            self.current = mainWindow.application.current
        else:
            self.application = None
            self.project = None
            self.current = None

        self.sampledDims = {}  #GWV: not sure what this is supposed to do
        self.ids = []  # list of currently displayed NmrAtom ids + <all>

        # main window

        # Frame-1: NmrAtoms list, hidden when table first opens - left hand side
        # GST NOT TRUE ANY MORE?
        self._nmrAtomListFrameWidth = 150
        self.nmrAtomListFrame = Frame(parent, grid=(0, 0), gridSpan=(1, 1), setLayout=True)

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

        # Frame-2: peaks
        self.frame2 = Frame(parent, grid=(0, 1), gridSpan=(1, 1), setLayout=True)
        self.peaksLabel = Label(self.frame2, 'Peaks assigned to NmrAtom(s):', bold=True,
                                grid=(0, 0), gridSpan=(1, 4))

        self.peaksLabel.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        self.frame2.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)

        # initialise the currently attached dataFrame
        self._hiddenColumns = ['Pid']
        self.dataFrameObject = None

        # initialise the table and put in the second column (first hidden at start)
        super().__init__(parent=self.frame2,
                         mainWindow=self.mainWindow,
                         dataFrameObject=None,
                         setLayout=True,
                         autoResize=True, multiSelect=True,
                         selectionCallback=selectionCallback,
                         actionCallback=actionCallback,
                         grid=(1, 0), gridSpan=(1, 4),
                         enableDelete=False, enableSearch=False
                         )
        self.moduleParent = moduleParent

        # self._registerNotifiers() - removed for testing

        self._peakList = None
        self._nmrResidues = []
        self._nmrAtoms = []

        # update if current.nmrResidue is defined
        if self.application.current.nmrResidue is not None and self.application.current.chemicalShiftList is not None:
            self._updateModuleCallback({NMRRESIDUES: [self.application.current.nmrResidue]},
                                       updateFromNmrResidues=True)

        # set the required table notifiers
        self.setTableNotifiers(tableClass=None,
                               rowClass=Peak,
                               cellClassNames=None,
                               tableName='peakList', rowName='peak',
                               changeFunc=self._refreshTable,
                               className=self.attributeName,
                               updateFunc=self._refreshTable,
                               tableSelection='_peakList',
                               pullDownWidget=None,
                               callBackClass=Peak,
                               selectCurrentCallBack=self._selectOnTableCurrentPeaksNotifierCallback,
                               moduleParent=moduleParent)

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
            with self.blockWidgetSignals(self.attachedNmrAtomsList):
                self.attachedNmrAtomsList.clear()

                # NOTE:ED - should we only display those nmrAtoms that are in the chemicalShift table? - similarly for peaks
                self.ids = [atm.id for atm in nmrAtoms if not (atm.isDeleted or atm._flaggedForDelete)]

                self.attachedNmrAtomsList.addItems(self.ids)
                self.attachedNmrAtomsList.selectItems(_select)

            # populate peak table with the correct peaks
            self._updatePeakTable(nmrAtoms, messageAll=updateFromNmrResidues)

        else:
            # data is None, so update from current settings - should only be called from _refreshNmrAtoms

            _select = self.attachedNmrAtomsList.getSelectedTexts()
            # there is currently a hidden list widget containing the nmrAtom ids
            with self.blockWidgetSignals(self.attachedNmrAtomsList):
                self.attachedNmrAtomsList.clear()
                _nmrAtoms = [atm for _nmrRes in self._nmrResidues if not (_nmrRes.isDeleted or _nmrRes._flaggedForDelete) for atm in _nmrRes.nmrAtoms if not (atm.isDeleted or atm._flaggedForDelete)]
                self.ids = [atm.id for atm in _nmrAtoms]
                self.attachedNmrAtomsList.addItems(self.ids)
                self.attachedNmrAtomsList.selectItems(_select)

            # populate peak table with the correct peaks
            self._updatePeakTable(_nmrAtoms, messageAll=updateFromNmrResidues)

    def _selectOnTableCurrentPeaksNotifierCallback(self, data):
        """
        Callback from a notifier to highlight the peaks on the peak table
        :param data:
        """

        currentPeaks = data['value']
        self._selectOnTableCurrentPeaks(currentPeaks)

    def _selectOnTableCurrentPeaks(self, currentPeaks):
        """
        Highlight the list of peaks on the table
        :param currentPeaks:
        """
        self.highlightObjects(currentPeaks)

    def _updatePeakTableCallback(self, data=None):
        """
        Update the peakTable using item.text (which contains a NmrAtom pid or <all>)
        """
        # something has been selected, so get all selected items
        numTexts = self.attachedNmrAtomsList.count()
        selectedTexts = self.attachedNmrAtomsList.getSelectedTexts()
        nmrAtoms = [self.project.getByPid('NA:' + id) for id in selectedTexts]

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
            nmrAtoms = [self.project.getByPid('NA:' + id) for id in self.attachedNmrAtomsList.getTexts()]

            # # populate the table with valid nmrAtoms
            # self._updatePeakTable([atm for atm in nmrAtoms if atm is not None], messageAll=True)

        self._peakList = _emptyObject()
        allPeaks = list(set([pk for nmrAtom in nmrAtoms if nmrAtom for pk in nmrAtom.assignedPeaks]))
        chemicalShiftList = self.moduleParent.chemicalShiftTable._chemicalShiftListPulldown.getSelectedObject()
        if chemicalShiftList:
            spectra = chemicalShiftList.spectra #show peaks only for spectra currently available for the selected CSL
            peaks = [peak for peak in allPeaks if peak.peakList.spectrum in spectra]
        else:
            peaks = []

        self._peakList.peaks = peaks

        self.populateTable(rowObjects=self._peakList.peaks,
                           columnDefs=self.getColumns())

        ids = [atm.id for atm in nmrAtoms]

        if messageAll:
            self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % ALL)
        else:
            atomList = ', '.join([str(id) for id in ids])
            self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % atomList)  # nmrAtom.id)

        # if messageAll:
        #
        #     self._peakList = _emptyObject()
        #     self._peakList.peaks = list(set([pk for nmrAtom in nmrAtoms for pk in nmrAtom.assignedPeaks]))
        #
        #     with self._projectBlanking():
        #         self._dataFrameObject = self.getDataFrameFromList(table=self,
        #                                                           buildList=self._peakList.peaks,
        #                                                           colDefs=self.getColumns(),
        #                                                           hiddenColumns=self._hiddenColumns)
        #
        #     # self.assignedPeaksTable.setObjects(peaks)
        #     # highlight current.nmrAtom in the list widget
        #     self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
        #     self.peaksLabel.setText('Peaks assigned to NmrAtom(s): %s' % ALL)
        # else:
        #     pid = 'NA:' + id
        #     nmrAtom = self.application.project.getByPid(pid)
        #
        #     if nmrAtom is not None:
        #         with self._projectBlanking():
        #             self._dataFrameObject = self.getDataFrameFromList(table=self,
        #                                                               buildList=nmrAtom.assignedPeaks,
        #                                                               colDefs=self.getColumns(),
        #                                                               hiddenColumns=self._hiddenColumns)
        #
        #         # self.assignedPeaksTable.setObjects(nmrAtom.assignedPeaks)
        #         # highlight current.nmrAtom in the list widget
        #
        #         self.attachedNmrAtomsList.setCurrentRow(self.ids.index(id))
        #         atomList = ', '.join([str(id) for id in self.ids if id != ALL])
        #         self.peaksLabel.setText('Peaks assigned to NmrAtom: %s' % id)  # nmrAtom.id)

    def getColumns(self):
        """get columns for initialisation of table"""
        columns = ColumnClass([('Peak', lambda pk: pk.id, '', None, None),
                               ('Pid', lambda pk: pk.pid, 'Pid of peak', None, None),
                               ('_object', lambda pk: pk, 'Object', None, None),
                               ('serial', lambda pk: pk.serial, '', None, None)])
        tipTexts = []
        # get the maxmimum number of dimensions from all spectra in the project
        numDim = max([sp.dimensionCount for sp in self.application.project.spectra] + [1])

        for i in range(numDim):
            j = i + 1
            c = Column('Assign F%d' % j,
                       lambda pk, dim=i: getPeakAnnotation(pk, dim),
                       'NmrAtom assignments of peak in dimension %d' % j,
                       None,
                       None)
            columns._columns.append(c)

            # columns.append(c)
            # tipTexts.append('NmrAtom assignments of peak in dimension %d' % j)

            sampledDim = self.sampledDims.get(i)
            if sampledDim:
                text = 'Sampled\n%s' % sampledDim.conditionVaried
                tipText = 'Value of sampled plane'
                unit = sampledDim

            else:
                text = 'Pos F%d' % j
                tipText = 'Peak position in dimension %d' % j
                unit = 'ppm'

            c = Column(text,
                       lambda pk, dim=i, unit=unit: getPeakPosition(pk, dim, unit),
                       tipText,
                       None,
                       '%0.3f')
            columns._columns.append(c)

            # columns.append(c)
            # tipTexts.append(tipText)
        columns._columns.append(Column('height', lambda pk: pk.height, '', None, None))
        columns._columns.append(Column('volume', lambda pk: pk.volume, '', None, None))
        return columns

    def _getPeakHeight(self, peak):
        """
        Returns the height of the specified peak as formatted string or 'None' if undefined
        """
        if peak.height:
            return '%7.2E' % float(peak.height)
        else:
            return '%s' % None

    def _refreshTable(self, *args):
        self.update()
