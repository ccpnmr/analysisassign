"""
AtomSelector module: 

used for backbone and sidechain assignments to set the nmrAtom of a (number of) selected
peaks to a specific nucleus. Thus far, it set the 'y-axis' dimension of the peak (to be
made flexible in the settings tab). Options other than 'protein' are not yet implemented.

Original by SS
First rework by GWV
Reworked by EJB
"""

#TODO Needs cleanup
"""
- buttons are destroyed (!?) and create with every refresh; 
- atom type prediction in sidechain mode (when type of nmrResidue is not yet known) gives
  the result of the last residue (i.e. Tyr), rather then an average
- Should retain predictions of selected peaks when switching backbone/sidechain
# - assignSelected should be refactored to be proper for all general cases. Done

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
__dateModified__ = "$dateModified: 2024-03-21 11:51:37 +0000 (Thu, March 21, 2024) $"
__version__ = "$Revision: 3.2.2 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import typing
import copy
from PyQt5 import QtCore, QtGui, QtWidgets
from contextlib import contextmanager
from ccpn.core.Peak import Peak
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.NmrAtom import UnknownIsotopeCode
from ccpn.core.lib import Pid
from ccpn.core.lib.AssignmentLib import isInterOnlyExpt, getNmrAtomPrediction, CCP_CODES
from ccpn.core.lib.AssignmentLib import peaksAreOnLine
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets.PulldownList import PulldownList
from ccpn.ui.gui.widgets.RadioButton import RadioButton
from ccpn.ui.gui.widgets.RadioButtons import RadioButtons
from ccpn.ui.gui.widgets.Widget import Widget
from ccpn.ui.gui.widgets.Spacer import Spacer
from ccpn.ui.gui.widgets.Frame import Frame, ScrollableFrame
from ccpn.ui.gui.widgets.MessageDialog import showWarning
from ccpn.ui.gui.widgets.PulldownListsForObjects import NmrResiduePulldown, NmrChainPulldown

from ccpn.ui.gui.lib.GuiNotifier import GuiNotifier
from ccpn.ui.gui.widgets.DropBase import DropBase

from ccpn.util.Common import makeIterableList, _truncateText
from ccpn.core.lib.AssignmentLib import PROTEIN_NEF_ATOM_NAMES, NEF_ATOM_NAMES_SORTED
from ccpn.util.Logging import getLogger
from ccpn.util.Constants import DNA_ATOMS, RNA_ATOMS, DNA_ATOM_NAMES, RNA_ATOM_NAMES, ALL_DNARNA_ATOMS_SORTED
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.lib.AssignmentLib import _assignNmrAtomsToPeaks


logger = getLogger()

# TODO:ED Add DNA, RNA structures to the list
# MOLECULE_TYPES = ['protein', 'DNA', 'RNA', 'carbohydrate', 'other']
MOLECULE_TYPES = ['protein']
BACKBONEATOMS = ['H', 'N', 'CA', 'CB', 'C', 'HA', 'HB']
ADDITIONALBACKBONEATOMS = ['H', 'N', 'C']

MSG = '<Not-defined. Select any to start>'
PROTEIN_MOLECULE = 'protein'
DNA_MOLECULE = 'DNA'
RNA_MOLECULE = 'RNA'

BUTTON_MINX = 80  # should be based on font-size
BUTTON_MINY = 24

DEFAULT_BUTTON = """QRadioButton { background-color: %s }
                   QRadioButton::hover { background-color: %s}""" % ('lightgrey', 'white')
DEFAULT_COLOURS = ('whitesmoke', 'lightgrey')
GREEN_COLOURS = ('palegreen', 'mediumseagreen')
ORANGE_COLOURS = ('gold', 'orange')
RED_COLOURS = ('lightpink', 'tomato')


class _RButton(RadioButton):
    _enterColour = 'white'
    _leaveColour = 'lightgrey'

    def mouseReleaseEvent(self, e: QtGui.QMouseEvent) -> None:
        self.group()._parent._nmrAtomButtonsCallback(self)

    def enterEvent(self, a0: QtCore.QEvent) -> None:
        if self._enterColour:
            self.setStyleSheet(f'QRadioButton {{ background-color: {self._enterColour} }}')

        super(_RButton, self).enterEvent(a0)

    def leaveEvent(self, a0: QtCore.QEvent) -> None:
        if self._leaveColour:
            self.setStyleSheet(f'QRadioButton {{ background-color: {self._leaveColour} }}')

        super(_RButton, self).leaveEvent(a0)

    def setBackgroundColours(self, enterColour, leaveColour):
        """Set the enter/leave-event colours
        """
        self._enterColour = enterColour
        self._leaveColour = leaveColour
        if self._leaveColour:
            self.setStyleSheet(f'QRadioButton {{ background-color: {self._leaveColour} }}')


#=========================================================================================
# NmrAtomAssignerModule
#=========================================================================================

class NmrAtomAssignerModule(CcpnModule):
    """
    Module to be used with PickAndAssignModule for prediction of nmrAtom names and assignment of nmrAtoms
    to peak dimensions
    Responds to current.nmrResidue and current.peaks; accepts nmrResidue pid drops
    """
    className = 'NmrAtomAssignerModule'

    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    defaultSettingsState = 0
    settingsPosition = 'top'

    def __init__(self, mainWindow=None, name='NmrAtomAssigner', nmrAtom=None):

        super().__init__(mainWindow=mainWindow, name=name)

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        if mainWindow:
            self.application = mainWindow.application
            self.project = mainWindow.application.project
            self.current = mainWindow.application.current

        # module attributes
        self._thisPeaks = self.current and self.current.peaks
        self._thisNmrChain = self.current and self.current.nmrResidue and self.current.nmrResidue.nmrChain
        self._thisNmrResidue = self.current and self.current.nmrResidue
        self._thisShift = 0.0
        self._thisDim = 0
        self._thisAxes = None

        self._setWidgets()

        self.buttons = {}
        self._registerNotifiers()
        self._updateWidget()

    def _setWidgets(self):
        """Set up the widgets
        """
        # Settings Widget
        self._ASwidget = Widget(self.settingsWidget, setLayout=True,
                                grid=(0, 0), vAlign='top', hAlign='left')
        row = 0
        self.selectionLabel = Label(self._ASwidget, 'Select NmrAtom by', grid=(row, 0))
        self.selectionRadioButtons = RadioButtons(self._ASwidget, texts=['AtomType', 'Axis'], selectedInd=1,
                                                  callback=self._selectionCallback, grid=(row, 1))
        self.selectAtomType, self.selectAxisCode = self.selectionRadioButtons.radioButtons
        # pulldown for Molecule type
        row += 1
        self.molTypeLabel = Label(self._ASwidget, 'Molecule Type', grid=(row, 0))
        self.molTypePulldown = PulldownList(self._ASwidget, grid=(row, 1), texts=MOLECULE_TYPES,
                                            callback=self._changeMoleculeType)
        row += 1
        self.modeTypeLabel = Label(self._ASwidget, 'Mode', grid=(row, 0))
        self.modeRadioButtons = RadioButtons(self._ASwidget, texts=['Backbone', 'All'], selectedInd=0,
                                             callback=self._createButtonsCallback, grid=(row, 1))
        self.selectBackboneButton = self.modeRadioButtons.getRadioButton('Backbone')
        self.selectAllButton = self.modeRadioButtons.getRadioButton('All')
        # modifiers for side-chain
        row += 1
        self.offsetLabel = Label(self._ASwidget, 'Offset', grid=(row, 0))
        self.offsetSelector = PulldownList(self._ASwidget, grid=(row, 1), texts=['0', '-1', '+1'],
                                           callback=self._offsetPullDownCallback)
        self._sidechainModifiers = [self.offsetLabel, self.offsetSelector]
        # set size policies to allow the main widget to overlap the settings, cleaner display
        self._ASwidget.setMinimumSize(self._ASwidget.sizeHint())
        self.settingsWidget.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.settingsWidget.setContentsMargins(5, 5, 5, 5)
        for w in self._sidechainModifiers:
            w.hide()
        self._residueFrame = ScrollableFrame(parent=self.mainWidget,
                                             showBorder=False, setLayout=True,
                                             acceptDrops=True, grid=(0, 0), gridSpan=(1, 1), spacing=(5, 5))
        self._scrollAreaWidget = self._residueFrame._scrollArea
        self._residueFrame.insertCornerWidget()
        # self._residueFrame = Frame(self.mainWidget, setLayout=True, acceptDrops=True, showBorder=False, spacing=(5, 5))
        self._residueFrame.setContentsMargins(5, 5, 5, 5)

        resRow = 0
        _f = Frame(self._residueFrame, setLayout=True, showBorder=False, grid=(resRow, 0), gridSpan=(1, 3))
        self._peaksLabel = Label(_f, 'Assigning Peak(s):', bold=True, grid=(0, 0), hPolicy='minimal')
        self.currentPeaksLabel = Label(_f, grid=(0, 1), gridSpan=(1, 2), hPolicy='minimal', hAlign='l')
        resRow += 1

        _f = Frame(self._residueFrame, setLayout=True, showBorder=False, grid=(resRow, 0), gridSpan=(1, 3))
        self._nmrChainPulldown = NmrChainPulldown(_f, mainWindow=self.mainWindow,
                                                  labelText='NmrChain:', showSelectName=True,
                                                  setCurrent=False,
                                                  callback=self._nmrChainPullDownCallback,
                                                  grid=(0, 0), hPolicy='minimal', minimumWidths=None,
                                                  sizeAdjustPolicy=QtWidgets.QComboBox.AdjustToContents)

        self._nmrResiduePulldown = NmrResiduePulldown(_f, mainWindow=self.mainWindow,
                                                      labelText='NmrResidue:', useIds=False, showSelectName=False,
                                                      setCurrent=True, followCurrent=True,
                                                      filterFunction=self._filterResidues,
                                                      grid=(0, 1), hPolicy='minimal', minimumWidths=None,
                                                      sizeAdjustPolicy=QtWidgets.QComboBox.AdjustToContents)

        self._newNmrResidueButton = Button(_f, text='New', grid=(0, 2), gridSpan=(1, 1),
                                           callback=self._newNmrResidueCallback, hPolicy='minimal')
        self._newNmrResidueButton.setToolTip('Create new nmrResidue in current chain')
        self._nmrResidueEditButton = Button(_f, text='Edit', grid=(0, 3), gridSpan=(1, 1),
                                            callback=self._nmrResidueEditCallback, hPolicy='minimal')
        self._nmrResidueEditButton.setToolTip('Edit current nmrResidue')
        Spacer(_f, 2, 2, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed,
               grid=(1, 4), gridSpan=(1, 1))

        resRow += 1
        self._labelFrame = Frame(self._residueFrame, setLayout=True, showBorder=False, grid=(resRow, 0),
                                 gridSpan=(1, 4))
        labRow = 0
        # modifier for atomCode
        self.axisCodeLabel = Label(self._labelFrame, 'Assign by axis', grid=(labRow, 0))
        self.axisCodeOptions = RadioButtons(self._labelFrame, selectedInd=0, texts=['C'],
                                            callback=self._changeAxisCode, grid=(labRow, 1),
                                            halign='l')
        # labRow += 1
        # modifier for atomType - overlay the above
        self.atomTypeLabel = Label(self._labelFrame, 'Assign by atomType', grid=(labRow, 0))
        self.atomTypeOptions = RadioButtons(self._labelFrame, selectedInd=1, texts=['H', 'C', 'N', 'Other'],
                                            callback=self._changeAtomType, grid=(labRow, 1),
                                            halign='l')

        Spacer(self._labelFrame, 16, 2, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed,
               grid=(labRow, 2), gridSpan=(1, 1))
        self.shiftlabel = Label(self._labelFrame, 'Shift/diverge', grid=(labRow, 3))

        labRow += 1
        Spacer(self._labelFrame, 2, 2, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed,
               grid=(labRow, 4), gridSpan=(1, 1))

        resRow += 1
        self._assignWidget = Frame(self._residueFrame, setLayout=True, showBorder=False, grid=(resRow, 0), spacing=(5, 5))
        self._scrollAreaWidget.setWidget(self._residueFrame)
        resRow += 1
        # add spacer to stop columns changing width
        Spacer(self._residueFrame, 2, 2, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding,
               grid=(resRow, 4), gridSpan=(1, 1))
        resRow += 1
        self.buttonGroup = QtWidgets.QButtonGroup()
        #self.buttonGroup.buttonClicked.connect(self._nmrAtomButtonsCallback)
        self.buttonGroup.setExclusive(False)
        self.buttonGroup._parent = self

    def _registerNotifiers(self):
        """Register notifiers for the module
        """
        # self.setNotifier(self.project, [Notifier.CHANGE, Notifier.CREATE, Notifier.DELETE],
        #                  NmrAtom.className, callback=self._nmrResidueCallBack, onceOnly=True)
        # self.setNotifier(self.project, [Notifier.CHANGE],
        #                  Peak.className, callback=self._nmrResidueCallBack, onceOnly=True)

        # update on current peak/nmrResidue change
        self.setCurrentNotifier(Peak._pluralLinkName, callback=self._currentPeaksCallback, onceOnly=True)
        self.setCurrentNotifier(NmrResidue._pluralLinkName, callback=self._currentNmrResiduesCallback, onceOnly=True)

        # notifiers for name-change, peak-update
        self.setNotifier(self.project, [Notifier.RENAME],
                         'NmrResidue', self._updateNmrResidue, onceOnly=True)
        self.setNotifier(self.project, [Notifier.CHANGE],
                         targetName=Peak.className,
                         callback=self._updateShiftFromPeaks,
                         onceOnly=True)

        # allow drop-event
        self.setGuiNotifier(self._residueFrame, [GuiNotifier.DROPEVENT],
                            [DropBase.PIDS], callback=self._handleNmrResidue)

    def _unRegisterNotifiers(self):
        """clean up the notifiers
        """
        # _closeModule() will do most of them
        self._nmrResiduePulldown.unRegister()
        self._nmrChainPulldown.unRegister()

    def _closeModule(self):
        self._unRegisterNotifiers()
        super()._closeModule()

    #================================================================================================================
    # callbacks and functionalities
    #================================================================================================================

    def _filterResidues(self, pids):
        """Filter function for the residue pulldown
        Add any nmrResidue defined by selected peaks at the top of the list
        """

        # first time hack, as during initialising this routine is called to populate
        # the pulldown; however _nmrResidue is not yet defined
        # if not hasattr(self, '_nmrResiduePulldown'):
        #     return pids

        if not (self._thisNmrChain and self._thisNmrChain):
            return pids

        # nmrChain = self._nmrChainPulldown.getSelectedObject()

        # For selected peaks: get the pids of nmrResidues of assigned nmrAtoms
        newPids = []
        for peak in self.current.peaks:
            if peak:
                for assignment in peak.assignments:
                    for nmrAtom in assignment:
                        if nmrAtom:
                            newPids.append(nmrAtom.nmrResidue.pid)
        newPids = sorted(set(newPids))

        def _isOk(pid):
            if self._thisNmrChain is None:
                # No filtering
                return True
            # nmrResidue = self._nmrResiduePulldown.value2object(pid)
            nmrRes = self.project.getByPid(pid)
            return (nmrRes is not None and nmrRes.nmrChain == self._thisNmrChain)

        newPids = newPids + [pid for pid in pids if _isOk(pid)]
        return newPids

    def _newNmrResidueCallback(self):
        """Callback to create a new nmrResidue and add to the current chain
        """
        # if nmrChain := self._nmrChainPulldown.getSelectedObject():
        #     # create a new nmrResidue as required
        #     nmrResidue = self._fetchNmrResidue(nmrChain)
        #     self.current.nmrResidue = nmrResidue

        if self._thisNmrChain:
            # create a new nmrResidue as required
            self._thisNmrResidue = self._fetchNmrResidue(self._thisNmrChain)
            self.current.nmrResidue = self._thisNmrResidue

    def _nmrResidueEditCallback(self, data):
        """Callback to edit the current nmrResidue
        """
        # call popup on current nmrResidue
        from ccpn.ui.gui.popups.NmrResiduePopup import NmrResidueEditPopup

        popup = NmrResidueEditPopup(parent=self.mainWindow, mainWindow=self.mainWindow,
                                    obj=self.current.nmrResidue)
        popup.exec_()

    def _nmrChainPullDownCallback(self, value):
        """Callback for the NmrChain selection
        """
        self._thisNmrChain = self.project.getByPid(value)

        self._nmrResiduePulldown.update()

    def _handleNmrResidue(self, dataDict):
        """drop event handler to accept NmrResidue pids
        """
        if pids := dataDict.get(DropBase.PIDS):
            objs = [self.project.getByPid(pid) for pid in pids]
            if nmrResidues := [obj for obj in objs if obj is not None and isinstance(obj, NmrResidue)]:
                self.current.nmrResidues = nmrResidues

    def _togglePressedButton(self, pressedButton=None):
        """Ensures only a button at the time is checked, yet allows to uncheck a radio button. If pressedButton is None: unchecks all
        """
        if pressedButton:
            pressedButton.setChecked(not pressedButton.isChecked())

        for button in self.buttonGroup.buttons():
            if button != pressedButton:
                button.setChecked(False)

    def _nmrAtomButtonsCallback(self, pressedButton):

        self._togglePressedButton(pressedButton)

        from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar

        with undoBlockWithoutSideBar():
            try:
                if pressedButton.isChecked():
                    self._assignSelected(atomName=pressedButton._atomName, offSet=pressedButton._offSet)
                else:
                    self._deassignSelected(atomName=pressedButton._atomName, offSet=pressedButton._offSet)

            except Exception as es:
                showWarning(str(self.windowTitle()), str(es))
                self._togglePressedButton()  # uncheck all if any error

    def _deassignSelected(self, atomName, offSet):
        nmrResidue = self._getCorrectResidue(self.current.nmrResidue, offSet, atomName)
        if nmrAtom := nmrResidue.getNmrAtom(atomName.translate(Pid.remapSeparators)):
            self.deassignAtomFromSelectedPeaks(self.current.peaks, nmrAtom)

    def _assignSelected(self, atomName, offSet):
        nmrResidue = self._getCorrectResidue(self.current.nmrResidue, offSet, atomName)
        if not nmrResidue:
            getLogger().warning('Error creating new nmrResidue')
            raise ValueError('Error creating new nmrResidue')

        if nmrAtom := nmrResidue.fetchNmrAtom(name=atomName):
            if self.selectAxisCode.isChecked():
                self.assignNmrAtomsToPeaks(nmrAtom=nmrAtom, peaks=self.current.peaks)
            else:
                _assignNmrAtomsToPeaks(nmrAtoms=[nmrAtom], peaks=self.current.peaks)

            self._setCheckedButtonOfAssignedAtoms(nmrResidue,
                                                  offSet=offSet)  # this so that only assigned atoms are checked.

    def _createButtonsCallback(self):
        self._updateWidget()
        return

    # def _nmrResidueCallBack(self, data=None):
    #     "Callback if current.nmrResidue changes"
    #
    #     if self.current.nmrResidue:
    #         self._updateWidget()
    #     else:
    #         self._assignWidgetHide()
    #         # self.currentNmrResidueLabel.setText(MSG)
    #         # self._nmrResiduePulldown.select(MSG)

    def _assignWidgetShow(self):
        # show the required widgets
        self._assignWidget.show()
        self._showSelectionButtons()

    def _assignWidgetHide(self):
        # hide all the widgets
        self._assignWidget.hide()
        self.atomTypeLabel.hide()
        self.atomTypeOptions.hide()
        self.axisCodeLabel.hide()
        self.axisCodeOptions.hide()
        self.shiftlabel.hide()

    def _offsetPullDownCallback(self, tmp=None):
        """Callback if offset pullDown changes
        """
        if self._thisNmrResidue:
            self._updateWidget()

    def _setPeaksLabel(self):
        """ update the peaks label from current-peaks
        """
        pks = self._thisPeaks
        if pks and None not in pks:
            splitter = ', '
            pText = _truncateText(splitter.join([p.id for p in pks]), splitter=splitter)
            self.currentPeaksLabel.setToolTip(splitter.join([p.id for p in pks]))
            self.currentPeaksLabel.setText(pText)
        else:
            self.currentPeaksLabel.setText(MSG)

    def _setPeakAxisCodes(self, peaks):
        if not peaks:
            return

        maxLen = 0
        refAxisCodes = None
        for peak in peaks:
            if len(peak.axisCodes) > maxLen:
                maxLen = len(peak.axisCodes)
                refAxisCodes = list(peak.axisCodes)

        if not maxLen:
            return

        axisLabels = [set() for _ in range(maxLen)]

        mappings = {}
        for peak in peaks:
            matchAxisCodes = peak.axisCodes

            self.axisCodeMapping(mappings, matchAxisCodes, refAxisCodes)
            self.axisCodeMapping(mappings, refAxisCodes, matchAxisCodes)

            # example of mappings dict - includes mapping from both sides
            # ('Hn', 'C', 'Nh')
            # {'Hn': {'Hn'}, 'Nh': {'Nh'}, 'C': {'C'}}
            # {'Hn': {'H', 'Hn'}, 'Nh': {'Nh'}, 'C': {'C'}}
            # {'CA': {'C'}, 'Hn': {'H', 'Hn'}, 'Nh': {'Nh'}, 'C': {'CA', 'C'}}
            # {'CA': {'C'}, 'Hn': {'H', 'Hn'}, 'Nh': {'Nh'}, 'C': {'CA', 'C'}}

        self.peakIndex = {}
        # go through the peaks
        for peak in peaks:
            self.peakIndex[peak] = [0 for _ in range(len(peak.axisCodes))]

            # get the peak dimension axisCode, nd see if is already there
            for peakDim, peakAxis in enumerate(peak.axisCodes):

                if peakAxis in refAxisCodes:
                    self.peakIndex[peak][peakDim] = refAxisCodes.index(peakAxis)
                    axisLabels[self.peakIndex[peak][peakDim]].add(peakAxis)

                else:
                    # if the axisCode is not in the reference list then find the mapping from the dict
                    for k, v in mappings.items():
                        if peakAxis in v:
                            # refAxisCodes[dim] = k
                            self.peakIndex[peak][peakDim] = refAxisCodes.index(k)
                            axisLabels[refAxisCodes.index(k)].add(peakAxis)

        # peakCodes = set()
        # for peak in peaks:
        #     # for code in peak.peakList.spectrum.isotopeCodes:
        #     for code in peak.axisCodes:
        #         peakCodes.add(code)
        # peakCodes = sorted(list(peakCodes), key=CcpnSorting.stringSortKey)
        #
        # # peakCodes = peaks[0].peakList.spectrum.spectrumDisplay.axisCodes
        # peakCodes = ['H', 'C', 'N', 'Other']

        axisLabels = [', '.join(ax) for ax in axisLabels]
        self.axisCodeOptions.setButtons(texts=axisLabels, tipTexts=axisLabels, silent=True)
        if not self.axisCodeOptions.getSelectedText():
            self.axisCodeOptions.setIndex(0, blockSignals=True)

    @staticmethod
    def axisCodeMapping(mappings, matchAxisCodes, refAxisCodes):
        from ccpn.core.lib.AxisCodeLib import getAxisCodeMatch

        mapping = getAxisCodeMatch(refAxisCodes, matchAxisCodes)
        for k, v in mapping.items():
            if v not in mappings:
                mappings[v] = {k}
            else:
                mappings[v].add(k)

    def _setPeakAtomCodes(self):
        atomCodes = ['H', 'C', 'N', 'Other']
        self.axisCodeOptions.setButtons(texts=list(atomCodes), tipTexts=list(atomCodes), silent=True)
        if not self.axisCodeOptions.getSelectedText():
            self.axisCodeOptions.setIndex(0, blockSignals=True)

    def _blockEvents(self):
        """Block all updates/signals/notifiers in the module.
        """
        # self.setUpdatesEnabled(False)
        self.blockSignals(True)

    def _unblockEvents(self):
        """Unblock all updates/signals/notifiers in the module.
        """
        self.blockSignals(False)
        # self.setUpdatesEnabled(True)

    @contextmanager
    def _moduleBlocking(self):
        """Context manager to handle blocking, unblocking of the module.
        """
        self._blockEvents()
        try:
            # pass control to the calling function
            yield

        except Exception as es:
            raise es
        finally:
            self._unblockEvents()

    def _updateNmrResidue(self, data):
        """Update the widget after renaming an nmrResidue
        """
        if (nmrRes := data[Notifier.OBJECT]) and nmrRes == self._thisNmrResidue:
            # re-populate the pulldown texts
            self._updateWidget()

    def _updateWidget(self, data=None):  # also used as notifier callback function
        """Update the widget to reflect the proper state
        """
        with self._moduleBlocking():
            # populate the label with the peak-pids
            self._setPeaksLabel()

            nmrRes = self._thisNmrResidue
            self._assignWidgetHide()
            if nmrRes is None:
                return

            self._nmrResiduePulldown.select(nmrRes.pid)

            pks = self._thisPeaks
            # if pks and None not in pks:
            #     self._setPeakAxisCodes(pks)

            if self.selectBackboneButton.isChecked():
                for w in self._sidechainModifiers:
                    w.hide()
                self._createBackBoneButtons()
            elif self.selectAllButton.isChecked():
                for w in self._sidechainModifiers:
                    w.show()
                self._createSideChainButtons()
            self._setCheckedButtonOfAssignedAtoms(nmrRes)

            # add a spacer to the radiobutton box - probably not a good thing to do here :|
            #   use new table-widget instead?
            # Spacer(self._assignWidget, 3, 3,
            #        QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding,
            #        grid=(30, 30), gridSpan=(1, 1))

            if pks and None not in pks:
                # temporarily restrict highlighting to backbone only, as the side-chain version needs improving
                if self.selectBackboneButton.isChecked():
                    self._predictHighlight(pks)
                elif self.selectAllButton.isChecked():
                    self._predictHighlight(pks)
                self._assignWidgetShow()

            # fill the shift-label
            if pks and (peak := pks[0]):
                if None not in pks:
                    self._setPeakAxisCodes(pks)
                self._updatePeakShiftLabel(peak)
                self._updatePeakPulldown(peak)

    def _updatePeakShiftLabel(self, peak):
        """Update the shift value in the shift-label
        """
        dim = self._getValidAxisCodeIndex()
        if self.selectAxisCode.isChecked():
            nmrs = makeIterableList(peak.dimensionNmrAtoms[dim]) if 0 <= dim < len(peak.dimensionNmrAtoms) else []
        else:
            nmrs = makeIterableList(peak.assignedNmrAtoms)

        label = self.axisCodeOptions.get()
        if nmrs:
            # if there are assignments then use the average of the chemical-shifts
            if chs := [shift.value for nmr in nmrs for shift in nmr.chemicalShifts if not shift.isDeleted and shift.value is not None]:
                shift = sum(chs) / len(chs)
                self.shiftlabel.setText(f'{label}: {shift:.3f}')
            else:
                self.shiftlabel.setText(f'{label}: None')

        else:
            # use the ppmPosition for the specified dimension
            ppm = peak.ppmPositions
            if dim < len(ppm) and ppm[dim] is not None:
                self.shiftlabel.setText(f'{label}: {ppm[dim]:.3f}')
            else:
                self.shiftlabel.setText(f'{label}:')

    def _updatePeakPulldown(self, peak):
        """Update the pulldowns from specified peak
        """
        dim = self._getValidAxisCodeIndex()
        if self.selectAxisCode.isChecked():
            nmrs = makeIterableList(peak.dimensionNmrAtoms[dim]) if 0 <= dim < len(peak.dimensionNmrAtoms) else []
        else:
            nmrs = makeIterableList(peak.assignedNmrAtoms)

        if nmrs:
            self._nmrResiduePulldown.select(nmrs[0].nmrResidue)
            self._nmrChainPulldown.select(nmrs[0].nmrResidue.nmrChain.pid)
            # self.current.nmrResidue = nmrs[0].nmrResidue

    @staticmethod
    def _removeOffsetFromButtonText(text: str):
        p = text.split(' ')
        return p[0] if p else text

    def _setCheckedButtonOfAssignedAtoms(self, nmrResidue, offSet='0'):
        """setChecked the radioButton Of Assigned Nmr Atoms.
        This makes sure that if a peak is selected and assigned to an nmrAtom, the relative button is checked """

        if not self.current.peak:
            return
        if not nmrResidue:
            return

        peaks = self.current.peaks
        currentDisplayedButtons = self.buttonGroup.buttons()
        buttonsToCheck = []

        currentAxis = self._getValidAxisCodeIndex()

        for peak in peaks:
            counts = set()
            if self.selectAxisCode.isChecked():
                peakList = makeIterableList(peak.dimensionNmrAtoms[currentAxis]) if currentAxis < len(peak.dimensionNmrAtoms) else []
            else:
                peakList = makeIterableList(peak.assignedNmrAtoms)

            for assignedNmrAtom in peakList:  #makeIterableList(peak.assignedNmrAtoms):
                if assignedNmrAtom in nmrResidue.nmrAtoms:
                    for button in currentDisplayedButtons:
                        if offSet == '0':
                            if assignedNmrAtom and assignedNmrAtom.name == button.getText():
                                counts.add(button)
                        elif button._offSet == offSet:
                            if assignedNmrAtom and assignedNmrAtom.name == button._atomName:
                                counts.add(button)

                else:  #Try to search in + and - 1 offset
                    for offset in ['-1', '+1']:
                        if r := self._getNmrResidue(nmrResidue.nmrChain, sequenceCode=nmrResidue.mainNmrResidue.sequenceCode + offset):
                            if assignedNmrAtom in r.nmrAtoms:
                                for button in currentDisplayedButtons:
                                    if assignedNmrAtom:
                                        btext = self.atomLabel(assignedNmrAtom.name, offset)
                                        if btext == button.getText():
                                            counts.add(button)
            buttonsToCheck.append(list(counts))

        buttonsToCheck = makeIterableList(buttonsToCheck)
        if len(buttonsToCheck) >= len(peaks):
            for b in currentDisplayedButtons:
                if b in buttonsToCheck:
                    b.setChecked(True)
                else:
                    b.setChecked(False)
        else:
            self._togglePressedButton()

    def _getValidAxisCodeIndex(self):
        return self.axisCodeOptions.getIndex() or 0

    def _getValidAxisCode(self, numChars=1):
        """Get the valid axis code from the buttons, numChars is included as this may be needed for DNA/RNA
        """
        code = self.axisCodeOptions.getSelectedText()
        return code[:numChars] if code else ''

        # if code:
        #     for cc in code:
        #         if cc.isalpha():
        #             return cc.upper()[0]
        # return '-'

    def _getValidAtomType(self):
        """Get the valid atom type from the buttons, numChars is included as this may be needed for DNA/RNA
        """
        return self.atomTypeOptions.getSelectedText()

    def _createBackBoneButtons(self):
        self._cleanupPickAndAssignWidget()
        for w in self._sidechainModifiers: w.hide()

        atoms = BACKBONEATOMS

        self.buttons = {}
        # seems to be deleting the same widgets as _clean
        for b in list(self.buttonGroup.buttons()):
            self.buttonGroup.removeButton(b)
            del b

        #     b.setParent(None)
        #     print('>>>deleting', b)
        #     del b

        # rowCount = self._assignWidget.layout().rowCount()
        # colCount = self._assignWidget.layout().columnCount()
        #
        # for r in range(1, rowCount):
        #   for m in range(colCount):
        #     item = self._assignWidget.layout().itemAtPosition(r, m)
        #     if item:
        #       if item.widget():
        #         item.widget().hide()
        #     self._assignWidget.layout().removeItem(item)

        if self.current.nmrResidue:
            rows = 0
            cols = 0
            validAxisCode = self._getValidAxisCode()
            validAtomType = self._getValidAtomType()

            if not validAxisCode:
                return 0, 0

            for ii, atom in enumerate(atoms):
                self.buttons[atom] = []

                if self.selectAxisCode.isChecked():
                    # display by axis codes
                    if not atom.startswith(validAxisCode):
                        continue

                elif validAtomType != 'Other' and not atom.startswith(validAtomType):
                    continue
                elif validAtomType == 'Other':
                    # display by atom types
                    if atom[0] in ['H', 'C', 'N']:
                        continue

                # # skip if startswith these atomTypes
                # if not self.cCheckBox.isChecked() and atom.startswith('C'):
                #     continue
                # if not self.hCheckBox.isChecked() and atom.startswith('H'):
                #     continue
                # if not self.nCheckBox.isChecked() and atom.startswith('N'):
                #     continue
                # if not self.otherCheckBox.isChecked() and not atom.startswith('C') \
                #         and not atom.startswith('H') \
                #         and not atom.startswith('N'):
                #     continue

                innerCols = 0
                for jj, offset in enumerate(['-1', '0', '+1']):
                    btext = self.atomLabel(atom, offset)
                    button = _RButton(self._assignWidget, text=btext, grid=(rows, jj),
                                      callback=None)  #partial(self.assignSelected, offset, atom))
                    button.setMinimumSize(BUTTON_MINX, BUTTON_MINY)
                    self.buttonGroup.addButton(button)
                    button._atomName = atom
                    button._offSet = offset
                    # button.clicked.connect(self._buttonCallback)

                    self.buttons[atom].append(button)

                    innerCols += 1
                rows += 1
                cols = max(cols, innerCols)

            # self._predictAssignments(self.current.peaks)
            return rows, cols

    def _createSideChainButtons(self):
        self._cleanupPickAndAssignWidget()
        for w in self._sidechainModifiers: w.show()

        ii, jj = self._updateChainLayout()
        # self._predictAssignments(self.current.peaks)
        return ii, jj

    def _changeAxisCode(self):
        self._toggleBox()

    def _changeAtomType(self):
        self._toggleBox()

    def _selectionCallback(self):
        self._showSelectionButtons()
        self._toggleBox()

    def _showSelectionButtons(self):
        if self.current.peaks:
            if self.selectAtomType.isChecked():
                self.axisCodeOptions.hide()
                self.axisCodeLabel.hide()
                self.atomTypeLabel.show()
                self.atomTypeOptions.show()
            else:
                self.axisCodeOptions.show()
                self.axisCodeLabel.show()
                self.atomTypeLabel.hide()
                self.atomTypeOptions.hide()

            self.shiftlabel.show()

    def _toggleBox(self):
        if self.selectBackboneButton.isChecked():
            for w in self._sidechainModifiers: w.hide()
        elif self.selectAllButton.isChecked():
            for w in self._sidechainModifiers: w.show()
        self._updateWidget()

    @staticmethod
    def _getAtomsForButtons(atomList, atomName):
        [atomList.remove(atom) for atom in sorted(atomList) if not atom.startswith(atomName)]

    @staticmethod
    def _removeAtomsForButtons(atomList, atomName):
        [atomList.remove(atom) for atom in sorted(atomList) if atom.startswith(atomName)]

    @staticmethod
    def _getAtomButtonList(residueType=None):

        additionalAtoms = list(ADDITIONALBACKBONEATOMS)
        alphaAtoms = list(NEF_ATOM_NAMES_SORTED['alphas'])
        betaAtoms = list(NEF_ATOM_NAMES_SORTED['betas'])
        gammaAtoms = list(NEF_ATOM_NAMES_SORTED['gammas'])
        moreGammaAtoms = list(NEF_ATOM_NAMES_SORTED['moreGammas'])
        deltaAtoms = list(NEF_ATOM_NAMES_SORTED['deltas'])
        moreDeltaAtoms = list(NEF_ATOM_NAMES_SORTED['moreDeltas'])
        epsilonAtoms = list(NEF_ATOM_NAMES_SORTED['epsilons'])
        moreEpsilonAtoms = list(NEF_ATOM_NAMES_SORTED['moreEpsilons'])
        zetaAtoms = list(NEF_ATOM_NAMES_SORTED['zetas'])
        etaAtoms = list(NEF_ATOM_NAMES_SORTED['etas'])
        moreEtaAtoms = list(NEF_ATOM_NAMES_SORTED['moreEtas'])

        atomButtonList = [additionalAtoms,
                          alphaAtoms, betaAtoms, gammaAtoms, moreGammaAtoms, deltaAtoms, moreDeltaAtoms,
                          epsilonAtoms, moreEpsilonAtoms, zetaAtoms, etaAtoms, moreEtaAtoms]

        if not residueType or not isinstance(residueType, str):
            return atomButtonList
        residueType = residueType.upper()
        if residueType in PROTEIN_NEF_ATOM_NAMES:
            residueAtoms = PROTEIN_NEF_ATOM_NAMES[residueType]
            residueAdditional = [atom for atom in additionalAtoms if atom in residueAtoms]
            residueAlphas = [atom for atom in alphaAtoms if atom in residueAtoms]
            residueBetas = [atom for atom in betaAtoms if atom in residueAtoms]
            residueGammas = [atom for atom in gammaAtoms if atom in residueAtoms]
            residueMoreGammas = [atom for atom in moreGammaAtoms if atom in residueAtoms]
            residueDeltas = [atom for atom in deltaAtoms if atom in residueAtoms]
            residueMoreDeltas = [atom for atom in moreDeltaAtoms if atom in residueAtoms]
            residueEpsilons = [atom for atom in epsilonAtoms if atom in residueAtoms]
            residueMoreEpsilons = [atom for atom in moreEpsilonAtoms if atom in residueAtoms]
            residueZetas = [atom for atom in zetaAtoms if atom in residueAtoms]
            residueEtas = [atom for atom in etaAtoms if atom in residueAtoms]
            residueMoreEtas = [atom for atom in moreEtaAtoms if atom in residueAtoms]
            return [residueAdditional, residueAlphas, residueBetas, residueGammas, residueMoreGammas, residueDeltas, residueMoreDeltas, residueEpsilons, residueMoreEpsilons, residueZetas, residueEtas, residueMoreEtas]

    @staticmethod
    def _getDnaRnaButtonList(atomList=None, residueType=None):
        residueAtomButtonList = copy.deepcopy(ALL_DNARNA_ATOMS_SORTED)

        if residueType and atomList:
            residueAtoms = atomList[residueType]

            for atomType in ALL_DNARNA_ATOMS_SORTED.keys():
                atomTypeList = ALL_DNARNA_ATOMS_SORTED[atomType]
                for atom in atomTypeList:
                    if atom not in residueAtoms:
                        residueAtomButtonList[atomType].remove(atom)

        return [residueAtomButtonList[atom] for atom in residueAtomButtonList.keys()]

    def _removeCodes(self, atomButtonList):
        if self.selectAxisCode.isChecked():
            # add atoms for the axisCode selected
            [self._getAtomsForButtons(atomList, self._getValidAxisCode()) for atomList in atomButtonList]

        else:
            # add atoms for atom type selected
            validAtomType = self._getValidAtomType()
            if validAtomType == 'C':
                [self._getAtomsForButtons(atomList, 'C') for atomList in atomButtonList]

            elif validAtomType == 'H':
                [self._getAtomsForButtons(atomList, 'H') for atomList in atomButtonList]

            elif validAtomType == 'N':
                [self._getAtomsForButtons(atomList, 'N') for atomList in atomButtonList]

            elif validAtomType == 'Other':
                for atomList in atomButtonList:
                    [self._removeAtomsForButtons(atomList, 'C') for atomList in atomButtonList]
                    [self._removeAtomsForButtons(atomList, 'H') for atomList in atomButtonList]
                    [self._removeAtomsForButtons(atomList, 'N') for atomList in atomButtonList]

    def _updateChainLayout(self):

        atomButtonList = []
        # needs more work to allow DNA/RNA molecules
        if self.molTypePulldown.currentText() == PROTEIN_MOLECULE:
            # group atoms in useful categories based on usage
            atomButtonList = self._getAtomButtonList()

        elif self.molTypePulldown.currentText() == DNA_MOLECULE:
            # testing DNA/RNA button-list
            atomButtonList = self._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DT')

        elif self.molTypePulldown.currentText() == RNA_MOLECULE:
            # testing DNA/RNA button-list
            atomButtonList = self._getDnaRnaButtonList(RNA_ATOM_NAMES, 'G')

        self._removeCodes(atomButtonList)

        # if self.selectAxisCode.isChecked():
        #     # add atoms for the axisCode selected
        #     [self._getAtomsForButtons(atomList, self._getValidAxisCode()) for atomList in atomButtonList]
        #
        # else:
        #     # add atoms for atom type selected
        #     validAtomType = self._getValidAtomType()
        #     if validAtomType == 'C':
        #         [self._getAtomsForButtons(atomList, 'C') for atomList in atomButtonList]
        #
        #     elif validAtomType == 'H':
        #         [self._getAtomsForButtons(atomList, 'H') for atomList in atomButtonList]
        #
        #     elif validAtomType == 'N':
        #         [self._getAtomsForButtons(atomList, 'N') for atomList in atomButtonList]
        #
        #     elif validAtomType == 'Other':
        #         for atomList in atomButtonList:
        #             [self._removeAtomsForButtons(atomList, 'C') for atomList in atomButtonList]
        #             [self._removeAtomsForButtons(atomList, 'H') for atomList in atomButtonList]
        #             [self._removeAtomsForButtons(atomList, 'N') for atomList in atomButtonList]
        # [atomList.remove(atom) for atom in sorted(atomList) if not atom.startswith('C') \
        #  and not atom.startswith('H') \
        #  and not atom.startswith('N')]

        # # Activate button for Carbons
        # if not self.cCheckBox.isChecked():
        #     [self._getAtomsForButtons(atomList, 'C') for atomList in atomButtonList]
        #
        # if not self.hCheckBox.isChecked():
        #     [self._getAtomsForButtons(atomList, 'H') for atomList in atomButtonList]
        #
        # if not self.nCheckBox.isChecked():
        #     [self._getAtomsForButtons(atomList, 'N') for atomList in atomButtonList]
        #
        # if not self.otherCheckBox.isChecked():
        #     for atomList in atomButtonList:
        #         [atomList.remove(atom) for atom in sorted(atomList) if not atom.startswith('C') \
        #          and not atom.startswith('H') \
        #          and not atom.startswith('N')]

        # rowCount = self._assignWidget.layout().rowCount()
        # colCount = self._assignWidget.layout().columnCount()
        #
        # for r in range(1, rowCount):
        #     for m in range(colCount):
        #         item = self._assignWidget.layout().itemAtPosition(r, m)
        #         if item:
        #             if item.widget():
        #                 item.widget().hide()
        #         self._assignWidget.layout().removeItem(item)

        rows = 0
        cols = 0
        if self.current.nmrResidue:
            # self.currentNmrResidueLabel.setText(self.current.nmrResidue.id)
            self._nmrResiduePulldown.select(self.current.nmrResidue.pid)

            self.buttons = {}
            # seems to be deleting the same widgets as _clean
            for b in list(self.buttonGroup.buttons()):
                self.buttonGroup.removeButton(b)
                del b

            if not self.current.nmrResidue.residueType:
                for ii, atomList in enumerate(atomButtonList):

                    for jj, atom in enumerate(atomList):
                        self.buttons[atom] = []
                        offset = self.offsetSelector.currentText()
                        bText = self.atomLabel(atom, offset)
                        button = _RButton(self._assignWidget, text=bText, grid=(rows, jj), hAlign='t', )
                        # callback=partial(self.assignSelected, offset, atom))
                        button._atomName = atom
                        button._offSet = offset
                        button.setMinimumSize(BUTTON_MINX, BUTTON_MINY)
                        self.buttonGroup.addButton(button)
                        self.buttons[atom].append(button)

                        cols = max(cols, jj + 1)

                    if atomList:
                        rows += 1
            else:
                if self.offsetSelector.currentText() == '-1':
                    nmrResidue = self.current.nmrResidue.previousNmrResidue
                elif self.offsetSelector.currentText() == '+1':
                    nmrResidue = self.current.nmrResidue.nextNmrResidue
                else:
                    nmrResidue = self.current.nmrResidue
                residueType = nmrResidue.residueType.upper() if nmrResidue else None

                atomButtonList2 = self._getAtomButtonList(residueType)
                self._removeCodes(atomButtonList2)

                for ii, atomList in enumerate(atomButtonList2):
                    for jj, atom in enumerate(atomList):
                        self.buttons[atom] = []

                        offset = self.offsetSelector.currentText()
                        bText = self.atomLabel(atom, offset)
                        button = _RButton(self._assignWidget, text=bText, grid=(rows, jj), hAlign='t', )
                        # callback=partial(self.assignSelected, self.offsetSelector.currentText(), atom))
                        # button = Button(self._assignWidget, text=atom, grid=(ii, jj), hAlign='t',
                        #         callback=partial(self.assignSelected, self.offsetSelector.currentText(), atom))
                        button._atomName = atom
                        button._offSet = offset
                        self.buttonGroup.addButton(button)
                        button.setMinimumSize(BUTTON_MINX, BUTTON_MINY)
                        # button.setAutoExclusive(True)

                        self.buttons[atom].append(button)

                        cols = max(cols, jj + 1)

                    if atomList:
                        rows += 1

        return cols, rows

    @staticmethod
    def _showMoreAtomButtons(buttons, moreButton):
        if moreButton.isChecked():
            [button.show() for button in buttons]
        else:
            [button.hide() for button in buttons]

    @staticmethod
    def _removeWidget(widget, removeTopWidget=False):
        """Destroy a widget and all it's contents
        """

        def deleteItems(layout):
            if layout is not None:
                while layout.count():
                    item = layout.takeAt(0)
                    widget = item.widget()
                    if widget is not None:
                        widget.setVisible(False)
                        widget.setParent(None)
                        del widget

        deleteItems(widget.getLayout())
        if removeTopWidget:
            del widget

    def _cleanupPickAndAssignWidget(self):
        self._removeWidget(self._assignWidget)

    @staticmethod
    def atomLabel(atom, offset, showAll=False):
        if showAll:
            return str(f'{atom} [i]' if offset == '0' else f'{atom} [i{offset}]')
        else:
            return str(atom if offset == '0' else f'{atom} [i{offset}]')

    # NOT NEEDED
    # def checkAssignedAtoms(self, nmrResidue, atoms, predictAtoms, checkMode='backbone'):
    #   """
    #   Check if the i-1, i, i+1 nmrAtoms for the current residue exist
    #   :param nmrResidue:
    #   :return foundAtoms - dict containing True for each found nmrAtom:
    #   """
    #   foundAtoms = {}
    #   if checkMode == 'backbone':
    #     # all residues are displayed so use the central mainNmrResidue
    #     nmrResidue = self._getMainNmrResidue(nmrResidue)
    #
    #   # atoms = ['H', 'N', 'CA', 'CB', 'CO', 'HA', 'HB']
    #   for ii, atom in enumerate(atoms):
    #     for jj, offset in enumerate(['-1', '0', '+1']):
    #       bText = self.atomLabel(atom, offset)
    #
    #       if self._checkAssignedAtom(nmrResidue, offset, atom):
    #         foundAtoms[bText] = True
    #         if atom in self.buttons.keys():
    #           for button in self.buttons[atom]:
    #             if button.getText() == bText:
    #
    #               # colour the button if the atom exists
    #               if bText in predictAtoms:
    #                 score = predictAtoms[bText]
    #
    #                 if score >= 85:
    #                   button.setStyleSheet('background-color: mediumseagreen')
    #                 elif 50 < score < 85:
    #                   button.setStyleSheet('background-color: lightsalmon')
    #                 if score < 50:
    #                   button.setStyleSheet('background-color: mediumvioletred')
    #
    #               # else: # Users don't need to know/ or be notified here that the nmrAtom exist in the selected nmrResidue.
    #               #   button.setStyleSheet('background-color: cornflowerblue')
    #
    #
    #   return foundAtoms

    def _getNmrResidue(self, nmrChain, sequenceCode: typing.Union[int, str] = None,
                       residueType: str = None) -> typing.Optional[NmrResidue]:
        partialId = f'{nmrChain.id}.{str(sequenceCode).translate(Pid.remapSeparators)}.'

        if ll := self.project.getObjectsByPartialId(className='NmrResidue', idStartsWith=partialId):
            return ll[0]
        else:
            return nmrChain.getNmrResidue(sequenceCode)

    def _fetchNmrResidue(self, nmrChain, sequenceCode: typing.Union[int, str] = None,
                         residueType: str = None) -> typing.Optional[NmrResidue]:
        partialId = f'{nmrChain.id}.{str(sequenceCode).translate(Pid.remapSeparators)}.'

        if ll := self.project.getObjectsByPartialId(className='NmrResidue', idStartsWith=partialId):
            return ll[0]
        else:
            return nmrChain.fetchNmrResidue(sequenceCode)

    def _getCorrectResidue(self, nmrResidue, offset: str, atomType: str):
        if offset == '-1' and '-1' not in nmrResidue.sequenceCode:
            if not (r := nmrResidue.previousNmrResidue):
                r = self._fetchNmrResidue(nmrResidue.nmrChain, sequenceCode=f'{nmrResidue.sequenceCode}-1')

        elif offset == '+1' and '+1' not in nmrResidue.sequenceCode:
            if not (r := nmrResidue.nextNmrResidue):
                r = self._fetchNmrResidue(nmrResidue.nmrChain, sequenceCode=f'{nmrResidue.sequenceCode}+1')

        else:
            r = nmrResidue

        return r

    # NOT NEEDED
    # def _checkAssignedAtom(self, nmrResidue, offset:int, atomType:str):
    #   '''
    #   This checks only if an NMRAtom exists not if is assigned to a peak!
    #   :param nmrResidue:
    #   :param offset:
    #   :param atomType:
    #   :return:
    #   '''
    #   r = self._getCorrectResidue(nmrResidue=nmrResidue, offset=offset, atomType=atomType)
    #   if r:
    #     atom = r.getNmrAtom(atomType.translate(Pid.remapSeparators))
    #     if atom and not atom.isDeleted:
    #       return r
    #     else:
    #       return None
    #   else:
    #     return None

    # def _getMainNmrResidue(self, nmrResidue):
    #   if nmrResidue:
    #     if nmrResidue.relativeOffset and nmrResidue.relativeOffset != 0:
    #       return nmrResidue.mainNmrResidue
    #
    #   return nmrResidue

    def _assignDimension(self):
        """
        update the peak assignment and create a event to update the module
        """
        pass

    def assignNmrAtomsToPeaks(self, nmrAtom, peaks):
        """Assign the nmrAtom to the dimension
        """
        if not peaks: return
        if not nmrAtom: return

        index = self._getValidAxisCodeIndex()
        for peak in peaks:
            for ii, axisCode in enumerate(peak.axisCodes):
                if self.peakIndex[peak][ii] == index:
                    peak.assignDimension(axisCode, nmrAtom)
                    isotopeCode = peak.peakList.spectrum.getByAxisCodes('isotopeCodes', [axisCode], exactMatch=True)[-1]
                    if nmrAtom.isotopeCode in [UnknownIsotopeCode, None]:
                        nmrAtom._setIsotopeCode(isotopeCode)

                    # update the label after the assignment has completed, otherwise notifier is too early for shift update?
                    self._updatePeakShiftLabel(peak)

    def deassignAtomFromSelectedPeaks(self, peaks, nmrAtom):
        """Deassign the nmrAtom from the dimension
        """
        if not peaks: return
        if not nmrAtom: return

        # newAssignedAtoms = ()
        index = self._getValidAxisCodeIndex()
        for peak in peaks:

            if self.selectAxisCode.isChecked():
                # de-assign by axis code dimension
                for ii, axisCode in enumerate(peak.axisCodes):
                    if self.peakIndex[peak][ii] == index:
                        peak.assignDimension(axisCode, None)

            else:
                # de-assign by atom types
                peakDimNmrAtoms = list(peak.dimensionNmrAtoms)
                for dim, dimNmrAtoms in enumerate(peakDimNmrAtoms):
                    if nmrAtom in dimNmrAtoms:
                        dimNmrAtoms = list(dimNmrAtoms)
                        dimNmrAtoms.remove(nmrAtom)

                        peakDimNmrAtoms[dim] = dimNmrAtoms

                peak.dimensionNmrAtoms = peakDimNmrAtoms

                # for subTuple in peak.assignedNmrAtoms:
                #     a = tuple(None if na is nmrAtom else na for na in subTuple)
                #     newAssignedAtoms += (a,)

                # try:
                #     peak.assignedNmrAtoms = newAssignedAtoms
                # except Exception as es:
                #     pass

    def _returnButtonsToNormal(self):
        """
        Returns all buttons in Atom Selector to original colours and style.
        """
        self._assignWidget.setStyleSheet(DEFAULT_BUTTON)
        for btnList in self.buttons.values():
            for btn in btnList:
                btn.setBackgroundColours(*DEFAULT_COLOURS)

    def _currentNmrResiduesCallback(self, data):
        """Callback for the nmrResidues notifier
        """
        # set to the first current nmrResidue
        if (curRess := data[Notifier.VALUE]):
            self._thisNmrResidue = curRess[0]
        else:
            self._thisNmrResidue = None

        self._updateWidget()

    def _currentPeaksCallback(self, data):
        """Callback for the peaks notifier
        """
        self._thisPeaks = data[Notifier.VALUE]

        self._nmrResiduePulldown.update()
        self._updateWidget()

    # def _predictAssignments(self, peaks: typing.List[Peak]):
    #     """
    #     Predicts atom type for selected peaks and highlights the relevant buttons with confidence of
    #     that assignment prediction, green is very confident, orange is less confident.
    #     """
    #     self._assignWidgetHide()
    #
    #     if self.current.nmrResidue and self.current.peaks:
    #         self._predictHighlight(peaks)
    #
    #     self._assignWidgetShow()

    def _updateShiftFromPeaks(self, data):
        """Update the chemical-shift if notified of a peak-change
        """
        if (peak := data[Notifier.OBJECT]) and peak == self.current.peak:
            # update the shift
            self._updatePeakShiftLabel(peak)

    def _predictHighlight(self, peaks: typing.List[Peak]):
        """Highlight the predictions in the atomName table
        """
        self._returnButtonsToNormal()
        # if self.current.nmrResidue is None or len(peaks) == 0:
        #     self._assignWidgetHide()
        #     return
        #
        # # make sure that the widget is visible
        # self._assignWidgetShow()

        # make sure that you have buttons!
        if not self.buttonGroup.buttons():
            return

        # check if peaks coincide
        for dim in range(peaks[0].spectrum.dimensionCount):
            if not peaksAreOnLine(peaks, dim):
                logger.debug('dimension %s: peaksAreonLine=False' % dim)
                return

        types = {peak.spectrum.experimentType for peak in peaks}
        anyInterOnlyExperiments = any(isInterOnlyExpt(x) for x in types)

        logger.debug('peaks=%s' % (peaks,))
        logger.debug('types=%s, anyInterOnlyExperiments=%s' % (types, anyInterOnlyExperiments))

        peak = peaks[0]
        peakListViews = [peakListView for peakListView in self.project.peakListViews if
                         peakListView.peakList == peak.peakList]

        if peakListViews:
            spectrumIndices = peakListViews[0].spectrumView.dimensionIndices

            # for the 1D case, this is (0, None)
            if len(spectrumIndices) < 2 or spectrumIndices[1] is None:
                return

            isotopeCode = peak.spectrum.isotopeCodes[spectrumIndices[1]]

            # backbone
            if self.selectBackboneButton.isChecked():
                self._highlightBackboneAtomNames(anyInterOnlyExperiments, isotopeCode, peak, spectrumIndices)

            # side-chain is checked
            elif self.selectAllButton.isChecked():

                if self.current.nmrResidue.residueType == '':
                    # In this case, we loop over all CCP_CODES (i.e. residue types)
                    predictedAtomTypes = []
                    for residueType in CCP_CODES:
                        for tp, score in getNmrAtomPrediction(residueType, peak.position[spectrumIndices[1]],
                                                              isotopeCode):
                            if len(tp) > 0 and score > 50:
                                predictedAtomTypes.append((tp, score))

                else:
                    # bigger offsets?
                    if self.offsetSelector.currentText() == '-1':
                        nmrResidue = self.current.nmrResidue.previousNmrResidue
                    elif self.offsetSelector.currentText() == '+1':
                        nmrResidue = self.current.nmrResidue.nextNmrResidue
                    else:
                        nmrResidue = self.current.nmrResidue

                    # don't need to predict if there is no nmrResidue
                    if not nmrResidue:
                        return

                    predictedAtomTypes = getNmrAtomPrediction(nmrResidue.residueType.title(),
                                                              peak.position[spectrumIndices[1]], isotopeCode)

                # print('>predictAtomTypes>', predictedAtomTypes)
                # find the maximum of each atomType
                predictedDict = {}
                for tp, score in predictedAtomTypes:
                    # if tp[1] not in predictedDict:
                    #     predictedDict[tp[1]] = [tp[0], score]
                    # else:
                    #     if score > predictedDict[tp[1]][1]:
                    #         predictedDict[tp[1]] = [tp[0], score]
                    pd = predictedDict.setdefault(tp[1], [tp[0], score])
                    pd[1] = max(pd[1], score)

                # print ('>>>predictedDict', predictedDict)

                foundPredictList = {}
                currentOffset = self.offsetSelector.currentText()
                for atomDictType in predictedDict.keys():
                    bText = self.atomLabel(atomDictType, currentOffset)
                    for atomType, buttons in self.buttons.items():  # get the correct button list
                        if atomDictType == atomType:
                            for button in buttons:
                                if bText == button.getText():

                                    if atomType in ['CA', 'CB']:

                                        # match the colouring above
                                        if (currentOffset == '-1' and anyInterOnlyExperiments) or \
                                                (currentOffset == '0' and not anyInterOnlyExperiments):

                                            # print('>type[1], atomType, button>', atomDictType, bText)
                                            foundPredictList[self.atomLabel(atomDictType, currentOffset)] = score

                                            if score >= 85:
                                                button.setBackgroundColours(*GREEN_COLOURS)
                                            elif 50 < score < 85:
                                                button.setBackgroundColours(*ORANGE_COLOURS)
                                            if score < 50:
                                                button.setBackgroundColours(*RED_COLOURS)

                                    else:

                                        # print('>type[1], atomType, button>', atomDictType, bText)
                                        foundPredictList[self.atomLabel(atomDictType, currentOffset)] = score

                                        if score >= 85:
                                            button.setBackgroundColours(*GREEN_COLOURS)
                                        elif 50 < score < 85:
                                            button.setBackgroundColours(*ORANGE_COLOURS)
                                        if score < 50:
                                            button.setBackgroundColours(*RED_COLOURS)

                # new routine to colour any existing atoms
                # atomButtonList = self._getAtomButtonList()
                # atomButtonList = [x for i in atomButtonList for x in i]
                # foundAtoms = self.checkAssignedAtoms(self.current.nmrResidue, atomButtonList,
                #                                      foundPredictList, 'sideChain')

    def _highlightBackboneAtomNames(self, anyInterOnlyExperiments, isotopeCode, peak, spectrumIndices):
        # get the predictions aligned to the y-axis of the spectrumDisplay
        predictedAtomTypes = [
            getNmrAtomPrediction(ccpCode, peak.position[spectrumIndices[1]], isotopeCode, strict=True)
            for ccpCode in CCP_CODES]
        refinedPreds = [(tp[0][0][1], tp[0][1]) for tp in predictedAtomTypes if len(tp) > 0]
        atomPredictions = {atomPred for atomPred, score in refinedPreds if score >= 85}

        # list containing those atoms that exist - used for colouring in 'checkAssignedAtoms'
        foundPredictList = {}
        for atomPred in atomPredictions:
            for atm in ['CA', 'CB']:
                self._colourAtomName(anyInterOnlyExperiments, atm, atomPred, foundPredictList)

        # new routine to colour any existing atoms
        # foundAtoms = self.checkAssignedAtoms(self.current.nmrResidue, ATOM_TYPES,
        #                                      foundPredictList, 'backbone')

    def _colourAtomName(self, anyInterOnlyExperiments, atm, atomPred, foundPredictList):
        """Colour an atomName in the table
        """
        if atomPred == atm and self.buttons[atm]:
            if anyInterOnlyExperiments:
                self.buttons[atm][0].setBackgroundColours(*GREEN_COLOURS)
                foundPredictList[self.atomLabel(atm, '-1')] = 100
            else:
                self.buttons[atm][0].setBackgroundColours(*GREEN_COLOURS)
                self.buttons[atm][1].setBackgroundColours(*GREEN_COLOURS)
                foundPredictList[self.atomLabel(atm, '-1')] = 100
                foundPredictList[self.atomLabel(atm, '0')] = 100

    def _changeMoleculeType(self, data):
        """
        change the  available atomList depending on the moleculeType
        :param data: str from pullDown
        """
        pass

    @staticmethod
    def getResidueTypes(moleculeType: str = 'protein'):
        """
        return a list of residue types assiciated with the moleculeType
        :param moleculeType - str ['protein', 'DNA', 'RNA', 'carbohydrate', 'other']
        :return list of str:
        """
        if moleculeType not in MOLECULE_TYPES:
            return None
        if moleculeType == 'protein':
            return list(PROTEIN_NEF_ATOM_NAMES.keys())


if __name__ == '__main__':
    from ccpn.ui.gui.widgets.Application import TestApplication
    from ccpn.ui.gui.widgets.TextEditor import TextEditor


    app = TestApplication()

    popup = NmrAtomAssignerModule()

    textBox = TextEditor(popup.mainWidget, grid=(3, 0), gridSpan=(1, 1))

    all_atoms = {atom: [] for atom in ['O', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'P']}

    atomList = [DNA_ATOMS, RNA_ATOMS]
    startAtoms = ['DA', 'DC', 'DG', 'DT', 'A', 'G', 'C', 'U']

    for atomText in atomList:
        atoms = {}
        atomText = atomText.split('\n')
        for line in atomText:
            if ll := line.split():
                if ll[0] in startAtoms:
                    if ll[0] in atoms:
                        atoms[ll[0]].append(ll[1])
                    else:
                        atoms[ll[0]] = [ll[1]]

                    if ll[1] != 'P' and ll[1] not in all_atoms[ll[1][1]]:
                        all_atoms[ll[1][1]].append(ll[1])
                        all_atoms[ll[1][1]].sort()

        textBox.append(str(atoms))
    all_atoms['P'] = ['P']
    textBox.append(str(all_atoms))

    popup._nmrResidue.setText('NmrResidue here')

    print(popup._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DT'))
    print(popup._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DC'))
    print(popup._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DA'))
    print(popup._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DG'))

    print(popup._getDnaRnaButtonList(RNA_ATOM_NAMES, 'G'))
    print(popup._getDnaRnaButtonList(RNA_ATOM_NAMES, 'U'))
    print(popup._getDnaRnaButtonList(RNA_ATOM_NAMES, 'A'))
    print(popup._getDnaRnaButtonList(RNA_ATOM_NAMES, 'C'))

    popup.show()
    popup.raise_()
    app.start()
