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
__modifiedBy__ = "$modifiedBy: Geerten Vuister $"
__dateModified__ = "$dateModified: 2021-12-09 16:21:40 +0000 (Thu, December 09, 2021) $"
__version__ = "$Revision: 3.0.4 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

# import os
import typing
import copy
from functools import partial
from collections import OrderedDict
from PyQt5 import QtCore, QtGui, QtWidgets
from contextlib import contextmanager
from ccpn.core.Peak import Peak
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.NmrAtom import NmrAtom
from ccpn.core.lib import Pid
from ccpn.core.lib.AssignmentLib import isInterOnlyExpt, getNmrAtomPrediction, CCP_CODES
from ccpn.core.lib.AssignmentLib import peaksAreOnLine
from ccpn.ui.gui.modules.CcpnModule import CcpnModule
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.CheckBox import CheckBox
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
from ccpnmodel.ccpncore.lib.assignment.ChemicalShift import PROTEIN_ATOM_NAMES, ALL_ATOMS_SORTED
from ccpn.core.lib.AssignmentLib import PROTEIN_NEF_ATOM_NAMES, NEF_ATOM_NAMES_SORTED
from ccpn.util.Logging import getLogger
from ccpn.core.lib.Notifiers import Notifier
from ccpn.core.lib.AssignmentLib import _assignNmrAtomsToPeaks
from ccpn.ui.gui.widgets.ScrollArea import ScrollArea
from ccpn.ui.gui.guiSettings import BORDERNOFOCUS_COLOUR
from ccpn.core.lib import CcpnSorting


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

BUTTON_MINX = 70
BUTTON_MINY = 24

DEFAULT_BUTTON = """QRadioButton { background-color: %s }
                   QRadioButton::hover { background-color: %s}""" % ('lightgrey', 'white')
GREEN_BUTTON = """QRadioButton { background-color: %s }
                   QRadioButton::hover { background-color: %s}""" % ('mediumseagreen', 'palegreen')
ORANGE_BUTTON = """QRadioButton { background-color: %s }
                   QRadioButton::hover { background-color: %s}""" % ('orange', 'gold')
RED_BUTTON = """QRadioButton { background-color: %s }
                   QRadioButton::hover { background-color: %s}""" % ('tomato', 'lightpink')


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

        # modifiers for sidechain
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
        self._nmrChain = NmrChainPulldown(_f, mainWindow=self.mainWindow,
                                          labelText='NmrChain:', showSelectName=True,
                                          setCurrent=False,
                                          callback=self._nmrChainPullDownCallback,
                                          grid=(0, 0), hPolicy='minimal', minimumWidths=None,
                                          sizeAdjustPolicy=QtWidgets.QComboBox.AdjustToContents)
        self._nmrResidue = NmrResiduePulldown(_f, mainWindow=self.mainWindow,
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
                                            callback=self._changeAxisCode, grid=(labRow, 1))
        labRow += 1

        # modifier for atomType
        self.atomTypeLabel = Label(self._labelFrame, 'Assign by atomType', grid=(labRow, 0))
        self.atomTypeOptions = RadioButtons(self._labelFrame, selectedInd=1, texts=['H', 'C', 'N', 'Other'],
                                            callback=self._changeAtomType, grid=(labRow, 1))
        labRow += 1

        Spacer(self._labelFrame, 2, 2, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed,
               grid=(labRow, 2), gridSpan=(1, 1))
        resRow += 1

        self._assignWidget = Frame(self._residueFrame, setLayout=True, showBorder=False, grid=(resRow, 0), spacing=(5, 5))
        self._scrollAreaWidget.setWidget(self._residueFrame)
        resRow += 1

        # add spacer to stop columns changing width
        Spacer(self._residueFrame, 2, 2, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding,
               grid=(resRow, 4), gridSpan=(1, 1))
        resRow += 1

        self.buttonGroup = QtWidgets.QButtonGroup()
        self.buttonGroup.buttonClicked.connect(self._nmrAtomButtonsCallback)
        self.buttonGroup.setExclusive(False)

        self.buttons = {}
        self._registerNotifiers()
        self._updateWidget()

    def _registerNotifiers(self):
        """Register notifiers for the module
        """
        # self.setNotifier(self.project, [Notifier.CHANGE, Notifier.CREATE, Notifier.DELETE],
        #                  NmrAtom.className, callback=self._nmrResidueCallBack, onceOnly=True)
        # self.setNotifier(self.project, [Notifier.CHANGE],
        #                  Peak.className, callback=self._nmrResidueCallBack, onceOnly=True)
        self.setNotifier(self.current, [Notifier.CURRENT],
                         Peak._pluralLinkName, callback=self._currentPeaksCallback, onceOnly=True)
        self.setNotifier(self.current, [Notifier.CURRENT],
                         NmrResidue._pluralLinkName, callback=self._currentNmrResiduesCallback, onceOnly=True)
        self.setGuiNotifier(self._residueFrame, [GuiNotifier.DROPEVENT],
                            [DropBase.PIDS], callback=self._handleNmrResidue)
        self.setNotifier(self.project, [Notifier.RENAME],
                         'NmrResidue', self._updateNmrResidue, onceOnly=True)

    def _unRegisterNotifiers(self):
        """clean up the notifiers
        """
        # _closeModule() will do most of them
        self._nmrResidue.unRegister()
        self._nmrChain.unRegister()

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
        # the pulldown; however _nmrResidue is not yet defined then
        if not hasattr(self, '_nmrResidue'): return pids

        nmrChain = self._nmrChain.getSelectedObject()
        #print('>>> filtering pids on:', nmrChain)

        # For selected peaks: get the pids of nmrResidues of assigned nmrAtoms
        newPids = []
        for peak in self.current.peaks:
            if peak:
                for assignment in peak.assignments:
                    for nmrAtom in assignment:
                        if nmrAtom:
                            newPids.append(nmrAtom.nmrResidue.pid)
        newPids = list(set(newPids))
        newPids.sort()

        def _isOk(pid):
            if nmrChain is None:
                # No filtering
                return True
            nmrResidue = self._nmrResidue.value2object(pid)
            if nmrResidue is not None and nmrResidue.nmrChain == nmrChain:
                # nmrResidue is part of the filtered nmrChain
                return True
            return False

        newPids = newPids + [pid for pid in pids if _isOk(pid)]
        return newPids

    def _newNmrResidueCallback(self):
        """Callback to create a new nmrResidue and add to the current chain
        """
        nmrChain = self._nmrChain.getSelectedObject()
        if nmrChain:
            nmrResidue = self._fetchNmrResidue(nmrChain)
            self.current.nmrResidue = nmrResidue

    def _nmrResidueEditCallback(self, data):
        """Callback to edit the current nmrResidue
        """
        # call popup on current nmrResidue
        from ccpn.ui.gui.popups.NmrResiduePopup import NmrResidueEditPopup

        popup = NmrResidueEditPopup(parent=self.mainWindow, mainWindow=self.mainWindow,
                                    obj=self.current.nmrResidue)
        popup.exec_()

    def _nmrChainPullDownCallback(self, value):
        "Callback for the NmrChain selection"
        self._nmrResidue.update()

    def _handleNmrResidue(self, dataDict):
        """drop event handler to accept NmrResidue pids
        """
        pids = dataDict.get(DropBase.PIDS)
        if pids:
            objs = [self.project.getByPid(pid) for pid in pids]
            nmrResidues = [obj for obj in objs if (not obj is None) and isinstance(obj, NmrResidue)]
            if nmrResidues:
                self.current.nmrResidues = nmrResidues

    def _togglePressedButton(self, pressedButton=None):
        '''Ensures only a button at the time is checked, yet allows to uncheck a radio button. If pressedButton is None: unchecks all'''
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
        nmrAtom = nmrResidue.getNmrAtom(atomName.translate(Pid.remapSeparators))
        if nmrAtom:
            self.deassignAtomFromSelectedPeaks(self.current.peaks, nmrAtom)

    def _assignSelected(self, atomName, offSet):
        nmrResidue = self._getCorrectResidue(self.current.nmrResidue, offSet, atomName)
        if not nmrResidue:
            getLogger().warning('Error creating new nmrResidue')
            raise ValueError('Error creating new nmrResidue')

        nmrAtom = nmrResidue.fetchNmrAtom(name=atomName)
        if nmrAtom:
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
    #         # self._nmrResidue.select(MSG)

    def _assignWidgetShow(self):
        self._assignWidget.show()
        self._showSelectionButtons()

    def _assignWidgetHide(self):
        self._assignWidget.hide()
        self.atomTypeLabel.hide()
        self.atomTypeOptions.hide()
        self.axisCodeLabel.hide()
        self.axisCodeOptions.hide()

    def _offsetPullDownCallback(self, tmp=None):
        "Callback if offset pullDown changes"
        if self.current.nmrResidue:
            self._updateWidget()
            # self._predictAssignments(self.current.peaks)

    def _setPeaksLabel(self):
        " update the peaks label from current.peaks"
        if self.current.peaks and None not in self.current.peaks:
            splitter = ', '
            pText = _truncateText(splitter.join([p.id for p in self.current.peaks]), splitter=splitter)
            self.currentPeaksLabel.setToolTip(splitter.join([p.id for p in self.current.peaks]))
            self.currentPeaksLabel.setText(pText)
        else:
            self.currentPeaksLabel.setText(MSG)

    def _setPeakAxisCodes(self, peaks):

        from ccpn.core.lib.AxisCodeLib import getAxisCodeMatch

        if peaks:

            maxLen = 0
            refAxisCodes = None
            for peak in peaks:
                if len(peak.axisCodes) > maxLen:
                    maxLen = len(peak.axisCodes)
                    refAxisCodes = list(peak.axisCodes)

            if not maxLen:
                return

            axisLabels = [set() for ii in range(maxLen)]

            mappings = {}
            for peak in peaks:
                matchAxisCodes = peak.axisCodes

                mapping = getAxisCodeMatch(refAxisCodes, matchAxisCodes)
                for k, v in mapping.items():
                    if v not in mappings:
                        mappings[v] = set([k])
                    else:
                        mappings[v].add(k)

                mapping = getAxisCodeMatch(matchAxisCodes, refAxisCodes)
                for k, v in mapping.items():
                    if v not in mappings:
                        mappings[v] = set([k])
                    else:
                        mappings[v].add(k)

                # example of mappings dict - includes mapping from both sides
                # ('Hn', 'C', 'Nh')
                # {'Hn': {'Hn'}, 'Nh': {'Nh'}, 'C': {'C'}}
                # {'Hn': {'H', 'Hn'}, 'Nh': {'Nh'}, 'C': {'C'}}
                # {'CA': {'C'}, 'Hn': {'H', 'Hn'}, 'Nh': {'Nh'}, 'C': {'CA', 'C'}}
                # {'CA': {'C'}, 'Hn': {'H', 'Hn'}, 'Nh': {'Nh'}, 'C': {'CA', 'C'}}

                self.peakIndex = {}
                # go through the peaks
                for peak in peaks:
                    self.peakIndex[peak] = [0 for ii in range(len(peak.axisCodes))]

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

    def _setPeakAtomCodes(self):
        atomCodes = ['H', 'C', 'N', 'Other']
        self.axisCodeOptions.setButtons(texts=list(atomCodes), tipTexts=list(atomCodes), silent=True)

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
        nmrResidue = data[Notifier.OBJECT]
        if nmrResidue and nmrResidue == self.current.nmrResidue:
            self._updateWidget()

    def _updateWidget(self, dataDict=None):  # also used as notifier callback function
        """Update the widget to reflect the proper state"""
        # try:
        ii = jj = 0
        with self._moduleBlocking():
            self._setPeaksLabel()
            if self.current.nmrResidue is not None:
                # self.currentNmrResidueLabel.setText(self.current.nmrResidue.id)
                self._nmrResidue.select(self.current.nmrResidue.pid)
                self._assignWidgetHide()

                # dimensionalities = set([len(peak.position) for peak in self.current.peaks])
                # if len(dimensionalities) > 1:
                #     self.currentPeaksLabel.setText(MSG)
                #     getLogger().warning('Not all peaks have the same number of dimensions.')
                #     return False

                if self.current.peaks and None not in self.current.peaks:
                    self._setPeakAxisCodes(self.current.peaks)

                if self.selectBackboneButton.isChecked():
                    for w in self._sidechainModifiers:
                        w.hide()
                    ii, jj = self._createBackBoneButtons()
                elif self.selectAllButton.isChecked():
                    for w in self._sidechainModifiers:
                        w.show()
                    ii, jj = self._createSideChainButtons()
                self._setCheckedButtonOfAssignedAtoms(self.current.nmrResidue)

                # add a spacer to the radiobutton box
                Spacer(self._assignWidget, 3, 3,
                       QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding,
                       grid=(30, 30), gridSpan=(1, 1))

                if self.current.peaks and None not in self.current.peaks:
                    # temporarily restrict highlighting to backbone only, as the side-chain version needs improving
                    if self.selectBackboneButton.isChecked():
                        self._predictHighlight(self.current.peaks)
                    self._assignWidgetShow()
            else:
                self._assignWidgetHide()

    def _removeOffsetFromButtonText(self, text: str):
        p = text.split(' ')
        if len(p) > 0:
            return p[0]
        else:
            return text

    def _setCheckedButtonOfAssignedAtoms(self, nmrResidue, offSet='0'):
        """setChecked the radioButton Of Assigned Nmr Atoms.
        This makes sure that if a peak is selected and and assigned to an nmrAtom, the relative button is checked """

        if not self.current.peak: return
        if not nmrResidue: return

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
                        if assignedNmrAtom:
                            if offSet == '0':
                                if assignedNmrAtom.name == button.getText():
                                    counts.add(button)
                            elif button._offSet == offSet:
                                if assignedNmrAtom.name == button._atomName:
                                    counts.add(button)

                else:  #Try to search in + and - 1 offset
                    for offset in ['-1', '+1']:
                        r = self._getNmrResidue(nmrResidue.nmrChain,
                                                sequenceCode=nmrResidue.mainNmrResidue.sequenceCode + offset)
                        if r:
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
        return self.axisCodeOptions.getIndex()

    def _getValidAxisCode(self, numChars=1):
        """Get the valid axis code from the buttons, numChars is included as this may be needed for DNA/RNA
        """
        code = self.axisCodeOptions.getSelectedText()
        return code[0:numChars]

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

        # wb104 27 Jun 2017: changed _cleanupPickAndAssignWidget so removes widgets in reverse order
        # not sure if there was anything else leading to this cludge (and one below) though
        # cludge: don't know why I have to do this for the button to appear:
        ###_Label = Label(self, text='')
        ###self._assignWidget.layout().addWidget(_Label, 0, 0)
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

                else:
                    # display by atom types
                    if validAtomType != 'Other' and not atom.startswith(validAtomType):
                        continue
                    elif validAtomType == 'Other':
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
                    button = RadioButton(self._assignWidget, text=btext, grid=(rows, jj),
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

        # see comment about cludge above
        # cludge: don't know why I have to do this for the button to appear: TODO: fix this
        ###_label= Label(self._assignWidget, '',  hAlign='l')
        ###self._assignWidget.layout().addWidget(_label, 0, 0, QtCore.Qt.AlignRight)

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

    def _toggleBox(self):
        if self.selectBackboneButton.isChecked():
            for w in self._sidechainModifiers: w.hide()
        elif self.selectAllButton.isChecked():
            for w in self._sidechainModifiers: w.show()
        self._updateWidget()

    def _getAtomsForButtons(self, atomList, atomName):
        [atomList.remove(atom) for atom in sorted(atomList) if not atom.startswith(atomName)]

    def _removeAtomsForButtons(self, atomList, atomName):
        [atomList.remove(atom) for atom in sorted(atomList) if atom.startswith(atomName)]

    def _getAtomButtonList(self, residueType=None):

        additionalAtoms = list(ADDITIONALBACKBONEATOMS)
        alphaAtoms = [x for x in NEF_ATOM_NAMES_SORTED['alphas']]
        betaAtoms = [x for x in NEF_ATOM_NAMES_SORTED['betas']]
        gammaAtoms = [x for x in NEF_ATOM_NAMES_SORTED['gammas']]
        moreGammaAtoms = [x for x in NEF_ATOM_NAMES_SORTED['moreGammas']]
        deltaAtoms = [x for x in NEF_ATOM_NAMES_SORTED['deltas']]
        moreDeltaAtoms = [x for x in NEF_ATOM_NAMES_SORTED['moreDeltas']]
        epsilonAtoms = [x for x in NEF_ATOM_NAMES_SORTED['epsilons']]
        moreEpsilonAtoms = [x for x in NEF_ATOM_NAMES_SORTED['moreEpsilons']]
        zetaAtoms = [x for x in NEF_ATOM_NAMES_SORTED['zetas']]
        etaAtoms = [x for x in NEF_ATOM_NAMES_SORTED['etas']]
        moreEtaAtoms = [x for x in NEF_ATOM_NAMES_SORTED['moreEtas']]

        atomButtonList = [additionalAtoms,
                          alphaAtoms, betaAtoms, gammaAtoms, moreGammaAtoms, deltaAtoms, moreDeltaAtoms,
                          epsilonAtoms, moreEpsilonAtoms, zetaAtoms, etaAtoms, moreEtaAtoms]

        if residueType and isinstance(residueType, str):
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
                residueAtomButtonList = [residueAdditional,
                                         residueAlphas, residueBetas, residueGammas, residueMoreGammas,
                                         residueDeltas, residueMoreDeltas, residueEpsilons,
                                         residueMoreEpsilons, residueZetas, residueEtas, residueMoreEtas]
                return residueAtomButtonList

        else:
            return atomButtonList

    def _getDnaRnaButtonList(self, atomList=None, residueType=None):
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

        # needs more work to allow DNA/RNA molecules
        if self.molTypePulldown.currentText() == PROTEIN_MOLECULE:
            # group atoms in useful categories based on usage
            atomButtonList = self._getAtomButtonList()

        elif self.molTypePulldown.currentText() == DNA_MOLECULE:
            # testing DNA/RNA buttonlist
            atomButtonList = self._getDnaRnaButtonList(DNA_ATOM_NAMES, 'DT')

        elif self.molTypePulldown.currentText() == RNA_MOLECULE:
            # testing DNA/RNA buttonlist
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
            self._nmrResidue.select(self.current.nmrResidue.pid)

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
                        button = RadioButton(self._assignWidget, text=bText, grid=(rows, jj), hAlign='t', )
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
                        button = RadioButton(self._assignWidget, text=bText, grid=(rows, jj), hAlign='t', )
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

    def _showMoreAtomButtons(self, buttons, moreButton):
        if moreButton.isChecked():
            [button.show() for button in buttons]
        else:
            [button.hide() for button in buttons]

    def _removeWidget(self, widget, removeTopWidget=False):
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

    def atomLabel(self, atom, offset, showAll=False):
        if showAll:
            return str(atom + ' [i]' if offset == '0' else atom + ' [i' + offset + ']')
        else:
            return str(atom if offset == '0' else atom + ' [i' + offset + ']')

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
        partialId = '%s.%s.' % (nmrChain.id, str(sequenceCode).translate(Pid.remapSeparators))
        ll = self.project.getObjectsByPartialId(className='NmrResidue', idStartsWith=partialId)
        if ll:
            return ll[0]
        else:
            return nmrChain.getNmrResidue(sequenceCode)

    def _fetchNmrResidue(self, nmrChain, sequenceCode: typing.Union[int, str] = None,
                         residueType: str = None) -> typing.Optional[NmrResidue]:
        partialId = '%s.%s.' % (nmrChain.id, str(sequenceCode).translate(Pid.remapSeparators))
        ll = self.project.getObjectsByPartialId(className='NmrResidue', idStartsWith=partialId)
        if ll:
            return ll[0]
        else:
            return nmrChain.fetchNmrResidue(sequenceCode)

    def _getCorrectResidue(self, nmrResidue, offset: str, atomType: str):
        name = atomType
        r = None
        if offset == '-1' and '-1' not in nmrResidue.sequenceCode:
            r = nmrResidue.previousNmrResidue
            if not r:
                r = self._fetchNmrResidue(nmrResidue.nmrChain, sequenceCode=nmrResidue.sequenceCode + '-1')
        elif offset == '+1' and '+1' not in nmrResidue.sequenceCode:
            r = nmrResidue.nextNmrResidue
            if not r:
                r = self._fetchNmrResidue(nmrResidue.nmrChain, sequenceCode=nmrResidue.sequenceCode + '+1')
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

    def deassignAtomFromSelectedPeaks(self, peaks, nmrAtom):
        """Deassign the nmrAtom from the dimension
        """
        if not peaks: return
        if not nmrAtom: return

        newAssignedAtoms = ()
        index = self._getValidAxisCodeIndex()
        for peak in peaks:

            if self.selectAxisCode.isChecked():
                # deassign by axis code dimension
                for ii, axisCode in enumerate(peak.axisCodes):
                    if self.peakIndex[peak][ii] == index:
                        peak.assignDimension(axisCode, None)

            else:
                # deassign by atom types
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

    def _currentNmrResiduesCallback(self, data):
        "Callback for the nmrResidues notifier"
        self._updateWidget()

    def _currentPeaksCallback(self, data):
        "Callback for the peaks notifier"
        peaks = data[Notifier.VALUE]
        # self._setPeaksLabel()
        self._updateWidget()
        self._nmrResidue.update()

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

    def _predictHighlight(self, peaks: typing.List[Peak]):

        self._returnButtonsToNormal()
        # if self.current.nmrResidue is None or len(peaks) == 0:
        #     self._assignWidgetHide()
        #     return
        #
        # # make sure that the widget is visible
        # self._assignWidgetShow()

        # make sure that you have buttons!
        if len(self.buttonGroup.buttons()) == 0:
            return

        # check if peaks coincide
        for dim in range(peaks[0].peakList.spectrum.dimensionCount):
            if not peaksAreOnLine(peaks, dim):
                logger.debug('dimension %s: peaksAreonLine=False' % dim)
                return

        types = set(peak.peakList.spectrum.experimentType for peak in peaks)
        anyInterOnlyExperiments = any(isInterOnlyExpt(x) for x in types)

        logger.debug('peaks=%s' % (peaks,))
        logger.debug('types=%s, anyInterOnlyExperiments=%s' % (types, anyInterOnlyExperiments))

        peak = peaks[0]
        peakListViews = [peakListView for peakListView in self.project.peakListViews if
                         peakListView.peakList == peak.peakList]

        if peakListViews:
            spectrumIndices = peakListViews[0].spectrumView.axes

            # for the 1D case, this is (0, None)
            if spectrumIndices[1] is None:
                return

            isotopeCode = peak.peakList.spectrum.isotopeCodes[spectrumIndices[1]]

            # backbone
            if self.selectBackboneButton.isChecked():
                predictedAtomTypes = [
                    getNmrAtomPrediction(ccpCode, peak.position[spectrumIndices[1]], isotopeCode, strict=True)
                    for ccpCode in CCP_CODES]
                refinedPreds = [(type[0][0][1], type[0][1]) for type in predictedAtomTypes if len(type) > 0]
                atomPredictions = set()
                for atomPred, score in refinedPreds:
                    if score > 90:
                        atomPredictions.add(atomPred)

                # list containing those atoms that exist - used for colouring in 'checkAssignedAtoms'
                foundPredictList = {}
                for atomPred in atomPredictions:
                    if atomPred == 'CB' and self.buttons['CB']:
                        if anyInterOnlyExperiments:
                            self.buttons['CB'][0].setStyleSheet(GREEN_BUTTON)
                            foundPredictList[self.atomLabel('CB', '-1')] = 100
                        else:
                            self.buttons['CB'][0].setStyleSheet(GREEN_BUTTON)
                            self.buttons['CB'][1].setStyleSheet(GREEN_BUTTON)
                            foundPredictList[self.atomLabel('CB', '-1')] = 100
                            foundPredictList[self.atomLabel('CB', '0')] = 100
                    if atomPred == 'CA' and self.buttons['CA']:
                        if anyInterOnlyExperiments:
                            self.buttons['CA'][0].setStyleSheet(GREEN_BUTTON)
                            foundPredictList[self.atomLabel('CA', '-1')] = 100
                        else:
                            self.buttons['CA'][0].setStyleSheet(GREEN_BUTTON)
                            self.buttons['CA'][1].setStyleSheet(GREEN_BUTTON)
                            foundPredictList[self.atomLabel('CA', '-1')] = 100
                            foundPredictList[self.atomLabel('CA', '0')] = 100

                # new routine to colour any existing atoms
                # foundAtoms = self.checkAssignedAtoms(self.current.nmrResidue, ATOM_TYPES,
                #                                      foundPredictList, 'backbone')

            # sidechain is checked
            elif self.selectAllButton.isChecked():
                foundPredictList = {}

                if self.current.nmrResidue.residueType == '':
                    # In this case, we loop over all CCP_CODES (i.e. residue types)
                    predictedAtomTypes = []
                    for residueType in CCP_CODES:
                        for type, score in getNmrAtomPrediction(residueType, peak.position[spectrumIndices[1]],
                                                                isotopeCode):
                            if len(type) > 0 and score > 50:
                                predictedAtomTypes.append((type, score))

                else:
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
                for type, score in predictedAtomTypes:
                    if type[1] not in predictedDict:
                        predictedDict[type[1]] = (type[0], score)
                    else:
                        if score > predictedDict[type[1]][1]:
                            predictedDict[type[1]] = (type[0], score)
                # print ('>>>predictedDict', predictedDict)

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
                                                button.setStyleSheet(GREEN_BUTTON)
                                            elif 50 < score < 85:
                                                button.setStyleSheet(ORANGE_BUTTON)
                                            if score < 50:
                                                button.setStyleSheet(RED_BUTTON)

                                    else:

                                        # print('>type[1], atomType, button>', atomDictType, bText)
                                        foundPredictList[self.atomLabel(atomDictType, currentOffset)] = score

                                        if score >= 85:
                                            button.setStyleSheet(GREEN_BUTTON)
                                        elif 50 < score < 85:
                                            button.setStyleSheet(ORANGE_BUTTON)
                                        if score < 50:
                                            button.setStyleSheet(RED_BUTTON)

                # new routine to colour any existing atoms
                # atomButtonList = self._getAtomButtonList()
                # atomButtonList = [x for i in atomButtonList for x in i]
                # foundAtoms = self.checkAssignedAtoms(self.current.nmrResidue, atomButtonList,
                #                                      foundPredictList, 'sideChain')

    def _changeMoleculeType(self, data):
        """
        change the  available atomList depending on the moleculeType
        :param data - str from pullDown:
        """
        pass

    def getResidueTypes(self, moleculeType: str = 'protein'):
        """
        return a list of residue types assiciated with the moleculeType
        :param moleculeType - str ['protein', 'DNA', 'RNA', 'carbohydrate', 'other']
        :return list of str:
        """
        if moleculeType in MOLECULE_TYPES:
            if moleculeType == 'protein':
                return [atomName for atomName in PROTEIN_NEF_ATOM_NAMES.keys()]
        else:
            return None


#TODO: clean this up to a proper place
DNA_ATOMS = """
Res    Name     Atom     Count        Min.        Max.       Avg.     Std Dev

DA          H2       H        950          2.60         8.64       7.49       0.80       
DA          H61      H        145          2.56        13.40       6.83       1.16       
DA          H62      H        144          5.31        12.90       6.92       1.11       
DA          H8       H        1097         2.35         9.10       8.00       0.83       
DA          H1'      H        1092         4.19         8.50       6.04       0.33       
DA          H2'      H        1083         0.84         4.82       2.66       0.37       
DA          H2''     H        1065         0.56         6.89       2.79       0.37       
DA          H3'      H        1031         4.00         6.04       4.94       0.19       
DA          H4'      H        819          2.08         9.13       4.34       0.34       
DA          H5'      H        449          2.79         5.95       4.06       0.30       
DA          H5''     H        406          2.96         8.41       4.07       0.45       
DA          C2       C        134        151.00       159.33     154.69       1.25       
DA          C4       C        2          150.30       150.60     150.45       0.21       
DA          C5       C        2          119.80       119.90     119.85       0.07       
DA          C6       C        2          157.00       157.40     157.20       0.28       
DA          C8       C        145        138.00       142.79     141.33       0.94       
DA          C1'      C        138         81.50        91.90      85.19       1.26       
DA          C2'      C        123         36.23        42.93      40.58       1.05       
DA          C3'      C        99          72.71        83.26      78.46       1.70       
DA          C4'      C        68          84.43        91.33      87.21       1.19       
DA          C5'      C        45          47.77        70.96      67.18       3.42       
DA          N1       N        3          223.40       227.00     225.20       1.80       
DA          N3       N        3          214.50       216.40     215.73       1.07       
DA          N6       N        9           77.40        81.70      80.22       1.29       
DA          N7       N        2          233.50       234.10     233.80       0.42       
DA          N9       N        1          170.20       170.20     170.20       0.00       
DA          P        P        239        -45.29        35.48      -2.32       5.04       

DC          H41      H        593         -0.26        12.08       7.40       1.15       
DC          H42      H        576         -0.26        10.62       7.48       1.07       
DC          H5       H        1263        -6.22         8.25       5.40       0.89       
DC          H6       H        1277       -14.72         8.31       7.36       0.98       
DC          H1'      H        1269         2.17         6.60       5.72       0.69       
DC          H2'      H        1251       -29.45         6.93       2.11       1.42       
DC          H2''     H        1235       -20.82         8.48       2.39       1.08       
DC          H3'      H        1177        -7.37         8.62       4.72       0.65       
DC          H4'      H        923        -38.59        28.54       4.18       2.40       
DC          H5'      H        458        -18.35        11.19       4.17       1.81       
DC          H5''     H        381        -41.21        15.96       3.98       3.50       
DC          C2       C        2          158.80       159.20     159.00       0.28       
DC          C4       C        2          167.80       168.00     167.90       0.14       
DC          C5       C        125         84.55        99.90      98.12       2.25       
DC          C6       C        129        133.73       144.70     142.76       1.63       
DC          C1'      C        134         77.90        88.53      86.94       1.65       
DC          C2'      C        105         30.62        42.27      40.34       1.45       
DC          C3'      C        101         66.41        80.30      76.57       2.75       
DC          C4'      C        52          77.80        88.39      85.90       1.93       
DC          C5'      C        33          59.30        71.72      66.28       2.63       
DC          N1       N        1          150.70       150.70     150.70       0.00       
DC          N3       N        1          196.30       196.30     196.30       0.00       
DC          N4       N        16          95.10       100.41      98.29       1.33       
DC          P        P        332        -14.20        52.72      -1.65       6.05       

DG          H1       H        1299        -2.91        13.94      11.62       2.06       
DG          H21      H        322        -26.86        18.44       7.46       3.00       
DG          H22      H        309        -26.86        18.44       6.36       3.02       
DG          H8       H        1975         6.07         9.31       7.78       0.25       
DG          H1'      H        1979         2.53        12.45       5.90       0.45       
DG          H2'      H        1932         1.43         8.23       2.67       0.41       
DG          H2''     H        1902         1.16        11.50       2.72       0.42       
DG          H3'      H        1863         0.00        13.24       4.97       0.68       
DG          H4'      H        1547         0.00        13.26       4.46       0.98       
DG          H5'      H        885          0.00         7.82       4.12       0.36       
DG          H5''     H        761          0.00         7.76       4.10       0.35       
DG          C2       C        2          156.50       156.70     156.60       0.14       
DG          C4       C        2          153.10       154.00     153.55       0.64       
DG          C5       C        2          117.40       117.60     117.50       0.14       
DG          C6       C        2          161.00       161.40     161.20       0.28       
DG          C8       C        182        134.38       142.93     138.20       1.62       
DG          C1'      C        166         74.90        91.82      85.07       2.02       
DG          C2'      C        128         32.64        49.00      40.23       2.00       
DG          C3'      C        124         67.72        81.41      77.69       2.27       
DG          C4'      C        77          75.27        90.64      86.92       2.19       
DG          C5'      C        57          58.50        71.42      67.06       2.21       
DG          N1       N        33         142.74       147.70     145.74       1.96       
DG          N2       N        5           75.10        76.10      75.52       0.38       
DG          N7       N        4          236.90       238.30     237.28       0.68       
DG          N9       N        1          168.50       168.50     168.50       0.00       
DG          P        P        406        -20.92        47.10      -1.39       6.34       

DT          H3       H        575          2.12        14.53      12.85       1.97       
DT          H6       H        1420         1.61         8.10       7.36       0.32       
DT          H7       H        8            1.46         1.79       1.65       0.13       
DT          H71      H        1249         0.18         7.53       1.60       0.41       
DT          H72      H        1248         0.18         7.53       1.60       0.41       
DT          H73      H        1248         0.18         7.53       1.60       0.41       
DT          H1'      H        1412         1.65         6.72       5.88       0.54       
DT          H2'      H        1402         0.50         6.21       2.16       0.47       
DT          H2''     H        1383         0.78         4.87       2.39       0.36       
DT          H3'      H        1337         3.79        14.11       4.80       0.47       
DT          H4'      H        997          2.00        14.15       4.16       0.55       
DT          H5'      H        608          2.54        14.20       4.20       1.56       
DT          H5''     H        547          0.98         7.19       3.97       0.40       
DT          C2       C        2          152.90       153.40     153.15       0.35       
DT          C4       C        2          168.50       169.10     168.80       0.42       
DT          C5       C        9           14.30       113.60      36.43      43.64       
DT          C6       C        229        135.80       143.60     139.24       0.92       
DT          C7       C        59          11.96        14.90      14.19       0.62       
DT          C1'      C        215         80.60        92.22      86.88       1.37       
DT          C2'      C        195         34.38        42.96      39.99       0.88       
DT          C3'      C        176         71.58        80.10      77.49       1.64       
DT          C4'      C        64          78.80        89.31      86.04       1.46       
DT          C5'      C        41          59.40        68.40      66.36       1.91       
DT          N1       N        1          142.90       142.90     142.90       0.00       
DT          N3       N        41         156.40       160.50     159.27       0.71       
DT          P        P        261         -5.65        47.79      -1.36       7.03       
"""

RNA_ATOMS = """
Res    Name     Atom     Count        Min.        Max.       Avg.     Std Dev

A          H2       H        1617         0.00         9.39       7.63       0.47       
A          H61      H        298          0.00        10.85       7.41       1.18       
A          H62      H        265          0.00         9.98       6.78       1.11       
A          H8       H        1592         7.00         9.20       8.01       0.27       
A          H1'      H        1616         4.56         7.16       5.88       0.21       
A          HO2'     H        44           4.27         7.75       5.75       1.18       
A          H2'      H        1414         3.61         5.91       4.59       0.22       
A          H3'      H        1213         3.91         5.45       4.63       0.21       
A          H4'      H        1042         0.00         5.29       4.46       0.25       
A          H5'      H        803          0.00         5.03       4.30       0.27       
A          H5''     H        786          0.00       153.90       4.37       5.35       
A          C2       C        1015         0.00       167.71     151.71       8.98       
A          C4       C        35           0.00       151.89     131.04      47.78       
A          C5       C        35           0.00       153.63     104.96      39.80       
A          C6       C        37           0.00       159.60     139.05      49.16       
A          C8       C        993          0.00       145.80     138.18       9.80       
A          C1'      C        944         84.23       130.39      91.64       3.57       
A          C2'      C        703          4.75        99.00      75.25       5.81       
A          C3'      C        640          0.00       101.40      73.68       4.97       
A          C4'      C        629          0.00        92.60      82.53       4.09       
A          C5'      C        540          0.00       109.20      66.22       6.79       
A          N1       N        253          0.00       229.70     216.11      30.52       
A          N3       N        203          0.00       231.08     210.51      30.33       
A          N6       N        157          0.00        97.01      79.33      13.30       
A          N7       N        151          0.00       237.10     224.40      37.21       
A          N9       N        189          0.00       176.80     166.48      24.64       
A          P        P        176         -5.96         2.43      -2.97       1.64       

C          H41      H        973          5.64        10.38       7.91       0.72       
C          H42      H        951          5.07         9.06       7.25       0.69       
C          H5       H        1865         4.20         7.09       5.48       0.29       
C          H6       H        1915         5.62         8.46       7.68       0.21       
C          H1'      H        1871         3.42         6.42       5.55       0.24       
C          HO2'     H        66           3.96         9.70       5.82       1.49       
C          H2'      H        1578         0.00         5.07       4.32       0.24       
C          H3'      H        1389         0.00         5.63       4.41       0.23       
C          H4'      H        1135         0.00        12.58       4.33       0.37       
C          H5'      H        850          0.00         7.17       4.26       0.38       
C          H5''     H        871          0.00         5.04       4.10       0.46       
C          C2       C        65           0.00       189.36     136.58      59.06       
C          C4       C        78           0.00       169.52     147.08      53.62       
C          C5       C        1101         0.00       141.05      97.28       4.33       
C          C6       C        1132        93.40       145.02     139.71       7.35       
C          C1'      C        1076        84.00       118.95      92.68       2.10       
C          C2'      C        839          0.00        99.40      75.76       5.32       
C          C3'      C        780          0.00       104.90      73.18       6.41       
C          C4'      C        735          0.00        92.60      81.94       5.86       
C          C5'      C        609          0.00       110.20      64.72       9.47       
C          N1       N        229          0.00       178.83     145.28      30.98       
C          N3       N        181          0.00       204.30     184.53      43.16       
C          N4       N        390          0.00       104.50      96.67      10.21       
C          P        P        186         -5.13         0.62      -3.15       1.65       

G          H1       H        1527         6.09        14.32      12.35       0.98       
G          H21      H        444          0.00         9.39       7.37       1.42       
G          H22      H        402          0.00         9.06       6.31       1.22       
G          H8       H        2293         0.00         8.63       7.60       0.39       
G          H1'      H        2257         3.41         7.62       5.66       0.33       
G          HO2'     H        53           4.40         7.10       5.77       1.16       
G          H2'      H        1931         3.26         6.30       4.58       0.25       
G          H3'      H        1724         3.80         5.78       4.56       0.26       
G          H4'      H        1405         0.00         5.11       4.43       0.22       
G          H5'      H        1187         2.89         5.42       4.28       0.23       
G          H5''     H        1157         2.55         5.11       4.18       0.23       
G          C2       C        59           0.00       162.50     131.48      56.36       
G          C4       C        44           0.00       158.62     120.73      61.96       
G          C5       C        78           0.00       167.24     106.95      40.80       
G          C6       C        86           0.00       162.74     142.40      49.29       
G          C8       C        1429         0.00       146.73     135.38       8.27       
G          C1'      C        1244        79.80        95.05      91.50       1.97       
G          C2'      C        975          0.00        99.40      75.27       4.88       
G          C3'      C        887          0.00       102.20      73.85       5.27       
G          C4'      C        855         72.08        96.10      82.65       2.63       
G          C5'      C        758         50.36       108.20      66.41       6.34       
G          N1       N        891          0.00       166.07     145.58      13.10       
G          N2       N        151          0.00        83.38      71.00      17.00       
G          N3       N        25           0.00       234.10      97.88      80.79       
G          N7       N        184          0.00       240.99     221.87      50.58       
G          N9       N        280          0.00       176.40     164.01      30.01       
G          P        P        208         -6.00         0.56      -2.86       1.59       

U          H3       H        898          1.20       154.54      13.24       4.89       
U          H5       H        1524         4.20         7.83       5.47       0.31       
U          H6       H        1590         5.81         8.52       7.75       0.20       
U          H1'      H        1559         3.69         6.59       5.61       0.26       
U          HO2'     H        45           4.08         9.74       6.06       1.51       
U          H2'      H        1357         0.00         6.63       4.37       0.27       
U          H3'      H        1117         0.00         7.83       4.48       0.25       
U          H4'      H        949          0.00         4.82       4.36       0.23       
U          H5'      H        697          0.00         4.85       4.24       0.30       
U          H5''     H        713          0.00         4.86       4.15       0.29       
U          C2       C        80           0.00       183.20     144.31      38.05       
U          C4       C        86           0.00       169.70     157.82      39.47       
U          C5       C        944          0.00       154.13     102.90       6.62       
U          C6       C        975          0.00       169.95     140.14       8.53       
U          C1'      C        990         80.62        96.30      91.90       2.14       
U          C2'      C        737          0.00        99.80      75.04       4.49       
U          C3'      C        656          0.00       102.50      73.74       5.43       
U          C4'      C        643          0.00        94.42      82.70       4.14       
U          C5'      C        525          0.00       106.60      65.07       5.67       
U          N1       N        209          0.00       192.14     147.96      25.60       
U          N3       N        553          0.00       167.30     158.99      16.67       
U          O4       O        5            0.00         0.00       0.00       0.00       
U          P        P        149         -5.30         1.58      -2.99       1.78       
"""

DNA_ATOM_NAMES = {
    'DT': ['H3', 'H6', 'H7', 'H71', 'H72', 'H73', "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''",
           'C2', 'C4', 'C5', 'C6', 'C7', "C1'", "C2'", "C3'", "C4'", "C5'",
           'N1', 'N3',
           'P'],
    'DC': ['H41', 'H42', 'H5', 'H6', "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''",
           'C2', 'C4', 'C5', 'C6', "C1'", "C2'", "C3'", "C4'", "C5'",
           'N1', 'N3', 'N4',
           'P'],
    'DA': ['H2', 'H61', 'H62', 'H8', "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''",
           'C2', 'C4', 'C5', 'C6', 'C8', "C1'", "C2'", "C3'", "C4'", "C5'",
           'N1', 'N3', 'N6', 'N7', 'N9', 'P'],
    'DG': ['H1', 'H21', 'H22', 'H8', "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''",
           'C2', 'C4', 'C5', 'C6', 'C8', "C1'", "C2'", "C3'", "C4'", "C5'",
           'N1', 'N2', 'N7', 'N9',
           'P']
    }
RNA_ATOM_NAMES = {
    'G': ['H1', 'H21', 'H22', 'H8', "H1'", "HO2'", "H2'", "H3'", "H4'", "H5'", "H5''",
          'C2', 'C4', 'C5', 'C6', 'C8', "C1'", "C2'", "C3'", "C4'", "C5'",
          'N1', 'N2', 'N3', 'N7', 'N9',
          'P'],
    'U': ['H3', 'H5', 'H6', "H1'", "HO2'", "H2'", "H3'", "H4'", "H5'", "H5''",
          'C2', 'C4', 'C5', 'C6', "C1'", "C2'", "C3'", "C4'", "C5'",
          'N1', 'N3', 'O4',
          'P'],
    'A': ['H2', 'H61', 'H62', 'H8', "H1'", "HO2'", "H2'", "H3'", "H4'", "H5'", "H5''",
          'C2', 'C4', 'C5', 'C6', 'C8', "C1'", "C2'", "C3'", "C4'", "C5'",
          'N1', 'N3', 'N6', 'N7', 'N9',
          'P'],
    'C': ['H41', 'H42', 'H5', 'H6', "H1'", "HO2'", "H2'", "H3'", "H4'", "H5'", "H5''",
          'C2', 'C4', 'C5', 'C6', "C1'", "C2'", "C3'", "C4'", "C5'",
          'N1', 'N3', 'N4',
          'P']
    }

ALL_DNARNA_ATOMS_SORTED = OrderedDict([
    ('O', ["HO2'"]),
    ('1', ["C1'", 'H1', "H1'", 'N1']),
    ('2', ['C2', "C2'", 'H2', "H2'", "H2''", 'H21', 'H22', 'N2']),
    ('3', ["C3'", 'H3', "H3'", 'N3']),
    ('4', ['C4', "C4'", "H4'", 'H41', 'H42', 'N4', 'O4']),
    ('5', ['C5', "C5'", 'H5', "H5'", "H5''"]),
    ('6', ['C6', 'H6', 'H61', 'H62', 'N6']),
    ('7', ['C7', 'H7', 'H71', 'H72', 'H73', 'N7']),
    ('8', ['C8', 'H8']),
    ('9', ['N9']),
    ('P', ['P'])
    ])

if __name__ == '__main__':
    from ccpn.ui.gui.widgets.Application import TestApplication
    from ccpn.ui.gui.widgets.TextEditor import TextEditor


    app = TestApplication()

    popup = NmrAtomAssignerModule()

    textBox = TextEditor(popup.mainWidget, grid=(3, 0), gridSpan=(1, 1))

    all_atoms = dict()
    for atom in ['O', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'P']:
        all_atoms[atom] = []

    atomList = [DNA_ATOMS, RNA_ATOMS]
    startAtoms = ['DA', 'DC', 'DG', 'DT', 'A', 'G', 'C', 'U']

    for atomText in atomList:
        atoms = {}
        atomText = atomText.split('\n')
        for line in atomText:
            ll = line.split()
            if ll:
                if ll[0] in startAtoms:
                    if ll[0] in atoms.keys():
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
