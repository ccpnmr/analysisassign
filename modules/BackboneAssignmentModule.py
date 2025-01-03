"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (https://www.ccpn.ac.uk) 2014 - 2025"
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
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2025-01-03 18:50:14 +0000 (Fri, January 03, 2025) $"
__version__ = "$Revision: 3.2.11 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:40 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import typing
from collections import OrderedDict
from PyQt5 import QtWidgets, QtCore
from ccpn.AnalysisAssign.lib.scoring import getNmrResidueMatches
from ccpn.core.ChemicalShift import ChemicalShift
from ccpn.core.NmrResidue import NmrResidue
from ccpn.core.NmrChain import NmrChain
from ccpn.core.lib.ContextManagers import undoBlockWithoutSideBar
from ccpn.ui.gui.guiSettings import getColours, DIVIDER
from ccpn.ui.gui.lib.SpectrumDisplayLib import makeStripPlot
from ccpn.ui.gui.lib.StripLib import matchAxesAndNmrAtoms, markNmrAtoms
from ccpn.ui.gui.lib.StripLib import navigateToNmrResidueInDisplay
from ccpn.ui.gui.lib.alignWidgets import alignWidgets
from ccpn.ui.gui.modules.NmrResidueTable import NmrResidueTableModule, LINKTOPULLDOWNCLASS
from ccpn.ui.gui.widgets.CheckBox import CheckBox
from ccpn.ui.gui.widgets.CompoundWidgets import PulldownListCompoundWidget, CheckBoxCompoundWidget, \
    SpinBoxCompoundWidget
from ccpn.ui.gui.widgets.MessageDialog import showWarning, progressManager, showYesNo
from ccpn.ui.gui.widgets.PulldownListsForObjects import ChemicalShiftListPulldown
from ccpn.ui.gui.widgets.DropBase import DropBase
from ccpn.ui.gui.widgets.Font import getTextDimensionsFromFont
from ccpn.ui.gui.widgets.PlaneToolbar import STRIPLABEL_CONNECTDIR, STRIPLABEL_CONNECTNONE, \
    STRIPCONNECT_LEFT, STRIPCONNECT_RIGHT
from ccpn.ui.gui.widgets.Tabs import Tabs
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.ui.gui.widgets.HLine import LabeledHLine, HLine
from ccpn.util.decorators import logCommand
from ccpn.util.Logging import getLogger


ALL = '<Use all>'
MINMATCHES = 1
MAXMATCHES = 20
DEFAULTMATCHES = 3
STRIPBACKBONE = 'backboneAssignment'
MARKCONNECTED = False
EXTRAWIDTH = 150
EXTRAOFFSET = 100


class BackboneAssignmentModule(NmrResidueTableModule):
    """Class implementing the module
    """
    className = 'BackboneAssignmentModule'

    includeSettingsWidget = True
    maxSettingsState = 2  # states are defined as: 0: invisible, 1: both visible, 2: only settings visible
    settingsPosition = 'left'
    settingsMinimumSizes = (500, 200)

    includeDisplaySettings = True
    activePulldownClass = NmrChain
    activePulldownInitialState = True

    registeredExtensions = set()

    def __init__(self, mainWindow=None, name='Backbone Assignment'):
        """Initialise the module widgets
        """
        super().__init__(mainWindow=mainWindow, name=name, selectFirstItem=True)

        # Derive application, project, and current from mainWindow
        self.mainWindow = mainWindow
        if mainWindow:
            self.application = mainWindow.application
            self.project = mainWindow.application.project
            self.current = mainWindow.application.current
            self.nmrChains = self.application.project.nmrChains
        else:
            self.application = self.project = self.current = self.nmrChains = None

        # add a new checkbox to the header in the main-widget area
        self.matchCheckBoxWidget = CheckBox(self.tableFrame, grid=(1, 2), checked=True, text='Find matches')
        self.tableFrame.addWidgetToPos(self.matchCheckBoxWidget, row=0, col=2)

        self._createSettingsWidgets()
        self._stripNotifiers = []  # list to store GuiNotifiers for strips

        ## main table options
        self._tableWidget.multiSelect = True
        self._tableWidget.setSelectionMode(self._tableWidget.SingleSelection)
        self._tableWidget.setActionCallback(self.navigateToNmrResidueCallBack)
        self.mainWidget.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
        self.layout.setContentsMargins(0, 1, 0, 0)

    def _createSettingsWidgets(self):
        self.settingsWidget.setContentsMargins(5, 5, 5, 5)
        self.settingsTabWidget = Tabs(self.settingsWidget, setLayout=True, grid=(0, 0))
        ##  General Tab
        self.generalSettingsFrame = Frame(self.settingsWidget, setLayout=True)
        self.generalSettingsFrame.getLayout().setAlignment(QtCore.Qt.AlignTop)
        self.settingsTabWidget.addTab(self.generalSettingsFrame, 'General')
        self.generalSettingsFrame.getLayout().addWidget(self.nmrResidueTableSettings)
        self._setupGeneralSettings()

        ##  Extensions Tab
        self.extensionsSettingsFrame = Frame(self.settingsWidget, setLayout=True)
        self.extensionsSettingsFrame.getLayout().setAlignment(QtCore.Qt.AlignTop)
        self.settingsTabWidget.addTab(self.extensionsSettingsFrame, 'Extensions')
        self._addExtensionsToSettings()

    def _setupGeneralSettings(self):
        """add to the layout of the general settings widgets"""
        ### Settings ###

        # change defaults setting inherited from NmrResidueTableModule
        self.nmrResidueTableSettings.sequentialStripsWidget.checkBox.setChecked(True)
        if self.nmrResidueTableSettings.displaysWidget:
            self.nmrResidueTableSettings.displaysWidget.addPulldownItem(0)

        # colWidth0 = 180
        texts = ['i-1 Matches to show:',
                 'i+1 Matches to show:',
                 'Match SpectrumDisplay',
                 'Search SpectrumDisplay',
                 'Match CA NmrAtoms',
                 'Match CB NmrAtoms',
                 'Match C NmrAtoms',
                 'ChemicalShiftList']
        _, maxDim = getTextDimensionsFromFont(textList=texts)
        colWidth0 = maxDim.width()

        focusRow = row = self.nmrResidueTableSettings.maxRows  ## Number of widgets of NmrResidueTable - add extra widgets below
        col = 0

        row += 1
        HLine(parent=self.nmrResidueTableSettings, grid=(row, 0), gridSpan=(1, 2), colour=getColours()[DIVIDER],
              height=15)

        # new match module pulldown list
        row += 1
        self.matchWidget = PulldownListCompoundWidget(self.nmrResidueTableSettings, labelText=texts[2],
                                                      fixedWidths=(colWidth0, colWidth0, None), grid=(row, col),
                                                      gridSpan=(1, 2),
                                                      vAlign='top', hAlign='left',
                                                      )
        self.matchWidget.setPreSelect(self._fillMatchWidget)
        self._fillMatchWidget()
        self.matchWidget.pulldownList.setIndex(0)

        # Number of matches to show
        row += 1
        self.numberOfMinusMatchesWidget = SpinBoxCompoundWidget(self.nmrResidueTableSettings,
                                                                grid=(row, col), gridSpan=(1, 2),
                                                                vAlign='top', hAlign='left',
                                                                fixedWidths=(colWidth0, colWidth0 // 3, None),
                                                                labelText=texts[0],
                                                                minimum=1, maximum=MAXMATCHES,
                                                                value=DEFAULTMATCHES
                                                                )
        row += 1
        self.numberOfPlusMatchesWidget = SpinBoxCompoundWidget(self.nmrResidueTableSettings,
                                                               grid=(row, col), gridSpan=(1, 2),
                                                               vAlign='top', hAlign='left',
                                                               fixedWidths=(colWidth0, colWidth0 // 3, None),
                                                               labelText=texts[1],
                                                               minimum=1, maximum=MAXMATCHES,
                                                               value=DEFAULTMATCHES
                                                               )

        row += 1
        self.showSearchInMatch = CheckBoxCompoundWidget(self.nmrResidueTableSettings,
                                                        grid=(row, col), gridSpan=(1, 2), vAlign='top', hAlign='left',
                                                        fixedWidths=(colWidth0, None),
                                                        orientation='left',
                                                        labelText='Show Search Strip in Match Module',
                                                        checked=False
                                                        )

        # new search module pulldown list
        row += 1
        self.targetWidget = PulldownListCompoundWidget(self.nmrResidueTableSettings, labelText=texts[3],
                                                       fixedWidths=(colWidth0, colWidth0, None), grid=(row, col),
                                                       gridSpan=(1, 2),
                                                       vAlign='top', hAlign='left',
                                                       )
        self.targetWidget.setPreSelect(self._fillTargetWidget)
        self._fillTargetWidget()
        self.targetWidget.pulldownList.setIndex(0)

        row += 1
        self.focusYAxis = CheckBoxCompoundWidget(self.nmrResidueTableSettings,
                                                 grid=(row, col), gridSpan=(1, 2), vAlign='top', hAlign='left',
                                                 fixedWidths=(colWidth0, None),
                                                 orientation='left',
                                                 labelText='Focus Y-Axis',
                                                 checked=True
                                                 )

        # re-order, move the sequential checkbox to here, swap with focusYAxis checkbox - not nice method :|
        row += 1
        self.nmrResidueTableSettings.layout().addWidget(self.nmrResidueTableSettings.sequentialStripsWidget, row, col,
                                                        1, 2)
        self.nmrResidueTableSettings.layout().addWidget(self.focusYAxis, focusRow, col, 1, 2)

        # Select which NmrAtoms to match
        # VAH: This could probably be improved by putting the check-boxes into a group "NmrAtoms to Match"
        # and then automatically using the label text 'CA', 'CB' etc. in the _setNmrAtomsToMatch function.
        row += 1
        self.matchCA = CheckBoxCompoundWidget(self.nmrResidueTableSettings,
                                              grid=(row, col), gridSpan=(1, 2), vAlign='top', hAlign='left',
                                              fixedWidths=(colWidth0, None),
                                              orientation='left',
                                              labelText='Match CA NmrAtoms',
                                              callback=self._setNmrAtomsToMatch,
                                              checked=True
                                              )
        row += 1
        self.matchCB = CheckBoxCompoundWidget(self.nmrResidueTableSettings,
                                              grid=(row, col), gridSpan=(1, 2), vAlign='top', hAlign='left',
                                              fixedWidths=(colWidth0, None),
                                              orientation='left',
                                              labelText='Match CB NmrAtoms',
                                              callback=self._setNmrAtomsToMatch,
                                              checked=True
                                              )
        row += 1
        self.matchC = CheckBoxCompoundWidget(self.nmrResidueTableSettings,
                                             grid=(row, col), gridSpan=(1, 2), vAlign='top', hAlign='left',
                                             fixedWidths=(colWidth0, None),
                                             orientation='left',
                                             labelText='Match C NmrAtoms',
                                             callback=self._setNmrAtomsToMatch,
                                             checked=False
                                             )
        self._setNmrAtomsToMatch()

        row += 1
        HLine(parent=self.nmrResidueTableSettings, grid=(row, 0), gridSpan=(1, 2), colour=getColours()[DIVIDER],
              height=15)

        # Chemical shift list selection
        row += 1
        self.shiftListWidget = ChemicalShiftListPulldown(self.nmrResidueTableSettings, self.mainWindow,
                                                         grid=(row, col), gridSpan=(1, 2),
                                                         vAlign='top', hAlign='left',
                                                         fixedWidths=(colWidth0, colWidth0, None),
                                                         callback=self._setupShiftDicts, default=None
                                                         )
        self._setupShiftDicts()

        self._activeLinkCheckbox = self.activePulldownClass and getattr(self.nmrResidueTableSettings,
                                                                        LINKTOPULLDOWNCLASS, None)

        # align the widgets in the settings-widget
        alignWidgets(self.nmrResidueTableSettings)

    @staticmethod
    def registerExtension(cls, extension):
        from ccpn.AnalysisAssign.modules.backboneExtensions.BackboneAssignmentExtensionABC import \
            BackboneAssignmentExtensionFrame

        if issubclass(extension, BackboneAssignmentExtensionFrame):
            cls.registeredExtensions.add(extension)
        else:
            getLogger().warning('Cannot register Extension for this module. Ensure the format is correct')

    def _addExtensionsToSettings(self):
        """ Add registered extensions to the Settings Panel. """
        from ccpn.AnalysisAssign.modules.backboneExtensions import _loadAssignExtensions

        try:
            _loadAssignExtensions()

            registeredExtensions = self.registeredExtensions
            for extensionFrameObj in registeredExtensions:
                extensionFrame = extensionFrameObj(guiModule=self)
                hLine = LabeledHLine(self, text=extensionFrame.NAME)
                self.extensionsSettingsFrame.getLayout().addWidget(hLine)
                self.extensionsSettingsFrame.getLayout().addWidget(extensionFrame)
        except Exception as err:
            getLogger().warning(f"Some Extensions failed to load {err}")

    def _fillMatchWidget(self):
        ll = ['> select-to-add <'] + [display.pid for display in self.mainWindow.spectrumDisplays]
        thisText = self.matchWidget.getText()
        self.matchWidget.pulldownList.setData(texts=ll)
        if thisText:
            self.matchWidget.select(thisText, True)

    def _fillTargetWidget(self):
        ll = ['> select-to-add <'] + [display.pid for display in self.mainWindow.spectrumDisplays]
        thisText = self.targetWidget.getText()
        self.targetWidget.pulldownList.setData(texts=ll)
        if thisText:
            self.targetWidget.select(thisText, True)

    def _getDisplays(self):
        """return list of displays to navigate"""
        displays = []

        if self.nmrResidueTableSettings and self.nmrResidueTableSettings.displaysWidget:
            dGids = self.nmrResidueTableSettings.displaysWidget.getTexts()  # gids of displays
            if len(dGids) == 0: return displays

            matchGids = self.matchWidget.getText()  # gid of the match module
            targetGids = self.targetWidget.getText()  # gid of the target module - don't discard for the minute

            if ALL in dGids:
                displays = [dp for dp in self.application.ui.mainWindow.spectrumDisplays if
                            dp.pid not in (matchGids, targetGids)]
            else:
                displays = [self.application.getByGid(gid) for gid in dGids if gid not in (ALL, matchGids, targetGids)]

            displays = [display for display in displays if display is not None]

        return displays

    def _getMatchDisplays(self):
        """return list of displays to display matches
        """
        mGids = self.matchWidget.getText()  # gid of the match displays
        displays = [self.application.getByGid(gid) for gid in (mGids,)]
        displays = [display for display in displays if display is not None]
        return displays

    def _getTargetDisplays(self):
        """return list of displays to display targets
        """
        mGids = self.targetWidget.getText()  # gid of the match displays
        displays = [self.application.getByGid(gid) for gid in (mGids,)]
        displays = [display for display in displays if display is not None]
        return displays

    def navigateToNmrResidueCallBack(self, selection, lastItem):
        """Navigate in selected displays to nmrResidue; skip if none defined
        """
        try:
            if not (objs := list(lastItem[self._tableWidget._OBJECT])):
                return
        except Exception as es:
            getLogger().debug2(f'{self.__class__.__name__}.navigateToNmrResidueCallBack: No selection\n{es}')
            return

        nmrResidue = objs[0] if isinstance(objs, (tuple, list)) else objs
        self.navigateToNmrResidue(nmrResidue)

    @logCommand(get='self')
    def navigateToNmrResidue(self, nmrResidue):
        """Navigate in selected displays to nmrResidue; skip if no displays defined
        If matchCheckbox is checked, also call findAndDisplayMatches
        """
        displays = self._getDisplays()
        # if len(displays) == 0 and self.nmrResidueTableSettings.displaysWidget:
        #     getLogger().warning('Undefined display module(s); select in settings first')
        #     showWarning('startAssignment', 'Undefined display module(s);\nselect in settings first')
        #     return

        matchIndex = self.matchWidget.getIndex()
        targetIndex = self.targetWidget.getIndex()
        if self.matchCheckBoxWidget.isChecked() and matchIndex == 0:
            getLogger().warning('Undefined Match module; select Match module in Settings or unselect "Find matches"')
            showWarning('startAssignment',
                        'Undefined Match module;\nSelect your Match module in the Backbone Assignment Settings panel'
                        'or unselect "Find matches"')
            return

        if self.matchCheckBoxWidget.isChecked() and targetIndex == 0 and not self.showSearchInMatch.isChecked():
            getLogger().warning(
                    'Undefined Search module; select Search module in Settings, unselect "Find matches" or select "Show Search Strip in Match Module"')
            showWarning('startAssignment',
                        'Undefined Search module;\nSelect your Search module in the Backbone Assignment Settings panel, '
                        'unselect "Find matches" or\nselect "Show Search Strip in Match Module" in the Settings')
            return

        if (matchIndex == targetIndex) and matchIndex != 0:
            getLogger().warning('Match module and Search module cannot be the same')
            showWarning('startAssignment', 'Match module and Search module cannot be the same')
            return

        with undoBlockWithoutSideBar():

            # optionally clear the marks
            if self.nmrResidueTableSettings.autoClearMarksWidget.checkBox.isChecked():
                self.mainWindow.clearMarks()

            # clear any notifiers of previous strips
            for notifier in self._stripNotifiers:
                notifier.unRegister()
                del (notifier)
            self._stripNotifiers = []

            nr = nmrResidue.mainNmrResidue
            targetDisplays = self._getTargetDisplays()
            matchDisplays = self._getMatchDisplays()

            # navigate to the other displays - not matchDisplay
            for display in displays:

                # display.showAllStripHeaders()  # tag all headers with backboneAssignment module as handler

                if len(display.strips) > 0:

                    # if contains 2D's (e.g. a hsqc) then keep zoom
                    if display.spectrumViews[0].spectrum.dimensionCount <= 2:
                        newWidths = []  #_getCurrentZoomRatio(display.strips[0].viewBox.viewRange())
                    else:
                        # set the width in case of nD (n>2)
                        _widths = {'H': 2.5, 'C': 1.0, 'N': 1.0}
                        _ac = display.strips[0].axisCodes[0]
                        _w = _widths.setdefault(_ac[0], 1.0)
                        newWidths = [_w, 'full']

                    strips = navigateToNmrResidueInDisplay(nr, display, stripIndex=0,
                                                           widths=newWidths,
                                                           showSequentialResidues=(len(display.axisCodes) > 2) and
                                                                                  self.nmrResidueTableSettings.sequentialStripsWidget.checkBox.isChecked(),
                                                           markPositions=False,
                                                           #self.nmrResidueTableSettings.markPositionsWidget.checkBox.isChecked()
                                                           showDropHeaders=display in targetDisplays,
                                                           )

                    # for st, strip in enumerate(strips):
                    #     if strip is not None:
                    #         strip.header.handle = STRIPBACKBONE
                    #         strip.header.headerVisible = True
                    strips[0].spectrumDisplay.setColumnStretches(True)
                    # need better way to make sure that the floating axis updates
                    strips[0]._CcpnGLWidget.emitYAxisChanged(allStrips=True)

            # navigate to the targetDisplays
            for display in targetDisplays:

                display.showAllStripHeaders()  # tag all headers with backboneAssignment module as handler

                if len(display.strips) > 0:

                    # add a mask to ignore the position for the Y-axis
                    axisMask = [True] * len(display.axisCodes)
                    axisMask[1] = False

                    newWidths = []  #_getCurrentZoomRatio(display.strips[0].viewBox.viewRange())
                    strips = navigateToNmrResidueInDisplay(nr, display, stripIndex=0,
                                                           widths=newWidths,
                                                           showSequentialResidues=(len(display.axisCodes) > 2) and
                                                                                  self.nmrResidueTableSettings.sequentialStripsWidget.checkBox.isChecked(),
                                                           markPositions=False,
                                                           #self.nmrResidueTableSettings.markPositionsWidget.checkBox.isChecked()
                                                           showDropHeaders=True,
                                                           axisMask=axisMask
                                                           )
                    for strip in strips:
                        if strip is not None:
                            strip.header.handle = STRIPBACKBONE
                            strip.header.headerVisible = True
                    strips[0].spectrumDisplay.setColumnStretches(True)
                    strips[0]._CcpnGLWidget.emitYAxisChanged(allStrips=True)

                    # ejb
                    # if 'i-1' residue, take CA CB, and take H, N from the 'i' residue (.mainNmrResidue)
                    # check if contains '-1' in pid, is this robust? no :)
                    #
                    # VAH:
                    # Changed, so marks are drawn for the C atoms that are being matched and the base
                    # atoms specified here. Relies on use of NEF atom names, but makes it easier to
                    # make more generic at a later stage.
                    baseNmrAtoms = ['H', 'N']
                    if self.nmrResidueTableSettings.markPositionsWidget.checkBox.isChecked():
                        if nmrResidue.relativeOffset is not None and nmrResidue.relativeOffset != 0:
                            # offset residue (not necessarily i-1!) so need to split the match nmrAtoms
                            # (e.g. CA/CB) from the base nmrAtoms (e.g. N, H)
                            nmrAtomsOffset = nmrAtomsFromResidue(nmrResidue)
                            nmrAtomsCentre = nmrAtomsFromResidue(nmrResidue.mainNmrResidue)

                            nmrAtoms = [naOffset for naOffset in nmrAtomsOffset if
                                        naOffset.name in self.nmrAtomsToMatch]
                            nmrAtoms.extend(naCentre for naCentre in nmrAtomsCentre if naCentre.name in baseNmrAtoms)

                        elif MARKCONNECTED:
                            nmrAtoms = [na for na in nmrAtomsFromResidue(nmrResidue.mainNmrResidue)
                                        if na.name in self.nmrAtomsToMatch
                                        or na.name in baseNmrAtoms]
                        else:
                            nmrAtoms = [na for na in nmrResidue.mainNmrResidue.nmrAtoms
                                        if na.name in self.nmrAtomsToMatch
                                        or na.name in baseNmrAtoms]
                        markNmrAtoms(mainWindow=self.mainWindow, nmrAtoms=nmrAtoms, guiTarget=strips[0])

            if self.matchCheckBoxWidget.isChecked():
                self.findAndDisplayMatches(nmrResidue)

            # select the order for copying YAxis values
            if self.focusYAxis.isChecked():
                # align the target modules to the display module
                self._setDisplayPosWidth(matchDisplays, targetDisplays)
            else:
                # align the display modules to the target module
                self._setDisplayPosWidth(targetDisplays, matchDisplays)

        # update current to trigger other modules
        if self._activeLinkCheckbox and self._activeLinkCheckbox.isChecked():
            self.current.nmrChain = nmrResidue.nmrChain
        self.current.nmrResidue = nmrResidue

    @staticmethod
    def _setDisplayPosWidth(matchDisplays, targetDisplays):

        if matchDisplays and matchDisplays[0].strips:
            # get the current position/width of the first match display
            matchPos = matchDisplays[0].strips[0].getAxisPosition(axisIndex=1)
            matchWidth = matchDisplays[0].strips[0].getAxisWidth(axisIndex=1)

            for target in targetDisplays:
                # align to the match module
                for tgStrip in target.strips:
                    tgStrip.setAxisPosition(axisIndex=1, position=matchPos, rescale=False, update=False)
                    tgStrip.setAxisWidth(axisIndex=1, width=matchWidth, rescale=True,
                                         update=True)  #(tgStrip == target.strips[-1]))
                target.strips[0]._CcpnGLWidget.emitYAxisChanged(allStrips=True)

    def findAndDisplayMatches(self, nmrResidue):
        """Find and displays the matches to nmrResidue"""

        # If NmrResidue is a -1 offset NmrResidue, set queryShifts as value from self.interShifts dictionary
        # Set matchShifts as self.intraShifts
        if nmrResidue.relativeOffset == -1:
            if nmrResidue not in self.interShifts:
                queryShifts = []
            else:
                queryShifts = [shift for shift in self.interShifts[nmrResidue]
                               if (shift and not shift.isDeleted) and shift.nmrAtom and (
                                       shift.nmrAtom.name in self.nmrAtomsToMatch)]
            matchShifts = self.intraShifts

        # If NmrResidue is not an offset NmrResidue, set queryShifts as value from self.intraShifts dictionary
        # Set matchShifts as self.interShifts
        elif nmrResidue.relativeOffset == 0 or nmrResidue.relativeOffset is None:
            if nmrResidue not in self.intraShifts:
                queryShifts = []
            else:
                queryShifts = [shift for shift in self.intraShifts[nmrResidue]
                               if (shift and not shift.isDeleted) and shift.nmrAtom.name in self.nmrAtomsToMatch]
            matchShifts = self.interShifts

        # If NmrResidue has offset other than -1 or 0/None, tell user that we are not able to match
        else:
            getLogger().warning(
                    f"Assignment matching not supported for NmrResidue offset {nmrResidue.relativeOffset}. Matching display skipped")
            return

        assignMatrix = getNmrResidueMatches(queryShifts, matchShifts, 'averageQScore')

        # some old code which will match the query shifts against ALL shifts
        #    queryShifts = [shift for shift in self.allShifts[nmrResidue]
        #                   if shift.nmrAtom.isotopeCode == '13C']
        #    assignMatrix = getNmrResidueMatches(queryShifts, self.allShifts, 'averageQScore')

        if not assignMatrix.values():
            getLogger().info(f'No matches found for NmrResidue: {nmrResidue.pid} - no matching isotopeCodes?')
            return
        self._createMatchStrips(nmrResidue, assignMatrix)

    def _processDroppedNmrResidrueLabel(self, data, toLabel=None, plusChain=None):
        if toLabel and toLabel.obj:
            self._processDroppedNmrResidue(data, toLabel.obj, plusChain)

    def _processDroppedNmrResidue(self, data, nmrResidue, plusChain=None):
        """Process the dropped NmrResidue id"""

        droppedNmrResidue = None
        if DropBase.TEXT in data and len(data[DropBase.TEXT]) > 0:
            droppedNmrResidue = self.application.project.getByPid(data[DropBase.TEXT])
        if droppedNmrResidue is None:
            showWarning(str(self.windowTitle()), 'Backbone assignment: invalid dropped item')
            getLogger().warning('Backbone assignment: invalid "pid" of dropped item')
            return

        if not isinstance(droppedNmrResidue, NmrResidue):
            showWarning(str(self.windowTitle()), 'Backbone assignment: item is not an nmrResidue')
            getLogger().warning('Backbone assignment: item is not an nmrResidue')
            return

        getLogger().debug(f'nmrResidue: {nmrResidue}, droppedNmrResidue: {droppedNmrResidue}')
        if droppedNmrResidue == nmrResidue:
            return

        allNmrResidues = nmrResidue._getAllConnectedList()

        isPlus = data[STRIPLABEL_CONNECTDIR] if STRIPLABEL_CONNECTDIR in data else STRIPLABEL_CONNECTNONE
        if isPlus == STRIPCONNECT_RIGHT:
            data['shiftLeftMouse'] = False
        elif isPlus == STRIPCONNECT_LEFT:
            data['shiftLeftMouse'] = True
        else:
            data['shiftLeftMouse'] = None

        index = allNmrResidues.index(nmrResidue)
        lenNmr = len(allNmrResidues) - 1

        if data['shiftLeftMouse'] and not plusChain and index == 0:
            # okay to connect to left
            okay = True

        elif not data['shiftLeftMouse'] and plusChain and index == lenNmr:
            # okay to connect to right
            okay = True

        elif data['shiftLeftMouse'] and plusChain and index == lenNmr:
            # check connecting i-1 nmrResidue to the right
            yesNo = showYesNo(str(self.windowTitle()), "Trying to connect 'i-1' nmrResidue to end of chain.\n\n"
                                                       "Do you want to continue?")
            getLogger().warning("Trying to connect 'i-1' nmrResidue to end of chain")
            if not yesNo:
                return

            # force the connection to the start of chain
            data['shiftLeftMouse'] = False
            okay = True

        elif not data['shiftLeftMouse'] and not plusChain and index == 0:
            # check connecting i+1 nmrResidue to the left
            yesNo = showYesNo(str(self.windowTitle()), "Trying to connect 'i+1' nmrResidue to start of chain.\n\n"
                                                       "Do you want to continue?")
            getLogger().warning("Trying to connect 'i+1' nmrResidue to start of chain")
            if not yesNo:
                return

            # force the connection to the end of chain
            data['shiftLeftMouse'] = True
            okay = True

        else:
            # connecting to the middle of a stretch - may do disconnect later
            showWarning(str(self.windowTitle()), "Illegal connection, cannot connect to the middle of a chain")
            getLogger().warning("Illegal connection, cannot connect to the middle of a chain")
            return

        # silence the update of the nmrResidueTable as we will to an explicit update later
        # put in try/finally block because otherwise if exception thrown in the following code
        # (which can happen) then you no longer get updates of the NmrResidue table

        if data['shiftLeftMouse']:
            progressText = f"connecting  {droppedNmrResidue.pid}  >  {nmrResidue.pid}"
        else:
            progressText = f"connecting  {nmrResidue.pid}  <  {droppedNmrResidue.pid}"

        with undoBlockWithoutSideBar():
            with progressManager(self.mainWindow, progressText):

                try:
                    matchNmrResidue = None
                    try:  # display popup warning
                        if data['shiftLeftMouse']:
                            # leftShift drag; connect to previous

                            if not nmrResidue.residue and not droppedNmrResidue.residue:
                                nmrResidue.connectPrevious(droppedNmrResidue)

                            elif nmrResidue.residue and not droppedNmrResidue.residue:
                                # connected an unassigned nmrChain to the current assigned chain
                                if droppedNmrResidue.nmrChain.id == '@-':
                                    # assume that it is the only one
                                    droppedNmrResidue.nmrChain.assignSingleResidue(droppedNmrResidue,
                                                                                   nmrResidue.residue.previousResidue)
                                else:
                                    nRes = nmrResidue.residue
                                    for _ii in range(len(droppedNmrResidue.nmrChain.mainNmrResidues)):
                                        nRes = nRes.previousResidue
                                    droppedNmrResidue.nmrChain.assignConnectedResidues(nRes)

                            elif not nmrResidue.residue:
                                # connected an assigned chain to the current unassigned nmrChain

                                if nmrResidue.nmrChain.id == '@-':
                                    # assume that it is the only one
                                    nmrResidue.nmrChain.assignSingleResidue(nmrResidue,
                                                                            droppedNmrResidue.residue.nextResidue)
                                else:
                                    nmrResidue.nmrChain.assignConnectedResidues(droppedNmrResidue.residue.nextResidue)

                            matchNmrResidue = droppedNmrResidue.getOffsetNmrResidue(offset=-1)
                            if matchNmrResidue is None:
                                # Non -1 residue - stay with current
                                getLogger().info(f"NmrResidue {droppedNmrResidue} has no i-1 residue to display")
                                matchNmrResidue = nmrResidue

                        else:
                            if not nmrResidue.residue and not droppedNmrResidue.residue:
                                nmrResidue.connectNext(droppedNmrResidue)

                            elif nmrResidue.residue and not droppedNmrResidue.residue:
                                # connected an unassigned nmrChain to the current assigned chain
                                if droppedNmrResidue.nmrChain.id == '@-':
                                    # assume that it is the only one
                                    droppedNmrResidue.nmrChain.assignSingleResidue(droppedNmrResidue,
                                                                                   nmrResidue.residue.nextResidue)
                                else:
                                    droppedNmrResidue.nmrChain.assignConnectedResidues(nmrResidue.residue.nextResidue)

                            elif not nmrResidue.residue:
                                # connected an assigned chain to the current unassigned nmrChain

                                if nmrResidue.nmrChain.id == '@-':
                                    # assume that it is the only one
                                    nmrResidue.nmrChain.assignSingleResidue(nmrResidue,
                                                                            droppedNmrResidue.residue.previousResidue)
                                else:
                                    dropRes = droppedNmrResidue.residue
                                    for _ in range(len(nmrResidue.nmrChain.mainNmrResidues)):
                                        dropRes = dropRes.previousResidue
                                    nmrResidue.nmrChain.assignConnectedResidues(dropRes)

                            matchNmrResidue = droppedNmrResidue

                    except Exception as es:
                        showWarning('Connect NmrResidue', str(es))
                    finally:
                        if matchNmrResidue:
                            self.navigateToNmrResidue(matchNmrResidue)

                            # # update the NmrResidueTable
                            # getLogger().info('>>>DISPLAYTABLE', droppedNmrResidue.nmrChain, self.project.nmrChains)
                            # self.nmrResidueTable.displayTableForNmrChain(droppedNmrResidue.nmrChain)

                            # from ccpn.ui.gui.lib.OpenGL.CcpnOpenGL import GLNotifier
                            #
                            # GLSignals = GLNotifier(parent=self)
                            # GLSignals.emitEvent(triggers=[GLNotifier.GLMARKS])

                except Exception as es:
                    getLogger().warning(str(es))

        # # update the NmrResidueTable - outside of the undoBlock for notifiers to catch up
        # print(f'   dropped {droppedNmrResidue.nmrChain.pid}')
        # self.tableFrame._modulePulldown.select(droppedNmrResidue.nmrChain.pid)
        # self._tableWidget._update(useSelected=True)  # droppedNmrResidue.nmrChain)

        from ccpn.ui.gui.lib.OpenGL.CcpnOpenGL import GLNotifier

        GLSignals = GLNotifier(parent=self)
        GLSignals.emitEvent(triggers=[GLNotifier.GLMARKS])

    # def _centreStripForNmrResidue(self, nmrResidue, strip):
    #     """
    #     Centre y-axis of strip based on chemical shifts of from NmrResidue.nmrAtoms
    #     """
    #     if not nmrResidue:
    #         getLogger().warning('No NmrResidue specified')
    #         return
    #
    #     if not strip:
    #         getLogger().warning('No Strip specified')
    #         return
    #
    #     yShifts = matchAxesAndNmrAtoms(strip, nmrResidue.nmrAtoms)[strip.axisOrder[1]]
    #     yShiftValues = [x.value for x in yShifts]
    #     if yShiftValues:
    #         yPosition = (max(yShiftValues) + min(yShiftValues)) / 2
    #         yWidth = max(yShiftValues) - min(yShiftValues) + 10
    #         strip.orderedAxes[1].position = yPosition
    #         if strip._CcpnGLWidget.aspectRatioMode == 0:
    #             strip.orderedAxes[1].width = yWidth
    #
    #         try:
    #             axisCode = strip.axisCodes[1]
    #             if strip._CcpnGLWidget.aspectRatioMode == 0:
    #                 strip._CcpnGLWidget.setAxisPosition(axisCode=axisCode, position=yPosition, update=False)
    #                 strip._CcpnGLWidget.setAxisWidth(axisCode=axisCode, width=yWidth, update=False)
    #                 strip._CcpnGLWidget._rescaleAllAxis()
    #             else:
    #                 strip._CcpnGLWidget.setAxisPosition(axisCode=axisCode, position=yPosition, update=True)
    #
    #         except Exception as es:
    #             getLogger().debugGL('OpenGL widget not instantiated')

    @staticmethod
    def _centreCcpnStripsForNmrResidue(nmrResidue, strips):
        """
        Centre y-axis of strip based on chemical shifts of from NmrResidue.nmrAtoms
        """
        if not nmrResidue:
            getLogger().warning('No NmrResidue specified')
            return

        if not strips:
            getLogger().warning('No Strip specified')
            return

        yShifts = matchAxesAndNmrAtoms(strips[0], nmrResidue.nmrAtoms)[strips[0].axisOrder[1]]

        if yShiftValues := [x.value for x in yShifts]:
            _minPpmWidths = {'H': 0.5, 'C': 8.0, 'N': 2.0}  # based on standard ratios

            yPosition = (max(yShiftValues) + min(yShiftValues)) / 2
            yWidth = max(yShiftValues) - min(yShiftValues)

            axisCode = strips[0].axisCodes[1]
            yHeight = strips[0]._CcpnGLWidget.mainViewHeight() or 1

            minPpm = 1.0 if axisCode[0] not in _minPpmWidths else _minPpmWidths[axisCode[0]]

            # add increase of EXTRAWIDTH pixels
            dY = max(yWidth, minPpm) / max(1, (yHeight - EXTRAWIDTH))

            # add offset for the top, and extra height
            yPos = yPosition - (EXTRAOFFSET - (EXTRAWIDTH / 2)) * dY
            yW = max(yWidth, minPpm) + EXTRAWIDTH * dY

            # this should rescale all in spectrumDisplay
            strips[0].setAxisPosition(axisIndex=1, position=yPos, rescale=False, update=False)
            strips[0].setAxisWidth(axisIndex=1, width=yW, rescale=True, update=True)

    def _setupShiftDicts(self, *args):
        """
        Creates three ordered dictionaries containing a) intra-residue, b) -1 offset (inter) and c) all
        shifts for all NmrResidues in the project.
        """
        self.intraShifts = OrderedDict()
        self.interShifts = OrderedDict()
        self.allShifts = OrderedDict()

        if chemicalShiftList := self.application.project.getByPid(self.shiftListWidget.pulldownList.currentText()):
            for nmrResidue in self.application.project.nmrResidues:
                nmrAtoms = list(nmrResidue.nmrAtoms)
                shifts = [chemicalShiftList.getChemicalShift(atom) for atom in nmrAtoms
                          if chemicalShiftList.getChemicalShift(atom) is not None]
                if nmrResidue.relativeOffset == -1:
                    self.interShifts[nmrResidue] = shifts
                elif nmrResidue.relativeOffset == 0 or nmrResidue.relativeOffset is None:
                    self.intraShifts[nmrResidue] = shifts
                self.allShifts[nmrResidue] = shifts

    def _setNmrAtomsToMatch(self):
        self.nmrAtomsToMatch = []
        if self.matchCA.isChecked():
            self.nmrAtomsToMatch.append('CA')
        if self.matchCB.isChecked():
            self.nmrAtomsToMatch.append('CB')
        if self.matchC.isChecked():
            self.nmrAtomsToMatch.append('C')

    def _createMatchStrips(self, nmrResidue, assignMatrix: typing.Tuple[
        typing.Dict[NmrResidue, typing.List[ChemicalShift]], typing.List[float]]):
        """
        Creates strips in match module corresponding to the best assignment possibilities
        in the assignMatrix.
        """
        if not assignMatrix:
            getLogger().warning('No assignment matrix specified')
            return

        # Assignment score has the format {score: nmrResidue} where score is a float
        # assignMatrix[0] is a dict {score: nmrResidue} assignMatrix[1] is a concurrent list of scores
        # numberOfMatches = int(self.numberOfMatchesWidget.getText())
        assignmentScores = ([-1] if self.showSearchInMatch.isChecked() else []) + sorted(list(assignMatrix.keys()))[
                                                                                  :MAXMATCHES]
        scoreAssignment = [''] if self.showSearchInMatch.isChecked() else []
        scoreLabelling = [''] if self.showSearchInMatch.isChecked() else []
        nmrAtomPairs = [(None, None)] if self.showSearchInMatch.isChecked() else []

        matchDirection = 0
        for assignmentScore in assignmentScores:
            if assignmentScore == -1:
                continue

            matchResidue = assignMatrix[assignmentScore]
            if matchResidue.sequenceCode.endswith('-1'):
                iNmrResidue = matchResidue.mainNmrResidue

                # this is where the nmrResidue can be dropped in the existing nmrChain
                scoreLabelling.append('[ i+1 ]')
                matchDirection = 1
            else:
                iNmrResidue = matchResidue
                scoreLabelling.append('[ i-1 ]')
                matchDirection = -1

            scoreAssignment.append('[ %i' % int(100 - min(1000 * assignmentScore, 100)) + '% ]')

            nmrAtomPairs.append((iNmrResidue.fetchNmrAtom(name='N', isotopeCode='N'),
                                 iNmrResidue.fetchNmrAtom(name='H', isotopeCode='H')))

        numberOfMatches = 1 if self.showSearchInMatch.isChecked() else 0
        if matchDirection == 1:
            # numberOfMatches = int(self.numberOfPlusMatchesWidget.getText())
            numberOfMatches += self.numberOfPlusMatchesWidget.getValue()
        else:
            # numberOfMatches = int(self.numberOfMinusMatchesWidget.getText())
            numberOfMatches += self.numberOfMinusMatchesWidget.getValue()

        nmrAtomPairs = nmrAtomPairs[:numberOfMatches]
        scoreAssignment = scoreAssignment[:numberOfMatches]
        scoreLabelling = scoreLabelling[:numberOfMatches]

        for module in self._getMatchDisplays():

            # skip of the module if not defined - possibly in the case that spectrumDisplays have been closed
            if not module:
                continue

            # set the first strip to the search-strip, not nice here
            if self.showSearchInMatch.isChecked():
                # add a mask to ignore the position for the Y-axis
                axisMask = [True] * len(module.axisCodes)
                axisMask[1] = False

                newWidths = []  #_getCurrentZoomRatio(display.strips[0].viewBox.viewRange())
                navigateToNmrResidueInDisplay(nmrResidue, module, stripIndex=0,
                                              widths=newWidths,
                                              showSequentialResidues=False,
                                              markPositions=False,
                                              showDropHeaders=True,
                                              axisMask=axisMask,
                                              keepExistingStrips=True,
                                              )
                if (strip := module.strips[0]) is not None:
                    strip.header.handle = STRIPBACKBONE
                    strip.header.headerVisible = True

            makeStripPlot(module, nmrAtomPairs)

            # make the matching strips
            for ii, strip in enumerate(module.strips):
                if None in nmrAtomPairs[ii]:
                    # skip the dummy atom-pairs
                    continue

                nmrResiduePid = nmrAtomPairs[ii][0].nmrResidue.pid

                # strip.setStripLabelText(nmrResiduePid)
                # strip.showStripLabel()
                # strip.setStripLabelisPlus(True if scoreLabelling[ii].startswith('i+1') else False)
                # strip.setStripResidueIdText(scoreLabelling[ii])
                # strip.showStripResidueId()
                # strip.setStripResidueDirText(scoreAssignment[ii])
                # strip.showStripResidueDir()

                strip.header.reset()
                strip.header.setLabelText(position='l', text=scoreLabelling[ii])
                strip.header.setLabelText(position='c', text=nmrResiduePid)

                # TODO:ED need to improve this
                # strip.header.setLabelConnectDir(position='c', connectDir=STRIPCONNECT_LEFT if scoreLabelling[ii].startswith('i-1') else STRIPCONNECT_RIGHT)
                strip.header.setLabelConnectDir(position='c', connectDir=STRIPCONNECT_LEFT if 'i-1' in scoreLabelling[
                    ii] else STRIPCONNECT_RIGHT)
                strip.header.setLabelText(position='r', text=scoreAssignment[ii])

                # disable dropping onto these labels
                strip.header.setLabelObject(position='l', obj=None)
                strip.header.setLabelObject(position='c', obj=None)
                strip.header.setLabelObject(position='r', obj=None)

                strip.header.handle = STRIPBACKBONE
                strip.header.headerVisible = True

            # self._centreStripForNmrResidue(assignMatrix[assignmentScores[0]], module.strips[0])
            self._centreCcpnStripsForNmrResidue(
                    assignMatrix[assignmentScores[1 if self.showSearchInMatch.isChecked() else 0]], module.strips)
            module.setColumnStretches(stretchValue=True)

            # this forces a refresh/rescale of all strips in the spectrumDisplay
            module.strips[0]._CcpnGLWidget.emitYAxisChanged(allStrips=True)

    def _closeModule(self):
        """
        Re-implementation of the closeModule method of the CcpnModule class required
        """
        for display in self._getDisplays() + self._getMatchDisplays():
            if display:
                display.hideAllStripHeaders(handle=STRIPBACKBONE)
        super()._closeModule()


def nmrAtomsFromResidue(nmrResidue):
    """
    Retrieve a list of nmrAtoms from nmrResidue
    """
    # nmrResidue = nmrResidue.mainNmrResidue
    nmrResidues = []
    if previousNmrResidue := nmrResidue.previousNmrResidue:
        nmrResidues.append(previousNmrResidue)
    nmrResidues.append(nmrResidue)
    if nextNmrResidue := nmrResidue.nextNmrResidue:
        nmrResidues.append(nextNmrResidue)

    nmrAtoms = []
    for nr in nmrResidues:
        nmrAtoms.extend(nr.nmrAtoms)

    return nmrAtoms


def nmrAtomsFromOffsets(nmrResidue):
    """
    Retrieve a list of nmrAtoms from nmrResidue
    """
    # nmrResidue = nmrResidue.mainNmrResidue
    nmrResidues = [nmrResidue]
    if nmrResidue.offsetNmrResidues:
        nmrResidues.extend(nmrResidue.offsetNmrResidues)

    nmrAtoms = []
    for nr in nmrResidues:
        nmrAtoms.extend(nr.nmrAtoms)

    return nmrAtoms


#=====  Just some code to 'save' =====
# def hasNmrResidue(nmrChain, residueCode):
#     "Simple function to check if sequenCode is found within the nmrResidues of nmrChain"
#     resCodes = [res.sequenceCode for res in nmrChain.nmrResidues]
#     return (residueCode in resCodes)
#
#
# def endOfchain(nmrResidue):
#     # changes to end of connected chain; not a good idea
#     if nmrResidue.nmrChain.isConnected:
#         if nmrResidue.sequenceCode.endswith('-1'):
#             nmrResidue = nmrResidue.nmrChain.mainNmrResidues[0].getOffsetNmrResidue(-1)
#         else:
#             nmrResidue = nmrResidue.nmrChain.mainNmrResidues[-1]
#     return nmrResidue
#
#
# def getPids(fromObject, attributeName):
#     "Get a list of pids fromObject.attributeName or None on error"
#     if not hasattr(fromObject, attributeName): return None
#     return [obj.pid for obj in getattr(fromObject, attributeName)]
#
#
#===== end code save =====


if __name__ == '__main__':
    from ccpn.ui.gui.widgets.Application import TestApplication


    app = TestApplication()

    popup = BackboneAssignmentModule()

    popup.show()
    popup.raise_()
    app.start()
