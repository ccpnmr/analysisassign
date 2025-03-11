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
__modifiedBy__ = "$modifiedBy: Luca Mureddu $"
__dateModified__ = "$dateModified: 2025-03-11 19:28:27 +0000 (Tue, March 11, 2025) $"
__version__ = "$Revision: 3.3.1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2016-05-23 10:02:47 +0100 (Mon, 23 May 2016) $"
#=========================================================================================
# Start of code
#=========================================================================================

import math


def qScore(value1: float, value2: float):
    if value1 + value2 == 0: return 0 # otherwise is a ZeroDivision error.
    return math.sqrt(((value1 - value2) ** 2) / ((value1 + value2) ** 2))



def averageQScore(valueLists):
    if len(valueLists[0]) == 0: return 0 # otherwise is a ZeroDivision error.
    score = sum([qScore(valueList[0], valueList[1]) for valueList in valueLists]) / len(valueLists[0])
    return score



def euclidean(valueList):
    score = sum([(scoringValue[0] - scoringValue[1]) ** 2 for scoringValue in valueList])
    return math.sqrt(score)


functionDict = {
    'averageQScore': averageQScore,
    'euclidean': euclidean,
    }


def getNmrResidueMatches(queryShifts, matchNmrResiduesDict, scoringMethod, isotopeCode='13C'):
    """Method to calculate scores for matching nmrAtoms/isotopeCodes from chemicalShifts
    """
    scoringMatrix = {}

    for res, mShifts in matchNmrResiduesDict.items():

        if res.isDeleted:
            continue

        # get the matching isotopeCodes
        mShifts2 = [shift for shift in mShifts
                    if (shift and shift.nmrAtom) and shift.nmrAtom.isotopeCode == isotopeCode]
        qShifts2 = [shift for shift in queryShifts
                    if (shift and shift.nmrAtom) and shift.nmrAtom.isotopeCode == isotopeCode]

        # get the matching atomNames
        scoringValues = [(mShift.value, qShift.value)
                         for mShift in mShifts2
                         for qShift in qShifts2
                         if (mShift != qShift) and (mShift.value is not None) and (qShift.value is not None) and (mShift.nmrAtom.name == qShift.nmrAtom.name)]

        if scoringValues and len(scoringValues) == len(qShifts2):  # or should this be len(queryShifts)?
            score = functionDict[scoringMethod](scoringValues)
            scoringMatrix[score] = res

    return scoringMatrix
