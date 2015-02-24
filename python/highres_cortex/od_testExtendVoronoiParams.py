#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright CEA (2014).
# Copyright Université Paris XI (2014).
#
# Contributor: Olga Domanova <olga.domanova@cea.fr>.
#
# This file is part of highres-cortex, a collection of software designed
# to process high-resolution magnetic resonance images of the cerebral
# cortex.
#
# This software is governed by the CeCILL licence under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/or redistribute the software under the terms of the CeCILL
# licence as circulated by CEA, CNRS and INRIA at the following URL:
# <http://www.cecill.info/>.
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the licence, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of scientific
# software, that may mean that it is complicated to manipulate, and that
# also therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL licence and that you accept its terms.
#
#
# this scripts compares results of heat calculation for several cases:
# the original is the heat calculated in the whole cortex
# these temperature values are compared to the heat calculated in some "cut out" region
# e.g. in Heschl gyrus and in Planum Temporale. 
# These cut-out regions were extended. The question is to find the required "extension size"
# to keep values in the ROI similar to those in the original volume.
# The maximal difference should not be > 1%.

# this function will take the original volume, the mask of ROIs, and a list of volumes that were
# created at various extensions.
# the maximal difference is calculated and reported.

# example how to run this file:
# python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_testExtendVoronoiParams.py -p lg140146 -d /neurospin/lnao/dysbrain/optimiseParmaters_ExtendVoronoi/ -s R

import random
from soma import aims, aimsalgo
import subprocess
from optparse import OptionParser
from scipy.stats import mode
import sys, glob, os, os.path, subprocess, sys, time, timeit
import numpy as np
import highres_cortex.od_cutOutRois
from soma.aims import volumetools


def findDiffInROI(vol1, vol2, mask):
    """
    this function takes original volume, the volume to be compared
    and the mask of ROIs
    It finds the maximal difference between the volumes in the ROI
    """
    arr1 = np.array(vol1, copy = False)
    arr2 = np.array(vol2, copy = False)
    arrMask = np.array(mask, copy = False)

    arr1cut = arr1[arrMask != 0]
    arr2cut = arr2[arrMask != 0]
    
    diff = np.max(np.abs(arr1cut - arr2cut))
    return(diff)
    
    
    
if __name__ == '__main__':    
    realPatientID = None
    directory = None
    realSide = 'L'

    parser = OptionParser('Find differences between the volumes in a ROI')    
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    parser.add_option('-d', dest='directory', help='directory')
    parser.add_option('-s', dest='realSide', help='Hemisphere to be processed: L or R. L is default')   
    options, args = parser.parse_args(sys.argv)
    print options
    print args   
    
    if options.directory is None:
        print >> sys.stderr, 'New: exit. no directory given'
        sys.exit(1)
    else:
        directory = options.directory           

    if options.realPatientID is None:
        print >> sys.stderr, 'New: exit. no patient ID given'
        sys.exit(1)
    else:
        realPatientID = options.realPatientID     
  
    if options.realSide is not None:
        realSide = options.realSide

    pathToFullCortex = directory + '%s/%s_T1inT2_ColumnsNew/' %(realPatientID, realPatientID)
    pathToCutCortex = directory + '%s/%s_T1inT2_ColumnsCutNew' %(realPatientID, realPatientID)     # add 3It, 5It, ..
    pathToROIMask = ''
    f = open(directory + '%s/%s_extendVoronoiParams.txt' %(realPatientID, realPatientID), "w")

    volFull = aims.read(pathToFullCortex + 'heat/heat_%s_%s_noSulci.nii.gz' %(realPatientID, realSide)) 
    f.write('Differences to the file ' + pathToFullCortex + 'heat/heat_%s_%s_noSulci.nii.gz' %(realPatientID, realSide) + '\n')
    mask = aims.read(directory + '%s/%s_T1inT2_ColumnsCutNew20It/voronoiCorr_%s_%s_cut_noSulci.nii.gz' %(realPatientID, realPatientID, realPatientID, realSide))

    # get a list of cut cortex volumes
    pathToVolsCut = glob.glob(pathToCutCortex + '[0-9]*It/heat/heat_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide))

    for i in range(len(pathToVolsCut)):
        volCut = aims.read(pathToVolsCut[i])
        d = findDiffInROI(volFull, volCut, mask)
        print 'for the file ', pathToVolsCut[i], ' the maxDiff to fullCortexHeat = ', str(d)
        f.write('for the file ' + pathToVolsCut[i] + ' the maxDiff to fullCortexHeat = ' + str(d) + '\n')
        
    f.close()    

    
    
    
    
    
    