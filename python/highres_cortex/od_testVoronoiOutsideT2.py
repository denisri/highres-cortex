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
#
# this function will verify whether the Voronoi segmentation of the ROIs is 
# outside of the T2 acquired images# further, it will notify if voronoi touches the lower or the upper 4 layers of T2

# example how to run this file:
# with cortical columns
#python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_extractProfiles.py -p ad140157 -s R -d /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/ad140157/ -c 3

# or without cortical columns
#python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_extractProfiles.py -p ad140157 -s R -d /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/ad140157/

import random
from soma import aims, aimsalgo
import subprocess
from optparse import OptionParser
from scipy.stats import mode
import sys, glob, os, os.path, subprocess, sys, time, timeit
import numpy as np
import highres_cortex.od_cutOutRois
from soma.aims import volumetools
import matplotlib.pyplot as plt
import math
from sklearn import preprocessing

if __name__ == '__main__':
    
    realPatientID = None
    directory = None
    realSide = 'L'
    columnDiameter = None
    workOnLaptop = False
    reExtractProfilesCutBorders = False
    pathToNobiasT2 = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S4/'
    pathToNobiasT2_new = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S16/'
    minColumnSize = 10 # depends on the diameter
    heights = []
    minColumnSizes = []
    ignoreCommunColums = True
    parser = OptionParser('Extract profiles from T2 nobias data using cortex-density-coordinates in ROIs')    
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    parser.add_option('-s', dest='realSide', help='Hemisphere to be processed: L or R. L is default')   
    parser.add_option('-d', dest='directory', help='directory')
    parser.add_option('-c', dest='columnDiameter', help='columnDiameter to work with')
    parser.add_option('-g', dest='ignoreCommunColums', action = 'store_false', help='Select if want to INCLUDE into calculation cortical colums found in several ROIs (Excluding them is default')
    parser.add_option('-l', dest='workOnLaptop', action = 'store_true', help='Select if working on laptop (neurospin DB location is different. False is default') 
    options, args = parser.parse_args(sys.argv)
    print options
    print args   
    
    if options.directory is None:
        print >> sys.stderr, 'New: exit. no directory given'
        sys.exit(0)
    else:
        directory = options.directory           

    if options.realPatientID is None:
        print >> sys.stderr, 'New: exit. no patient ID given'
        sys.exit(0)
    else:
        realPatientID = options.realPatientID     

    if options.realSide is not None:
        realSide = options.realSide
        
    if options.columnDiameter is not None:
        columnDiameter = int(options.columnDiameter)
        
    if options.ignoreCommunColums is not None:
        ignoreCommunColums = options.ignoreCommunColums

    if options.workOnLaptop is not None:
            workOnLaptop = options.workOnLaptop      
            # if true, then processes are run on the laptop. Change locations of neurospin DBs
            #pathToNobiasT2 = pathToNobiasT2.replace('/neurospin/lnao/', '/nfs/neurospin/lnao/')
            #pathToNobiasT2_new = pathToNobiasT2_new.replace('/neurospin/lnao/', '/nfs/neurospin/lnao/')  
            
            # saving the data locally on the laptop
            pathToNobiasT2 = pathToNobiasT2.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')
            pathToNobiasT2_new = pathToNobiasT2_new.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')           
        
    pathToCoord = directory + '%s_T1inT2_ColumnsCutNew20It/isovolume/' %(realPatientID)
    pathToMask = directory + '%s_T1inT2_ColumnsCutNew20It/' %(realPatientID)
    volsCoord = glob.glob(pathToCoord + 'pial-volume-fraction_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide))
    volValue = aims.read(pathToNobiasT2 + '%s_NewNobiasT2_cropped.nii.gz' %(realPatientID))
    volValue2 = aims.read(pathToNobiasT2_new + '%s_NewT2_cropped.nii.gz' %(realPatientID))
    volsMask = glob.glob(pathToMask + 'voronoiCorr_%s_%s_cut_noSulci.nii.gz' %(realPatientID, realSide))
    print 'volsCoord = ' , pathToCoord + 'pial-volume-fraction_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)
    print 'volsMask = ', pathToMask + 'voronoiCorr_%s_%s_cut_noSulci.nii.gz' %(realPatientID, realSide)
    print 'volValue2 = ', pathToNobiasT2_new + '%s_NewT2_cropped.nii.gz' %(realPatientID)     
    volClassifNoBorders = directory + '%s_T1inT2_ColumnsCutNew20It/GWsegm_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realPatientID, realSide)
















