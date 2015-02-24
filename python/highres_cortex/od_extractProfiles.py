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
# this function will take two volumes: then one with "x coordinates" of voxels
# and the other with values
# the values are colelcted and stored with their respective coordinates
# plots are created
# profiles can be extracted in certain ROIs defined by a mask


# example how to run this file:
#python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_testExtendVoronoiParams.py -p he140338 -d /neurospin/lnao/dysbrain/optimiseParmaters_ExtendVoronoi/
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


def extractProfiles(volCoord, volValue, volMask = None):
    arrCoord = np.asarray(volCoord)
    arrValue = np.asarray(volValue)

    # apply the mask if it was given
    if volMask is not None:
        mask = np.asarray(volMask)        
        # get these ROIs
        arrCoord1 = arrCoord[mask != 0]
        arrValue1 = arrValue[mask != 0]
        
        coords = arrCoord1[arrCoord1 != 0]
        values = arrValue1[arrCoord1 != 0]
    else :
        coords = arrCoord[arrCoord != 0]
        values = arrValue[arrCoord != 0]
        
    # 2 arrays of coordinates and values. Plot them
    res = []
    res.append(coords)
    res.append(values)
    return(res)       
    
    
if __name__ == '__main__':
    
    realPatientID = None
    directory = None
    realSide = 'L'

    parser = OptionParser('Extract profiles from T2 nobias data using cortex-density-coordinates in ROIs')    
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    parser.add_option('-s', dest='realSide', help='Hemisphere to be processed: L or R. L is default')   
    parser.add_option('-d', dest='directory', help='directory')
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

    pathToCoord = directory + '%s/%s_T1inT2_ColumnsCutNew20It/isovolume/' %(realPatientID, realPatientID)
    pathToNobiasT2 = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S4/'
    pathToNobiasT2_new = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S16/'
    pathToMask = directory + '%s/%s_T1inT2_ColumnsCutNew20It/' %(realPatientID, realPatientID)
    
    volCoord = aims.read(pathToCoord + 'pial-volume-fraction_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide))
    volValue = aims.read(pathToNobiasT2 + '%s_NewNobiasT2_cropped.nii.gz' %(realPatientID))
    volValue2 = aims.read(pathToNobiasT2_new + '%s_NewT2_cropped.nii.gz' %(realPatientID))
    volMask = aims.read(pathToMask + 'voronoiCorr_%s_%s_cut_noSulci.nii.gz' %(realPatientID, realSide))

    result = extractProfiles(volCoord, volValue, volMask)
    coordinates = result[0]
    intensities = result[1]

    # plot the data
    #plt.plot(coordinates, 'bo')
    #plt.plot(intensities, 'bo')
    plt.plot(coordinates, intensities, '.', c = 'b')
    plt.title('Profile in ROI')   # subplot 211 title
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')
    plt.savefig(directory + '%s/%s_%s_It20_nobiasT2vsCorticalDepthROI.png' %(realPatientID, realPatientID, realSide))

    
    # repeat for the NEW nobias images!
    result2 = extractProfiles(volCoord, volValue2, volMask)
    coordinates2 = result2[0]
    intensities2 = result2[1]

    # plot the data
    #plt.plot(coordinates, 'bo')
    #plt.plot(intensities, 'bo')
    plt.plot(coordinates2, intensities2, '.', c = 'r')
    plt.title('Profile in ROI')   # subplot 211 title
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')
    plt.savefig(directory + '%s/%s_%s_It20_2nobiasT2vsCorticalDepthROI.png' %(realPatientID, realPatientID, realSide))
    
    plt.clf()
    plt.close()
    
    plt.plot(coordinates2, intensities2, '.', c = 'r')
    plt.title('Profile in ROI')   # subplot 211 title
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')
    plt.savefig(directory + '%s/%s_%s_It20_newNobiasT2vsCorticalDepthROI.png' %(realPatientID, realPatientID, realSide))

 
    
    
    
    