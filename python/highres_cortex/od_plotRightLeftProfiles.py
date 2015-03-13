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
# this function will plot already extracted profiles from 2 hemispheres on one plot


# example how to run this file:
#python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_extractProfiles.py -p ac140159 -s L -d /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/
#python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_extractProfiles.py -p ad140157 -s L -d /neurospin/lnao/dysbrain/testNewLittleRegion/ad140157_2/

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
    
if __name__ == '__main__':
    
    realPatientID = None
    directory = None
    realSide = 'L'

    parser = OptionParser('Extract profiles from T2 nobias data using cortex-density-coordinates in ROIs')    
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    #parser.add_option('-s', dest='realSide', help='Hemisphere to be processed: L or R. L is default')   
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

    pathToProfL = directory + '%s_L_profiles2.txt' %(realPatientID)
    pathToProfR = directory + '%s_R_profiles2.txt' %(realPatientID)
    
    # check if both these profiles exist
    profL = glob.glob(pathToProfL)
    profR = glob.glob(pathToProfR)
    
    if len(profL) != 1 or len(profR) != 1:
        # abort the calculation, as too many or not a single texture file was found
        f = open(directory + '%s_compareProfilesStat.txt' %(realPatientID), "w")
        print 'abort the calculation, as too many or not a single profL or R file was found'
        f.write('abort the calculation, as ' + str(len(profL)) + ' profL and ' + str(len(profR)) + ' profR profile files were found' + '\n')
        f.close()
        sys.exit(0)
 
    
    
    numbL, coordL, valueL = np.loadtxt(pathToProfL, skiprows = 1, unpack = True)
    numbR, coordR, valueR = np.loadtxt(pathToProfR, skiprows = 1, unpack = True)
       
    # plot the data
    plt.plot(coordL, valueL, '.', c = 'b')
    plt.title('Profile in ROI')   # subplot 211 title
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')

    plt.plot(coordR, valueR, '.', c = 'r')
    plt.title('Profile in ROI')   # subplot 211 title
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')
    plt.savefig(directory + '%s_%s_LvsR_2nobiasT2vsCorticalDepthROI.png' %(realPatientID, realSide))
    
    plt.clf()
    plt.close()
      
    ### now plot and save the data for individual ROIs
    #iDs = result2[2]
    #listOfCoords = result2[3]
    #listOfValues = result2[4]
    #for i in range(len(iDs)):
        #print 'i = ', i, ' work with id', iDs[i]
        #currCoords = listOfCoords[i]
        #currValues = listOfValues[i]
        #plt.plot(currCoords, currValues, '.', c = 'b')
        #plt.title('Profile in ROI')   # subplot 211 title
        #plt.xlabel('Cortical depth')
        #plt.ylabel('T2-nobias intensity')
        #plt.savefig(directory + '%s_%s_It20_nobiasT2vsCorticalDepth_ROI_' %(realPatientID, realSide) + str(iDs[i]) + '.png')
        #plt.clf()
        #plt.close()
        
    

        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    