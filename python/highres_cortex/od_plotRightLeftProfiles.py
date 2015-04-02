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
#python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_plotRightLeftProfiles.py -p ad140157 -d /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/ad140157/

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
    #realSide = 'L'

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
    plt.plot(coordL, valueL, '.', c = 'b', label = 'L')
    plt.title('Profile in all ROIs')   # subplot 211 title
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')
    plt.plot(coordR, valueR, '.', c = 'r', label = 'R')
    plt.legend(loc='upper right', numpoints = 1)
    plt.savefig(directory + '%s_LvsR_2nobiasT2.png' %(realPatientID))    
    
    # TODO: delete later if no need: save profiles also to the outer folder!
    #plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/' + '%s_LvsR_2nobiasT2.png' %(realPatientID))    
    plt.clf()
    plt.close()
    
    # now plot L vs R in various ROIs
    # af140169_R_profiles2_ROI_21.txt, af140169_R_profiles2_ROI_11.txt and the same with L
    pathToROIsProfL = directory + '%s_L_profiles2_ROI_[0-9]*.txt' %(realPatientID)
    #print 'pathToROIsProfL'
    #print pathToROIsProfL

    pathToROIsProfR = directory + '%s_R_profiles2_ROI_[0-9]*.txt' %(realPatientID)
    
    # check if both these profiles exist. and how many are there
    profROIsL = glob.glob(pathToROIsProfL)
    profROIsR = glob.glob(pathToROIsProfR)
    
    #print 'profROIsL'
    #print profROIsL
    
    for i in range(len(profROIsL)):
        # read in the respective files
        numbROIsL, coordROIsL, valueROIsL = np.loadtxt(profROIsL[i], skiprows = 1, unpack = True)
        numbROIsR, coordROIsR, valueROIsR = np.loadtxt(profROIsR[i], skiprows = 1, unpack = True)
        
        # get the ROI ID
        iD = profROIsL[i].split('_L_profiles2_ROI_')[1]
        iD = iD.split('.txt')[0]

        # plot the data
        plt.plot(coordROIsL, valueROIsL, '.', c = 'b', label = 'L')
        plt.title('Profile in ROI %s ' %(iD))   # subplot 211 title
        plt.xlabel('Cortical depth')
        plt.ylabel('T2-nobias intensity')
        plt.plot(coordROIsR, valueROIsR, '.', c = 'r', label = 'R')
        plt.legend(loc='upper right', numpoints = 1)
        plt.savefig(directory + '%s_LvsR_2nobiasT2_ROI_%s.png' %(realPatientID, iD))   
        
        # TODO: delete later if no need: save profiles also to the outer folder!
        #plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/' + '%s_LvsR_2nobiasT2_ROI_%s.png' %(realPatientID, iD))    
        plt.clf()
        plt.close() 
      
            
       
    # plot these plots into 1 image
    fig = plt.figure(figsize=(21, 6)) #, dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(131)
    ax1.plot(coordL, valueL, '.', c = 'b', label = 'L')
    plt.title('Profile in all ROIs')   # subplot 211 title
    #ax1.xlabel('Cortical depth')
    #ax1.ylabel('T2-nobias intensity')
    ax1.plot(coordR, valueR, '.', c = 'r', label = 'R')
    plt.legend(loc='upper right', numpoints = 1)
   
    for i in range(len(profROIsL)):
        # read in the respective files
        numbROIsL, coordROIsL, valueROIsL = np.loadtxt(profROIsL[i], skiprows = 1, unpack = True)
        numbROIsR, coordROIsR, valueROIsR = np.loadtxt(profROIsR[i], skiprows = 1, unpack = True)
        
        # get the ROI ID
        iD = profROIsL[i].split('_L_profiles2_ROI_')[1]
        iD = iD.split('.txt')[0]

        ax2 = fig.add_subplot(1,3, 2 + i)
        ax2.plot(coordROIsL, valueROIsL, '.', c = 'b', label = 'L')
        plt.title('Profile in ROI %s' %(iD))   # subplot 211 title
        #ax2.xlabel('Cortical depth')
        #ax2.ylabel('T2-nobias intensity')
        ax2.plot(coordROIsR, valueROIsR, '.', c = 'r', label = 'R')
        plt.legend(loc='upper right', numpoints = 1)

    #plt.show()
    plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/' + '%s_LvsR_allVsROIs.png' %(realPatientID), bbox_inches='tight')    
    plt.clf()
    plt.close()
    # plot it for sufficiently large testBatchColumnsExtrProfiles

      
  
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    