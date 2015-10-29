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
    columnDiameter = None
    heatCaluclation = None      # version of the heat volume (and so the corresponding columns) to use

    parser = OptionParser('Extract profiles from T2 nobias data using cortex-density-coordinates in ROIs')    
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    parser.add_option('-c', dest='columnDiameter', help='columnDiameter to work with')
    #parser.add_option('-s', dest='realSide', help='Hemisphere to be processed: L or R. L is default')   
    parser.add_option('-d', dest='directory', help='directory')
    parser.add_option('-j', dest='heatCaluclation', help='Version of the heat calculation program: old or new') 
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
        
    if options.columnDiameter is not None:
        columnDiameter = int(options.columnDiameter)    
        
    if options.heatCaluclation is None:
        print >> sys.stderr, 'New: exit. No version for the heat calculation was given'
        sys.exit(1)
    else:
        heatCaluclation = options.heatCaluclation

    pathToColumnResults = ''
    # set it depending on the version of the heat equation that has to be used
    if heatCaluclation == 'old':
        print 'selected the old version for the heat calculation'
        pathToColumnResults = directory + '%s_T1inT2_ColumnsCutNew20It_NewDB/' %(realPatientID)        
    elif heatCaluclation == 'new':
        print 'selected the new version for the heat calculation'
        pathToColumnResults = directory + '%s_T1inT2_ColumnsCutNew20It_NewDB_NewHeat/' %(realPatientID)        
    else:
        print 'selected the wrong keyword for the version for the heat calculation. Exit'
        sys.exit(0)    

    pathToProfL = pathToColumnResults + '%s_L_profiles2.txt' %(realPatientID)
    pathToProfR = pathToColumnResults + '%s_R_profiles2.txt' %(realPatientID)
    
    # check if both these profiles exist
    profL = glob.glob(pathToProfL)
    profR = glob.glob(pathToProfR)
    
    if len(profL) != 1 or len(profR) != 1:
        # abort the calculation, as too many or not a single texture file was found
        f = open(pathToColumnResults + '%s_compareProfilesStat.txt' %(realPatientID), "w")
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
    plt.savefig(pathToColumnResults + '%s_LvsR_2nobiasT2_%sHeat.png' %(realPatientID, heatCaluclation))  
    
    # also save this plot to a higher laval folder for comparison
    higherDirLvsR = directory.split(realPatientID)[0] + 'LvsR/'       # create this dir if it does not exist yet
    if not os.path.exists(higherDirLvsR):
        os.makedirs(higherDirLvsR)
        print 'created directory ', higherDirLvsR
    plt.savefig(higherDirLvsR + '%s_LvsR_2nobiasT2_%sHeat.png' %(realPatientID, heatCaluclation))  
    plt.clf()
    plt.close()
        
    
    # now plot L vs R in various ROIs
    # af140169_R_profiles2_ROI_21.txt, af140169_R_profiles2_ROI_11.txt and the same with L
    pathToROIsProfL = pathToColumnResults + '%s_L_profiles2_scaledNoOutlAddOutl_ROI_[0-9]*.txt' %(realPatientID)
    print 'pathToROIsProfL'
    print pathToROIsProfL
    pathToROIsProfR = pathToColumnResults + '%s_R_profiles2_scaledNoOutlAddOutl_ROI_[0-9]*.txt' %(realPatientID)
    
    # check if both these profiles exist. and how many are there
    profROIsL = glob.glob(pathToROIsProfL)
    profROIsR = glob.glob(pathToROIsProfR)
    
    print 'profROIsL'
    print profROIsL
    print 'profROIsR'
    print profROIsR
    numOfCommonRLregions = 0

    for i in range(len(profROIsL)):
        # check if this L profile is from the same mask ROI as the right one:
        # check if the corresponding R - file exists
        profROIsRcorrespond = glob.glob(profROIsL[i].replace('_L_', '_R_'))
        if len(profROIsRcorrespond) == 1:      
            numOfCommonRLregions += 1
            print 'corresponding files ', profROIsL[i], ' and ', profROIsRcorrespond[0]
            # read in the respective files
            numbROIsL, coordROIsL, valueROIsL = np.loadtxt(profROIsL[i], skiprows = 1, unpack = True)
            numbROIsR, coordROIsR, valueROIsR = np.loadtxt(profROIsRcorrespond[0], skiprows = 1, unpack = True)
            
            # get the ROI ID
            iD = (profROIsL[i].split('_ROI_')[1]).split('.txt')[0] #/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles_NewDB/ad140157/ad140157_T1inT2_ColumnsCutNew20It_NewDB/ad140157_L_profiles2_scaledNoOutlAddOutl_ROI_11.txt
            
            # plot the data
            plt.plot(coordROIsL, valueROIsL, '.', c = 'b', label = 'L')
            plt.title('Profile in ROI %s ' %(iD))   # subplot 211 title
            plt.xlabel('Cortical depth')
            plt.ylabel('T2-nobias intensity')
            plt.plot(coordROIsR, valueROIsR, '.', c = 'r', label = 'R')
            plt.legend(loc='upper right', numpoints = 1)
            plt.savefig(pathToColumnResults + '%s_LvsR_2nobiasT2_%sHeat_ROI_%s.png' %(realPatientID, heatCaluclation, iD))               
            # save profiles also to the outer folder!    
            plt.savefig(higherDirLvsR + '%s_LvsR_2nobiasT2_%sHeat_ROI_%s.png' %(realPatientID, heatCaluclation, iD))  
            plt.clf()
            plt.close() 
    print 'found ', numOfCommonRLregions, ' numOfCommonRLregions '       
            
       
    # for better comparison: plot the same plots but given max min scaled intensity ranges (-6, 10)





    # plot these plots into 1 image
    fig = plt.figure(figsize=(21, 6)) #, dpi=80, facecolor='w', edgecolor='k')
    numOfCommonRLregions += 1
    ax1 = fig.add_subplot(1,numOfCommonRLregions,1)
    ax1.plot(coordL, valueL, '.', c = 'b', label = 'L')
    ax1.set_title('Profile in all ROIs')   # subplot 211 title
    ax1.set_xlabel('Cortical depth')
    ax1.set_ylabel('T2-nobias intensity')
    ax1.plot(coordR, valueR, '.', c = 'r', label = 'R')
    #added
    ax1.set_ylim([-6, 10])
    ax1.legend(loc='upper right', numpoints = 1)
   
    for i in range(len(profROIsL)):
        profROIsRcorrespond = glob.glob(profROIsL[i].replace('_L_', '_R_'))
        if len(profROIsRcorrespond) == 1: 
            # read in the respective files
            numbROIsL, coordROIsL, valueROIsL = np.loadtxt(profROIsL[i], skiprows = 1, unpack = True)
            numbROIsR, coordROIsR, valueROIsR = np.loadtxt(profROIsRcorrespond[0], skiprows = 1, unpack = True)
            
            # get the ROI ID
            iD = (profROIsL[i].split('_ROI_')[1]).split('.txt')[0]
            ax2 = fig.add_subplot(1,numOfCommonRLregions, 2 + i)
            ax2.plot(coordROIsL, valueROIsL, '.', c = 'b', label = 'L')
            ax2.set_title('Profile in ROI %s' %(iD))   # subplot 211 title
            ax2.set_xlabel('Cortical depth')
            ax2.set_ylabel('T2-nobias intensity')
            ax2.plot(coordROIsR, valueROIsR, '.', c = 'r', label = 'R')
            #added
            ax2.set_ylim([-6, 10])
            ax2.legend(loc='upper right', numpoints = 1)

    #plt.show()
    print 'save the image /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/%s_LvsR_allVsROIs.png' %(realPatientID)
    plt.savefig(pathToColumnResults + '%s_LvsR_%sHeat_allVsROIs_uniqueY.png' %(realPatientID, heatCaluclation), bbox_inches='tight')               
    # save also to the outer folder!    
    plt.savefig(higherDirLvsR + '%s_LvsR_%sHeat_allVsROIs_uniqueY.png' %(realPatientID, heatCaluclation), bbox_inches='tight')              
    plt.clf()
    plt.close()
    print 'saved the image /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/%s_LvsR_allVsROIs.png' %(realPatientID)
    
    
    # plot it for sufficiently large testBatchColumnsExtrProfiles
    # for a given diameter
    # find all thresholds
    # plot e.g. ROI 21 LvsR for diam = 5 size > 392
    
    ## find those thresholds where non zero number of columns is available. for left
    ##largeL = glob.glob(directory + '/diam%s/%s_L_IDs_diam%s_over_[0-9]*.txt' %(columnDiameter, realPatientID, columnDiameter))
    #largeL = glob.glob(directory + '%s_L_nobiasT2_ROIs_11_21_diam%s_over[0-9]*_exclCommun.png' %(realPatientID, columnDiameter))
    ##ag140439_R_nobiasT2_ROIs_11_21_diam3_over197_exclCommun.png
    #print 'found largeL'
    #print largeL
    ## and for right realSide
    ##largeR = glob.glob(directory + '/diam%s/%s_R_IDs_diam%s_over_[0-9]*.txt' %(columnDiameter, realPatientID, columnDiameter))
    #largeR = glob.glob(directory + '%s_R_nobiasT2_ROIs_11_21_diam%s_over[0-9]*_exclCommun.png' %(realPatientID, columnDiameter))
    #print 'found largeR'
    #print largeR
    ## find common elements
    ##largeLnum = [int((x.split('_over_')[1]).split('.txt')[0]) for x in largeL]
    #largeLnum = [int((x.split('_over')[1]).split('_exclCommun.png')[0]) for x in largeL]
    #print largeLnum
    ##largeRnum = [int((x.split('_over_')[1]).split('.txt')[0]) for x in largeR]
    #largeRnum = [int((x.split('_over')[1]).split('_exclCommun.png')[0]) for x in largeR]
    #print largeRnum
    
    ## find common elelemnts
    #common = set(largeLnum) & set(largeRnum)
    #print common
    
    ## plot LvsR for various ROIs for these thresholds
    #for i in common:
        #print i
        ## find all points from the L hemisphere in ROI 11 that are in columns larger i
        ## open file _L_ColumnInfo.txt - see which columns to take into which ROIs
        
    
    
    

    

      
  
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    