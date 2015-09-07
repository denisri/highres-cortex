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
# this function calculate statistics from the profiles 


# example how to run this file on a laptop:
# python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_analyseProfiles.py -p at140353 -s L -c 3 -d /volatile/od243208/neurospin/testBatchColumnsExtrProfiles_NewDB/at140353/ -l -j new

# on a local machine:
# python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_analyseProfiles.py -p at140353 -s L -c 3 -d /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles_NewDB/at140353/ -j new

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
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
from sklearn import datasets
from pylab import plot,show
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances
from sklearn import preprocessing
from sklearn.decomposition import PCA
from operator import add

if __name__ == '__main__':
    
    healthyList = ['fg140290', 'af140169', 'ml140175', 'ac140159', 'md140208', 'at140353'] #, 'js140266', 'lg140146', 'he140338', 'cb140330']
    dyslexicList = ['js140311', 'ad140157', 'ag140439', 'sg140335']
    realPatientID = None
    directory = None
    threshold = None
    #realSide = 'L'
    columnDiameter = None
    workOnLaptop = False
    numOfLayersToUse = 6
    heatCaluclation = None      # version of the heat volume (and so the corresponding columns) to use
    pathToNobiasT2_newCroppedDB = '/neurospin/lnao/dysbrain/brainvisa_db_T2_cropped/dysbrain/'   
    pathToNew_T1inT2DB = '/neurospin/lnao/dysbrain/brainvisa_db_T1_in_T2_space/dysbrain/'
    pathToNobiasT2_newCroppedDB = '/neurospin/lnao/dysbrain/brainvisa_db_T2_cropped/dysbrain/'   


#    corticalIntervals = [0, 0.1, 0.2, 0.5, 0.62, 0.82, 1.0]
#    corticalIntervals = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.62, 0.72, 0.82, 0.91, 1.0]
#    corticalIntervals = [0, 0.1, 0.2, 0.35, 0.5, 0.62, 0.72, 0.82, 1.0]
    #corticalIntervals = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    #addedName = '_1PpDE'         # 1 point per depth entity
    addedName = '_1PpL'        # 1 point per Layer
    corticalIntervals = [0, 0.1, 0.2, 0.5, 0.6, 0.8, 1.0]

    parser = OptionParser('Analyze profiles from T2 nobias data using cortex-density-coordinates in ROIs')    
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    parser.add_option('-s', dest='realSide', help='Hemisphere to be processed: L or R. L is default')   
    parser.add_option('-c', dest='columnDiameter', help='columnDiameter to work with')
    parser.add_option('-n', dest='numOfLayersToUse', help='numOfLayersToUse for the cortex analysis (6 is default)')
    #parser.add_option('-t', dest='threshold', help='threshold to work with')
    parser.add_option('-d', dest='directory', help='directory')
    parser.add_option('-l', dest='workOnLaptop', action = 'store_true', help='Select if working on laptop (neurospin DB location is different. False is default') 
    parser.add_option('-j', dest='heatCaluclation', help='Version of the heat calculation program: old or new') 
    options, args = parser.parse_args(sys.argv)
    print options
    print args   
    
    #realPatientID = 'at140353'
    #realSide = 'L'
    #columnDiameter = 3
    #directory = '/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles_NewDB/at140353/'
    #heatCaluclation = 'new'
    
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
        
    if options.realSide is None:
        print >> sys.stderr, 'New: exit. no realSide given'
        sys.exit(1)
    else:
        realSide = options.realSide  
        
    if options.numOfLayersToUse is not None:
        print >> sys.stderr, 'New: given numOfLayersToUse.'
        numOfLayersToUse = options.numOfLayersToUse
    else:
        print >> sys.stderr, 'Keep numOfLayersToUse = 6'  
        
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
     	
    if options.workOnLaptop is not None:
	workOnLaptop = options.workOnLaptop      
	# if true, then processes are run on the laptop. Change locations of neurospin DBs
	# TODO! check! which one to use!
	#directory = directory.replace('/neurospin/lnao/', '/nfs/neurospin/lnao/')        
	directory = directory.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')  
        pathToNobiasT2_newCroppedDB = pathToNobiasT2_newCroppedDB.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')
        pathToNew_T1inT2DB = pathToNew_T1inT2DB.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')	
        pathToNobiasT2_newCroppedDB = pathToNobiasT2_newCroppedDB.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')

    addedColumnsDiamName = ''
    pathForFiles = pathToColumnResults
    
    if options.columnDiameter is not None:
        columnDiameter = int(options.columnDiameter)
        addedColumnsDiamName = '_diam%s' %(columnDiameter)
        # if the result was for the columns, find a folder for it
        pathForFiles = pathForFiles + 'diam%s/'%(columnDiameter)
        print '############################################ found profile directory = ', pathForFiles

    sizeSorted = pathForFiles + 'sortedBySize/'
    divGradnSorted = pathForFiles + 'sortedByDivGradn/'
    corticalLayers = pathForFiles + 'corticalLayers/'
    classification = pathForFiles + 'classification/'

    if not os.path.exists(sizeSorted):
        os.makedirs(sizeSorted)
    if not os.path.exists(divGradnSorted):
        os.makedirs(divGradnSorted)
    if not os.path.exists(corticalLayers):
        os.makedirs(corticalLayers)
    if not os.path.exists(classification):
        os.makedirs(classification)
        print 'created classification = ', classification
    else:
        print 'existed classification = ', classification
    healthyDir = directory.split(realPatientID)[0] + 'healthy/'
    dyslexicDir = directory.split(realPatientID)[0] + 'dyslexic/'
    if not os.path.exists(healthyDir):
        os.makedirs(healthyDir)
    if not os.path.exists(dyslexicDir):
        os.makedirs(dyslexicDir)
    healthyDir_NewHeat = directory.split(realPatientID)[0] + 'healthy_NewHeat/'
    dyslexicDir_NewHeat = directory.split(realPatientID)[0] + 'dyslexic_NewHeat/'
    if not os.path.exists(healthyDir_NewHeat):
        os.makedirs(healthyDir_NewHeat)
    if not os.path.exists(dyslexicDir_NewHeat):
        os.makedirs(dyslexicDir_NewHeat)
       
    outerDir = ''
    # decide - current subject is - so find the folder to write out files
    if realPatientID in healthyList:
        if heatCaluclation == 'new':        
            outerDir = healthyDir_NewHeat
        else:
            outerDir = healthyDir
    elif realPatientID in dyslexicList:
        if heatCaluclation == 'new':
            outerDir = dyslexicDir_NewHeat
        else:
            outerDir = dyslexicDir
    
    #new folder for avg ROI profiles for cortical layers
    outerDirAvgProfCortLayers = outerDir + 'avgProfCortLayers/'
    if not os.path.exists(outerDirAvgProfCortLayers):
        os.makedirs(outerDirAvgProfCortLayers)
    
    #if options.threshold is None:
        #print >> sys.stderr, 'New: exit. no threshold given'
        #sys.exit(1)
    #else:
        #threshold = options.threshold  
        
    # TODO: for test. maybe delete it later
    # plot histograms for the intensities in CSF, cortex and WM
    pathToClassifNoBorders = pathToColumnResults + 'GWsegm_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)
    volGWNoborders = aims.read(pathToClassifNoBorders)
    arrGWNoborders = np.array(volGWNoborders)
    pathToClassifWithBorders = pathToColumnResults + '/dist/classif_with_outer_boundaries_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)
    volGWborders = aims.read(pathToClassifWithBorders)
    arrGWborders = np.array(volGWborders)
    pathToNobiasT2_new = pathToNobiasT2_newCroppedDB + '%s/t1mri/t2_resamp/%s.nii.gz' %(realPatientID, realPatientID)
    volT2 = aims.read(pathToNobiasT2_new)
    arrT2 = np.array(volT2)
    statFileName = classification + "statFile_analyzeProf_%s_%s.txt" %(realPatientID, realSide)
    print 'statFileName = ', statFileName
    f = open(statFileName, "w")
    f.write('pathToClassifNoBorders\t' + pathToColumnResults + 'GWsegm_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide) + '\n')
    f.write('pathToClassifWithBorders\t' + pathToColumnResults + '/dist/classif_with_outer_boundaries_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide) + '\n')
    f.write('pathToNobiasT2_new\t' + pathToNobiasT2_newCroppedDB + '%s/t1mri/t2_resamp/%s.nii.gz' %(realPatientID, realPatientID) + '\n')    








    

  ###############################################################################################
    # TODO: plot e.g. avg for PT for all healthy subjects on one plot. To compare their vaiance
    meanListOfLists = []
    stdvListOfLists = []
    namesListOfLists = []
    fig = plt.figure(figsize=(3 * 7, 6 * 2))
    # for the beginning do it for old heat, olny, only L hemisphere
    # iterate for all ROIs, PT, Heschl
    roiS = ['allROIs', 'ROI_11', 'ROI_21']
    for j in range(len(roiS)):	
	print '=============================== work with ROI  ' + roiS[j] + ' =================================='
	ax = fig.add_subplot(2,3,(j + 1))
	meanList = []
	stdvList = []
	namesList = []
	for realPatientIDcurr in healthyList:
	    # read in files and plot
	    dirToRead = directory.split(realPatientID)[0] + realPatientIDcurr + '/'
	    print dirToRead, '-----------------------------'
	    filesToRead = glob.glob(dirToRead + '%s_T1inT2_ColumnsCutNew20It_NewDB/%s_%s_CortLayersROI_%s.txt' %(realPatientIDcurr, realPatientIDcurr, realSide, roiS[j]))
	    print filesToRead		    #/volatile/od243208/neurospin/testBatchColumnsExtrProfiles_NewDB/ad140157/ad140157_T1inT2_ColumnsCutNew20It_NewDB/ad140157_L_CortLayersROI_ROI_11.txt
	    
	    # read in these files and collect the info
	    if len(filesToRead) != 1:
		print 'too much or no info'
	    else:
		currCoords, currMeans, currStdvs  = np.loadtxt(filesToRead[0], skiprows = 1, usecols = (1, 2, 3), unpack=True)
		meanList.append(currMeans)
		stdvList.append(currStdvs)
		currName = realPatientIDcurr + '_' + roiS[j]
		namesList.append(currName)
		print 'plot the errorbar'
		ax.errorbar(currCoords, currMeans, currStdvs, linestyle='solid', marker='^', label = currName)
			
	# collected the data on healthy subjects. plot on 1 plot
	ax.set_title('Mean profiles for healthy subjects in ' + roiS[j])
	ax.set_xlabel('Cortical depth') 
	ax.set_ylabel('T2-nobias intensity')                 
	ax.legend(loc='upper right', numpoints = 1)  

	ax2 = fig.add_subplot(2,3,(j + 4))
	meanList2 = []
	stdvList2 = []
	namesList2 = []
	for realPatientIDcurr in dyslexicList:
	    # read in files and plot
	    dirToRead = directory.split(realPatientID)[0] + realPatientIDcurr + '/'
	    print dirToRead, '-----------------------------'
	    filesToRead = glob.glob(dirToRead + '%s_T1inT2_ColumnsCutNew20It_NewDB/%s_%s_CortLayersROI_%s.txt' %(realPatientIDcurr, realPatientIDcurr, realSide, roiS[j]))
	    print filesToRead		    #/volatile/od243208/neurospin/testBatchColumnsExtrProfiles_NewDB/ad140157/ad140157_T1inT2_ColumnsCutNew20It_NewDB/ad140157_L_CortLayersROI_ROI_11.txt
	    
	    # read in these files and collect the info
	    if len(filesToRead) != 1:
		print 'too much or no info'
	    else:
		currCoords, currMeans, currStdvs  = np.loadtxt(filesToRead[0], skiprows = 1, usecols = (1, 2, 3), unpack=True)
		meanList.append(currMeans)
		stdvList.append(currStdvs)
		currName = realPatientIDcurr + '_' + roiS[j]
		namesList.append(currName)
		print 'plot the errorbar'
		ax2.errorbar(currCoords, currMeans, currStdvs, linestyle='solid', marker='^', label = currName)
			
	# collected the data on healthy subjects. plot on 1 plot
	ax2.set_title('Mean profiles for dyslexic subjects in ' + roiS[j])
	ax2.set_xlabel('Cortical depth') 
	ax2.set_ylabel('T2-nobias intensity')                 
	ax2.legend(loc='upper right', numpoints = 1)  
	
    plt.savefig(directory.split(realPatientID)[0] + 'healthyVsDyslexic_MeanStdv.png')
    plt.clf()
    plt.close()
	    
      
    
    sys.exit(0)











    #collect individual values for the separate layers. For all ROIs, for PT, for Heschl. To compare then variance between healthy/dyslexic subjects
    #at140353_L_profiles2, at140353_L_profiles2_scaledNoOutlAddOutl_ROI_11, at140353_L_profiles2_scaledNoOutlAddOutl_ROI_21
    allROIs = pathToColumnResults + '%s_%s_profiles2.txt' %(realPatientID, realSide)    
    separateROIs = pathToColumnResults + '%s_%s_profiles2_scaledNoOutlAddOutl_ROI_*.txt' %(realPatientID, realSide)    
    allROIsFiles = glob.glob(allROIs)
    # check how many files there are. only if one -> no ambiguity
    if len(allROIsFiles) != 1:
        print 'abort the calculation, as too many or not a single allROIsFile was found'
        f.write('abort the calculation, as ' + str(len(allROIsFiles)) + ' allROIsFile were found' + '\n')
        f.close()
        sys.exit(0)    
        
    allROIsFile = allROIsFiles[0]  
    # find files for separate ROIs
    separeteROIsFiles = glob.glob(separateROIs) # should be 11 and 21
    # extract ROI names
    roiNames = [(i.split('.txt')[0]).split('_ROI_')[1] for i in separeteROIsFiles]
    allFiles = [allROIsFile]
    allFiles.extend(separeteROIsFiles)
    names = ['allROIs']
    names.extend([('ROI_' + i) for i in roiNames])    
    fig = plt.figure(figsize=(1 * 7, 6)) #, dpi=80, facecolor='w', edgecolor='k')
    for fileROI, name, num in zip(allFiles, names, range(len(names))):
        # read in the file
        print 'read in the file ', fileROI
        #/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles_NewDB/at140353/at140353_T1inT2_ColumnsCutNew20It_NewDB_NewHeat/at140353_L_profiles2_scaledNoOutlAddOutl_ROI_21.txt
        currCoords, currValues  = np.loadtxt(fileROI, skiprows = 1, usecols = (1, 2), unpack=True)
        #print zip(currCoords[0:5], currValues[0:5])
        means = []
        stdvs = []
        xCoords = []
        
        # extract individual layers, calculate mean, stdv, plot
        for c in range(len(corticalIntervals) - 1):
            start = corticalIntervals[c]
            stop = corticalIntervals[c + 1]
            # find all points where the depth coordinate is between these thresholds
            thisLayerValues = currValues[np.where((currCoords >= start) & (currCoords < stop))]                
            # save their means and stdvs
            means.append(np.mean(thisLayerValues))
            stdvs.append(np.std(thisLayerValues))
            xCoords.append((start + stop) / 2.0)                        
           
        # collected the data for all layers. now plot for this particular ROI  
        #plt.errorbar(xCoords, means, stdvs, linestyle='solid', marker='^', label = name)
        # to plot it a bit shifted:
        xCoords_new = [xxx + ((len(names)/2 - num))*0.005 for xxx in xCoords]
        plt.errorbar(xCoords_new, means, stdvs, linestyle='solid', marker='^', label = name)

        dataCort = open(pathToColumnResults + '%s_%s_CortLayersROI_%s.txt' %(realPatientID, realSide, name), "w")
        dataCort.write('CorticalLayer\tAvgCoord\tMeanValue\tStdValue\n')            

        for j in range(len(xCoords)):
            dataCort.write(str(j + 1) + '\t' + str(xCoords[j]) + '\t' + str("%.4f" % means[j]) + '\t' + str("%.4f" % stdvs[j]) + '\n')    
        dataCort.close()
    
    plt.title('Mean and stdv values in cortical layers in ROIs') 
    plt.xlabel('Cortical depth') 
    plt.ylabel('T2-nobias intensity')                 
    plt.legend(loc='upper right', numpoints = 1)   
    plt.savefig(pathToColumnResults + '%s_%s_cortLayersInROIs.png' %(realPatientID, realSide))
    #save also to outer folder
    plt.savefig(outerDirAvgProfCortLayers + '%s_%s_cortLayersInROIs.png' %(realPatientID, realSide))
    plt.clf()
    plt.close()




    sys.exit(0)
#################################################################################################################################################################
################################################################ analyze individual cortical columns. Classify them #############################################
#################################################################################################################################################################


    arr150 = arrT2[arrGWborders == 150]
    arr50 = arrT2[arrGWborders == 50]
    arr100 = arrT2[arrGWborders == 100]
    arr200 = arrT2[arrGWborders == 200]    
    
    plt.hist(arr100, label = '100 (cortex)', alpha = 0.3)
    plt.hist(arr150, label = '150 (border to WM)', alpha = 0.7)
    plt.hist(arr50, label = '50 (border to CSF)', alpha = 0.7)
    plt.title('Histogram of intensities')   # subplot 211 title
    plt.xlabel('T2 nobias intensities')
    plt.ylabel('Number of voxels')
    plt.legend(loc='upper right')
    #plt.show()
    plt.savefig(directory + '%s_%s_histOfT2nobiasIntensities.png' %(realPatientID, realSide))    
    plt.clf()
    plt.close()
            
    # read in the columns info file
    toLookFor = pathForFiles + '%s_%s_ColumnInfo_*' %(realPatientID, realSide)    
    print 'columns info file toLookFor = ', toLookFor
    infoFileNames = glob.glob(toLookFor)                              
    # check how many files there are. only if one -> no ambiguity
    if len(infoFileNames) != 1:
        print 'abort the calculation, as too many or not a single column info file was found'
        f.write('abort the calculation, as ' + str(len(infoFileNames)) + ' column info files were found' + '\n')
        f.close()
        sys.exit(0)    

    infoFileName = infoFileNames[0]    
    # extract the keyword
    addedColumnsDiamName = (infoFileName.split('_ColumnInfo_')[1]).split('.txt')[0]
    print 'extracted addedColumnsDiamName = ', addedColumnsDiamName    
    #ad140157_L_ColumnInfo_diam3_scaledNoOutlAddOutl.txt                              
    #ad140157_L_profiles2_diam3_scaledNoOutlAddOutl_ROI_511
    columnID, size, avgX, avgY, avgZ, avgDivGradn  = np.loadtxt(infoFileName, skiprows = 1, usecols = (0, 1, 2, 3, 4, 5), unpack=True)
    #print zip(columnID, size, avgX, avgY, avgZ, avgDivGradn)
        
    # calculatethe size threshold for columns - (volume of a cylinder, diam is given, resolution 0.5, height - 2mm min, and divided by 4 for conical cases, or too narrow cases)
    sizeThreshold = np.round(columnDiameter * columnDiameter * 4.0 / 4.0 * np.pi * 4.0 / 5.0) 
    print 'sizeThreshold = ', sizeThreshold
    f.write('sizeThreshold\t' + str(sizeThreshold) + '\n')    
    
    # find columns larger than this threshold
    largeIDs = columnID[np.where(size >= sizeThreshold)]
    print 'largeIDs = ', largeIDs
    f.write('percentage of largeIDs is \t' + str((len(largeIDs) * 1.0/len(columnID)*100.0)) + '% \n') 
    print 'percentage of largeIDs is \t' + str((len(largeIDs) * 1.0/len(columnID)*100.0))
    f.write('largeIDs\t' + str(largeIDs) + '\n') 
    completeN = 0
    incompleteN = 0
    complete5LayersN = 0    
    
    # create data structure for the collected data! the first column: ROI ID !!!!!!!!!!
    arrProfiles6Layers = []
    arrProfilesID6Layers = []
    arrVar6Layers = []
    arrFeatures6Layers = []
    
    arrProfiles5Layers = []
    arrProfilesID5Layers = []
    arrVar5Layers = []
    arrFeatures5Layers = []
    
    for i in range(len(largeIDs)):
        # iterate through these profiles, calculate means, stdv in layers, plot figures, save files
        # find profiles of these large columns: ad140157_L_profiles2_diam3_ROI_259.txt - find all profile files like this one
        currID = int(largeIDs[i])
        #print 'currID = ', currID                        
        #ad140157_L_profiles2_diam3_scaledNoOutlAddOutl_ROI_478.txt        - new names after scaling and adding scaled outliers!
        profileFile = glob.glob(pathForFiles + '%s_%s_profiles2_%s_ROI_%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(currID)))[0]
        #print 'work with file  ', profileFile
        currCoords, currValues  = np.loadtxt(profileFile, skiprows = 1, usecols = (1, 2), unpack=True)
        #print zip(depthCoord, value)
        
        fig = plt.figure(figsize=(2 * 7, 6)) #, dpi=80, facecolor='w', edgecolor='k')
        axPoints = fig.add_subplot(1,2,1)
        axMeans = fig.add_subplot(1,2,2)
        axPoints.plot(currCoords, currValues, '.', c = 'b') #, label = 'ROI '+ str(maskROIids[j]))
        axPoints.set_title('Profile in ROI %s' %(str(currID)))
        axPoints.set_xlabel('Cortical depth') 
        axPoints.set_ylabel('T2-nobias intensity')                 
        axPoints.legend(loc='upper right', numpoints = 1)   
        
        # do it only for columns, where data in EACH cortical layer is available
        means = []
        stdvs = []
        xCoords = []
        complete = True     # shows whether points in every cortical layer were found        
        # do the same for the case when we 'ignore' the layer I
        complete5Layers = True     # shows whether points in 5 cortical layers (II - VI) were found
        
        for c in range(len(corticalIntervals) - 1):
            start = corticalIntervals[c]
            stop = corticalIntervals[c + 1]
            # find all points where the depth coordinate is between these thresholds
            #thisLayerCoords = currCoords[np.where((currCoords >= start) & (currCoords < stop)]
            thisLayerValues = currValues[np.where((currCoords >= start) & (currCoords < stop))]
            # check the length: if zero: ignore this column
            if (complete & (len(thisLayerValues) == 0)):
                complete = False                          
            if (complete5Layers & (c != 0) & (len(thisLayerValues) == 0)):
                complete5Layers = False
                
            # save their means and stdvs
            means.append(np.mean(thisLayerValues))
            stdvs.append(np.std(thisLayerValues))
            xCoords.append((start + stop) / 2.0)            
            #print ' cortLayer = ', c, ' len(thisLayerValues)=  ', len(thisLayerValues), ' complete5Layers= ',  complete5Layers, ' complete = ', complete
            # save the individual values
            indValuesInLayers.append(thisLayerValues)
            
        # collected the data for all layers. now plot for this particular column  
        axMeans.set_title('Mean and stdv values in cortical layers in maskROI %s' %(str(currID)))  
        axMeans.set_xlabel('Cortical depth') 
        axMeans.set_ylabel('T2-nobias intensity')                 
        axMeans.legend(loc='upper right', numpoints = 1)  
        
        if complete:   # plot the profile for all 6 cortical layers
            axMeans.errorbar(xCoords, means, stdvs, linestyle='solid', marker='^', ecolor = 'r')
            plt.savefig(corticalLayers + '%s_%s_completeCortLay%s_size%s_ROI_%s%s.png' %(realPatientID, realSide, addedColumnsDiamName, str(len(currCoords)), str(currID), addedName))                
            # save txt files with means and stdvs !!
            dataCort = open(corticalLayers + '%s_%s_CortLay%s_ROI_%s%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(currID), addedName), "w")
            dataCort.write('CorticalLayer\tAvgCoord\tMeanValue\tStdValue\n')            
    
            for j in range(len(xCoords)):
                dataCort.write(str(j + 1) + '\t' + str(xCoords[j]) + '\t' + str("%.4f" % means[j]) + '\t' + str("%.4f" % stdvs[j]) + '\n')    
            dataCort.close()  
            completeN += 1
            
            # also plot the profile for the 5 cortical layers II - VI
            axMeans.clear()
            axMeans.errorbar(xCoords[1:], means[1:], stdvs[1:], linestyle='solid', marker='^', ecolor = 'r')
            plt.savefig(corticalLayers + '%s_%s_5CortLay%s_size%s_ROI_%s%s.png' %(realPatientID, realSide, addedColumnsDiamName, str(len(currCoords)), str(currID), addedName))   
                        
            # add the infor from this profile to the data
            arrProfilesID6Layers.append(currID)
            arrProfiles6Layers.append(means)
            arrVar6Layers.append(stdvs)                  
            
            # create 'full' feature vectors, including profiles, avg coordinates of columns, sizes, avgDivGrads.... for 6 layers
            currFeatureVector = []
            currFeatureVector.extend(means)
            currFeatureVector.extend(stdvs)     
            # find the info on this column in the columnInfo file. Line numbers are different!!
            infoFileID = np.where(columnID == currID)[0]
            #print 'size ', size, 'type ', type(size)
            #print size[infoFileID], 'type ', type(size[infoFileID]) 
            #print (size[infoFileID])[0], 'type ', type((size[infoFileID])[0])
            featuresFromInfoFile = [(size[infoFileID])[0], (avgX[infoFileID])[0], (avgY[infoFileID])[0], (avgZ[infoFileID])[0], (avgDivGradn[infoFileID])[0]]
            #print 'featuresFromInfoFile = ', featuresFromInfoFile
            #print 'currRealROIid ', currID, ' line number in columnInfofile', infoFileID, ' columnId from infoFile ', columnID[infoFileID], ' wrong Id ', columnID[i]  
            currFeatureVector.extend(featuresFromInfoFile)      # avg column coordinates in VOXELS!!! # TODO? need to transform into mm ?? maybe not, as we do scaling
            arrFeatures6Layers.append(currFeatureVector)
            
            # add the infor from this profile to the data for 5 layers
            arrProfilesID5Layers.append(currID)
            arrProfiles5Layers.append(means[1:])  
            arrVar5Layers.append(stdvs[1:])          
            
            # create 'full' feature vectors for 5 layers
            currFeatureVector5 = []
            currFeatureVector5.extend(means[1:])
            currFeatureVector5.extend(stdvs[1:])
            currFeatureVector5.extend(featuresFromInfoFile)     # avg column coordinates in VOXELS!!! # TODO? need to transform into mm ?? maybe not, as we do scaling
            arrFeatures5Layers.append(currFeatureVector5)            
        else :
            # now check, whether this profile is complete for 5 layers (II - VI)
            if complete5Layers == False:
                # this profile is incomplete even for 5 layers
                #print 'incomplete ID , even for 5 layers ', str(currID)
                axMeans.errorbar(xCoords, means, stdvs, linestyle='None', marker='^', ecolor = 'r')
                plt.savefig(corticalLayers + '%s_%s_incompleteCortLay%s_size%s_ROI_%s%s.png' %(realPatientID, realSide, addedColumnsDiamName, str(len(currCoords)), str(currID), addedName))  
                # do NOT save separate txt files, as the data is incomplete and can not be used for clustering
                incompleteN += 1
            else:   
                # this profile is incomplete for 6 layers, but is complete for 5 layers
                axMeans.errorbar(xCoords[1:], means[1:], stdvs[1:], linestyle='solid', marker='^', ecolor = 'r') # plot means and stdvs for layers II - VI
                plt.savefig(corticalLayers + '%s_%s_5CortLay%s_size%s_ROI_%s%s.png' %(realPatientID, realSide, addedColumnsDiamName, str(len(currCoords)), str(currID), addedName))                
                # save txt files with means and stdvs !!
                dataCort = open(corticalLayers + '%s_%s_CortLay%s_ROI_%s%s_5Layers.txt' %(realPatientID, realSide, addedColumnsDiamName, str(currID), addedName), "w")
                dataCort.write('CorticalLayer\tAvgCoord\tMeanValue\tStdValue\n')            
        
                for j in range(len(xCoords) - 1):
                    dataCort.write(str(j + 2) + '\t' + str(xCoords[j + 1]) + '\t' + str("%.4f" % means[j + 1]) + '\t' + str("%.4f" % stdvs[j + 1]) + '\n')    
                dataCort.close()  
                complete5LayersN += 1    
                
                # save the data for analysis
                arrProfilesID5Layers.append(currID)
                arrProfiles5Layers.append(means[1:])
                arrVar5Layers.append(stdvs[1:])
                
                # create 'full' feature vectors for 5 layers!!!!!!!!!!!
                currFeatureVector5 = []
                currFeatureVector5.extend(means[1:])
                currFeatureVector5.extend(stdvs[1:])
                infoFileID = np.where(columnID == currID)
                featuresFromInfoFile = [size[infoFileID], avgX[infoFileID], avgY[infoFileID], avgZ[infoFileID], avgDivGradn[infoFileID]]
                currFeatureVector5.extend(featuresFromInfoFile) # avg column coordinates in VOXELS!!! # TODO? need to transform into mm ?? maybe not, as we do scaling
                arrFeatures5Layers.append(currFeatureVector5)
#            plt.savefig(pathForFiles + '%s_%s_nobiasT2%s_ROI_' %(realPatientID, realSide, addedColumnsDiamName) + str(iDs[i]) + '.png')        
            # write out the same info, but sorted by sizes
#            plt.savefig(sizeSorted + '%s_%s_nobiasT2%s_size%s_ROI_' %(realPatientID, realSide, addedColumnsDiamName, str(len(currCoords))) + str(iDs[i]) + '.png')
            
            ## now save the same plots but sorted by their avg gradients              
#            plt.savefig(divGradnSorted + '%s_%s_nobiasT2%s_avgDivGradn%s_ROI_%s.png' %(realPatientID, realSide, addedColumnsDiamName, strCurrGradnStr, str(iDs[i])))
        plt.clf()
        plt.close()     
    
    print '*********************************************** complete = ', completeN, ', incomplete = ', incompleteN, ', complete5LayersN = ', complete5LayersN, ' from total of ', len(columnID)
    # save these id-s of the columns
    f.write('arrProfilesID6Layers\t' + str(arrProfilesID6Layers) + '\n') 
    arrProfilesID6LayersFile = open(classification + 'arrProfilesID6Layers_%s_%s.txt' %(realPatientID, realSide), 'w')
    arrProfilesID6LayersFile.write('arrProfilesID6Layers\t' + str(arrProfilesID6Layers) + '\n')
    arrProfilesID6LayersFile.close()
    # save 5 layers largeIDs
    arrProfilesID5LayersFile = open(classification + 'arrProfilesID5Layers_%s_%s.txt' %(realPatientID, realSide), 'w')
    arrProfilesID5LayersFile.write('arrProfilesID5Layers\t' + str(arrProfilesID5Layers) + '\n')
    arrProfilesID5LayersFile.close()    

    # save extracted feature vectors
    np.save(classification + 'arrFeatures6Layers_%s_%s.npy' %(realPatientID, realSide), arrFeatures6Layers)
    np.save(classification + 'arrFeatures5Layers_%s_%s.npy' %(realPatientID, realSide), arrFeatures5Layers)  
    #print arrFeatures6Layers

#################################################################################################################################################################
#######################################################  STATISTICAL ANALYSIS  ##################################################################################
#################################################################################################################################################################

    # work with 5 or 6 layers?
    data = []
    arrIDs = []
    if numOfLayersToUse == 6:
        data = arrFeatures6Layers
        arrIDs = arrProfilesID6Layers
    elif numOfLayersToUse == 5:
        data = arrFeatures5Layers
        arrIDs = arrProfilesID5Layers
    ############################ PCA, correlation analysis! cluster using these full feature vectors 

    #X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    #pca = PCA(n_components=2)
    #pca.fit(X)
    #PCA(copy=True, n_components=2, whiten=False)
    #print(pca.explained_variance_ratio_) 
    #[ 0.99244...  0.00755...]

    pca = PCA(n_components = 2)
    pca.fit(data)
    print(pca.explained_variance_ratio_) 
    f.write('PCA(n_components = 2), explained_variance_ratio_\t' + str(pca.explained_variance_ratio_) + '\n') 
    #[  7.94485291e-01   1.46662429e-01   3.95866252e-02   1.79609243e-02
    #5.54970106e-04   2.63572886e-04   1.88952565e-04   8.36929421e-05]   - if n_components = 8
    print 'total covered variance is ', np.sum(pca.explained_variance_ratio_)
    f.write('total covered variance is \t' + str(np.sum(pca.explained_variance_ratio_)) + '\n') 

    # looks like 2-3 components should be enough. print them and save them
    data_PCA2 = pca.fit_transform(data)
    print 'done PCA transform for 2 components. Save the result at ', classification + 'arrFeatures%sLayers_%s_%s_PCA2.npy' %(str(numOfLayersToUse), realPatientID, realSide)
    np.save(classification + 'arrFeatures%sLayers_%s_%s_PCA2.npy' %(str(numOfLayersToUse), realPatientID, realSide), data_PCA2)
    print data_PCA2
    
    # can now use this data for clustering
    
    
    
    
    
    # linear model of PCA components? - need a 'response variable'. - this will be the z-score
    



    # 1. Cluster the columns using vectors of their means and stdvs in the 6 'cortical layers'
    
    ## data generation
    #data = vstack((rand(150,2) + array([.5,.5]),rand(150,2)))

    ## computing K-Means with K = 2 (2 clusters)
    #centroids,_ = kmeans(data,2)
    ## assign each sample to a cluster
    #idx,_ = vq(data,centroids)

    ## some plotting using numpy's logical indexing
    #plot(data[idx==0,0],data[idx==0,1],'ob',
        #data[idx==1,0],data[idx==1,1],'or')
    #plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
    #show()





    
    # data scaling TODO!
    data_scaled = preprocessing.scale(data)
    data_PCA2_scaled = preprocessing.scale(data_PCA2)
    #print zip(data, data_scaled)    
    
    # TODO! try data normalization?
    
    
    ## computing K-Means with K = 2 (2 clusters)
    #centroids,_ = kmeans(data, 2)
    
    ## assign each sample to a cluster
    #idx,_ = vq(data,centroids)

    ## some plotting using numpy's logical indexing
    #plot(data[idx==0,0],data[idx==0,1],'ob',
        #data[idx==1,0],data[idx==1,1],'or')
    #plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
    #show()

    #print 





    # Gael Varoquaux
    
    np.random.seed(5)

#    centers = [[1, 1], [-1, -1], [1, -1]]
#    iris = datasets.load_iris()
#    X = iris.data
#    y = iris.target

    #estimators = {'k_means_iris_3': KMeans(n_clusters=3),
                #'k_means_iris_8': KMeans(n_clusters=8),
                #'k_means_iris_bad_init': KMeans(n_clusters=3, n_init=1,
                                                #init='random')}
    estimators = {
		  'k_means_2': KMeans(n_clusters = 2),
                  'k_means_3': KMeans(n_clusters = 3),
                  'k_means_4': KMeans(n_clusters = 4)
                  #'agglomerative_2_eucl_avg' : AgglomerativeClustering(n_clusters = 2, linkage="average", affinity = "euclidean"),
                  #'agglomerative_3_eucl_avg' : AgglomerativeClustering(n_clusters = 3, linkage="average", affinity = "euclidean"),
		  #'agglomerative_2_cityblock_avg' : AgglomerativeClustering(n_clusters = 2, linkage="average", affinity = "cityblock"),
		  #'agglomerative_3_cityblock_avg' : AgglomerativeClustering(n_clusters = 3, linkage="average", affinity = "cityblock")
                  }
    
    # options!    
     
    
    ##2. scaled data. ALWAYS use scaled data! as intensities depend on pads, ... vary among patients    
    datasToUse = [data_scaled, data_PCA2_scaled]
    testNames = ['_dataScaled_%d_featureVect' %(numOfLayersToUse), '_dataPCA2Scaled_%d_featureVect' %(numOfLayersToUse)]
      
    
    

    #testName = 'initKMeans'    
    #for index, metric in enumerate(["cosine", "euclidean", "cityblock"]):
    #model = AgglomerativeClustering(n_clusters=n_clusters,
                                    #linkage="average", affinity=metric)
    #model.fit(X)
    
    for dataToUse, testName in zip(datasToUse, testNames):       
        labelsVariousK = []
        # take the initial volume with columns and colour it according to the labels
        inertiaVariousK = []
        f.write('testName\t' + testName + '\n') 

        for name, est in estimators.items():
            print 'start model ', name
            f.write('model\t' + name + '\n') 
            #est.fit(arrProfiles6Layers)
            est.fit(dataToUse)
            labels = est.labels_
            print labels
            labelsVariousK.append(labels)
            # save the labels as a file
            labelsFile = open(classification + 'labels_%s_%s_%s%s.txt' %(realPatientID, realSide, name, testName), 'w')
            labelsFile.write('labels\t' + str(labels) + '\n')
            labelsFile.close()
            # TODO! write out other parameters of the estimator!
            
            #inertia = est.inertia_
            #inertiaVariousK.append(inertia)
            
            # take the initial volume with columns and colour it according to the labels
            volColumns = aims.read(pathToColumnResults + '/column_regions/traverses_cortex_only_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide))        
            arrColumns = np.array(volColumns, copy = False)
            maxID = np.max(arrColumns)
            # compare the lengths
            print len(labels), len(data)
            print 'num of unique labels ', len(np.unique(arrColumns))
            
            # colour the initial columns volume into labels, and zero
            # to avoid the same colours:
            labelsUpdated = [((i + 1) * 10 + maxID) for i in labels]       
            #print zip(arrIDs, labels, labelsUpdated)
            
            
            # TODO : write out average profiles for each of the found clusters
            meanClusters = [[0 for x in range(numOfLayersToUse)] for x in range(est.n_clusters)]   # list of average profiles for columns belonging to a particular cluster
            stdvClusters = [[0 for x in range(numOfLayersToUse)] for x in range(est.n_clusters)]
            numColumnsInClusters = [0 for x in range(est.n_clusters)]
            numOfColumnsProcessed = 0
            
            meanClustersIndProf = [[0 for x in range(numOfLayersToUse)] for x in range(est.n_clusters)]   # list of average profiles for columns belonging to a particular cluster
            stdvClustersIndProf = [[0 for x in range(numOfLayersToUse)] for x in range(est.n_clusters)]
            numProfilesInClustersIndProf = [0 for x in range(est.n_clusters)]
            
            for iD, newLabel, newLabelUpdated in zip(arrIDs, labels, labelsUpdated):
                #numToReplace = len(np.where(arrColumns == iD)[0])
                #print iD, len(arrColumns[arrColumns == iD]), 'replaced by ',  newLabel, ' newLabelUpdated ', newLabelUpdated, ' numOfPoints ', numToReplace
                arrColumns[arrColumns == iD] = newLabelUpdated
                #print 'num of unique labels ', len(np.unique(arrColumns))
                
                # add the mean of the particular profile to the cluster mean, according to the label
                #print 'add to cluster ', newLabel, ' using ID ', iD, ' mean vector ', data[iD][0:numOfLayersToUse]
                #print 'add to cluster ', newLabel
                
                # find this column iD in the list of columns: arrIDs
                idInList = np.where(np.array(arrIDs) == iD)[0][0]   
                #print ' using ID ', iD, ' itsiDinList is ', idInList
                meanCurrProfile = data[idInList][0:numOfLayersToUse]
                stdvCurrProfile = data[idInList][numOfLayersToUse:2*numOfLayersToUse]
                
                #print ' mean vector to add ', meanCurrProfile
                numColumnsInClusters[newLabel] = numColumnsInClusters[newLabel] + 1
                meanClusters[newLabel] = map(add, meanClusters[newLabel], meanCurrProfile)
                stdvClusters[newLabel] = map(add, stdvClusters[newLabel],stdvCurrProfile)                
                #print 'updated meanClusters is ' , meanClusters
                numOfColumnsProcessed = numOfColumnsProcessed + 1
                #print 'numOfColumnsProcessed ', numOfColumnsProcessed
                
                ## TODO: need it ?? now calculate other cluster mean: from individual profiles
                #sumOfIndProfiles = 
                #meanClustersIndProf[newLabel] = map(add, meanClustersIndProf[newLabel], sumOfIndProfiles)
                
                
            meanClusters = [np.array(meanClusters)[xx]/numColumnsInClusters[xx] for xx in range(est.n_clusters)]
            stdvClusters = [np.array(stdvClusters)[xx]/numColumnsInClusters[xx] for xx in range(est.n_clusters)]
            # plot and save
            for xx in range(est.n_clusters):
                xCoords_new = [xxx + ((est.n_clusters/2 - xx))*0.005 for xxx in xCoords]
                #plt.errorbar(xCoords, meanClusters[xx], stdvClusters[xx], linestyle='solid', marker='^') #, ecolor = 'r', color = 'r')
                plt.errorbar(xCoords_new, meanClusters[xx], stdvClusters[xx], linestyle='solid', marker='^') # to shift points, to see the stdv
                            
            plt.title('Averages for clusters of average profile in cortical columns')   # subplot 211 title
            plt.xlabel('Cortical depth')
            plt.ylabel('T2-nobias intensity')
            plt.savefig(classification + 'AvgInClusters_AvgProfiles_%s_%s_%s%s.png' %(realPatientID, realSide, name, testName))  
            # save also to outer healthy/dyslexics folder
            plt.savefig(outerDir + 'AvgInClusters_AvgProfiles_%s_%s_%s%s.png' %(realPatientID, realSide, name, testName))  
            plt.clf()
            plt.close()
                
            # colour all columns that have not been labeled as zero
            arrColumns[arrColumns <= maxID] = 0        
            print 'num of unique labels ', len(np.unique(arrColumns))
            # save this new 'colouredVolume'        
            aims.write(volColumns, classification + 'corticalColumns_labeled_%s_%s_%s%s.nii.gz' %(realPatientID, realSide, name, testName))        
            
            
            
        
        #print labelsVariousK
        #print inertiaVariousK
        #[401337.82453829556, 306421.28899320849, 268000.21592927747]
        
        
    f.close()
    
    #sys.exit(0)   
      
      
      
      
    ###################################################################################################
      
      
    #fignum = 1
    #for name, est in estimators.items():
        #fig = plt.figure(fignum, figsize=(4, 3))
        #plt.clf()
        #ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

        #plt.cla()
        #est.fit(X)
        #labels = est.labels_

        #ax.scatter(X[:, 3], X[:, 0], X[:, 2], c=labels.astype(np.float))

        #ax.w_xaxis.set_ticklabels([])
        #ax.w_yaxis.set_ticklabels([])
        #ax.w_zaxis.set_ticklabels([])
        #ax.set_xlabel('Petal width')
        #ax.set_ylabel('Sepal length')
        #ax.set_zlabel('Petal length')
        #fignum = fignum + 1

    ## Plot the ground truth
    #fig = plt.figure(fignum, figsize=(4, 3))
    #plt.clf()
    #ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

    #plt.cla()

    #for name, label in [('Setosa', 0),
                        #('Versicolour', 1),
                        #('Virginica', 2)]:
        #ax.text3D(X[y == label, 3].mean(),
                #X[y == label, 0].mean() + 1.5,
                #X[y == label, 2].mean(), name,
                #horizontalalignment='center',
                #bbox=dict(alpha=.5, edgecolor='w', facecolor='w'))
    ## Reorder the labels to have colors matching the cluster results
    #y = np.choose(y, [1, 2, 0]).astype(np.float)
    #ax.scatter(X[:, 3], X[:, 0], X[:, 2], c=y)

    #ax.w_xaxis.set_ticklabels([])
    #ax.w_yaxis.set_ticklabels([])
    #ax.w_zaxis.set_ticklabels([])
    #ax.set_xlabel('Petal width')
    #ax.set_ylabel('Sepal length')
    #ax.set_zlabel('Petal length')
    #plt.show()




    
    
     
    
    
    
    
    
    ## which test?
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ## write out plots for healthy subjects
    ## for various mask ROIs
    #maskROIs = ['11', '21']
    #colours = ['b', 'g'] # , 'r', 'c', 'm', 'y', 'b']
    #for m, col in zip(maskROIs, colours):
        ## for R and L 
        #for realSide in ['L', 'R']:
            #healthyL_ROI11 = plt.figure(figsize=(28, 18)) #, dpi=80, facecolor='w', edgecolor='k')
            ##healthyR_ROI11 = plt.figure(figsize=(28, 18)) #, dpi=80, facecolor='w', edgecolor='k')
            #num = 1
            #for realPatientID in healthyList:
                #print '----------------- subject  ', realPatientID
                #pathToProfL11 = glob.glob(directory + '%s/diam%s/%s_%s_MaskROI%s_profiles_diam%s_over_%s.txt' %(realPatientID, str(columnDiameter), realPatientID, realSide, m, (columnDiameter), str(threshold)))            
                ##read in these files
                #numbROIsL, coordROIsL, valueROIsL = np.loadtxt(pathToProfL11[0], skiprows = 1, unpack = True)
                ## L
                #ax1 = healthyL_ROI11.add_subplot(3,3,num)
                #ax1.set_title('Profile in all %s mask ROI %s - %s' %(realPatientID, m, realSide))   # subplot 211 title
                #ax1.set_xlabel('Cortical depth')
                #ax1.set_ylabel('T2-nobias intensity')
                #ax1.plot(coordROIsL, valueROIsL, '.', c = col, label = realSide)
                ##ax1.legend(loc='upper right', numpoints = 1)
                #num += 1        
            #plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/Data_diam%s_over%s/healthy_%s_maskROI%s.png' %(str(columnDiameter), str(threshold), realSide, m), bbox_inches='tight')    
            #plt.clf()
            #plt.close()
    
            ## the same for dyslexics
            #healthyL_ROI11 = plt.figure(figsize=(28, 18)) #, dpi=80, facecolor='w', edgecolor='k')
            ##healthyR_ROI11 = plt.figure(figsize=(28, 18)) #, dpi=80, facecolor='w', edgecolor='k')
            #num = 1
            #for realPatientID in dyslexicList:
                #print '----------------- subject  ', realPatientID
                #pathToProfL11 = glob.glob(directory + '%s/diam%s/%s_%s_MaskROI%s_profiles_diam%s_over_%s.txt' %(realPatientID, str(columnDiameter), realPatientID, realSide, m, (columnDiameter), str(threshold)))            
                ##read in these files
                #numbROIsL, coordROIsL, valueROIsL = np.loadtxt(pathToProfL11[0], skiprows = 1, unpack = True)
                ## L
                #ax1 = healthyL_ROI11.add_subplot(3,3,num)
                #ax1.set_title('Profile in all %s mask ROI %s - %s' %(realPatientID, m, realSide))   # subplot 211 title
                #ax1.set_xlabel('Cortical depth')
                #ax1.set_ylabel('T2-nobias intensity')
                #ax1.plot(coordROIsL, valueROIsL, '.', c = col, label = realSide)
                ##ax1.legend(loc='upper right', numpoints = 1)
                #num += 1        
            #plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/Data_diam%s_over%s/dyslexic_%s_maskROI%s.png' %(str(columnDiameter), str(threshold), realSide, m), bbox_inches='tight')    
            #plt.clf()
            #plt.close()
    

    
    ###### try to change the order!!
    ## write out plots for healthy subjects
    ## for various mask ROIs
    #maskROIs = ['11', '21']
    #patientsLists = [healthyList, dyslexicList]
    #keywords = ['healthy', 'dyslexic']
    
    #colours = ['b', 'g'] # , 'r', 'c', 'm', 'y', 'b']    
    ## for R and L 
    #for realSide in ['L', 'R']:
        #for listt, keyword in zip(patientsLists, keywords):            
            #healthyL_allROIs = plt.figure(figsize=(28, 18)) #, dpi=80, facecolor='w', edgecolor='k')
            #num = 1
            #for realPatientID in listt:
                #print '----------------- subject  ', realPatientID
                #ax1 = healthyL_allROIs.add_subplot(3,3,num)
                #ax1.set_title('Profile in %s all mask ROIs - %s' %(realPatientID, realSide))   # subplot 211 title
                #ax1.set_xlabel('Cortical depth')
                #ax1.set_ylabel('T2-nobias intensity')
                #for m, col in zip(maskROIs, colours):
                    #pathToProfL11 = glob.glob(directory + '%s/diam%s/%s_%s_MaskROI%s_profiles_diam%s_over_%s.txt' %(realPatientID, str(columnDiameter), realPatientID, realSide, m, (columnDiameter), str(threshold)))            
                    ##read in these files
                    #numbROIsL, coordROIsL, valueROIsL = np.loadtxt(pathToProfL11[0], skiprows = 1, unpack = True)                    
                    #ax1.plot(coordROIsL, valueROIsL, '.', c = col, label = realSide)
                    ##ax1.legend(loc='upper right', numpoints = 1)
                #num += 1        
            #plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/Data_diam%s_over%s/%s_%s_allMaskROIs.png' %(str(columnDiameter), str(threshold), keyword, realSide), bbox_inches='tight')    
            #plt.clf()
            #plt.close()
                
         
         
         



                             
        
        
        
        
        
        
    #pathToProfL = directory + '%s_L_profiles2.txt' %(realPatientID)
    #pathToProfR = directory + '%s_R_profiles2.txt' %(realPatientID)
    
    ## check if both these profiles exist
    #profL = glob.glob(pathToProfL)
    #profR = glob.glob(pathToProfR)
    
    #if len(profL) != 1 or len(profR) != 1:
        ## abort the calculation, as too many or not a single texture file was found
        #f = open(directory + '%s_compareProfilesStat.txt' %(realPatientID), "w")
        #print 'abort the calculation, as too many or not a single profL or R file was found'
        #f.write('abort the calculation, as ' + str(len(profL)) + ' profL and ' + str(len(profR)) + ' profR profile files were found' + '\n')
        #f.close()
        #sys.exit(0) 
    
    
    #numbL, coordL, valueL = np.loadtxt(pathToProfL, skiprows = 1, unpack = True)
    #numbR, coordR, valueR = np.loadtxt(pathToProfR, skiprows = 1, unpack = True)
       
    ## plot the data
    #plt.plot(coordL, valueL, '.', c = 'b', label = 'L')
    #plt.title('Profile in all ROIs')   # subplot 211 title
    #plt.xlabel('Cortical depth')
    #plt.ylabel('T2-nobias intensity')
    #plt.plot(coordR, valueR, '.', c = 'r', label = 'R')
    #plt.legend(loc='upper right', numpoints = 1)
    #plt.savefig(directory + '%s_LvsR_2nobiasT2.png' %(realPatientID))    
    
    ## TODO: delete later if no need: save profiles also to the outer folder!
    ##plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/' + '%s_LvsR_2nobiasT2.png' %(realPatientID))    
    #plt.clf()
    #plt.close()
    
    ## now plot L vs R in various ROIs
    ## af140169_R_profiles2_ROI_21.txt, af140169_R_profiles2_ROI_11.txt and the same with L
    #pathToROIsProfL = directory + '%s_L_profiles2_ROI_[0-9]*.txt' %(realPatientID)
    ##print 'pathToROIsProfL'
    ##print pathToROIsProfL

    #pathToROIsProfR = directory + '%s_R_profiles2_ROI_[0-9]*.txt' %(realPatientID)
    
    ## check if both these profiles exist. and how many are there
    #profROIsL = glob.glob(pathToROIsProfL)
    #profROIsR = glob.glob(pathToROIsProfR)
    
    #print 'profROIsL'
    #print profROIsL
    #print 'profROIsR'
    #print profROIsR
    #numOfCommonRLregions = 0
    
    #for i in range(len(profROIsL)):
        ## check if this L profile is from the same mask ROI as the right one:
        ## check if the corresponding R - file exists
        #profROIsRcorrespond = glob.glob(profROIsL[i].replace('_L_', '_R_'))
        #if len(profROIsRcorrespond) == 1:      
            #numOfCommonRLregions += 1
            #print 'corresponding files ', profROIsL[i], ' and ', profROIsRcorrespond[0]
            ## read in the respective files
            #numbROIsL, coordROIsL, valueROIsL = np.loadtxt(profROIsL[i], skiprows = 1, unpack = True)
            #numbROIsR, coordROIsR, valueROIsR = np.loadtxt(profROIsRcorrespond[0], skiprows = 1, unpack = True)
            
            ## get the ROI ID
            #iD = (profROIsL[i].split('_L_profiles2_ROI_')[1]).split('.txt')[0] # '/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/ml140175/ml140175_R_profiles2_ROI_11.txt'
            ## plot the data
            #plt.plot(coordROIsL, valueROIsL, '.', c = 'b', label = 'L')
            #plt.title('Profile in ROI %s ' %(iD))   # subplot 211 title
            #plt.xlabel('Cortical depth')
            #plt.ylabel('T2-nobias intensity')
            #plt.plot(coordROIsR, valueROIsR, '.', c = 'r', label = 'R')
            #plt.legend(loc='upper right', numpoints = 1)
            #plt.savefig(directory + '%s_LvsR_2nobiasT2_ROI_%s.png' %(realPatientID, iD))   
            
            ## TODO: delete later if no need: save profiles also to the outer folder!
            ##plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/' + '%s_LvsR_2nobiasT2_ROI_%s.png' %(realPatientID, iD))    
            #plt.clf()
            #plt.close() 
    #print 'found ', numOfCommonRLregions, ' numOfCommonRLregions '       
            
       
    ## plot these plots into 1 image
    #fig = plt.figure(figsize=(21, 6)) #, dpi=80, facecolor='w', edgecolor='k')
    #numOfCommonRLregions += 1
    #ax1 = fig.add_subplot(1,numOfCommonRLregions,1)
    #ax1.plot(coordL, valueL, '.', c = 'b', label = 'L')
    #ax1.set_title('Profile in all ROIs')   # subplot 211 title
    #ax1.set_xlabel('Cortical depth')
    #ax1.set_ylabel('T2-nobias intensity')
    #ax1.plot(coordR, valueR, '.', c = 'r', label = 'R')
    #ax1.legend(loc='upper right', numpoints = 1)
   
    #for i in range(len(profROIsL)):
        #profROIsRcorrespond = glob.glob(profROIsL[i].replace('_L_', '_R_'))
        #if len(profROIsRcorrespond) == 1: 
            ## read in the respective files
            #numbROIsL, coordROIsL, valueROIsL = np.loadtxt(profROIsL[i], skiprows = 1, unpack = True)
            #numbROIsR, coordROIsR, valueROIsR = np.loadtxt(profROIsRcorrespond[0], skiprows = 1, unpack = True)
            
            ## get the ROI ID
            #iD = (profROIsL[i].split('_L_profiles2_ROI_')[1]).split('.txt')[0]
            #ax2 = fig.add_subplot(1,numOfCommonRLregions, 2 + i)
            #ax2.plot(coordROIsL, valueROIsL, '.', c = 'b', label = 'L')
            #ax2.set_title('Profile in ROI %s' %(iD))   # subplot 211 title
            #ax2.set_xlabel('Cortical depth')
            #ax2.set_ylabel('T2-nobias intensity')
            #ax2.plot(coordROIsR, valueROIsR, '.', c = 'r', label = 'R')
            #ax2.legend(loc='upper right', numpoints = 1)

    ##plt.show()
    #print 'save the image /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/%s_LvsR_allVsROIs.png' %(realPatientID)
    #plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/' + '%s_LvsR_allVsROIs.png' %(realPatientID), bbox_inches='tight')    
    #plt.clf()
    #plt.close()
    #print 'saved the image /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/%s_LvsR_allVsROIs.png' %(realPatientID)
    
    
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
        
    
    
    

    

      
  
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    