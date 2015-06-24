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


# example how to run this file:
#python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_plotHealthy_DyslexicProfiles.py -c 5 -t 314 -d /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/
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

if __name__ == '__main__':
    
    healthyList = ['fg140290', 'af140169', 'ml140175', 'ac140159', 'md140208', 'at140353'] #, 'js140266', 'lg140146', 'he140338', 'cb140330']
    dyslexicList = ['js140311', 'ad140157', 'ag140439', 'sg140335']
    realPatientID = None
    directory = None
    threshold = None
    #realSide = 'L'
    columnDiameter = None
    

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
    #parser.add_option('-t', dest='threshold', help='threshold to work with')
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
        
    if options.realSide is None:
        print >> sys.stderr, 'New: exit. no realSide given'
        sys.exit(1)
    else:
        realSide = options.realSide  
        
        
    addedColumnsDiamName = ''
    pathForFiles = directory
    
    if options.columnDiameter is not None:
        columnDiameter = int(options.columnDiameter)
        addedColumnsDiamName = '_diam%s' %(columnDiameter)
        # if the result was for the columns, find a folder for it
        pathForFiles = pathForFiles + 'diam%s/'%(columnDiameter)
        print '############################################ found profile directory = ', pathForFiles


    sizeSorted = pathForFiles + 'sortedBySize/'
    divGradnSorted = pathForFiles + 'sortedByDivGradn/'
    divGradnClasses = directory + 'divGradnClasses/'
    corticalLayers = pathForFiles + 'corticalLayers/'

    if not os.path.exists(sizeSorted):
        os.makedirs(sizeSorted)
    if not os.path.exists(divGradnSorted):
        os.makedirs(divGradnSorted)
    if not os.path.exists(divGradnClasses):
        os.makedirs(divGradnClasses)
    if not os.path.exists(corticalLayers):
        os.makedirs(corticalLayers)

    #if options.threshold is None:
        #print >> sys.stderr, 'New: exit. no threshold given'
        #sys.exit(1)
    #else:
        #threshold = options.threshold  
        
    # TODO: for test. maybe delete it later
    # plot histograms for the intensities in CSF, cortex and WM
    pathToClassifNoBorders = directory + '%s_T1inT2_ColumnsCutNew20It/GWsegm_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realPatientID, realSide)
    volGWNoborders = aims.read(pathToClassifNoBorders)
    arrGWNoborders = np.array(volGWNoborders)
    pathToClassifWithBorders = directory + '%s_T1inT2_ColumnsCutNew20It/dist/classif_with_outer_boundaries_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realPatientID, realSide)
    volGWborders = aims.read(pathToClassifWithBorders)
    arrGWborders = np.array(volGWborders)
    pathToNobiasT2_new = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S16/%s_NewT2_cropped.nii.gz' %(realPatientID)
    volT2 = aims.read(pathToNobiasT2_new)
    arrT2 = np.array(volT2)
    
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
    infoFileName = pathForFiles + '%s_%s_ColumnInfo%s.txt' %(realPatientID, realSide, addedColumnsDiamName)
    columnID, size  = np.loadtxt(infoFileName, skiprows = 1, usecols = (0, 1), unpack=True)
    print zip(columnID, size)
    
    # calculatethe size threshold for the columns
    sizeThreshold = np.round(columnDiameter * columnDiameter * 4.0 / 4.0 * np.pi * 4.0 / 4.0) #(volume of a cylinder, diameter is given, resolution 0.5, height - 2mm minimum, and divided by 4 for conical cases, or too narrow cases)
    print 'sizeThreshold = ', sizeThreshold
    
    # find columns larger than this threshold
    largeIDs = columnID[np.where(size >= sizeThreshold)]
    print 'largeIDs = ', largeIDs
    completeN = 0
    incompleteN = 0
    complete5LayersN = 0    
    
    # create data structure for the collected data! the first column: ROI ID !!!!!!!!!!
    arrProfiles6Layers = []
    arrProfilesID6Layers = []
    arrVar6Layers = []
    arrProfiles5Layers = []
    arrProfilesID5Layers = []
    arrVar5Layers = []
    
    for i in range(len(largeIDs)):
#    for i in range(1):
        # iterate through these profiles, calculate means, stdv in layers, plot figures, save files
        # find profiles of these large columns: ad140157_L_profiles2_diam3_ROI_259.txt - find all profile files like this one
        currID = int(largeIDs[i])
        #print 'currID = ', currID
        profileFile = pathForFiles + '%s_%s_profiles2%s_ROI_%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(currID))
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


        # collected the data for all layers. now plot for this particular column  
        axMeans.set_title('Mean and stdv values in cortical layers in maskROI %s' %(str(currID)))  
        axMeans.set_xlabel('Cortical depth') 
        axMeans.set_ylabel('T2-nobias intensity')                 
        axMeans.legend(loc='upper right', numpoints = 1)  
        
        if complete:
            # plot the profile for all 6 cortical layers
            axMeans.errorbar(xCoords, means, stdvs, linestyle='solid', marker='^', ecolor = 'r')
            plt.savefig(corticalLayers + '%s_%s_completeCortLay%s_size%s_ROI_%s%s.png' %(realPatientID, realSide, addedColumnsDiamName, str(len(currCoords)), str(currID), addedName))                
#                plt.savefig(sizeSorted + 'complete%s_size%s_ROI_' %(addedColumnsDiamName, str(len(currCoords))) + str(iDs[i]) + '.png')
            # save txt files with means and stdvs !!
            dataCort = open(corticalLayers + '%s_%s_CortLay%s_ROI_%s%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(currID), addedName), "w")
#                dataCort = open(pathForFiles + '%s_%s_CorticalLayers%s_ROI_%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(iDs[i])), "w")
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
            arrProfilesID5Layers.append(currID)
            arrProfiles5Layers.append(means[1:])  
            arrVar6Layers.append(stdvs)   
            arrVar5Layers.append(stdvs) 
        else :
            # now check, whether this profile is complete for 5 layers (II - VI)
            if complete5Layers == False:
                # this profile is incomplete even for 5 layers
                #print 'incomplete ID , even for 5 layers ', str(currID)
                axMeans.errorbar(xCoords, means, stdvs, linestyle='None', marker='^', ecolor = 'r')
                plt.savefig(corticalLayers + '%s_%s_incompleteCortLay%s_size%s_ROI_%s%s.png' %(realPatientID, realSide, addedColumnsDiamName, str(len(currCoords)), str(currID), addedName))  
                # do NOT save separate txt files, as the data is incomplete and can not be used for clustering
                incompleteN += 1
            else:   # this profile is incomplete for many points, but is complete for 5 layers
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
#            plt.savefig(pathForFiles + '%s_%s_nobiasT2%s_ROI_' %(realPatientID, realSide, addedColumnsDiamName) + str(iDs[i]) + '.png')        
            # write out the same info, but sorted by sizes
#            plt.savefig(sizeSorted + '%s_%s_nobiasT2%s_size%s_ROI_' %(realPatientID, realSide, addedColumnsDiamName, str(len(currCoords))) + str(iDs[i]) + '.png')
            
            ## now save the same plots but sorted by their avg gradients              
#            plt.savefig(divGradnSorted + '%s_%s_nobiasT2%s_avgDivGradn%s_ROI_%s.png' %(realPatientID, realSide, addedColumnsDiamName, strCurrGradnStr, str(iDs[i])))
            
        plt.clf()
        plt.close()     
    print '*********************************************** complete = ', completeN, ', incomplete = ', incompleteN, ', complete5LayersN = ', complete5LayersN 
  
    
######################################################## STATISTICAL ANALYSIS ####################################################################################

    # 1. Cluster the columns using vectors of their means and stdvs in the 6 'cortical layers'
    print (arrProfiles5Layers[0:10][:])   
    
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


    # data generation
    data = arrProfiles6Layers
    
    # computing K-Means with K = 2 (2 clusters)
    centroids,_ = kmeans(data, 2)
    
    # assign each sample to a cluster
    idx,_ = vq(data,centroids)

    # some plotting using numpy's logical indexing
    plot(data[idx==0,0],data[idx==0,1],'ob',
        data[idx==1,0],data[idx==1,1],'or')
    plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
    show()

    print 
    
    sys.exit()






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
    estimators = {'k_means_2': KMeans(n_clusters = 2),
                  'k_means_3': KMeans(n_clusters = 3),
                  'k_means_4': KMeans(n_clusters = 4)}
    
    labelsVariousK = []
    # take the initial volume with columns and colour it according to the labels
    
    
    for name, est in estimators.items():
        est.fit(arrProfiles6Layers)
        labels = est.labels_
        print labels
        labelsVariousK.append(labels)
        
        # take the initial volume with columns and colour it according to the labels
        volColumns = aims.read('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/at140353/at140353_T1inT2_ColumnsCutNew20It/traverses_cortex_only.nii.gz')
        arrColumns = np.array(volColumns, copy = False)
        # compare the lengths
        print len(labels)
        print len(arrProfilesID6Layers)
        
        
        for iD, newLabel in zip(arrProfilesID6Layers, labels):
            print iD, len(arrColumns[arrColumns == iD]), 'replaced by ',  newLabel
            arrColumns[arrColumns == iD] = newLabel
        # save this new 'colouredVolume'
        aims.write(volColumns, '/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/at140353/at140353_T1inT2_ColumnsCutNew20It/corticalColumns_labeled_%s.nii.gz' %(name))
        
        
        
        
      
      
      
    fignum = 1
    for name, est in estimators.items():
        fig = plt.figure(fignum, figsize=(4, 3))
        plt.clf()
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

        plt.cla()
        est.fit(X)
        labels = est.labels_

        ax.scatter(X[:, 3], X[:, 0], X[:, 2], c=labels.astype(np.float))

        ax.w_xaxis.set_ticklabels([])
        ax.w_yaxis.set_ticklabels([])
        ax.w_zaxis.set_ticklabels([])
        ax.set_xlabel('Petal width')
        ax.set_ylabel('Sepal length')
        ax.set_zlabel('Petal length')
        fignum = fignum + 1

    # Plot the ground truth
    fig = plt.figure(fignum, figsize=(4, 3))
    plt.clf()
    ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

    plt.cla()

    for name, label in [('Setosa', 0),
                        ('Versicolour', 1),
                        ('Virginica', 2)]:
        ax.text3D(X[y == label, 3].mean(),
                X[y == label, 0].mean() + 1.5,
                X[y == label, 2].mean(), name,
                horizontalalignment='center',
                bbox=dict(alpha=.5, edgecolor='w', facecolor='w'))
    # Reorder the labels to have colors matching the cluster results
    y = np.choose(y, [1, 2, 0]).astype(np.float)
    ax.scatter(X[:, 3], X[:, 0], X[:, 2], c=y)

    ax.w_xaxis.set_ticklabels([])
    ax.w_yaxis.set_ticklabels([])
    ax.w_zaxis.set_ticklabels([])
    ax.set_xlabel('Petal width')
    ax.set_ylabel('Sepal length')
    ax.set_zlabel('Petal length')
    plt.show()




    
    
     
    
    
    
    
    
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
        
    
    
    

    

      
  
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    