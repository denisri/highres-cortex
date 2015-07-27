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

def extractProfiles(volCoord, volValue, volMask = None):
    """
    this function takes 2 volumes, one with values (T2) and another with corresponding coordinates (cortex
    depth measure from Yann's scripts.) It also optionally takes a mask in which the profiles will be extracted.
    If no mask was given, profiles are extracted in the whole volume.
    Forfiles are extracted as a list of values and a list of coordinates.    
    """
    #print volCoord.header()
    
    arrCoord = np.asarray(volCoord)
    arrValue = np.asarray(volValue)

    # apply the mask if it was given
    listOfSeparateCoords = []
    listOfSeparateValues = []
   
    if volMask is not None:
        mask = np.asarray(volMask)        
        # get these ROIs
        arrCoord1 = arrCoord[mask != 0]
        arrValue1 = arrValue[mask != 0]
        
        coords = arrCoord1[arrCoord1 != 0]
        values = arrValue1[arrCoord1 != 0]
        
        ids = np.where(arrCoord1 != 0)
        print len(ids), ' len(ids) ', len(coords), ' len(coords) '
        # got the ids. now need to get x, y coordinates
        #print volCoord.getX()[1:10]
        
        ########################## now need to extract profiles in ROIs separately
        # get unique values for the mask
        roiIds = np.unique(mask[np.where(mask > 0)])
        for i in roiIds:
            if i != 0:
                print 'work with ROI ', i
                arrCoord1i = arrCoord[mask == i]
                arrValue1i = arrValue[mask == i]
                
                coordsi = arrCoord1i[arrCoord1i != 0]
                valuesi = arrValue1i[arrCoord1i != 0]
                
                listOfSeparateCoords.append(coordsi)
                listOfSeparateValues.append(valuesi)      
    else :
        coords = arrCoord[arrCoord != 0]
        values = arrValue[arrCoord != 0]
        
    # 2 arrays of coordinates and values. Plot them
    res = []
    res.append(coords)
    res.append(values)
    res.append(roiIds)
    res.append(listOfSeparateCoords)
    res.append(listOfSeparateValues)
    return(res)       
    

    
def extractProfilesInColumns(volCoord, volValue, volColumns, volDivGradn, divGrThr, volMask = None):
    """
    this function takes 2 volumes, one with values (T2) and another with corresponding coordinates (cortex
    depth measure from Yann's scripts.) It also takes a volume with "cortical columns" and will extract 
    lists of profiles for each of the columns separately.
    It also optionally takes a mask in which the profiles will be extracted.
    If no mask was given, profiles are extracted in the whole volume.
    Forfiles are extracted as a list of values and a list of coordinates. 
    
    modified: now also require the div_gradn image, to check whether a particular column is
    situated on a flat PT region (mean value of the div_gradn should be around zero)
    or on a curved region (mean value of the div_gradn should >> 0 or << 0)
    """
    #print volCoord.header()
    
    #TODO: check if this is correct. Attention: modifying the original data!!
    #arrColouredVol = np.array(volColumns, copy = True)
    #print 'len(arrColouredVol)= ', len(arrColouredVol)
    #arrColouredVol1 = arrColouredVol[mask != 0]
    #arrColouredVol11 = arrColouredVol1[arrCoord1 != 0]
    volColoured = aims.Volume(volColumns.getSizeX(), volColumns.getSizeY(), volColumns.getSizeZ(), volColumns.getSizeT(), 'S16')
    volColoured.header().update(volColumns.header())
    volColoured.fill(0)
    print 'volColoured.header() = ', volColoured.header()
    arrColouredVol = np.array(volColoured, copy = False)
        
    arrCoord = np.asarray(volCoord)
    arrValue = np.asarray(volValue)
    arrColumns = np.asarray(volColumns)
    arrDivGradn = np.asarray(volDivGradn)

    # apply the mask if it was given
    listOfSeparateCoords = []
    listOfSeparateValues = []
    # get columns in certain mask regions    
    listOfSeparateCoordsMask = []
    listOfSeparateValuesMask = []
    listOfROIsInMask = []
    listOfSizes = []
    listOfSeparateAVGCoords = []        # AVG coordinates of individual columns (to be used in clustering)
    strMeansInfo = ''
   
    if volMask is not None:
        mask = np.asarray(volMask)    
        
        # get the part (ROIs) of the volColumns where the mask is not zero
        roiColumns = arrColumns[mask != 0]
        
        # get these ROIs
        arrCoord1 = arrCoord[mask != 0]
        arrValue1 = arrValue[mask != 0]
        
        ## study T2 intensities for data cleaning        
        #arrValue1WM = arrValue[mask == 200] # get values in the WM
        #arrValue1Cortex = arrValue[mask == 100] # get values in the cortex
        #arrValue1CSF = arrValue[mask == 0] # get values in the CSF
        #fig = plt.figure(figsize=(3 * 7, 6)) #, dpi=80, facecolor='w', edgecolor='k')
        ## and now plot them
        #axWM = fig.add_subplot(1,3,1)
        #axWM.hist(arrValue1WM, bins = 50)
        #axWM.title('Histogram of T2 intensities in WM')   # subplot 211 title
        #axCortex = fig.add_subplot(1,3,2)
        #axCortex.hist(arrValue1Cortex, bins = 50)
        #axCortex.title('Histogram of T2 intensities in cortex')   # subplot 211 title
        #axCSF = fig.add_subplot(1,3,3)
        #axCSF.hist(arrValue1CSF, bins = 50)
        #axCSF.title('Histogram of T2 intensities in CSF')   # subplot 211 title
        ##plt.savefig(pathToColumnResults + '%s_%s_histOfColumnSizes%s.png' %(realPatientID, realSide, addedColumnsDiamName))    
        #plt.savefig('/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/ac140159/.png' %(realPatientID, realSide, addedColumnsDiamName))    
        #plt.clf()
        #plt.close()
        
        
        
        
        
        # TODO! before all this: ensure to work on zero-free data!!!!!!!!!!!!!!!
        # TODO! scale the values! due to different acquisition settings, can not compare profiles among subjects! -> need to scale
        # important! need to transform them to float before!!!
        strMeansInfo = 'Min, Mean, SD, Median, Max of unscaled T2 data: \n'
        print 'Min, Mean, SD, Median, Max of unscaled T2 data: \n'
        arrValue1 = arrValue1.astype('float')
        strMeansInfo = strMeansInfo + str(np.min(arrValue1)) + '\t' + str(np.mean(arrValue1)) + '\t' + str(np.std(arrValue1)) + '\t' + str(np.median(arrValue1)) + '\t' + str(np.max(arrValue1)) + '\n'
        print '-----------------------------------------', strMeansInfo
        
        # TODO! before scaling: eliminate 'outliers'
        # consider values lower than 1st percentile and greater than 99th percentile: outliers!
        # exclude them and perform scaling without them
        percentile1 = np.percentile(arrValue1, 1)
        percentile99 = np.percentile(arrValue1, 99)
        print 'computed percentile1 = ', percentile1, ' percentile99 = ', percentile99
        arrValue1OutliersSmall = arrValue1[arrValue1 <= percentile1]
        arrValue1OutliersBig = arrValue1[arrValue1 >= percentile99]
        arrValue1NoOutliers = arrValue1[(arrValue1 > percentile1) & (arrValue1 < percentile99)]
        print 'num of Voxelx changed from ', len(arrValue1), ' to ', len(arrValue1NoOutliers), ' eliminated ', len(arrValue1OutliersSmall), ' and ', len(arrValue1OutliersBig)
        
        # scale the T2 data without outliers
        arrValue1NoOutliersScaled = preprocessing.scale(arrValue1NoOutliers)
        
        # for comparison also scale the data WITH outliers
        arrValue1Scaled = preprocessing.scale(arrValue1)  
        strMeansInfo = strMeansInfo + '\n Min, Mean, SD, Median, Max of (T2 scaled with outliers) data: \n'
        strMeansInfo = strMeansInfo + str(np.min(arrValue1Scaled)) + '\t' + str(np.mean(arrValue1Scaled)) + '\t' + str(np.std(arrValue1Scaled)) + '\t' + str(np.median(arrValue1Scaled)) + '\t' + str(np.max(arrValue1Scaled)) + '\n'
        print strMeansInfo
        ###-------------------
        strMeansInfo = strMeansInfo + '\n Min, Mean, SD, Median, Max of (T2 scaled without outliers) data: \n'
        strMeansInfo = strMeansInfo + str(np.min(arrValue1NoOutliersScaled)) + '\t' + str(np.mean(arrValue1NoOutliersScaled)) + '\t' + str(np.std(arrValue1NoOutliersScaled)) + '\t' + str(np.median(arrValue1NoOutliersScaled)) + '\t' + str(np.max(arrValue1NoOutliersScaled)) + '\n'
        print strMeansInfo
         
        # for further analysis ADD the eliminated outliers, applied a correction to them        
        scaler = preprocessing.StandardScaler().fit(arrValue1NoOutliers)
        print scaler
        print scaler.mean_                                      
        print scaler.std_    
        arrValue1NoOutliersScaled2 = scaler.transform(arrValue1NoOutliers)
        # apply this transformation to the outliers        
        arrValue1OutliersSmallScaled = scaler.transform(arrValue1OutliersSmall)   
        arrValue1OutliersBigScaled = scaler.transform(arrValue1OutliersBig)   
        # add this data to the array
        idxSmall = [arrValue1 <= percentile1]
        idxOK = [(arrValue1 > percentile1) & (arrValue1 < percentile99)]
        idxBig = [arrValue1 >= percentile99]
        arrValue1[idxSmall] = arrValue1OutliersSmallScaled
        arrValue1[idxBig] = arrValue1OutliersBigScaled
        arrValue1[idxOK] = arrValue1NoOutliersScaled2
        # now arrValue1 contains all the original data, scaled!!! using only non-outliers data
        strMeansInfo = strMeansInfo + '\n Min, Mean, SD, Median, Max of (T2 scaled without outliers, added scaled outliers) data: \n'
        strMeansInfo = strMeansInfo + str(np.min(arrValue1)) + '\t' + str(np.mean(arrValue1)) + '\t' + str(np.std(arrValue1)) + '\t' + str(np.median(arrValue1)) + '\t' + str(np.max(arrValue1)) + '\n'
        print strMeansInfo
        #print arrValue1NoOutliersScaled.mean(axis=0)
        #print arrValue1NoOutliersScaled.mean()
        #print arrValue1NoOutliersScaled.std(axis=0)
        #print arrValue1NoOutliersScaled.std()
        #sys.exit(0)
        
        
        #arrValue1 = preprocessing.scale(arrValue1)              
        #print '-----------------------------------------', np.min(arrValue1), np.mean(arrValue1), np.std(arrValue1), np.median(arrValue1), np.max(arrValue1)
        ##################################################################################
        arrDivGradn1 = arrDivGradn[mask != 0]
        arrColouredVol1 = np.where(mask != 0)
        print 'arrColouredVol:', arrColouredVol.shape
        print 'mask:', np.where(mask != 0)
       
        coords = arrCoord1[arrCoord1 != 0]
        values = arrValue1[arrCoord1 != 0]
        divGradns = arrDivGradn1[arrCoord1 != 0]               
        ids = np.where(arrCoord1 != 0)
        print len(ids), ' len(ids) ', len(coords), ' len(coords) '
        # got the ids. now need to get x, y coordinates
        #print volCoord.getX()[1:10]
        
        
       ########################## now need to extract profiles in ROIs separately
        # get unique values for the mask
        #roiIds = np.unique(roiColumns[np.where(roiColumns > 0)])
        roiIds = np.unique(roiColumns[np.where(roiColumns > 0)])  # rois from the cortical columns
        roisMask = np.unique(mask[np.where(mask > 0)])  # rois from the mask
        listOfROIsInMask = [[] for x in range(len(roisMask) + 1)]
        print 'listOfROIsInMask = '
        print listOfROIsInMask
        
        #roiIds = np.unique(arrColumns)
        print 'found the following ROIs. their number  ', len(roiIds), ' from total number : ', len(np.unique(arrColumns))
        #print roiIds
        listOfAvgDivGradn = []
        listOfAvgDivGradnCateg = []     # it will be a list of 'classes' for the value of the average div gradn for a partucular column
        print 'len(arrColouredVol)= ', len(arrColouredVol)
        counters=[0,0,0]
        countVoxels=[0,0,0]

        for i in roiIds:
            #print 'roi ', i, len(arrCoord), len(roiColumns)
            # TODO! check if this is OK: I take the whole column even if only 1 voxel is inside the voronoiCorr. Cut it?
            arrCoord1i = arrCoord1[roiColumns == i]
            arrValue1i = arrValue1[roiColumns == i]  
            arrDivGradn1i = arrDivGradn1[roiColumns == i]  
            
            coordsi = arrCoord1i[arrCoord1i != 0]
            valuesi = arrValue1i[arrCoord1i != 0]
            divGradnsi = arrDivGradn1i[arrCoord1i != 0]
            
            avgDiv = np.mean(divGradnsi)
            #print '---------- work with ROI ', i, ' # of points= ', len(coordsi), 'realColumnSize= ', len(np.where(arrColumns == i)[0]), ' average divGradn = ', avgDiv, ' -------'
            #if math.isnan(avgDiv):
                #print len(coordsi), len(valuesi), len(divGradnsi), ' the whole digGradnList',  divGradnsi               
                
            listOfSeparateCoords.append(coordsi)
            listOfSeparateValues.append(valuesi)   
            listOfSizes.append(len(coordsi))
            listOfAvgDivGradn.append(avgDiv)
            
            
            # TODO! get AVG X, Y coordinates per column (might use it for clustering)
            #coordsThisROI = np.where(arrCoord1i != 0)
            coordsThisROI = np.where(arrColumns == i)   # coordinates in 'voxels'. need to transform into mm!
            #print 'coordsThisROI[:] = ', coordsThisROI[:]
            #print 'coordsThisROI[0] = ', coordsThisROI[0]
            xCoords, yCoords, zCoords, tCoords = np.where(arrColumns == i)
            meanROIVoxelCoords = [np.mean(xCoords), np.mean(yCoords), np.mean(zCoords)]
            #print 'xCoords ', xCoords
            #print 'yCoords ', yCoords
            #print 'zCoords ', zCoords
            #print 'ROI ', i, meanROIVoxelCoords, ' size = ', len(xCoords)
            listOfSeparateAVGCoords.append(meanROIVoxelCoords)

            #sys.exit()
                   
            # TODO: delete it later!!! colour the initial image into 3 colours:
            # if avg < - 0.1 blue ,  if 0.1 >= avg >= - 0.1 yellow,   if avg > 0.1 red     
            #print 'len(arrColouredVol) = ', len(arrColouredVol), len(arrValue), len(arrValue1), len(arrValue1i), 'len(arrColouredVol1) = ', len(arrColouredVol1), ' len(roiColumns)= ', len(roiColumns), ' i= ', i
            
            #print 'arrColouredVol1:', arrColouredVol1
            if avgDiv < divGrThr[0]:   # colour this column blue
                #print len(np.where(roiColumns == i))
                #print len(arrColouredVol1[roiColumns == i])
                #print 'colour 10'
                arrColouredVol[arrColouredVol1[0][roiColumns == i], 
                               arrColouredVol1[1][roiColumns == i],
                               arrColouredVol1[2][roiColumns == i],
                               arrColouredVol1[3][roiColumns == i]] = 10   
                counters[0] += 1
                countVoxels[0] += len(coordsi)
                listOfAvgDivGradnCateg.append(10)
            elif avgDiv < divGrThr[1]:  # colour this column yellow
                arrColouredVol[arrColouredVol1[0][roiColumns == i], 
                               arrColouredVol1[1][roiColumns == i],
                               arrColouredVol1[2][roiColumns == i],
                               arrColouredVol1[3][roiColumns == i]] = 20  
                counters[1] += 1
                countVoxels[1] += len(coordsi)
                listOfAvgDivGradnCateg.append(20)
            else:       # colour this column red
                arrColouredVol[arrColouredVol1[0][roiColumns == i], 
                               arrColouredVol1[1][roiColumns == i],
                               arrColouredVol1[2][roiColumns == i],
                               arrColouredVol1[3][roiColumns == i]] = 30
                #print 'colour 30' 
                counters[2] += 1
                countVoxels[2] += len(coordsi)
                listOfAvgDivGradnCateg.append(30)
        
        # save this new colouring scheme
        newVol = aims.Volume(arrColouredVol)
        print 'lower ', divGrThr[0], ' : ', counters[0]/float(len(roiIds)), ', lower ', divGrThr[1], ' : ', counters[1]/float(len(roiIds)), ', over ', divGrThr[1], ' : ', counters[2]/float(len(roiIds))   
        print 'lower ', divGrThr[0], ' : ', countVoxels[0]/float(sum(countVoxels)), ', lower ', divGrThr[1], ' : ', countVoxels[1]/float(sum(countVoxels)), ', over ', divGrThr[1], ' : ', countVoxels[2]/float(sum(countVoxels))
        
        # get IDs of columns in various mask regions
        #print 'len(roisMask) ', len(roisMask)
        for j in range(len(roisMask)):
            #print 'roisMask[j] ', roisMask[j]
            partOfColumns = arrColumns[mask == roisMask[j]]
            #print 'np.unique(partOfColumns) ', np.unique(partOfColumns)
            #print 'listOfROIsInMask[j] ', listOfROIsInMask[j]
            listOfROIsInMask[j] = np.unique(partOfColumns)
            print 'for mask ROI ', roisMask[j], ' the column IDs are ', len(listOfROIsInMask[j]), ' a list: ', str(listOfROIsInMask[j])
            
        # check which columns are in both ROIs
        # TODO!! decide what to do with these ROIs?? NOW - ignore them - 
        for i in range(len(roiIds)):
            # TODO! change it! need to ignore a column if it is present in any 2 mask ROIs, not necessarily 1st and 2nd
            if roiIds[i] != 0:                  #and roiIds[i] in listOfROIsInMask[0] and roiIds[i] in listOfROIsInMask[1]:
                includedIn = [int(roiIds[i] in listOfROIsInMask[w]) for w in range(len(listOfROIsInMask) - 1)]
                #print 'current ROI ', roiIds[i], ' is included in the following MASK rois ', includedIn
                
                if np.sum(includedIn) > 1:
                    print 'this ROI ', roiIds[i], ' is included in ', includedIn                
                    #print 'Column with ID ', roiIds[i], ' is in both regions'
                    # check what percentage of this column is in each mask ROI
                    numbers = []
                    for k in roisMask:                
                        w = arrColumns[(arrColumns == roiIds[i]) & (mask == k)] # get voxels in the respective mask ROI and the respective cortical column
                        numbers.append(len(w))                
                    print 'Column with ID ', roiIds[i], ' is in both regions. % in mask ROIs ', ' is ', numbers
                    
                    # analyze these numbers! TODO: Denis : leave it like this!! eliminate columns that are in both ROIs!!! as Voronoi can also contain errors!
                    # if one of the rois contains 10 times more voxels than all the rest ROIs, than we can accept this column and say that it belongs to this roi
                    # e.g. if column has 1200 voxels in ROI 11 and 13 voxels in ROI 21  - should it be dismissed??
                    
                    maxN = np.max(numbers)
                    ratio = maxN / (np.sum(numbers) - maxN)
                    if ratio >= 10:
                        # accept this column to the ROI: where its most voxels are
                        maxId = numbers.index(maxN) # maxId gives the ROI in the mask to which this column should be attributed
                        print 'ratio = ', ratio, '  maxId = ', maxId, ' mask ROI ', roisMask[maxId]
                        
                        # do not put onto ignore list! delete this column id from the lists of other mask ROIs
                        for k in range(len(roisMask)):
                            if roisMask[k] != roisMask[maxId]:
                                # for the case if there are 3 or more regions: need to check whether to delete from the third region
                                if roiIds[i] in listOfROIsInMask[k]:
                                    print 'delete the column number ', roiIds[i], ' from mask ROI ', roisMask[k]
                                    #print 'exclude from ', listOfROIsInMask[k]
                                    listOfROIsInMask[k] = np.delete(listOfROIsInMask[k], list(listOfROIsInMask[k]).index(roiIds[i]))   
                                    #print 'new list is ', listOfROIsInMask[k]
                    else:
                        # put it into ignore list
                        listOfROIsInMask[len(roisMask)].extend([roiIds[i]])               # a list of ROIs to ignore
                        print 'ignore the column number ', roiIds[i]     
                        # need to delete these columns from all other lists????
                #else:
                    #print 'Column with ID ', roiIds[i], ' is in one region - maskROI ', includedIn
        print 'final ignore list = ', listOfROIsInMask[len(roisMask)]                
                
            
        
    #else :   TODO!!!!!!!!!!!!!!!!!!!!!!!!!!
        #coords = arrCoord[arrCoord != 0]
        #values = arrValue[arrCoord != 0]
        
    # 2 arrays of coordinates and values. 
    res = []
    res.append(coords)
    res.append(values)
    res.append(roiIds)  # ids of the 'cortical' columns
    res.append(listOfSeparateCoords)
    res.append(listOfSeparateValues)
    res.append(listOfROIsInMask)
    res.append(listOfSizes)
    # add the info about div gradn averages
    res.append(listOfAvgDivGradn)
    res.append(volColoured)
    res.append(listOfAvgDivGradnCateg)
    res.append(listOfSeparateAVGCoords)
    res.append(strMeansInfo)
    return(res)       
    
    
    
# this function calculates 7-number statistics, prints it out and writes into a given file
def computeRanges(arr, fileF, keyWord):
    meanV = np.mean(arr)
    sdV = np.std(arr)
    lowerOutlierBorder = np.round(meanV - 3 * sdV)     
    upperOutlierBorder = np.round(meanV + 3 * sdV)     
    # lower and upper outlier borders - correspond to which percentiles?
    sizeN = len(arr)
    lowerPercentage = len(np.where(arr < lowerOutlierBorder)[0]) * 100.0 / sizeN
    upperPercentage = len(np.where(arr < upperOutlierBorder)[0]) * 100.0 / sizeN
    print 'for *%s* mean = ' %(keyWord), meanV, ' sd= ', sdV, ' mean +- 3sd = (', lowerOutlierBorder, ' , ', upperOutlierBorder, '), percentiles : (' , lowerPercentage, upperPercentage, ')'
    fT2intensInfo.write('for *' + keyWord + '* mean= ' + str(np.round(meanV)) + ' sd= ' + str(np.round(sdV)) + ' mean +- 3sd = (' + str(lowerOutlierBorder) + ' , '+ str(upperOutlierBorder) + '), percentiles : (' + str(lowerPercentage) + ' , ' + str(upperPercentage) + ')\n')
    fileF.write('------------------------------------------------------------------- \n')
    return(0)


# this function will verify whether the Voronoi segmentation of the ROIs is 
# outside of the T2 acquired images# further, it will notify if voronoi touches the lower or the upper 4 layers of T2
def testVoronoiOutsideT2(voronoiVol, T2Vol, numOfLayerToAvoid, fileStat = None):
    print '----------------------- start testVoronoiOutsideT2 --------------------------------------'
    arrVoronoi = np.asarray(voronoiVol)
    arrT2 = np.asarray(T2Vol)
    #print 'arrVoronoi ', np.shape(arrVoronoi),', arrT2 ', np.shape(arrT2)   
   
    # find z layers of T2 that are NOT completely zeros
    xCoordsAllT2, yCoordsAllT2, zCoordsAllT2, tCoordsAllT2 = np.where(arrT2 > 0)
    lowerNonZeroLayerT2 = np.max(zCoordsAllT2)       
    upperNonZeroLayerT2 = np.min(zCoordsAllT2)       
    print 'non zero z limits of T2 ', upperNonZeroLayerT2, lowerNonZeroLayerT2
          
    # compare it to the top and bottom z coordinates of the Voronoi
    xCoordsVor, yCoordsVor, zCoordsVor, tCoordsVor = np.where(arrVoronoi != 0)
    lowerNonZeroLayerVor = np.max(zCoordsVor)       
    upperNonZeroLayerVor = np.min(zCoordsVor) 
    print 'non zero z limits of Voronoi ', upperNonZeroLayerVor, lowerNonZeroLayerVor
    
    # write out info
    if fileStat is not None:
        fileStat.write('non zero z limits of T2 ' + str(upperNonZeroLayerT2) + ' , ' + str(lowerNonZeroLayerT2) + '\n')
        fileStat.write('non zero z limits of Voronoi ' + str(upperNonZeroLayerVor) + ' , ' + str(lowerNonZeroLayerVor) + '\n')
 
    if (upperNonZeroLayerVor < upperNonZeroLayerT2):
        overlapTop = upperNonZeroLayerVor - upperNonZeroLayerT2
        print 'Problem! Voronoi in zero T2 layers from top! Overlap of %s layers!' %(overlapTop)
        if fileStat is not None:
            fileStat.write('Problem! Voronoi in zero T2 layers from top! Overlap of %s layers!\n' %(str(overlapTop)))            
    # even if voronoi is NOT in the zero layer of T2 ,it might be in dark layers. Consider 4 top and bottom layers - too dark!
    if (upperNonZeroLayerVor < (upperNonZeroLayerT2 + numOfLayerToAvoid)):
        overlapTopDark = upperNonZeroLayerVor - (upperNonZeroLayerT2 + numOfLayerToAvoid)
        print 'Problem! Voronoi in %s dark T2 layers from top! Overlap of %s layers!' %(numOfLayerToAvoid, overlapTopDark)
        if fileStat is not None:
            fileStat.write('Problem! Voronoi in %s dark T2 layers from top! Overlap of %s layers!\n' %(str(numOfLayerToAvoid), str(overlapTopDark)))         
    if (lowerNonZeroLayerVor > lowerNonZeroLayerT2):
        overlapBottom = lowerNonZeroLayerVor - lowerNonZeroLayerT2
        print 'Problem! Voronoi in zero T2 layers from bottom! Overlap of %s layers!' %(overlapBottom)    
        if fileStat is not None:
            fileStat.write('Problem! Voronoi in zero T2 layers from bottom! Overlap of %s layers!\n' %(str(overlapBottom)))         
    # even if voronoi is NOT in the zero layer of T2 ,it might be in dark layers. Consider 4 top and bottom layers - too dark!
    if (lowerNonZeroLayerVor > (lowerNonZeroLayerT2 - numOfLayerToAvoid)):
        overlapBottomDark = lowerNonZeroLayerVor - (lowerNonZeroLayerT2 - numOfLayerToAvoid)
        print 'Problem! Voronoi in %s dark T2 layers from bottom! Overlap of %s layers!' %(numOfLayerToAvoid, overlapBottomDark)
        if fileStat is not None:
            fileStat.write('Problem! Voronoi in %s dark T2 layers from bottom! Overlap of %s layers!\n' %(str(numOfLayerToAvoid), str(overlapBottomDark))) 

    #print 'number of voxels in mask that are nonzero ', len(np.where(arrVoronoi != 0)[0])
    minValArrT2 = np.min(arrT2)
    maxValArrT2 = np.max(arrT2)
    #print minValArrT2, maxValArrT2
    
    arrT2[arrVoronoi == 0] = (maxValArrT2 + 10)
    arrT2[arrT2 > 0] = (maxValArrT2 + 10)
    # how many zero voxels are left
    xCoords, yCoords, zCoords, tCoords = np.where(arrT2 == 0)
    print 'number of T2 zero in nonZero Voronoi mask ', len(xCoords)
 
    #print zip(xCoords, yCoords, zCoords)
    #print 'sorted x ', np.sort(xCoords)
    #print 'sorted y ', np.sort(yCoords)
    #print 'sorted z ', np.sort(zCoords)
    #plt.hist(zCoords, bins = 100)
    unique, counts = np.unique(zCoords, return_counts = True)
    print 'z coordinates of zeros and their counts ', zip(unique, counts)
    if fileStat is not None:
        fileStat.write('number of T2 zero in nonZero Voronoi mask %s \n' %(len(xCoords))) 
        fileStat.write('z coordinates of zeros and their counts %s \n' %(zip(unique, counts)))
    
    #plt.show()

    #sys.exit(0)
    
    return(0)
    
    
    
    
    
if __name__ == '__main__':
    
    realPatientID = None
    directory = None
    realSide = 'L'
    columnDiameter = None
    workOnLaptop = False
    reExtractProfilesCutBorders = False
    #pathToNobiasT2 = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S4/'
    #pathToNobiasT2_new = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S16/'
    pathToNobiasT2_newCroppedDB = '/neurospin/lnao/dysbrain/brainvisa_db_T2_cropped/dysbrain/'   
    pathToNew_T1inT2DB = '/neurospin/lnao/dysbrain/brainvisa_db_T1_in_T2_space/dysbrain/'
    minColumnSize = 10 # depends on the diameter
    heights = []
    minColumnSizes = []
    ignoreCommunColums = True
    #divGradnThresholds = [-0.1, 0.1]
    #divGradnThresholds = [-0.2, 0.2]
    #divGradnThresholds = [-0.2, 0.3]
    #divGradnThresholds = [-0.3, 0.3]
    divGradnThresholds = [-0.2, 0.35]
    #divGradnThresholds = [-0.25, 0.35]
    #divGradnThresholds = [-0.5, 0.5] #just for test
    heatCaluclation = None      # version of the heat volume (and so the corresponding columns) to use


    parser = OptionParser('Extract profiles from T2 nobias data using cortex-density-coordinates in ROIs')    
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    parser.add_option('-s', dest='realSide', help='Hemisphere to be processed: L or R. L is default')   
    parser.add_option('-d', dest='directory', help='directory')
    parser.add_option('-c', dest='columnDiameter', help='columnDiameter to work with')
    parser.add_option('-g', dest='ignoreCommunColums', action = 'store_false', help='Select if want to INCLUDE into calculation cortical colums found in several ROIs (Excluding them is default')
    parser.add_option('-l', dest='workOnLaptop', action = 'store_true', help='Select if working on laptop (neurospin DB location is different. False is default') 
    parser.add_option('-j', dest='heatCaluclation', help='Version of the heat calculation program: old or new') 

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

    if options.heatCaluclation is None:
        print >> sys.stderr, 'New: exit. No version for the heat calculation was given'
        sys.exit(1)
    else:
        heatCaluclation = options.heatCaluclation

    if options.workOnLaptop is not None:
	    workOnLaptop = options.workOnLaptop      
	    # if true, then processes are run on the laptop. Change locations of neurospin DBs
	    #pathToNobiasT2 = pathToNobiasT2.replace('/neurospin/lnao/', '/nfs/neurospin/lnao/')
	    #pathToNobiasT2_new = pathToNobiasT2_new.replace('/neurospin/lnao/', '/nfs/neurospin/lnao/')  
	    
	    # saving the data locally on the laptop
	    #pathToNobiasT2 = pathToNobiasT2.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')
	    #pathToNobiasT2_new = pathToNobiasT2_new.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')  	
	    pathToNobiasT2_newCroppedDB = pathToNobiasT2_newCroppedDB.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')
            pathToNew_T1inT2DB = pathToNew_T1inT2DB.replace('/neurospin/lnao/dysbrain/', '/volatile/od243208/neurospin/')
            
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
    
    pathToCoord = pathToColumnResults + 'isovolume/'
    volsCoord = glob.glob(pathToCoord + 'pial-volume-fraction_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide))
    #volValue = aims.read(pathToNobiasT2 + '%s_NewNobiasT2_cropped.nii.gz' %(realPatientID))
    #volValue2 = aims.read(pathToNobiasT2_new + '%s_NewT2_cropped.nii.gz' %(realPatientID))
    volValue2 = aims.read(pathToNobiasT2_newCroppedDB + '%s/t1mri/t2_resamp/%s.nii.gz' %(realPatientID, realPatientID))
    volsMask = glob.glob(pathToColumnResults + 'voronoiCorr_%s_%s_cut_noSulci.nii.gz' %(realPatientID, realSide))
    print 'volsCoord = ' , pathToCoord + 'pial-volume-fraction_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)
    print 'volsMask = ', pathToColumnResults + 'voronoiCorr_%s_%s_cut_noSulci.nii.gz' %(realPatientID, realSide)
    print 'volValue2 = ', pathToNobiasT2_newCroppedDB + '%s/t1mri/t2_resamp/%s.nii.gz' %(realPatientID, realPatientID)     
    volClassifNoBorders = pathToColumnResults + 'GWsegm_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)
    fT2intensInfo = open(pathToColumnResults + '%s_%s_statFileProfiles.txt' %(realPatientID, realSide), "w")

    # test if all data is available
    if len(volsCoord) != 1 or len(volsMask) != 1:
        # abort the calculation, as too many or not a single texture file was found
        print 'abort the calculation, as too many or not a single volsCoord and volsMask file was found'
        fT2intensInfo.write('abort the calculation, as ' + str(len(volsCoord)) + ' volsCoord and ' + str(len(volsMask)) + ' volsMask files were found' + '\n')
        fT2intensInfo.close()
        sys.exit(0)    
        
    volCoord = aims.read(volsCoord[0])  
    volMask = aims.read(volsMask[0])
    
    
    
    ################ TODO ! ########  repeat this analysis for the new subjects!!!!!!!!!!!!!!!! #######################################################
    ## test how many zeros in the data 
    #testVoronoiOutsideT2(volMask, volValue2, 4, fT2intensInfo)    
    #sys.exit(0)
    ###################################################################################################################################################    
    
    # study T2 intensities for data cleaning  
    maskGW = np.asarray(aims.read(volClassifNoBorders))
    print 'maskGW = ' , volClassifNoBorders
    volFullGW = pathToNew_T1inT2DB + '%s/t1mri/t1m_resamp/default_analysis/segmentation/%sgrey_white_%s.nii.gz' %(realPatientID, realSide, realPatientID)    
    print 'volFullGW = ', volFullGW
    maskGWFull = np.asarray(aims.read(volFullGW))# also compare to all data available in T2, not only ROIs
    arrValue = np.asarray(volValue2)
    arrValue1WM = arrValue[maskGW == 200] # get values in the WM
    arrValue1Cortex = arrValue[maskGW == 100] # get values in the cortex
    arrValue1CSF = arrValue[maskGW == 0] # get values in the CSF
    #fig = plt.figure(figsize=(3 * 7, 6)) #, dpi=80, facecolor='w', edgecolor='k')
    fig = plt.figure(figsize=(3 * 7, 2* 6)) #, dpi=80, facecolor='w', edgecolor='k')
    # and now plot them
    axWM = fig.add_subplot(2,3,1)
    axWM.hist(arrValue1WM, bins = 50)    
    axWM.set_title('Histogram of T2 intensities in WM')   # subplot 211 title
    axCortex = fig.add_subplot(2,3,2)
    axCortex.hist(arrValue1Cortex, bins = 50)
    axCortex.set_title('Histogram of T2 intensities in cortex')   # subplot 211 title
    axCSF = fig.add_subplot(2,3,3)
    axCSF.hist(arrValue1CSF, bins = 50)
    axCSF.set_title('Histogram of T2 intensities in CSF')   # subplot 211 title
    
    # just for test: the same but EXCLUDING zeros
    arrValue1CortexNoZeros = arrValue1Cortex[arrValue1Cortex > 0]
    arrValue1WMNoZeros = arrValue1WM[arrValue1WM > 0]
    arrValue1CSFNoZeros = arrValue1CSF[arrValue1CSF > 0]
    axWMNoZeros = fig.add_subplot(2,3,4)
    axWMNoZeros.hist(arrValue1WMNoZeros, bins = 50)    
    axWMNoZeros.set_title('Histogram of T2 intensities in WM NoZeros')   # subplot 211 title
    axCortexNoZeros = fig.add_subplot(2,3,5)
    axCortexNoZeros.hist(arrValue1CortexNoZeros, bins = 50)
    axCortexNoZeros.set_title('Histogram of T2 intensities in cortex NoZeros')   # subplot 211 title
    axCSFNoZeros = fig.add_subplot(2,3,6)
    axCSFNoZeros.hist(arrValue1CSFNoZeros, bins = 50)
    axCSFNoZeros.set_title('Histogram of T2 intensities in CSF NoZeros')   # subplot 211 title    
    #plt.savefig(pathToColumnResults + '%s_%s_histOfT2intensitiesInROIs.png' %(realPatientID, realSide))    
    plt.savefig(pathToColumnResults + '%s_%s_histOfT2intensInROIs.png' %(realPatientID, realSide))    
    plt.clf()
    plt.close()
 
   
    arrValue1WMFull = arrValue[maskGWFull == 200] # get values in the WM
    arrValue1CortexFull = arrValue[maskGWFull == 100] # get values in the cortex
    arrValue1CSFFull = arrValue[maskGWFull == 0] # get values in the CSF
    fig = plt.figure(figsize=(3 * 7, 6)) #, dpi=80, facecolor='w', edgecolor='k')
    # and now plot them
    axWM = fig.add_subplot(1,3,1)
    axWM.hist(arrValue1WMFull, bins = 50)    
    axWM.set_title('Histogram of T2 intensities in WM')   # subplot 211 title
    axCortex = fig.add_subplot(1,3,2)
    axCortex.hist(arrValue1CortexFull, bins = 50)
    axCortex.set_title('Histogram of T2 intensities in cortex')   # subplot 211 title
    axCSF = fig.add_subplot(1,3,3)
    axCSF.hist(arrValue1CSFFull, bins = 50)
    axCSF.set_title('Histogram of T2 intensities in CSF')   # subplot 211 title
    plt.savefig(pathToColumnResults + '%s_%s_histOfT2intensitiesAllData.png' %(realPatientID, realSide))    
    plt.clf()
    plt.close()
    
    
    # print out the percentiles
    percentiles = [1, 5, 10, 90, 95, 99]
    # means and sd
    #meanV = np.mean(arrValue1Cortex)
    #sdV = np.std(arrValue1Cortex)
    #lowerOutlierBorder = np.round(meanV - 3 * sdV)      # cortex
    #upperOutlierBorder = np.round(meanV + 3 * sdV)
    # lower and upper outlier borders - correspond to which percentiles?
    sizeNCortex = len(arrValue1Cortex)
    sizeNCortexNoZeros = len(arrValue1CortexNoZeros)   
    #lowerPercentage = len(np.where(arrValue1Cortex < lowerOutlierBorder)[0]) * 100.0 / sizeN
    #upperPercentage = len(np.where(arrValue1Cortex < upperOutlierBorder)[0]) * 100.0 / sizeN 
    #print 'mean arrValue1Cortex = ', meanV, ' sd= ', sdV, ' mean +- 3sd = (', lowerOutlierBorder, ' , ', upperOutlierBorder, '), percentiles : (' , lowerPercentage, upperPercentage, ')'
    #fT2intensInfo.write('mean arrValue1Cortex = ' + str(np.round(meanV)) + ' sd= ' + str(np.round(sdV)) + ' mean +- 3sd = (' + str(lowerOutlierBorder) + ' , '+ str(upperOutlierBorder) + '), percentiles : (' + str(lowerPercentage) + ' , ' + str(upperPercentage) + ')\n')
    
    
    # TODO: looks like the problem is here:
    print 'problem?  ', np.mean(arrValue1Cortex)
    
    
    
    
    computeRanges(arrValue1Cortex, fT2intensInfo, 'cortexROIs')
    

    ######################## simulate the same EXCLUDING zeros (supposably segmentation errors)
    numZeros = len(np.where(arrValue1Cortex == 0)[0])
    print 'number of Voxels in Cortex ROIs with 0 intensity = ', numZeros
    fT2intensInfo.write('number of Voxels in Cortex ROIs with 0 intensity = ' + str(numZeros) + '\n')
    #meanVNoZeros = np.mean(arrValue1CortexNoZeros)
    #sdVNoZeros = np.std(arrValue1CortexNoZeros)
    #lowerOutlierBorderNoZeros = np.round(meanVNoZeros - 3 * sdVNoZeros)      # cortex, excluding zeros
    #upperOutlierBorderNoZeros = np.round(meanVNoZeros + 3 * sdVNoZeros)
    ## lower and upper outlier borders - correspond to which percentiles?
    #sizeNNoZeros = len(arrValue1CortexNoZeros)
    #lowerPercentageNoZeros = len(np.where(arrValue1CortexNoZeros < lowerOutlierBorderNoZeros)[0]) * 100.0 / sizeNNoZeros
    #upperPercentageNoZeros = len(np.where(arrValue1CortexNoZeros < upperOutlierBorderNoZeros)[0]) * 100.0 / sizeNNoZeros 
    #print 'mean arrValue1CortexNoZeros = ', meanVNoZeros, ' sdNoZeros= ', sdVNoZeros, ' meanNoZeros +- 3sdNoZeros = (', lowerOutlierBorderNoZeros, ' , ', upperOutlierBorderNoZeros, '), percentiles : (' , lowerPercentageNoZeros, upperPercentageNoZeros, ')'
    #fT2intensInfo.write('############ data with excluded zero values ############### \n')
    #fT2intensInfo.write('mean arrValue1CortexNoZeros = ' + str(np.round(meanVNoZeros)) + ' sd= ' + str(np.round(sdVNoZeros)) + ' mean +- 3sd = (' + str(lowerOutlierBorderNoZeros) + ' , '+ str(upperOutlierBorderNoZeros) + '), percentiles : (' + str(lowerPercentageNoZeros) + ' , ' + str(upperPercentageNoZeros) + ')\n')
    computeRanges(arrValue1CortexNoZeros, fT2intensInfo, 'cortexROIs_NoZeros')
    
    
    
    
    # WM    
    computeRanges(arrValue1WM, fT2intensInfo, 'WMROIs')
    numZerosWM = len(np.where(arrValue1WM == 0)[0])
    print 'number of Voxels in WM ROIs with 0 intensity = ', numZerosWM
    fT2intensInfo.write('number of Voxels in WM ROIs with 0 intensity = ' + str(numZerosWM) + '\n')
    #meanVNoZerosWM = np.mean(arrValue1WMNoZeros)
    #sdVNoZerosWM = np.std(arrValue1WMNoZeros)
    #lowerOutlierBorderNoZerosWM = np.round(meanVNoZerosWM - 3 * sdVNoZerosWM)      # cortex, excluding zeros
    #upperOutlierBorderNoZerosWM = np.round(meanVNoZerosWM + 3 * sdVNoZerosWM)
    ## lower and upper outlier borders - correspond to which percentiles?
    #sizeNNoZerosWM = len(arrValue1WMNoZeros)
    #lowerPercentageNoZerosWM = len(np.where(arrValue1WMNoZeros < lowerOutlierBorderNoZerosWM)[0]) * 100.0 / sizeNNoZerosWM
    #upperPercentageNoZerosWM = len(np.where(arrValue1WMNoZeros < upperOutlierBorderNoZerosWM)[0]) * 100.0 / sizeNNoZerosWM
    #print 'mean arrValue1WMNoZeros = ', meanVNoZerosWM, ' sdNoZerosWM= ', sdVNoZerosWM, ' meanNoZerosWM +- 3sdNoZerosWM = (', lowerOutlierBorderNoZerosWM, ' , ', upperOutlierBorderNoZerosWM, '), percentiles : (' , lowerPercentageNoZerosWM, upperPercentageNoZerosWM, ')'
    #fT2intensInfo.write('############ WM ############### \n')
    #fT2intensInfo.write('mean arrValue1WMNoZeros = ' + str(np.round(meanVNoZerosWM)) + ' sdWM= ' + str(np.round(sdVNoZerosWM)) + ' meanWM +- 3sdWM = (' + str(lowerOutlierBorderNoZerosWM) + ' , '+ str(upperOutlierBorderNoZerosWM) + '), percentiles : (' + str(lowerPercentageNoZerosWM) + ' , ' + str(upperPercentageNoZerosWM) + ')\n')
    computeRanges(arrValue1WMNoZeros, fT2intensInfo, 'WMROIs_NoZeros')    


    for p in percentiles:
        print 'percentile ', p, 'arrValue1Cortex = ', np.percentile(arrValue1Cortex, p), 'arrValue1CortexFull = ', np.percentile(arrValue1CortexFull, p), 'arrValue1CortexNoZeros = ', np.percentile(arrValue1CortexNoZeros, p)
        fT2intensInfo.write('percentile ' + str(p) + ' arrValue1Cortex = ' + str(np.percentile(arrValue1Cortex, p)) + ' arrValue1CortexFull = ' + str(np.percentile(arrValue1CortexFull, p)) + ' arrValue1CortexNoZeros = ' + str(np.percentile(arrValue1CortexNoZeros, p)) + '\n')
        
    fT2intensInfo.write('\n')
    #number ( % ) of voxels smaller then for CORTEX!
    smallerThen = [50, 40, 30, 20, 10, 5]
    for p in smallerThen:
        num = len(np.where(arrValue1Cortex < p)[0])
        numNoZeros = len(np.where(arrValue1CortexNoZeros < p)[0])
        print 'smallerThen ', p, ' arrValue1Cortex = ', num, ' from ', sizeNCortex, ' % = ', np.round(num * 100.0 / sizeNCortex),  ' arrValue1CortexNoZeros = ', numNoZeros, ' from ', sizeNCortexNoZeros, ' % = ', np.round(numNoZeros * 100.0 / sizeNCortexNoZeros) 
        fT2intensInfo.write('smallerThen ' + str(p) + ' arrValue1Cortex = ' + str(num) + ' from ' + str(sizeNCortex) + ' % = ' + str(np.round(num * 100.0 / sizeNCortex)) + ' arrValue1CortexNoZeros = ' + str(numNoZeros) + ' from ' + str(sizeNCortexNoZeros) + ' % = ' + str(np.round(numNoZeros * 100.0 / sizeNCortexNoZeros)) + '\n')    
    #sys.exit(0)

    
    # if no columns diameter was given, then extract profiles in the mask
    if columnDiameter is None:
        print '******************* no columns diameter was given, then extract profiles in the mask ***************************'    
        #result = extractProfiles(volCoord, volValue, volMask)
        # repeat for the NEW nobias images!
        result2 = extractProfiles(volCoord, volValue2, volMask)
    else:   # we work with columns! then calculate the required size. 
        columnDiameter = int(options.columnDiameter)
        print '******************* columns diameter was given = ', str(columnDiameter) , ' extract profiles in the columns ***************************'
        # human cortical thickness 2 - 4mm
        vs = volValue2.header()['voxel_size'][0:3]
        voxelSizeMax = np.max(vs)
        voxelSizeMin = np.min(vs)
        voxelSizeAvg = np.average(vs)
        heightMax = int(np.round(4.0 / voxelSizeMin))
        heightMin = int(np.round(2.0 / voxelSizeMax))
        print 'heights of the cortical columns found heightMin = ', heightMin, ' heightMax = ', heightMax
        #heights = range(heightMin, heightMax + 1)
        # just for test: to see the influence of really large columns
        heights = range(heightMin - 1, heightMax)
        print ' ##################################################### heights = ', heights
       # using this info: calculate minimal cortical column sizes
        
        # average size.     
        minColumnSize = int(np.round(columnDiameter * columnDiameter / voxelSizeAvg / voxelSizeAvg / 4 * np.pi * np.average(heights))) 
        print '##################################################### calculated the avgColumnSize = ', minColumnSize

        # a list of sizes        
        minColumnSizes = [int(k * columnDiameter * columnDiameter / voxelSizeAvg / voxelSizeAvg / 4 * np.pi) for k in heights]
        # modified! Consider narrow columns, or conical ones! Divide by a factor!     
        minColumnSizes = [int(np.min(minColumnSizes) / 4.0), int(np.min(minColumnSizes) / 2.0), int(np.min(minColumnSizes) / 1.5)] + minColumnSizes
        print 'new list of column sizes Consider narrow columns, or conical one ', minColumnSizes
                
        # start the processing for columns
        print '******************* columns diameter was given, then extract profiles in the columns ***************************'    
        # read in the columns file
        volName = ''
        if columnDiameter == 1:
            volName = pathToColumnResults + 'column_regions/merged_randomized_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)
        else:
            volName = pathToColumnResults + 'column_regions/merged_randomized_%s_%s_cut_noSulci_extended_diam%s.nii.gz' %(realPatientID, realSide, str(columnDiameter))
            
        volsColumns = glob.glob(volName)

        if len(volsColumns) != 1:# abort the calculation, as too many or not a single columns file was found
            print 'abort the calculation, as too many or not a single volsColumns file was found'
            fT2intensInfo.write('abort the calculation, as ' + str(len(volsColumns)) + ' volsColumns files were found' + '\n')
            fT2intensInfo.close()
            sys.exit(0)
            
        volColumns = aims.read(volsColumns[0])  
        print 'volColumns = ', volsColumns[0]        
        # read in the div_gradn file (for measuring the flatness of the resoective column region)
        volDivGradn = aims.read(pathToColumnResults + 'heat/heat_div_gradn_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide))          
        print 'volDivGradn = ', pathToColumnResults + 'heat/heat_div_gradn_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)
        #result = extractProfilesInColumns(volCoord, volValue, volColumns, volMask)
        # repeat for the NEW nobias images!
        
        
        # Do NOT need to re-extract the profiles in the merged-randomized volume, which was cut!!! (so that regions 50 and 150, corresponding to borders, were eliminated, but write this volume out)
        print 'cut borders from the volume ', volsColumns[0]
        print 'using the original file ', volClassifNoBorders
        # cut the borders
        subprocess.check_call(['AimsMerge', '-m', 'oo', '-l', '0', '-v', '0', '-i', volsColumns[0], '-M', volClassifNoBorders, '-o', pathToColumnResults + 'column_regions/traverses_without_CSF_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)])
        
        subprocess.check_call(['AimsMerge', '-m', 'oo', '-l', '200', '-v', '0', '-i', pathToColumnResults + 'column_regions/traverses_without_CSF_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide), '-M', volClassifNoBorders, '-o', pathToColumnResults + 'column_regions/traverses_cortex_only_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)])
              
        # TODO : delete it!      
        #sys.exit()
        
        result2 = extractProfilesInColumns(volCoord, volValue2, volColumns, volDivGradn, divGradnThresholds, volMask)   
        
        
        
        
        

        
            
    # work now only with the new nobias T2!!!!! commented the work with the old corrected nobias T2  
    #coordinates = result[0]
    #intensities = result[1]

    ## plot the data
    #plt.plot(coordinates, intensities, '.', c = 'b')
    #plt.title('Profile in ROI')   # subplot 211 title
    #plt.xlabel('Cortical depth')
    #plt.ylabel('T2-nobias intensity')
    #plt.savefig(pathToColumnResults + '%s_%s_It20_nobiasT2vsCorticalDepthROI.png' %(realPatientID, realSide))

    coordinates2 = result2[0]
    intensities2 = result2[1]

    ## plot the data
    #plt.plot(coordinates2, intensities2, '.', c = 'r')
    #plt.title('Profile in ROI')   # subplot 211 title
    #plt.xlabel('Cortical depth')
    #plt.ylabel('T2-nobias intensity')
    #plt.savefig(pathToColumnResults + '%s_%s_It20_2nobiasT2vsCorticalDepthROI.png' %(realPatientID, realSide))    
    #plt.clf()
    #plt.close()
    
    plt.plot(coordinates2, intensities2, '.', c = 'r')
    plt.title('Profile in all ROIs')   # subplot 211 title
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')
    plt.savefig(pathToColumnResults + '%s_%s_newNobiasT2ROIs.png' %(realPatientID, realSide))    
    plt.clf()
    plt.close()
  
    # save the data for further processing. 
    #data1 = open(pathToColumnResults + '%s_%s_profiles.txt' %(realPatientID, realSide), "w")
    data2 = open(pathToColumnResults + '%s_%s_profiles2.txt' %(realPatientID, realSide), "w")
    headerLine = '\t' + 'DepthCoord' + '\t' + 'Value'
    #data1.write(headerLine + '\n')
    data2.write(headerLine + '\n')
    for i in range(len(coordinates2)):
        #data1.write(str(i) + '\t' + str(coordinates[i]) + '\t' + str(intensities[i]) + '\n')
        data2.write(str(i) + '\t' + str(coordinates2[i]) + '\t' + str(intensities2[i]) + '\n')
    #data1.close()
    data2.close()
    
    ## now plot and save the data for individual ROIs
    iDs = result2[2]
    listOfCoords = result2[3]
    listOfValues = result2[4]
    listOfAvgDivGrads = result2[7]
    volColouredByAvgDivGrads = result2[8]
    listOfAvgDivGradsCateg = result2[9]         # a list of 3 categories : 10, 29 and 30, given according to the 2 set thresholds
    listOfAvgVoxelCoords = result2[10]
    strMeansInfo = result2[11]  # add it to the file
    fT2intensInfo.write('\n' + strMeansInfo + '\n')
    fT2intensInfo.close()
    #sys.exit(0)
    addedColumnsDiamName = ''
    pathForFiles = pathToColumnResults
    print '############################################ columnDiameter = ', str(columnDiameter)
    if columnDiameter is not None:
        addedColumnsDiamName = '_diam%s' %(columnDiameter)
        # TODO: delete it after the tests!
        addedColumnsDiamName = addedColumnsDiamName + '_scaledNoOutlAddOutl'        

        # if the result was for the columns, create a separate folder for it
        pathForFiles = pathForFiles + 'diam%s/'%(columnDiameter)
        if not os.path.exists(pathForFiles):
            os.makedirs(pathForFiles)
        print '############################################ created directory = ', pathForFiles

        # read in the results for columns IDs in various regions
        iDsInMaskROIs = result2[5]
        
    # write out a file with the info: Column ID, size, maskROI (belongs to it or not), To be ignored
    maskArray = np.asarray(volMask)
    maskROIids = np.unique(maskArray[np.where(maskArray > 0)])            
    dataS = open(pathForFiles + '%s_%s_ColumnInfo%s.txt' %(realPatientID, realSide, addedColumnsDiamName), "w")
    headerInfo = '\tColumnID\tSize\tAvgX\tAvgY\tAvgZ\tAvgDivGradn'
    
    # complete the header with the mask ROI info. If columnDiameter was given
    if columnDiameter is not None:        
        for i in range(len(maskROIids)):
            headerInfo += '\tMaskROI_' + str(maskROIids[i]) + '\t'    
        headerInfo += 'ToIgnore\n'    
    else:
        headerInfo += '\n'    
        
    dataS.write(headerInfo)    
    sizeSorted = pathForFiles + 'sortedBySize/'
    divGradnSorted = pathForFiles + 'sortedByDivGradn/'
    divGradnClasses = pathToColumnResults + 'divGradnClasses/'
    corticalLayers = pathForFiles + 'corticalLayers/'

    if not os.path.exists(sizeSorted):
        os.makedirs(sizeSorted)
    if not os.path.exists(divGradnSorted):
        os.makedirs(divGradnSorted)
    if not os.path.exists(divGradnClasses):
        os.makedirs(divGradnClasses)
    if not os.path.exists(corticalLayers):
        os.makedirs(corticalLayers)
       
    #test whether the length of ids and avgdivgradns is the same
    if len(iDs) != len(listOfAvgDivGrads):
        print 'lengths are unequal!!! ', len(iDs), len(listOfAvgDivGrads), ' force exit'
        sys.exit(0)
    
    for i in range(len(iDs)):
        #print '------------------ i = ', i, ' work with id', iDs[i]
        currCoords = listOfCoords[i]
        currValues = listOfValues[i]
        currAvgVoxelCoords = listOfAvgVoxelCoords[i]
        
        # get the info on the gradient
        currGradn = listOfAvgDivGrads[i]
        strCurrGradn = "%.3f" % currGradn
        strCurrGradnStr = str(strCurrGradn).replace('.','')           
        currLine = '\t%s\t%s\t%s\t%s\t%s\t%s' %(str(iDs[i]), len(currCoords), str("%.3f" % currAvgVoxelCoords[0]), str("%.3f" % currAvgVoxelCoords[1]), str("%.3f" % currAvgVoxelCoords[2]), strCurrGradn)
        # add info on Mask ROIs. if diameter was given
        if columnDiameter is not None:
            for j in range(len(maskROIids)):
                #print 'j = ', j, 'iDsInMaskROIs[j] = ', iDsInMaskROIs[j]
                if iDs[i] in iDsInMaskROIs[j]:
                    currLine += '\t1'
                    #print iDs[i], ' is in the list'
                else:
                    currLine += '\t0'
                    #print iDs[i], ' is NOT in the list'                    
                    
            # check if this element is in a list to be ignored
            #print 'a list to be ignored: ', iDsInMaskROIs[j + 1]
            if iDs[i] in iDsInMaskROIs[j + 1]:
                currLine += '\t1'
                #print 'ignore ', i
            else:
                currLine += '\t0'  
       
        # do not plot if it is zero
        if len(currCoords) != 0:   
            data2i = open(pathForFiles + '%s_%s_profiles2%s_ROI_%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(iDs[i])), "w")
            data2i.write(headerLine + '\n')
            for j in range(len(currCoords)):
                data2i.write(str(j) + '\t' + str(currCoords[j]) + '\t' + str(currValues[j]) + '\n')    
            data2i.close()              
            
        plt.clf()
        plt.close()     
        currLine += '\n'
        dataS.write(currLine)
    dataS.close()    
        
    
    # plot the ROIs from the mask (columns of any size) on one plot
    roiNames = ''
    for i in range(len(iDs)):
        roiNames += '_' + str(iDs[i])
        currCoords = listOfCoords[i]
        currValues = listOfValues[i]
        plt.plot(currCoords, currValues, '.', label = 'ROI ' + str(iDs[i]))
        
    plt.title('Profile in mask ROIs')   
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')
        
    if columnDiameter is not None:
        roiNames = '_allColumns' 
    else:
        # plot a legend only for the case without columns. otherwise - too many regions
        plt.legend(loc='upper right', numpoints = 1)
        
    plt.savefig(pathForFiles + '%s_%s_nobiasT2%s_ROIs%s' %(realPatientID, realSide, addedColumnsDiamName, roiNames) + '.png')
    plt.clf()
    plt.close() 
     
    strT = (str(divGradnThresholds[0]) + '_' +  str(divGradnThresholds[1])).replace('.','')         # transform thresholds into the string            
    # get results for various sizes of columns
    if columnDiameter is not None:
        # write out the volume coloured according to the average divGradn Values        
        aims.write(volColouredByAvgDivGrads, divGradnClasses + '%s_%s_nobiasT2%s_divGradnColour%s.nii.gz' %(realPatientID, realSide, addedColumnsDiamName, strT))
        # plot the histogram of column sizes
        sizes = result2[6]
        plt.hist(sizes, bins = 50)
        plt.title('Histogram of columns sizes for diameter %s ' %(columnDiameter))   # subplot 211 title
        plt.xlabel('Cortical columns sizes')
        plt.ylabel('Percentage')
        plt.savefig(pathToColumnResults + '%s_%s_histOfColumnSizes%s.png' %(realPatientID, realSide, addedColumnsDiamName))    
        plt.clf()
        plt.close()
        
        # plot the "theoretical" histogram of column heights
        theorHeights = [ int (w * 4 * voxelSizeAvg * voxelSizeAvg / np.pi / columnDiameter / columnDiameter) for w in sizes]
        plt.hist(theorHeights, bins = 50)
        plt.title('Histogram of theoretical heights of columns for diameter %s ' %(columnDiameter))   # subplot 211 title
        plt.xlabel('Theoretical heights of columns')
        plt.ylabel('Percentage')
        plt.savefig(pathToColumnResults + '%s_%s_histOfTheorHeights%s.png' %(realPatientID, realSide, addedColumnsDiamName))    
        plt.clf()
        plt.close()
        
        listsOfCoordsForLargeColumns = [[] for x in range(len(minColumnSizes))]
        listsOfValuesForLargeColumns = [[] for x in range(len(minColumnSizes))] 
        listsOfIDsForLargeColumns = [[] for x in range(len(minColumnSizes))] 
        listsOfNumbersLargeColumns = [0] * len(minColumnSizes)
        for i in range(len(iDs)):
            currCoords = listOfCoords[i]
            currValues = listOfValues[i]
            #print 'iDs ', iDs[i], len(currCoords)
            for j in range(len(minColumnSizes)):
                #print 'work with minColumnSizes[j] ', minColumnSizes[j]
                if len(currCoords) > minColumnSizes[j]:
                    # add these data to the respective list element                    
                    #print 'this column size ', len(currCoords), ' is larger than ', minColumnSizes[j], ' old len= ', len(listsOfCoordsForLargeColumns[j])
                    listsOfCoordsForLargeColumns[j].extend(currCoords)
                    listsOfValuesForLargeColumns[j].extend(currValues)
                    listsOfIDsForLargeColumns[j].extend([iDs[i]])
                    listsOfNumbersLargeColumns[j] = listsOfNumbersLargeColumns[j] + 1
                    #print ' new len= ', len(listsOfCoordsForLargeColumns[j])

        # now we have info for each min column size. plot it
        for j in range(len(minColumnSizes)):
            print 'number of columns over ', minColumnSizes[j], ' is ', listsOfNumbersLargeColumns[j], ' from ', len(iDs), ' from ', len(np.unique(np.asarray(volColumns)))
            # if no such columns were found - do NOT plot!!!
            if len(listsOfCoordsForLargeColumns[j]) != 0:                
                plt.plot(listsOfCoordsForLargeColumns[j], listsOfValuesForLargeColumns[j], '.', c = 'r')
                plt.title('Profile in cortical columns larger %s' %(minColumnSizes[j]))   # subplot 211 title
                plt.xlabel('Cortical depth')
                plt.ylabel('T2-nobias intensity')
                plt.savefig(pathToColumnResults + '%s_%s_newNobiasT2%s_over%s.png' %(realPatientID, realSide, addedColumnsDiamName, minColumnSizes[j]))    
                plt.clf()
                plt.close()
            
            # plot the data from all cortical columns and from the large ones together
            plt.plot(coordinates2, intensities2, '.', c = 'b')
            plt.title('Profile in cortical columns larger %s vs. all points' %(minColumnSizes[j]))   # subplot 211 title
            plt.xlabel('Cortical depth')
            plt.ylabel('T2-nobias intensity')
            plt.plot(listsOfCoordsForLargeColumns[j], listsOfValuesForLargeColumns[j], '.', c = 'r')
            # TODO!!! this plot is not plotted!! ???
            plt.clf()
            plt.close()    
            
            # write out a file with IDs of columns that are larger than this threshold
            if len(listsOfIDsForLargeColumns[j]) != 0:
                dataID = open(pathForFiles + '%s_%s_IDs%s_over_%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(minColumnSizes[j])), "w")
                string = ''
                for l in listsOfIDsForLargeColumns[j]:
                    string += str(l) + '\t'             
                dataID.write(string)
                dataID.close()
            
        # now plot avg profiles for various mask regions
        listOfIDsToIgnore = iDsInMaskROIs[len(iDsInMaskROIs) - 1]
        listsOfCoordsForMask = [[] for x in range(len(iDsInMaskROIs) - 1)]
        listsOfValuesForMask = [[] for x in range(len(iDsInMaskROIs) - 1)] 
        
        # lists of ROIs of various sizes
        listsOfCoordsForMaskVariousThr = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)]
        listsOfValuesForMaskVariousThr = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)] 
        listsOfColumnIDsForMaskVariousThr = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)]
        
        # lists of ROIs of various sizes accepted depending on the AvgDivGradn (ADG) thresholds
        listsOfCoordsForMaskVariousThrADG = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)]
        listsOfValuesForMaskVariousThrADG = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)] 
        listsOfColumnIDsForMaskVariousThrADG = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)]        
        # lists of ROIs of various sizes rejected depending on the AvgDivGradn (ADG) thresholds (too small)
        listsOfCoordsForMaskVariousThrADGs = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)]
        listsOfValuesForMaskVariousThrADGs = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)] 
        listsOfColumnIDsForMaskVariousThrADGs = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)]
        # lists of ROIs of various sizes rejected depending on the AvgDivGradn (ADG) thresholds (too big)
        listsOfCoordsForMaskVariousThrADGb = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)]
        listsOfValuesForMaskVariousThrADGb = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)] 
        listsOfColumnIDsForMaskVariousThrADGb = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)]
        
        maskArray = np.asarray(volMask)
        maskROIids = np.unique(maskArray[np.where(maskArray > 0)])
        
        for j in range(len(maskROIids)):
            listOfIDs = iDsInMaskROIs[j]
            for k in listOfIDs:
                criterion = None
                # check if we need to ignore this data
                if ignoreCommunColums: 
                    #print 'ignoreCommunColums is True'
                    criterion = k not in listOfIDsToIgnore and k != 0
                else:
                    #print 'INCLUDE commun columns'
                    criterion = k != 0

                if criterion:     
                    # get the ID of this column in the original data
                    #print 'iDs =', iDs, ' k = ', k
                    indexOfThisID = np.where(iDs == k)[0][0]
                    #print 'indexOfThisID ', indexOfThisID
                    currCoords = listOfCoords[indexOfThisID]
                    currValues = listOfValues[indexOfThisID]                    
                    #get the avg div gradn for the current column
                    currAvgDivGradn = listOfAvgDivGradsCateg[indexOfThisID]
                    
                    for t in range(len(minColumnSizes)):
                        #print t, j
                        if len(currCoords) > minColumnSizes[t]:
                            listsOfCoordsForMaskVariousThr[j][t].extend(currCoords)
                            listsOfValuesForMaskVariousThr[j][t].extend(currValues)  
                            # create a list of length len(currCoords), full of the respective columnID (k)
                            repeatedID = [k] * len(currValues)                            
                            listsOfColumnIDsForMaskVariousThr[j][t].extend(repeatedID) 
                            
                            # check if the AvgDivGradn is between the thresholds, then add the coordinates/values to another list!!!
                            if currAvgDivGradn == 20:
                                listsOfCoordsForMaskVariousThrADG[j][t].extend(currCoords)
                                listsOfValuesForMaskVariousThrADG[j][t].extend(currValues)  
                                listsOfColumnIDsForMaskVariousThrADG[j][t].extend(repeatedID)  
                            elif currAvgDivGradn == 10:
                                listsOfCoordsForMaskVariousThrADGs[j][t].extend(currCoords)
                                listsOfValuesForMaskVariousThrADGs[j][t].extend(currValues)  
                                listsOfColumnIDsForMaskVariousThrADGs[j][t].extend(repeatedID)  
                            else:
                                listsOfCoordsForMaskVariousThrADGb[j][t].extend(currCoords)
                                listsOfValuesForMaskVariousThrADGb[j][t].extend(currValues)  
                                listsOfColumnIDsForMaskVariousThrADGb[j][t].extend(repeatedID)  
                            
                            
        # plot the data for mask rois of large columns (various size thresholds). For ROIs separately, and together
        colours = ['b', 'g', 'r', 'c', 'm', 'y', 'b']        
        for t in range(len(minColumnSizes)):            
            numOfRegions = len(maskROIids) + 1  # mask rois and one plot with all together   
            fig = plt.figure(figsize=(numOfRegions * 7, 6)) #, dpi=80, facecolor='w', edgecolor='k')
            # and now plot the same to the last plot: all mask ROIs together
            axAll = fig.add_subplot(1,numOfRegions,numOfRegions)
            #plt.plot()
            #plt.title('Profile in Mask ROIs of large cortical columns')   
            #plt.xlabel('Cortical depth')
            #plt.ylabel('T2-nobias intensity')
            roiNames = ''
            #print '-------------------------------------------------------- plot for minColumnSizes = ', minColumnSizes[t]
            numP = 0    # total number of points plotted

            for j in range(len(maskROIids)):
                roiNames += '_' + str(maskROIids[j])
                #plt.plot(listsOfCoordsForMaskVariousThr[j][t], listsOfValuesForMaskVariousThr[j][t], '.', label = 'ROI '+ str(maskROIids[j]))                
                ax1 = fig.add_subplot(1,numOfRegions,j + 1)
                ax1.plot(listsOfCoordsForMaskVariousThr[j][t], listsOfValuesForMaskVariousThr[j][t], '.', c = colours[j], label = 'ROI '+ str(maskROIids[j]))
                #plt.title('Profile in maskROI %s' %(maskROIids[j]))   # subplot 211 title
                ax1.set_title('Profile in maskROI %s' %(maskROIids[j]))
                ax1.set_xlabel('Cortical depth') 
                ax1.set_ylabel('T2-nobias intensity')                 
                ax1.legend(loc='upper right', numpoints = 1)               
                
                # and now plot the same to the last plot: all mask ROIs together
                axAll.plot(listsOfCoordsForMaskVariousThr[j][t], listsOfValuesForMaskVariousThr[j][t], '.', label = 'ROI '+ str(maskROIids[j]))
                #print '----------------------------------------------------- len(listsOfCoordsForMaskVariousThr[',j,'][',t,']) = ', len(listsOfCoordsForMaskVariousThr[j][t])
                numP += len(listsOfCoordsForMaskVariousThr[j][t])
                
                # write out a file with columnID, coord, value (corresponding to this mask ROI and size threshold)
                if len(listsOfCoordsForMaskVariousThr[j][t]) != len(listsOfColumnIDsForMaskVariousThr[j][t]):
                    print ' --------------------------- PROBLEM !'
                else:
                    # write out profiles if diameter is > 1
                    if columnDiameter > 1 and len(listsOfCoordsForMaskVariousThr[j][t]) > 0:                        
                        f = open(pathForFiles + '%s_%s_MaskROI%s_profiles%s_over_%s.txt' %(realPatientID, realSide, str(maskROIids[j]), addedColumnsDiamName, str(minColumnSizes[t])), "w")
                        #print 'write out a file %s_%s_MaskROI%s_profiles%s_over_%s.txt' %(realPatientID, realSide, str(maskROIids[j]), addedColumnsDiamName, str(minColumnSizes[t])), ' number of elements ', len(listsOfCoordsForMaskVariousThr[j][t])
                        s = 'ColumnID\tCoord\tValue\n'
                        f.write(s)
                        for z in range(len(listsOfCoordsForMaskVariousThr[j][t])):
                            s = str(listsOfColumnIDsForMaskVariousThr[j][t][z]) + '\t' + str(listsOfCoordsForMaskVariousThr[j][t][z]) + '\t' + str(listsOfValuesForMaskVariousThr[j][t][z]) + '\n'                   
                            f.write(s)        
                        f.close()                        
            
            axAll.set_title('Profile in all maskROIs')
            axAll.set_xlabel('Cortical depth') 
            axAll.set_ylabel('T2-nobias intensity') 
            axAll.legend(loc='upper right', numpoints = 1)            
                            
            # modify the name if cortical columns communfor ROIs were included/excluded
            nameInclExcl = ''
            if ignoreCommunColums:
                #print 'Commun columns are ignored'
                nameInclExcl = 'exclCommun'
            else:
                #print 'Commun columns are INCLUDED'
                nameInclExcl = 'inclCommun'
            # do not save the picture if it is empty
            if numP != 0:
                plt.savefig(pathToColumnResults + '%s_%s_nobiasT2_ROIs%s' %(realPatientID, realSide, roiNames) + '%s_over%s_%s.png' %(addedColumnsDiamName, minColumnSizes[t], nameInclExcl), bbox_inches='tight')            
            plt.clf()
            plt.close()            
        
        # plot the same but also: included and rejected regions, according to the div_gradn volume
        for t in range(len(minColumnSizes)):            
            numOfRegions = len(maskROIids) + 1  # mask rois and one plot with all together   
            fig = plt.figure(figsize=(numOfRegions * 7, 6 * 4)) #, dpi=80, facecolor='w', edgecolor='k')
            # and now plot the same to the last plot: all mask ROIs together
            axAll = fig.add_subplot(4,numOfRegions,numOfRegions)        
            axAllSelected = fig.add_subplot(4,numOfRegions,numOfRegions * 2)    # a plot for Roi 11 and 21, for the "selected" columns
            axAllRejectedS = fig.add_subplot(4,numOfRegions,numOfRegions * 3)    # a plot for Roi 11 and 21, for the "rejected" columns, too small AvgDivGradn
            axAllRejectedB = fig.add_subplot(4,numOfRegions,numOfRegions * 4)    # a plot for Roi 11 and 21, for the "rejected" columns, too high AvgDivGradn            
            
            #plt.plot()
            #plt.title('Profile in Mask ROIs of large cortical columns')   
            ##plt.xlabel('Cortical depth')
            #plt.ylabel('T2-nobias intensity')
            roiNames = ''
            print '-------------------------------------------------------- plot for minColumnSizes = ', minColumnSizes[t]
            numP = 0    # total number of points plotted

            for j in range(len(maskROIids)):
                roiNames += '_' + str(maskROIids[j])
                #plt.plot(listsOfCoordsForMaskVariousThr[j][t], listsOfValuesForMaskVariousThr[j][t], '.', label = 'ROI '+ str(maskROIids[j]))                
                ax1 = fig.add_subplot(4,numOfRegions,j + 1)
                ax1.plot(listsOfCoordsForMaskVariousThr[j][t], listsOfValuesForMaskVariousThr[j][t], '.', c = colours[j], label = 'ROI '+ str(maskROIids[j]))
                #plt.title('Profile in maskROI %s' %(maskROIids[j]))   # subplot 211 title
                ax1.set_title('Profile in maskROI %s' %(maskROIids[j]))
                ax1.set_xlabel('Cortical depth') 
                ax1.set_ylabel('T2-nobias intensity')                 
                ax1.legend(loc='upper right', numpoints = 1)               
                
                # and now plot the same to the last plot: all mask ROIs together
                axAll.plot(listsOfCoordsForMaskVariousThr[j][t], listsOfValuesForMaskVariousThr[j][t], '.', label = 'ROI '+ str(maskROIids[j]))
                print '----------------------------------------------------- len(listsOfCoordsForMaskVariousThr[',j,'][',t,']) = ', len(listsOfCoordsForMaskVariousThr[j][t])
                numP += len(listsOfCoordsForMaskVariousThr[j][t])
                
                # write out a file with columnID, coord, value (corresponding to this mask ROI and size threshold)
                if len(listsOfCoordsForMaskVariousThr[j][t]) != len(listsOfColumnIDsForMaskVariousThr[j][t]):
                    print ' --------------------------- PROBLEM !'
                else:
                    # write out profiles if diameter is > 1
                    if columnDiameter > 1 and len(listsOfCoordsForMaskVariousThr[j][t]) > 0:                        
                        f = open(pathForFiles + '%s_%s_MaskROI%s_profiles%s_over_%s.txt' %(realPatientID, realSide, str(maskROIids[j]), addedColumnsDiamName, str(minColumnSizes[t])), "w")
                        #print 'write out a file %s_%s_MaskROI%s_profiles%s_over_%s.txt' %(realPatientID, realSide, str(maskROIids[j]), addedColumnsDiamName, str(minColumnSizes[t])), ' number of elements ', len(listsOfCoordsForMaskVariousThr[j][t])
                        s = 'ColumnID\tCoord\tValue\n'
                        f.write(s)
                        for z in range(len(listsOfCoordsForMaskVariousThr[j][t])):
                            s = str(listsOfColumnIDsForMaskVariousThr[j][t][z]) + '\t' + str(listsOfCoordsForMaskVariousThr[j][t][z]) + '\t' + str(listsOfValuesForMaskVariousThr[j][t][z]) + '\n'                   
                            f.write(s)        
                        f.close()
                        
                
                # plot the point from the columns selected or rejected according to the AvgDivGradn
                ax1selected = fig.add_subplot(4,numOfRegions,j + 1 + numOfRegions)
                ax1RejectedS = fig.add_subplot(4,numOfRegions,j + 1 + numOfRegions * 2)
                ax1RejectedB = fig.add_subplot(4,numOfRegions,j + 1 + numOfRegions * 3)
                
                ax1selected.plot(listsOfCoordsForMaskVariousThrADG[j][t], listsOfValuesForMaskVariousThrADG[j][t], '.', c = colours[j], label = 'ROI '+ str(maskROIids[j]) + ' selected AvgDivGradn')
                ax1RejectedS.plot(listsOfCoordsForMaskVariousThrADGs[j][t], listsOfValuesForMaskVariousThrADGs[j][t], '.', c = colours[j], label = 'ROI '+ str(maskROIids[j]) + ' too small AvgDivGradn')
                ax1RejectedB.plot(listsOfCoordsForMaskVariousThrADGb[j][t], listsOfValuesForMaskVariousThrADGb[j][t], '.', c = colours[j], label = 'ROI '+ str(maskROIids[j]) + ' too high AvgDivGradn')
                
                ax1selected.set_title('Profile in maskROI %s, %s < AvgDivGradn < %s' %(maskROIids[j], divGradnThresholds[0], divGradnThresholds[1]))
                ax1selected.set_xlabel('Cortical depth') 
                ax1selected.set_ylabel('T2-nobias intensity')                 
                ax1selected.legend(loc='upper right', numpoints = 1)                   
                ax1RejectedS.set_title('Profile in maskROI %s, %s > AvgDivGradn' %(maskROIids[j], divGradnThresholds[0]))
                ax1RejectedS.set_xlabel('Cortical depth') 
                ax1RejectedS.set_ylabel('T2-nobias intensity')                 
                ax1RejectedS.legend(loc='upper right', numpoints = 1) 
                ax1RejectedB.set_title('Profile in maskROI %s, %s < AvgDivGradn' %(maskROIids[j], divGradnThresholds[1]))
                ax1RejectedB.set_xlabel('Cortical depth') 
                ax1RejectedB.set_ylabel('T2-nobias intensity')                 
                ax1RejectedB.legend(loc='upper right', numpoints = 1)                
                
                # and now plot the same to the last plot: all mask ROIs together
                axAllSelected.plot(listsOfCoordsForMaskVariousThrADG[j][t], listsOfValuesForMaskVariousThrADG[j][t], '.', label = 'ROI '+ str(maskROIids[j]) + ' selected AvgDivGradn')
                axAllRejectedS.plot(listsOfCoordsForMaskVariousThrADGs[j][t], listsOfValuesForMaskVariousThrADGs[j][t], '.', label = 'ROI '+ str(maskROIids[j]) + ' too small AvgDivGradn')
                axAllRejectedB.plot(listsOfCoordsForMaskVariousThrADGb[j][t], listsOfValuesForMaskVariousThrADGb[j][t], '.', label = 'ROI '+ str(maskROIids[j]) + ' too high AvgDivGradn')
                
                print '----------------------------------------------------- len(listsOfCoordsForMaskVariousThrADG[',j,'][',t,']) = ', len(listsOfCoordsForMaskVariousThrADG[j][t])
                numP += len(listsOfCoordsForMaskVariousThr[j][t])
                
                # write out a file with columnID, coord, value (corresponding to this mask ROI and size threshold)
                if len(listsOfCoordsForMaskVariousThrADG[j][t]) != len(listsOfColumnIDsForMaskVariousThrADG[j][t]):
                    print ' --------------------------- PROBLEM !'
                else:
                    # write out profiles if diameter is > 1
                    if columnDiameter > 1 and len(listsOfCoordsForMaskVariousThrADG[j][t]) > 0:                        
                        f = open(pathForFiles + '%s_%s_MaskROI%s_profiles%s_over_%s_avgDivGradn_%s.txt' %(realPatientID, realSide, str(maskROIids[j]), addedColumnsDiamName, str(minColumnSizes[t]), strT), "w")
                        #print 'write out a file %s_%s_MaskROI%s_profiles%s_over_%s.txt' %(realPatientID, realSide, str(maskROIids[j]), addedColumnsDiamName, str(minColumnSizes[t])), ' number of elements ', len(listsOfCoordsForMaskVariousThr[j][t])
                        s = 'ColumnID\tCoord\tValue\n'
                        f.write(s)
                        for z in range(len(listsOfCoordsForMaskVariousThrADG[j][t])):
                            s = str(listsOfColumnIDsForMaskVariousThrADG[j][t][z]) + '\t' + str(listsOfCoordsForMaskVariousThrADG[j][t][z]) + '\t' + str(listsOfValuesForMaskVariousThrADG[j][t][z]) + '\n'                   
                            f.write(s)        
                        f.close()
                        #TODO: need to write out profiles of columns with too low or too high AvgDivGradn??
                        
                        
                # TODO : delete it. only for verification
                # check if this particular cortical column's AvgDivGradn Value is between the thresholds, the plot it to the second row of the figure
                if len(listsOfColumnIDsForMaskVariousThrADG[j][t]) != 0:  
                    uniqueCols = np.unique(listsOfColumnIDsForMaskVariousThrADG[j][t])
                    for tt in range(len(uniqueCols)):
                        num = np.where(np.array(iDs) == uniqueCols[tt])[0][0]
                        avgDivGrandClass = listOfAvgDivGradsCateg[num]
                        #print 'num is ', num, ' roi ', uniqueCols[tt], ' AvgDivGradn ', listOfAvgDivGrads[num], ' divGradnClass ', listOfAvgDivGradsCateg[num]
                        #if (avgDivGrandClass == 20):         # this is the id of the current column. need to get the "class" value for the AvgDivGradn
                            #print ' should be accepted'
                            
            
            axAll.set_title('Profile in all maskROIs')
            axAll.set_xlabel('Cortical depth') 
            axAll.set_ylabel('T2-nobias intensity') 
            axAll.legend(loc='upper right', numpoints = 1)            
                            
            axAllSelected.set_title('Profile in all maskROIs, %s < AvgDivGradn < %s' %(divGradnThresholds[0], divGradnThresholds[1]))
            axAllSelected.set_xlabel('Cortical depth') 
            axAllSelected.set_ylabel('T2-nobias intensity') 
            axAllSelected.legend(loc='upper right', numpoints = 1)            
            axAllRejectedS.set_title('Profile in all maskROIs, %s > AvgDivGradn' %(divGradnThresholds[0]))
            axAllRejectedS.set_xlabel('Cortical depth') 
            axAllRejectedS.set_ylabel('T2-nobias intensity') 
            axAllRejectedS.legend(loc='upper right', numpoints = 1)            
            axAllRejectedB.set_title('Profile in all maskROIs, %s < AvgDivGradn' %(divGradnThresholds[1]))
            axAllRejectedB.set_xlabel('Cortical depth') 
            axAllRejectedB.set_ylabel('T2-nobias intensity') 
            axAllRejectedB.legend(loc='upper right', numpoints = 1)
            
            # modify the name if cortical columns communfor ROIs were included/excluded
            nameInclExcl = ''
            if ignoreCommunColums:
                #print 'Commun columns are ignored'
                nameInclExcl = 'exclCommun'
            else:
                #print 'Commun columns are INCLUDED'
                nameInclExcl = 'inclCommun'
            # do not save the picture if it is empty
            if numP != 0:
                plt.savefig(pathToColumnResults + '%s_%s_nobiasT2_ROIs%s' %(realPatientID, realSide, roiNames) + '%s_over%s_%s_avgDivGradn_%s.png' %(addedColumnsDiamName, minColumnSizes[t], nameInclExcl, strT), bbox_inches='tight')            
            plt.clf()
            plt.close()    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
