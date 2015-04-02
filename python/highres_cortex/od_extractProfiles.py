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
#python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_extractProfiles.py -p ad140157 -s R -d /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/ad140157/ -c 9

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
    

    
def extractProfilesInColumns(volCoord, volValue, volColumns, minColumnSize, volMask = None):
    """
    this function takes 2 volumes, one with values (T2) and another with corresponding coordinates (cortex
    depth measure from Yann's scripts.) It also takes a volume with "cortical columns" and will extract 
    lists of profiles for each of the columns separately.
    It also optionally takes a mask in which the profiles will be extracted.
    If no mask was given, profiles are extracted in the whole volume.
    Forfiles are extracted as a list of values and a list of coordinates.    
    """
    #print volCoord.header()
    
    arrCoord = np.asarray(volCoord)
    arrValue = np.asarray(volValue)
    arrColumns = np.asarray(volColumns)

    # apply the mask if it was given
    listOfSeparateCoords = []
    listOfSeparateValues = []
    # get columns in certain mask regions    
    listOfSeparateCoordsMask = []
    listOfSeparateValuesMask = []
    listOfROIsInMask = []
    listOfSizes = []
    
   
    if volMask is not None:
        mask = np.asarray(volMask)    
        
        # get the part (ROIs) of the volColumns where the mask is not zero
        roiColumns = arrColumns[mask != 0]
        
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
        #roiIds = np.unique(roiColumns[np.where(roiColumns > 0)])
        roiIds = np.unique(roiColumns[np.where(roiColumns > 0)])  # rois from the cortical columns
        roisMask = np.unique(mask[np.where(mask > 0)])  # rois from the mask
        listOfROIsInMask = [[] for x in range(len(roisMask) + 1)]
        print 'listOfROIsInMask = '
        print listOfROIsInMask
        
        #roiIds = np.unique(arrColumns)
        print 'found the following ROIs. their number  ', len(roiIds), ' from total number : ', len(np.unique(arrColumns))
        #print roiIds
        
        for i in roiIds:
            # print len(arrCoord), len(roiColumns)
            arrCoord1i = arrCoord1[roiColumns == i]
            arrValue1i = arrValue1[roiColumns == i]           
            coordsi = arrCoord1i[arrCoord1i != 0]
            valuesi = arrValue1i[arrCoord1i != 0]
            
            #print 'work with ROI ', i, ' # of points= ', len(coordsi)
            listOfSeparateCoords.append(coordsi)
            listOfSeparateValues.append(valuesi)   
            listOfSizes.append(len(coordsi))
                        
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
            if roiIds[i] != 0 and roiIds[i] in listOfROIsInMask[0] and roiIds[i] in listOfROIsInMask[1]:
                #print 'Column with ID ', roiIds[i], ' is in both regions'
                listOfROIsInMask[len(roisMask)].extend([roiIds[i]])               # a list of ROIs to ignore
                
                # check what percentage of this column is in each mask ROI
                numbers = []
                for k in roisMask:                
                    w = arrColumns[(arrColumns == roiIds[i]) & (mask == k)] # get voxels in the respective mask ROI and the respective cortical column
                    numbers.append(len(w))                
                print 'Column with ID ', roiIds[i], ' is in both regions. % in mask ROIs ', ' is ', numbers
                
                ## TODO! analyze these numbers
                # TODO: Denis : leave it like this!! eliminate columns that are in both ROIs!!! as Voronoi can also contain errors!
                ## compare the max. to the second largest number. e.g. If the max is > 4 times than the second max : use it for the respective ROI. else: split.
                #sortedNums = sorted(range(len(numbers)), key=lambda k: numbers[k])
                #maxPart = numbers[sortedNums[len(numbers) - 1]]
                #max2Part = numbers[sortedNums[len(numbers) - 2]]
                #ratio = maxPart/max2Part
                #print 'ratio = ', ratio
                
                ## e.g. the max is > 4 times than the second max : use it for the respective ROI. else: split.
                #if ratio > 5.0:
                    #print 'this column will belong to ROI # ', sortedNums[len(numbers) - 1]
                    ## need to delete this ROI from the list sortedNums[len(numbers) - 2]
                    #toDeleteFrom = sortedNums[len(numbers) - 2]
                    #listOfROIsInMask[toDeleteFrom].delete()
                    
                    ## or create a tuple : (column-ROI, mask-ROI)                    
                    
                #else:
                    #print 'need to split it'
                
                
                
                
            
        
    #else :   TODO!!!!!!!!!!!!!!!!!!!!!!!!!!
        #coords = arrCoord[arrCoord != 0]
        #values = arrValue[arrCoord != 0]
        
    # 2 arrays of coordinates and values. 
    res = []
    res.append(coords)
    res.append(values)
    res.append(roiIds)
    res.append(listOfSeparateCoords)
    res.append(listOfSeparateValues)
    res.append(listOfROIsInMask)
    res.append(listOfSizes)
    return(res)       
    

    
if __name__ == '__main__':
    
    realPatientID = None
    directory = None
    realSide = 'L'
    columnDiameter = None
    workOnLaptop = False
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
	    pathToNobiasT2 = pathToNobiasT2.replace('/neurospin/lnao/', '/nfs/neurospin/lnao/')
	    pathToNobiasT2_new = pathToNobiasT2_new.replace('/neurospin/lnao/', '/nfs/neurospin/lnao/')
  
    pathToCoord = directory + '%s_T1inT2_ColumnsCutNew20It/isovolume/' %(realPatientID)
    pathToMask = directory + '%s_T1inT2_ColumnsCutNew20It/' %(realPatientID)
    
    print 'volsCoord = ' , pathToCoord + 'pial-volume-fraction_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)
    print 'volsMask = ', pathToMask + 'voronoiCorr_%s_%s_cut_noSulci.nii.gz' %(realPatientID, realSide)
    
    volsCoord = glob.glob(pathToCoord + 'pial-volume-fraction_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide))
    volValue = aims.read(pathToNobiasT2 + '%s_NewNobiasT2_cropped.nii.gz' %(realPatientID))
    volValue2 = aims.read(pathToNobiasT2_new + '%s_NewT2_cropped.nii.gz' %(realPatientID))
    volsMask = glob.glob(pathToMask + 'voronoiCorr_%s_%s_cut_noSulci.nii.gz' %(realPatientID, realSide))

    # test if all data is available
    if len(volsCoord) != 1 or len(volsMask) != 1:
        # abort the calculation, as too many or not a single texture file was found
        print 'abort the calculation, as too many or not a single volsCoord and volsMask file was found'
        ss = directory + '%s_%s_statFileProfiles.txt' %(realPatientID, realSide)        
        print 'f = ', ss

        f = open(directory + '%s_%s_statFileProfiles.txt' %(realPatientID, realSide), "w")
        f.write('abort the calculation, as ' + str(len(volsCoord)) + ' volsCoord and ' + str(len(volsMask)) + ' volsMask files were found' + '\n')
        f.close()
        sys.exit(0)    
        f.close()        
        
    volCoord = aims.read(volsCoord[0])  
    volMask = aims.read(volsMask[0])
    
    
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
        vs = volValue.header()['voxel_size'][0:3]
        voxelSizeMax = np.max(vs)
        voxelSizeMin = np.min(vs)
        voxelSizeAvg = np.average(vs)
        heightMax = int(np.round(4.0 / voxelSizeMin))
        heightMin = int(np.round(2.0 / voxelSizeMax))
        print 'heights of the cortical columns found heightMin = ', heightMin, ' heightMax = ', heightMax
        #heights = range(heightMin, heightMax + 1)
        # just for test: to see the influence of really large columns
        heights = range(heightMin, heightMax + 8)
        print ' ##################################################### heights = ', heights
       # using this info: calculate minimal cortical column sizes
        
        # average size         
        minColumnSize = int(np.round(columnDiameter * columnDiameter / voxelSizeAvg / voxelSizeAvg / 4 * np.pi * np.average(heights))) 
        print '##################################################### calculated the avgColumnSize = ', minColumnSize

        # a list of sizes        
        minColumnSizes = [int(k * columnDiameter * columnDiameter / voxelSizeAvg / voxelSizeAvg / 4 * np.pi) for k in heights]
        print 'calculated a range of possible column sizes: '
        print minColumnSizes
        
        # start the processing for columns
        print '******************* columns diameter was given, then extract profiles in the columns ***************************'    
        # read in the columns file
        volName = ''
        if columnDiameter == 1:
            volName = directory + '%s_T1inT2_ColumnsCutNew20It/column_regions/' %(realPatientID) + 'merged_randomized_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide)
        else:
            volName = directory + '%s_T1inT2_ColumnsCutNew20It/column_regions/' %(realPatientID) + 'merged_randomized_%s_%s_cut_noSulci_extended_diam%s.nii.gz' %(realPatientID, realSide, str(columnDiameter))
            
        volsColumns = glob.glob(volName)

        if len(volsColumns) != 1:# abort the calculation, as too many or not a single columns file was found
            print 'abort the calculation, as too many or not a single volsColumns file was found'
            f = open(directory + '%s_%s_statFileProfiles.txt' %(realPatientID, realSide), "w")
            f.write('abort the calculation, as ' + str(len(volsColumns)) + ' volsColumns files were found' + '\n')
            f.close()
            sys.exit(0)
            f.close()

        volColumns = aims.read(volsColumns[0])  
        print 'volColumns = ', volsColumns[0]
        #result = extractProfilesInColumns(volCoord, volValue, volColumns, minColumnSize, volMask)
        # repeat for the NEW nobias images!
        result2 = extractProfilesInColumns(volCoord, volValue2, volColumns, minColumnSize, volMask)      
    
    # work now only with the new nobias T2!!!!! commented the work with the old corrected nobias T2  
    #coordinates = result[0]
    #intensities = result[1]

    ## plot the data
    #plt.plot(coordinates, intensities, '.', c = 'b')
    #plt.title('Profile in ROI')   # subplot 211 title
    #plt.xlabel('Cortical depth')
    #plt.ylabel('T2-nobias intensity')
    #plt.savefig(directory + '%s_%s_It20_nobiasT2vsCorticalDepthROI.png' %(realPatientID, realSide))

    coordinates2 = result2[0]
    intensities2 = result2[1]

    ## plot the data
    #plt.plot(coordinates2, intensities2, '.', c = 'r')
    #plt.title('Profile in ROI')   # subplot 211 title
    #plt.xlabel('Cortical depth')
    #plt.ylabel('T2-nobias intensity')
    #plt.savefig(directory + '%s_%s_It20_2nobiasT2vsCorticalDepthROI.png' %(realPatientID, realSide))
    
    #plt.clf()
    #plt.close()
    
    plt.plot(coordinates2, intensities2, '.', c = 'r')
    plt.title('Profile in all ROIs')   # subplot 211 title
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')
    plt.savefig(directory + '%s_%s_newNobiasT2ROIs.png' %(realPatientID, realSide))    
    plt.clf()
    plt.close()
  
    # save the data for further processing. TODO: find information about their coordinates!!
    #data1 = open(directory + '%s_%s_profiles.txt' %(realPatientID, realSide), "w")
    data2 = open(directory + '%s_%s_profiles2.txt' %(realPatientID, realSide), "w")
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
    
    addedColumnsDiamName = ''
    pathForFiles = directory
    print '############################################ columnDiameter = ', str(columnDiameter)
    if columnDiameter is not None:
        addedColumnsDiamName = '_diam%s' %(columnDiameter)
        # if the result was for the columns, create a separate folder for it
        pathForFiles = pathForFiles + 'diam%s/'%(columnDiameter)
        if not os.path.exists(pathForFiles):
            os.makedirs(pathForFiles)
        print '############################################ created directory = ', pathForFiles

        # read in the results for columns IDs in various regions
        iDsInMaskROIs = result2[5]
     
    # write out a file with all columns and their sizes
    dataS = open(pathForFiles + '%s_%s_ColumnSizes%s.txt' %(realPatientID, realSide, addedColumnsDiamName), "w")
    dataS.write('\tColumnID\tSize\n')
    sizeSorted = pathForFiles + 'sortedBySize/'
    if not os.path.exists(sizeSorted):
        os.makedirs(sizeSorted)
        
    for i in range(len(iDs)):
        #print 'i = ', i, ' work with id', iDs[i]
        currCoords = listOfCoords[i]
        currValues = listOfValues[i]
        plt.plot(currCoords, currValues, '.', c = 'b')
        plt.title('Profile in ROI %s' %(str(iDs[i])))   # subplot 211 title
        plt.xlabel('Cortical depth')
        plt.ylabel('T2-nobias intensity')
        
        dataS.write('\t%s\t%s\n' %(str(iDs[i]), len(currCoords)))

        # do not plot if it is zero
        if len(currCoords) != 0:
            plt.savefig(pathForFiles + '%s_%s_nobiasT2%s_ROI_' %(realPatientID, realSide, addedColumnsDiamName) + str(iDs[i]) + '.png')        
            # write out the same info, but sorted by sizes
            plt.savefig(sizeSorted + '%s_%s_nobiasT2%s_size%s_ROI_' %(realPatientID, realSide, addedColumnsDiamName, str(len(currCoords))) + str(iDs[i]) + '.png')
            data2i = open(pathForFiles + '%s_%s_profiles2%s_ROI_%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(iDs[i])), "w")
            data2i.write(headerLine + '\n')
            for j in range(len(currCoords)):
                data2i.write(str(j) + '\t' + str(currCoords[j]) + '\t' + str(currValues[j]) + '\n')    
            data2i.close()  
            
        plt.clf()
        plt.close()        
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
     
                
    # get results for various sizes of columns
    if columnDiameter is not None:
        # plot the histogram of column sizes
        sizes = result2[6]
        plt.hist(sizes, bins = 50)
        plt.title('Histogram of columns sizes for diameter %s ' %(columnDiameter))   # subplot 211 title
        plt.xlabel('Cortical columns sizes')
        plt.ylabel('Percentage')
        plt.savefig(directory + '%s_%s_histOfColumnSizes_diam%s.png' %(realPatientID, realSide, str(columnDiameter)))    
        plt.clf()
        plt.close()
        
        # plot the "theoretical" histogram of column heights
        theorHeights = [ int (w * 4 * voxelSizeAvg * voxelSizeAvg / np.pi / columnDiameter / columnDiameter) for w in sizes]
        plt.hist(theorHeights, bins = 50)
        plt.title('Histogram of theoretical heights of columns for diameter %s ' %(columnDiameter))   # subplot 211 title
        plt.xlabel('Theoretical heights of columns')
        plt.ylabel('Percentage')
        plt.savefig(directory + '%s_%s_histOfTheorHeights_diam%s.png' %(realPatientID, realSide, str(columnDiameter)))    
        plt.clf()
        plt.close()
        
        listsOfCoordsForLargeColumns = [[] for x in range(len(heights))]
        listsOfValuesForLargeColumns = [[] for x in range(len(heights))] 
        listsOfIDsForLargeColumns = [[] for x in range(len(heights))] 
        listsOfNumbersLargeColumns = [0] * len(heights)
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
                plt.savefig(directory + '%s_%s_newNobiasT2_diam%s_over%s.png' %(realPatientID, realSide, str(columnDiameter), minColumnSizes[j]))    
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
            dataID = open(pathForFiles + '%s_%s_IDs%s_over_%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(minColumnSizes[j])), "w")
            string = ''
            for l in listsOfIDsForLargeColumns[j]:
                string += str(l) + '\t' 
                
            if len(listsOfIDsForLargeColumns[j]) == 0:  # print 0
                string += '0\t'
            
            dataID.write(string)
            dataID.close()
            
        # now plot avg profiles for various mask regions
        listOfIDsToIgnore = iDsInMaskROIs[len(iDsInMaskROIs) - 1]
        listsOfCoordsForMask = [[] for x in range(len(iDsInMaskROIs) - 1)]
        listsOfValuesForMask = [[] for x in range(len(iDsInMaskROIs) - 1)] 
        
        # lists of ROIs of various sizes
        listsOfCoordsForMaskVariousThr = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)]
        listsOfValuesForMaskVariousThr = [[[] for y in range(len(minColumnSizes))] for x in range(len(iDsInMaskROIs) - 1)] 

        
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
                    #print 'iDs ='
                    #print iDs
                    #print 'k = ', k
                    indexOfThisID = np.where(iDs == k)[0][0]
                    #print 'indexOfThisID ', indexOfThisID
                    currCoords = listOfCoords[indexOfThisID]
                    currValues = listOfValues[indexOfThisID]                    
                    
                    for t in range(len(minColumnSizes)):
                        #print t, j
                        if len(currCoords) > minColumnSizes[t]:
                            listsOfCoordsForMaskVariousThr[j][t].extend(currCoords)
                            listsOfValuesForMaskVariousThr[j][t].extend(currValues)   
                            
                            
                            
        # plot the data for mask rois of large columns (various size thresholds). For ROIs separately, and together
        for t in range(len(minColumnSizes)):
            plt.plot()
            plt.title('Profile in Mask ROIs of large cortical columns')   
            plt.xlabel('Cortical depth')
            plt.ylabel('T2-nobias intensity')
            roiNames = ''

            for j in range(len(maskROIids)):
                roiNames += '_' + str(maskROIids[j])
                plt.plot(listsOfCoordsForMaskVariousThr[j][t], listsOfValuesForMaskVariousThr[j][t], '.', label = 'ROI '+ str(maskROIids[j]))
                #plt.plot(listsOfCoordsForMaskVariousThr[j][t], listsOfValuesForMaskVariousThr[j][t], '.', c = 'b') 
                #plt.title('Profile in Mask ROI of large cortical columns')   
                #plt.xlabel('Cortical depth')
                #plt.ylabel('T2-nobias intensity')
                #plt.savefig(directory + '%s_%s_nobiasT2_ROI_' %(realPatientID, realSide) + str(maskROIids[j]) + '%s_over%s.png' %(addedColumnsDiamName, minColumnSizes[t]))
                #plt.clf()
                #plt.close()
            
            # modify the name if cortical columns communfor ROIs were included/excluded
            nameInclExcl = ''
            if ignoreCommunColums:
                #print 'Commun columns are ignored'
                nameInclExcl = 'exclCommun'
            else:
                #print 'Commun columns are INCLUDED'
                nameInclExcl = 'inclCommun'
            # do not save the picture if it is empty
            if len(listsOfCoordsForMaskVariousThr[j][t]) != 0:
                plt.legend(loc='upper right', numpoints = 1)
                plt.savefig(directory + '%s_%s_nobiasT2_ROIs%s' %(realPatientID, realSide, roiNames) + '%s_over%s_%s.png' %(addedColumnsDiamName, minColumnSizes[t], nameInclExcl))
            
            plt.clf()
            plt.close()

    
                
    #if columnDiameter is not None:
        #plt.plot(coordsFromLargeColumns, valuesFromLargeColumns, '.', c = 'r')
        #plt.title('Profile in large cortical columns')   # subplot 211 title
        #plt.xlabel('Cortical depth')
        #plt.ylabel('T2-nobias intensity')
        #plt.savefig(directory + '%s_%s_newNobiasT2_ColumnsOver%s.png' %(realPatientID, realSide, minColumnSize))    
        #plt.clf()
        #plt.close()
        
        ## plot the data from all cortical columns and from the large ones together
        #plt.plot(coordinates2, intensities2, '.', c = 'b')
        #plt.title('Profile in ROI')   # subplot 211 title
        #plt.xlabel('Cortical depth')
        #plt.ylabel('T2-nobias intensity')
        #plt.plot(coordsFromLargeColumns, valuesFromLargeColumns, '.', c = 'r')
        #plt.savefig(directory + '%s_%s_newNobiasT2_AllvsColumnsOver%s.png' %(realPatientID, realSide, minColumnSize)) 
        #plt.clf()
        #plt.close()    
        #print 'number of large columns is ', len(iDsLarge), ' from ', len(iDs), ' from ', len(np.unique(np.asarray(volColumns)))
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
