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
    

    
def extractProfilesInColumns(volCoord, volValue, volColumns, volDivGradn, divGrThr, minColumnSize, volMask = None):
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
    
   
    if volMask is not None:
        mask = np.asarray(volMask)    
        
        # get the part (ROIs) of the volColumns where the mask is not zero
        roiColumns = arrColumns[mask != 0]
        
        # get these ROIs
        arrCoord1 = arrCoord[mask != 0]
        arrValue1 = arrValue[mask != 0]
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
            #print '------------------------- work with ROI ', i, ' # of points= ', len(coordsi), ' average divGradn = ', avgDiv, ' ----------------'
            #if math.isnan(avgDiv):
                #print len(coordsi), len(valuesi), len(divGradnsi), ' the whole digGradnList',  divGradnsi               
                
            listOfSeparateCoords.append(coordsi)
            listOfSeparateValues.append(valuesi)   
            listOfSizes.append(len(coordsi))
            listOfAvgDivGradn.append(avgDiv)
            
                   
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
        #print newVol.header()
        #aims.write(newVol, '/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/testColouredVol6.nii.gz')
        #print 'unique values:', np.unique(arrColouredVol)
        #aims.write(volColoured, '/neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/sg140335/divGradnClasses/sg140335_R_testColouredVol_01_01.nii.gz')
        print 'lower ', divGrThr[0], ' : ', counters[0]/float(len(roiIds)), ', lower ', divGrThr[1], ' : ', counters[1]/float(len(roiIds)), ', over ', divGrThr[1], ' : ', counters[2]/float(len(roiIds))   
        print 'lower ', divGrThr[0], ' : ', countVoxels[0]/float(sum(countVoxels)), ', lower ', divGrThr[1], ' : ', countVoxels[1]/float(sum(countVoxels)), ', over ', divGrThr[1], ' : ', countVoxels[2]/float(sum(countVoxels))
        #print 'forced exit'
        #sys.exit(0)
        
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
                    
                    ## TODO! analyze these numbers
                    # TODO: Denis : leave it like this!! eliminate columns that are in both ROIs!!! as Voronoi can also contain errors!
                    # if one of the rois contains 10 times more voxels than all the rest ROIs, than we can accept this column and say that it belongs to this roi
                    # e.g. if column has 1200 voxels in ROI 11 and 13 voxels in ROI 21  - should it be dismissed??
                    
                    maxN = np.max(numbers)
                    ratio = maxN / (np.sum(numbers) - maxN)
                    if ratio >= 10:
                        # accept this column to the ROI: where its most voxels are
                        maxId = numbers.index(maxN) # maxId gives the ROI in the mask to which this column should be attributed
                        print 'ratio = ', ratio, '  maxId = ', maxId, ' mask ROI ', roisMask[maxId]
                        
                        # do not put onto ignore list!
                        # TODO! need to delete this column id from the lists of other mask ROIs
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
    #divGradnThresholds = [-0.1, 0.1]
    #divGradnThresholds = [-0.2, 0.2]
    #divGradnThresholds = [-0.2, 0.3]
    #divGradnThresholds = [-0.3, 0.3]
    divGradnThresholds = [-0.2, 0.35]
    #divGradnThresholds = [-0.25, 0.35]
    #divGradnThresholds = [-0.5, 0.5] #just for test
    # TODO! find real intervals!! May be even different for PT and Heschl!!
#    corticalIntervals = [0, 0.1, 0.2, 0.5, 0.62, 0.82, 1.0]
#    corticalIntervals = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.62, 0.72, 0.82, 0.91, 1.0]
#    corticalIntervals = [0, 0.1, 0.2, 0.35, 0.5, 0.62, 0.72, 0.82, 1.0]



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
        heights = range(heightMin, heightMax + 5)
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
        
        # read in the div_gradn file (for measuring the flatness of the resoective column region)
        volDivGradn = aims.read(directory + '%s_T1inT2_ColumnsCutNew20It/heat/heat_div_gradn_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realPatientID, realSide))          
        print 'volDivGradn = ', directory + '%s_T1inT2_ColumnsCutNew20It/heat/heat_div_gradn_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realPatientID, realSide)
        #result = extractProfilesInColumns(volCoord, volValue, volColumns, minColumnSize, volMask)
        # repeat for the NEW nobias images!
        result2 = extractProfilesInColumns(volCoord, volValue2, volColumns, volDivGradn, divGradnThresholds, minColumnSize, volMask)   
    
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
    listOfAvgDivGrads = result2[7]
    volColouredByAvgDivGrads = result2[8]
    listOfAvgDivGradsCateg = result2[9]         # a list of 3 categories : 10, 29 and 30, given according to the 2 set thresholds

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
        
    # write out a file with the info: Column ID, size, maskROI (belongs to it or not), To be ignored
    maskArray = np.asarray(volMask)
    maskROIids = np.unique(maskArray[np.where(maskArray > 0)])            
    dataS = open(pathForFiles + '%s_%s_ColumnInfo%s.txt' %(realPatientID, realSide, addedColumnsDiamName), "w")
    headerInfo = '\tColumnID\tSize\tAvgDivGradn'
    
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
       
    #test whether the length of ids and avgdivgradns is the same
    if len(iDs) != len(listOfAvgDivGrads):
        print 'lengths are unequal!!! ', len(iDs), len(listOfAvgDivGrads), ' force exit'
        sys.exit(0)
    
    for i in range(len(iDs)):
        print '------------------ i = ', i, ' work with id', iDs[i]
        currCoords = listOfCoords[i]
        currValues = listOfValues[i]
        
        # get the info on the gradient
        currGradn = listOfAvgDivGrads[i]
        strCurrGradn = "%.3f" % currGradn
        strCurrGradnStr = str(strCurrGradn).replace('.','')   
        currLine = '\t%s\t%s\t%s' %(str(iDs[i]), len(currCoords), strCurrGradn)
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
        #print 'forced exit'     #TODO: delete when not needed!
        #sys.exit(0)

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
                    #print 'iDs =', iDs
                    #print 'k = ', k
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
                plt.savefig(directory + '%s_%s_nobiasT2_ROIs%s' %(realPatientID, realSide, roiNames) + '%s_over%s_%s.png' %(addedColumnsDiamName, minColumnSizes[t], nameInclExcl), bbox_inches='tight')            
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
                plt.savefig(directory + '%s_%s_nobiasT2_ROIs%s' %(realPatientID, realSide, roiNames) + '%s_over%s_%s_avgDivGradn_%s.png' %(addedColumnsDiamName, minColumnSizes[t], nameInclExcl, strT), bbox_inches='tight')            
            plt.clf()
            plt.close()    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
