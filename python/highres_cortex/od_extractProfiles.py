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
    

    
def extractProfilesInColumns(volCoord, volValue, volColumns, volMask = None):
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
        roiIds = np.unique(roiColumns)
        #roiIds = np.unique(arrColumns)
        print 'found the following ROIs. their number  ', len(roiIds), ' from total number : ', len(np.unique(arrColumns))
        #print roiIds
        
        for i in roiIds:
            # print len(arrCoord), len(roiColumns)
            #arrCoord1i = arrCoord[roiColumns == i]
            arrCoord1i = arrCoord1[roiColumns == i]
            #arrValue1i = arrValue[roiColumns == i]
            arrValue1i = arrValue1[roiColumns == i]
           
            coordsi = arrCoord1i[arrCoord1i != 0]
            valuesi = arrValue1i[arrCoord1i != 0]
            
            print 'work with ROI ', i, ' # of points= ', len(coordsi)
            listOfSeparateCoords.append(coordsi)
            listOfSeparateValues.append(valuesi)      
            
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
    return(res)       
    

    
if __name__ == '__main__':
    
    realPatientID = None
    directory = None
    realSide = 'L'
    columnDiameter = None
    workOnLaptop = False
    pathToNobiasT2 = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S4/'
    pathToNobiasT2_new = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S16/'


    parser = OptionParser('Extract profiles from T2 nobias data using cortex-density-coordinates in ROIs')    
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    parser.add_option('-s', dest='realSide', help='Hemisphere to be processed: L or R. L is default')   
    parser.add_option('-d', dest='directory', help='directory')
    parser.add_option('-c', dest='columnDiameter', help='columnDiameter to work with')
    parser.add_option('-l', dest='workOnLaptop', action = 'store_true', help='Select if working on laptop (neurospin DB location is different. False is default') options, args = parser.parse_args(sys.argv)
    options, args = parser.parse_args(sys.argv)print options
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
        
    if options.columnDiameter is not None:
        columnDiameter = options.columnDiameter

    if options.workOnLaptop is not None:
	workOnLaptop = options.workOnLaptop      
	# if true, then processes are run on the laptop. Change locations of neurospin DBs
	pathToNobiasT2 = pathToNobiasT2.replace('/neurospin/lnao/', '/nfs/neurospin/lnao/')
	pathToNobiasT2_new = pathToNobiasT2_new.replace('/neurospin/lnao/', '/nfs/neurospin/lnao/')
  
    pathToCoord = directory + '%s_T1inT2_ColumnsCutNew20It/isovolume/' %(realPatientID)
    pathToNobiasT2 = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S4/'
    pathToNobiasT2_new = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/T2_nobias_FR5S16/'
    pathToMask = directory + '%s_T1inT2_ColumnsCutNew20It/' %(realPatientID)
    
    volsCoord = glob.glob(pathToCoord + 'pial-volume-fraction_%s_%s_cut_noSulci_extended.nii.gz' %(realPatientID, realSide))
    volValue = aims.read(pathToNobiasT2 + '%s_NewNobiasT2_cropped.nii.gz' %(realPatientID))
    volValue2 = aims.read(pathToNobiasT2_new + '%s_NewT2_cropped.nii.gz' %(realPatientID))
    volsMask = glob.glob(pathToMask + 'voronoiCorr_%s_%s_cut_noSulci.nii.gz' %(realPatientID, realSide))

    # test if all data is available
    f = open(directory + '%s_%s_statFileProfiles.txt' %(realPatientID, realSide), "w")
    if len(volsCoord) != 1 or len(volsMask) != 1:
        # abort the calculation, as too many or not a single texture file was found
        print 'abort the calculation, as too many or not a single volsCoord and volsMask file was found'
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
    else:
        print '******************* columns diameter was given, then extract profiles in the columns ***************************'    
        # read in the columns file
        volsColumns = glob.glob(directory + '%s_T1inT2_ColumnsCutNew20It/column_regions/' %(realPatientID) + 'merged_randomized_%s_%s_cut_noSulci_extended_diam%s.nii.gz' %(realPatientID, realSide, str(columnDiameter)))

        if len(volsColumns) != 1:# abort the calculation, as too many or not a single columns file was found
            print 'abort the calculation, as too many or not a single volsColumns file was found'
            f.write('abort the calculation, as ' + str(len(volsColumns)) + ' volsColumns files were found' + '\n')
            f.close()
            sys.exit(0)

        volColumns = aims.read(volsColumns[0])  
        print 'volColumns = ', volsColumns[0]
        #result = extractProfilesInColumns(volCoord, volValue, volColumns, volMask)
        # repeat for the NEW nobias images!
        result2 = extractProfilesInColumns(volCoord, volValue2, volColumns, volMask)

    f.close()        
    
    # work now only with the new nobias T2!!!!! commented the work with the old corrected nobias T2  
    #coordinates = result[0]
    #intensities = result[1]

    ## plot the data
    #plt.plot(coordinates, intensities, '.', c = 'b')
    #plt.title('Profile in ROI')   # subplot 211 title
    #plt.xlabel('Cortical depth')
    #plt.ylabel('T2-nobias intensity')
    #plt.savefig(directory + '%s_%s_It20_nobiasT2vsCorticalDepthROI.png' %(realPatientID, realSide))

    
    # repeat for the NEW nobias images!
    result2 = extractProfiles(volCoord, volValue2, volMask)
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
    plt.title('Profile in ROI')   # subplot 211 title
    plt.xlabel('Cortical depth')
    plt.ylabel('T2-nobias intensity')
    plt.savefig(directory + '%s_%s_It20_newNobiasT2vsCorticalDepthROI.png' %(realPatientID, realSide))    
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
    if columnDiameter is not None:
        addedColumnsDiamName = '_diam%s_' %(columnDiameter)
        # if the result was for the columns, create a separate folder for it
        pathForFiles = pathForFiles + 'diam%s/'%(columnDiameter)
        if not os.path.exists(pathForFiles):
            os.makedirs(pathForFiles)
    
    for i in range(len(iDs)):
        #print 'i = ', i, ' work with id', iDs[i]
        currCoords = listOfCoords[i]
        currValues = listOfValues[i]
        plt.plot(currCoords, currValues, '.', c = 'b')
        plt.title('Profile in ROI')   # subplot 211 title
        plt.xlabel('Cortical depth')
        plt.ylabel('T2-nobias intensity')
        plt.savefig(pathForFiles + '%s_%s_It20_nobiasT2vsCorticalDepth%s_ROI_' %(realPatientID, realSide, addedColumnsDiamName) + str(iDs[i]) + '.png')
        plt.clf()
        plt.close()
        
        data2i = open(pathForFiles + '%s_%s_profiles2%s_ROI_%s.txt' %(realPatientID, realSide, addedColumnsDiamName, str(iDs[i])), "w")
        data2i.write(headerLine + '\n')
        for j in range(len(currCoords)):
            data2i.write(str(j) + '\t' + str(currCoords[j]) + '\t' + str(currValues[j]) + '\n')
            
        data2i.close()


        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    