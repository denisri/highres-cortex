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


# an example how to run this script
#od_mainCorticalColumns.py -p ad140157 -c -d /volatile/od243208/brainvisa_manual/ad140157_cut_T1inT2_Columns/ -r

# od_mainCorticalColumns.py -p ml140175 -s L -c True -i /volatile/od243208/brainvisa_db_morphologist/dysbrain/ml140175/t1mri/reversed_t1map_2/default_analysis/segmentation/Lgw_interface_ml140175.nii.gz -d /volatile/od243208/brainvisa_manual/ml140175_test/

# this is the main script to run on a classified GW volume
# it launches scripts by Yann Leprince: dist, heat, isovolume, column-regions to compute 'cortical columns'


# od_mainCorticalColumns.py -p md140208 -s L -d /volatile/od243208/brainvisa_manual/md140208_T1inT2_ColumnsNew/ -e -r

from soma import aims, aimsalgo
from scipy.stats import mode
import sys, glob, os, subprocess, sys, time, timeit
import numpy as np
from optparse import OptionParser
import highres_cortex.od_cutOutRois

#read in the path and the directory
brainvisa_db_neurospin = '/neurospin/lnao/dysbrain/brainvisa_db_morphologist/dysbrain/'
brainvisa_raw_niftis = '/neurospin/lnao/dysbrain/raw_niftis/'
pathToTextures = '/neurospin/lnao/dysbrain/randomized_flipped_data/manual_work/'
pathToTrm = '/neurospin/lnao/dysbrain/imagesInNewT2Space_LinearCropped10/'  
patientID = None              # subject000
realSide = 'L'
hemisphere = 'left'

realPatientID = None  # ac140155
pathToClassifFile = None #'/volatile/od243208/brainvisa_manual/%s/GWNoInnerSulciSkel_%s_%s.nii.gz' %(realPatientID, realPatientID, realSide)
pathToT2File = None
data_directory = None
keyWord = None
eliminateSulci = False
pathToSulciFile = None
cutOut = False           # perform Voronoi on the seeds from the labelled texture and apply Yann's methods on the cut out region
toT2 = False           # transform volumes to T2 space or not. Recently decided to perform this transformation after Voronoi and subtracting sulci skeletons
workOnLaptop = False
workOnT1inT2Space = False
numberOfIt = 3   # number of iterations to extend (dilate) the selected regions

parser = OptionParser('Calculate iso-volumes and cortical columns')
parser.add_option('-p', dest='realPatientID', help='Patient ID')   # if nothing is given: exit
parser.add_option('-s', dest='realSide', help='Hemisphere to be processed: L or R. L is default')   
parser.add_option('-c', dest='cutOut', action = 'store_true', help='Work with the selected cortical regions. False is default')  
#parser.add_option('-t', dest='toT2', action = 'store_true', help='Transform all the volumes into T2 space. False is default')    
# obsolete option: do not need to transform to T2 anymore, as we work with images already transformed to the new T2 space
parser.add_option('-i', dest='pathToClassifFile', help='Path to the volume with labeled cortex (100), and white matter (200)')   # if nothing is given: exit
#parser.add_option('-j', dest='pathToT2File', help='Path to the original T2 volume')   # if nothing is given: exit
parser.add_option('-d', dest='data_directory', help='directory for the results') # if nothing is given exit
parser.add_option('-e', dest='eliminateSulci', action = 'store_true', help='Eliminate sulci skeletons from the GW segmentation volume. False is default') 
parser.add_option('-k', dest='keyWord', help='KeyWord for the result files (the patient ID and the hemisphere are set by default)') 
parser.add_option('-l', dest='workOnLaptop', action = 'store_true', help='Select if working on laptop (neurospin DB location is different. False is default') 
parser.add_option('-r', dest='workOnT1inT2Space', action = 'store_true', help='Select if working on with T1 images that were resampled into T2 space and then processed by Morphologist. False is default') 


##################### for tests ############
#realPatientID = 'ml140175'
#realSide = 'L'
#cutOut = True
#pathToClassifFile = brainvisa_db_neurospin +  'ml140175/t1mri/reversed_t1map_2/default_analysis/segmentation/Lgw_interface_ml140175.nii.gz'
#data_directory = '/volatile/od243208/brainvisa_manual/ml140175_test/' 
#eliminateSulci = True
#keyWord = 'ml140175_L'
############################################


options, args = parser.parse_args(sys.argv)
print options
print args

if options.realPatientID is None:
    print >> sys.stderr, 'New: exit. No patient ID was given'
    sys.exit(1)
else:
    realPatientID = options.realPatientID

if options.data_directory is None:
    print >> sys.stderr, 'New: exit. No directory for results was given'
    sys.exit(1)
else:
    data_directory = options.data_directory
    
if options.realSide is not None:
    realSide = options.realSide

if options.workOnT1inT2Space is not None:
    workOnT1inT2Space = options.workOnT1inT2Space      
    # if true, then processes are run on (T1 resampled in T2 space). Change locations of neurospin DBs
    brainvisa_db_neurospin = '/neurospin/lnao/dysbrain/brainvisa_db_highresLinearCropped10/dysbrain/'   
    numberOfIt = 30   # number of iterations to extend (dilate) the selected regions
    # for high resolution images increased the number of iterations! for better avoidance of border effects
    print 'numberOfIt = ', numberOfIt
    
if options.workOnLaptop is not None:
    workOnLaptop = options.workOnLaptop      
    # if true, then processes are run on the laptop. Change locations of neurospin DBs
    brainvisa_db_neurospin = brainvisa_db_neurospin.replace('/neurospin/lnao/', '/volatile/od243208/neurospin/lnao/')
    brainvisa_raw_niftis = brainvisa_raw_niftis.replace('/neurospin/lnao/', '/volatile/od243208/neurospin/lnao/')
    pathToTextures = pathToTextures.replace('/neurospin/lnao/', '/volatile/od243208/neurospin/lnao/')
    pathToTrm = pathToTrm.replace('/neurospin/lnao/', '/volatile/od243208/neurospin/lnao/')
    #brainvisa_db_neurospin = '/volatile/od243208/neurospin/lnao/dysbrain/brainvisa_db_morphologist/dysbrain/'
    #brainvisa_raw_niftis = '/volatile/od243208/neurospin/lnao/dysbrain/raw_niftis/'
    #pathToTextures = '/volatile/od243208/neurospin/lnao/dysbrain/randomized_flipped_data/manual_work/'

if options.pathToClassifFile is None:   # take the 'standard file'
    pathToClassifFile = brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/default_analysis/segmentation/%sgrey_white_%s.nii.gz' %(realSide, realPatientID)    
    print 'Took the standard classification file: ', pathToClassifFile
else:
    pathToClassifFile = options.pathToClassifFile
 
#if options.pathToT2File is None:   # take the 'standard file'
    #finder2 = aims.Finder()
    #pathToT2File = glob.glob(pathToTrm + 'T2/%s_NewT2_cropped.nii.gz' %(realPatientID))[0]   # ok only if there is one t2-file, or the first one of the t2-files is correct
    #finder2.check(pathToT2File)
    #print 'Took the nobias-T2-cropped file: ', pathToT2File
#else:
    #pathToT2File = options.pathToT2File

if options.keyWord is None:
    keyWord = '%s_%s' %(realPatientID, realSide)
else:
    keyWord = options.keyWord
    
if options.eliminateSulci is not None:
    eliminateSulci = options.eliminateSulci    
    
if options.cutOut is not None:
    cutOut = options.cutOut    
    
#if options.toT2 is not None:
    #toT2 = options.toT2    
        
    
# in the given directory create the subdirectory for the results
if not os.path.exists(data_directory):
    os.makedirs(data_directory)

pathToClassifFileOriginal = pathToClassifFile
    



############# todo before this processing : ###############
###################### check: work in space T1 or T2 #######################################################################
###################### check: work with the whole cortex or with the 'cut' and 'dilated' region? ###########################
###################### according to this: launch the cutRoi script, adapt the keyWord and the destination data directory ###
############################################################################################################################
    

# decided to take T1 segmentation, (: the given pathToClassifFile)
# 1. either cut out the regions of interest or not,     -> update the keyWord
# 2. eliminate sulci skeletons, -> update the keyWord
# 3. transform it and resample into T2 space
# 4. apply Yann's scripts
    
    
############################# 1. either cut out the regions of interest or not, update the keyWord ##############################
print 'cutOut is ', str(cutOut), 'type(cutOut) is ', type(cutOut)
if cutOut is True:
    print '###################################### cutOut ROIs ###################################################################'
    keyWord += '_cut'
    print 'updated keyWord is : ', keyWord
    # take the seeds from the texture and perform the Voronoi classification of the voxels
    print pathToClassifFile
    volGWBorder = aims.read(pathToClassifFile, 1)
    # find the path to the texture file
    if realSide == 'L':
        hemisphere = 'left'
    else:
        hemisphere = 'right'
        
    filesTex = glob.glob(pathToTextures + '%s/%s/' %(realPatientID, hemisphere) + 'subject[0-9]*_side[0-1]_textureNEWLITTLEROI.gii')
#    filesTex = glob.glob(pathToTextures + '%s/%s/' %(realPatientID, hemisphere) + 'subject[0-9]*_side[0-1]_textureNEWLITTLEROI_NearOccipital.gii')
    if len(filesTex) != 1:
        # abort the calculation, as too many or not a single texture file was found
        print 'abort the calculation, as too many or not a single texture file was found'
        statFileName = data_directory + "statFile_%s" %(keyWord)

        # note if the processing was performed on laptop:
        if workOnLaptop:
            statFileName += '_laptop.txt'
        else:
            statFileName += '.txt'

        f = open(statFileName, "w")
        f.write('abort the calculation, as ' + str(len(filesTex)) + ' texture files were found' + '\n')
        f.close()
        sys.exit(0)
    
    fileTex = filesTex[0]
    print 'found the texture file : ', fileTex
    texture = aims.read(fileTex) #    subject012_side0_texture.gii    
    # find the hemisphere file
    fileHemi = pathToTrm + 'Hemi/%s_hemi_%s_T1inNewT2_cropped.gii' %(realPatientID, realSide)

    # brainvisa_db_neurospin + '%s/t1mri/reversed_t1map_2/default_analysis/segmentation/mesh/%s_%shemi.gii' %(realPatientID, realPatientID, realSide)        
    print 'found the hemisphere file : ', fileHemi
    volHemi = aims.read(fileHemi)
    
    #################### problem!!! Texture is still in the "old space". Need to transform it to the new space
    # read in the transformation file (from T1 into T2 space)
   # pathTot1_to_t2 = pathToTrm + realPatientID + '/%s_t1_to_t2.trm' % (realPatientID)
  #  transfT1toT2 = aims.read(pathTot1_to_t2)
    ###################### todo!!!!!!!!!!!!!!

    
    
    # perform the Voronoi classification in the given GW segmentation volume using the seeds from the texture    
    #print volGWBorder.header()
    print '######################### start Voronoi #################################'
    volVoronoi = highres_cortex.od_cutOutRois.voronoiFromTexture(volGWBorder, texture, volHemi, 0, data_directory, keyWord)
    # created the Voronoi classification that was  'cleaned' from the 'zero' value from texture
    
    # update the GW classif file, as now we will work with the selected regions
    volGW = aims.read(pathToClassifFile)
    volGWCut = highres_cortex.od_cutOutRois.excludeROI(volGW, volVoronoi, 100, 0)
    pathToClassifFile = data_directory + 'GWsegm_%s.nii.gz' %(keyWord)
    aims.write(volGWCut, pathToClassifFile)
    
    
############################# 2. eliminate sulci skeletons if requested . update the keyWord #################################
print 'eliminateSulci is ', eliminateSulci, 'type(eliminateSulci) is ', type(eliminateSulci)
if eliminateSulci is True:
    print '###################################### eliminate sulci skeletons #################################################'
    pathToSulciFile = brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/default_analysis/folds/3.1/default_session_auto/segmentation/%sSulci_%s_default_session_auto.nii.gz' %(realSide, realPatientID)
    print 'found the sulci skeletons file : ', pathToSulciFile
    volSulci = aims.read(pathToSulciFile)
    arrSulci = np.array(volSulci, copy = False)

    # if Voronoi was performed : eliminate sulci from there. Then correct the volume. Then project onto GW-classification
    if cutOut is True:
        print '############################ eliminate sulci skeletons from Voronoi classification ##########################'
        volVor = aims.read(data_directory +  'voronoi_%s.nii.gz' %(keyWord))
        arrVor = np.array(volVor, copy = False)        
        arrVor[arrSulci != 0] = 0
        keyWord += '_noSulci'
        print 'updated keyWord is : ', keyWord

        aims.write(volVor, data_directory +  'voronoi_%s.nii.gz' %(keyWord))        

        print '######################## correct voxel classification after sulci eliminated from Voronoi  #####################'
        volVor = aims.read(data_directory +  'voronoi_%s.nii.gz' %(keyWord))
        volVorCorr = highres_cortex.od_cutOutRois.correctVoxelLabels(volVor, data_directory +  'voronoi_%s.nii.gz' %(keyWord), data_directory, keyWord, 8, 6)        
        aims.write(volVorCorr, data_directory +  'voronoiCorr_%s.nii.gz' %(keyWord))              
        
        ################### just for test: study diff #####################
        #arrVor1 = np.array(aims.read(data_directory +  'voronoi_%s.nii.gz' %(keyWord)), copy = False)
        #arrVor2 = np.array(aims.read(data_directory +  'voronoiCorr_%s.nii.gz' %(keyWord)), copy = False)
        #diff = np.where(arrVor1 != arrVor2)
        #print '######################## diff between Voronoi and VoronoiCorr   #####################', len(diff[0])

        #for i in range(len(diff[0])):
            #x = diff[0][i]
            #y = diff[1][i]
            #z = diff[2][i]
            ##print 'diff values at ', x, y, z, 'Voronoi = ', volVor.value(x, y, z),  'VoronoiCorr = ', volVorCorr.value(x, y, z)
        ####################################################################
 
        # update the GW classif file, as now we will work with the selected regions that were corrected
        volGW = aims.read(pathToClassifFile)
        volGWCut = highres_cortex.od_cutOutRois.excludeROI(volGW, volVorCorr, 100, 0)
        pathToClassifFile = data_directory + 'GWsegm_%s.nii.gz' %(keyWord)
        aims.write(volGWCut, pathToClassifFile)    
        print 'volVorCorr ', data_directory +  'voronoiCorr_%s.nii.gz' %(keyWord)
        
        # take the Voronoi-volume and dilate it ############################
        # therefor need to take the original GW volume and subtract sulci from there.   
        volGW = aims.read(pathToClassifFileOriginal)
        arrGW = np.array(volGW, copy = False)
        arrGW[arrSulci != 0] = 0
        # need to change a keyWord to save this file. delete '_cut' from there
        k1 = keyWord.replace('_cut','')
        pathToClassifFileOriginal = data_directory + 'GWsegm_%s.nii.gz' %(k1)
        aims.write(volGW, pathToClassifFileOriginal)          
        
      ##  numberOfIt = 1,2,3          # not enough
        #numberOfIt = 4          # looks good!
        #volVorCorrDil = highres_cortex.od_cutOutRois.dilateRoiConnectivity26(volVorCorr, volGW, [11, 21], 0, 100, numberOfIt)
        #aims.write(volVorCorrDil, data_directory +  'voronoiCorrDil%s_%s.nii.gz' %(numberOfIt, keyWord))      
        
        ######################### just for test:  study the 1st step of dilation ###################################
        #volLabel11Conn6Step1 = highres_cortex.od_cutOutRois.labelDirect6Neighbors(volVorCorr, volGW, 11, 0, 100, 12)
        #aims.write(volLabel11Conn6Step1, data_directory +  'voronoiCorrDil6_label11_step1_%s.nii.gz' %(keyWord))      

        #volLabel11Conn26Step1 = highres_cortex.od_cutOutRois.labelDirect26Neighbors(volVorCorr, volGW, 11, 0, 100, 12)
        #aims.write(volLabel11Conn26Step1, data_directory +  'voronoiCorrDil26_label11_step1_%s.nii.gz' %(keyWord))      

        ## study the 2nd step
        #volLabel11Conn6Step2 = highres_cortex.od_cutOutRois.labelDirect6Neighbors(volLabel11Conn6Step1, volGW, 12, 0, 100, 12)
        #aims.write(volLabel11Conn6Step2, data_directory +  'voronoiCorrDil6_label11_step2_%s.nii.gz' %(keyWord))      

        #volLabel11Conn26Step2 = highres_cortex.od_cutOutRois.labelDirect26Neighbors(volLabel11Conn26Step1, volGW, 12, 0, 100, 12)
        #aims.write(volLabel11Conn26Step2, data_directory +  'voronoiCorrDil26_label11_step2_%s.nii.gz' %(keyWord))      
        ##########################################################################################################
        
        # looks like 6-connectivity is inappropriate: too slow, and such a precision is not needed
        # and 3 iterations seem to be sufficiet, too        
        #iters = range(1,6)
        #for numberOfIt in iters:                    
            #volVorCorrDil26 = highres_cortex.od_cutOutRois.dilateRoiConnectivity26(volVorCorr, volGW, [11, 21], 0, 100, numberOfIt)
            #aims.write(volVorCorrDil26, data_directory +  'voronoiCorrDil26_it%s_%s.nii.gz' %(numberOfIt, keyWord))      
            
            ##volVorCorrDil6 = highres_cortex.od_cutOutRois.dilateRoiConnectivity6(volVorCorr, volGW, [11, 21], 0, 100, numberOfIt)
            ##aims.write(volVorCorrDil6, data_directory +  'voronoiCorrDil6_it%s_%s.nii.gz' %(numberOfIt, keyWord))  
            
        ###################### just for test ################################
        #arrVorCorr = np.array(volVorCorr, copy = False)
        #fileTex = glob.glob(pathToTextures + '%s/%s/' %(realPatientID, hemisphere) + 'subject[0-9]*_side[0-1]_texture.gii')[0]
        #print fileTex
        #texture = aims.read(fileTex) #    subject012_side0_texture.gii 
        #arrTex = np.array(texture, copy = False)
        #print '##################### just for test ### ', 'arrTex', np.unique(arrTex), 'arrVorCorr', np.unique(arrVorCorr)

        arrVorCorr = np.array(volVorCorr, copy = False)
        labelsUniq = np.unique(arrVorCorr)      # a list of unique labels to process, except background
        labelsUniq = labelsUniq[labelsUniq != 0]
        volVorCorrDil26 = highres_cortex.od_cutOutRois.dilateRoiConnectivity26(volVorCorr, volGW, labelsUniq, 0, 100, numberOfIt)
        aims.write(volVorCorrDil26, data_directory +  'voronoiCorrDil26_it%s_%s.nii.gz' %(numberOfIt, keyWord))      
        
        # now we can cut the GW volume to the "extended_cut", as now we will work with the selected regions that were corrected and extended. update the keyword, too
        volGW = aims.read(pathToClassifFileOriginal)
        volGWCut = highres_cortex.od_cutOutRois.excludeROI(volGW, volVorCorrDil26, 100, 0)
        keyWord += '_extended'
        pathToClassifFile = data_directory + 'GWsegm_%s.nii.gz' %(keyWord)
        aims.write(volGWCut, pathToClassifFile)     

        
    else: # if no cutting out was requested: directly subtract sulci from GW-classification
        print '############################ eliminate sulci skeletons from full GW classification ##########################'    
        volGW = aims.read(pathToClassifFile)
        arrGW = np.array(volGW, copy = False)
        arrGW[arrSulci != 0] = 0
        keyWord += '_noSulci'
        print 'updated keyWord is : ', keyWord
        pathToClassifFile = data_directory + 'GWsegm_%s.nii.gz' %(keyWord)
        aims.write(volGW, pathToClassifFile)  
    

############################# 3. do T2 transformation if requested . update the keyWord. Take original coordinates from the T1 volume #######################
#print 'toT2 is ', toT2, 'type(toT2) is ', type(toT2)
#if toT2 is True:
    #volT2 = aims.read(pathToT2File)             # find the original T2 file to compute the transformation
    #volClassif = aims.read(pathToClassifFile)
    #arrClassif = np.array(volClassif, copy = False)
    #print 'cortex voxels in arrClassif: ', len(np.where(arrClassif == 100)[0])
    #volT1 = aims.read(brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/%s.nii.gz' %(realPatientID))
    #print '################################################   start transformation into T2 of the file ', pathToClassifFile
    #print '################################################   using the T2   ', pathToT2File   
    #print '################################################   using the T1   ', brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/%s.nii.gz' %(realPatientID)   
    #volClassifT2 = highres_cortex.od_cutOutRois.transformResampleVolFromVol(volClassif, volT2, volT1) 
    #pathCopy = pathToClassifFile
    #keyWord += '_T2'
    #pathToClassifFile = data_directory + 'GWsegm_%s.nii.gz' %(keyWord)
    #print '################################################   the updated classif file in T2 space is  ', pathToClassifFile
    #aims.write(volClassifT2, pathToClassifFile)
    ## just for test: 
    #arrClassifT2 = np.array(aims.read(pathToClassifFile), copy = False)
    #print 'cortex voxels in arrClassifT2: ', len(np.where(arrClassifT2 == 100)[0])  
    
    ################## just for test: needed to find problems with transformation done/or not donr in T1 and GW classification #########    
    ##volClassif = aims.read(pathCopy)
    ##volClassifT2_2 = highres_cortex.od_cutOutRois.transformResampleVolFromVol(volClassif, volT2) 
    ##aims.write(volClassifT2_2, data_directory + 'GWsegm_%s_notFromT1.nii.gz' %(keyWord))
    ### just for test: 
    ##arrClassifT2_2 = np.array(aims.read(data_directory + 'GWsegm_%s_notFromT1.nii.gz' %(keyWord)), copy = False)
    ##print 'cortex voxels in arrClassifT2_2: ', len(np.where(arrClassifT2_2 == 100)[0])  

    
    ############################ just for test. T2 transform NOT from T1 volume. Or the 'manual' one ####################
    ##volClassif = aims.read(pathCopy)
    ##print '########################### just for test ################ : manually transform into T2 '
    ##print 'pathToClassifFile ', pathCopy
    ##print 'T1 file ', brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/%s.nii.gz' %(realPatientID)
    ##finder = aims.Finder()
    ##finder.check(brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/%s.nii.gz' %(realPatientID))
    ##transformT1 = finder.header()['transformations'][0]
    ##affineTransformT1 = aims.AffineTransformation3d(transformT1)
    ### get T2 transformation without reading a file
    ##finder2 = aims.Finder()
    ##fileT2 = glob.glob(pathToT2File)[0]   # ok only if there is one t2-file, or the first one of the t2-files is correct
    ##print 'T2 file ', pathToT2File
    ##finder2.check(fileT2)
    ##transformT2 = finder2.header()['transformations'][0]
    ##dimsT2 = finder2.header()['volume_dimension']
    ##affineTransformT2 = aims.AffineTransformation3d(transformT2)
    ##t1_to_t2 = affineTransformT2.inverse() * affineTransformT1
    ##resamplerT1 = aims.ResamplerFactory_S16().getResampler(0)
    ##resamplerT1.setRef(volClassif)
    ##volumeT1_resamp = resamplerT1.doit(t1_to_t2, dimsT2[0], dimsT2[1], dimsT2[2], finder2.header()['voxel_size'])
    ##aims.write(volumeT1_resamp, data_directory + 'GWsegm_%s_manuallyFromT1.nii.gz' %(keyWord)) 
    
    ######################## just for test: how to improve skeletons in T2 ##############################################
    #if eliminateSulci is True:
        #print '1'
        #pathToSulciFile = brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/default_analysis/folds/3.1/default_session_auto/segmentation/%sSulci_%s_default_session_auto.nii.gz' %(realSide, realPatientID)
        #volSulci = aims.read(pathToSulciFile)
        #volSulciT2 = highres_cortex.od_cutOutRois.transformResampleVolFromVol(volSulci, volT2, volT1) 
        #aims.write(volSulciT2, data_directory + 'sulciSkel_%s_T2.nii.gz' %(keyWord)) 
        #print '2'

        ## read in sulci skeletons volume. dilate it
        #volSulciDil = aims.read(pathToSulciFile, 1)
        #volSulciDil = aimsalgo.AimsMorphoDilation(volSulciDil, 1)
        #print '3'
        #aims.write(volSulciDil, data_directory + 'sulciSkel_%s_dilated.nii.gz' %(keyWord)) 
        ## transform into T2 space. Better? Less holes?
        #volSulciDilT2 = highres_cortex.od_cutOutRois.transformResampleVolFromVol(volSulciDil, volT2, volT1) 
        #print '4'
        #aims.write(volSulciDilT2, data_directory + 'sulciSkel_%s_dilatedT2.nii.gz' %(keyWord)) 
        ## erode this volume. come to original one
        #volSulciDilT2 = aims.read(data_directory + 'sulciSkel_%s_dilatedT2.nii.gz' %(keyWord), 1)
        #volSulciDilT2Eroded = aimsalgo.AimsMorphoErosion(volSulciDilT2, 1)
        #print '5'
        #aims.write(volSulciDilT2Eroded, data_directory + 'sulciSkel_%s_dilatedT2eroded.nii.gz' %(keyWord))      

        ################### try another approach: tansform the vol1 into T2. and create skeletons there ##########################
        #fileCortex = brainvisa_db_neurospin + '%s/t1mri/reversed_t1map_2/default_analysis/segmentation/%scortex_%s.nii.gz' %(realPatientID, realSide, realPatientID)        
        #volCortex = aims.read(fileCortex)
        #volCortexT2 = highres_cortex.od_cutOutRois.transformResampleVolFromVol(volCortex, volT2, volT1)
        #aims.write(volCortexT2, data_directory + 'volCortexT2_%s.nii.gz' %(keyWord)) 

        ## read in T1 - with corrected grey level.
        #volT1nobias = aims.read(brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/default_analysis/nobias_%s.nii.gz' %(realPatientID))
        #volT1nobiasT2 = highres_cortex.od_cutOutRois.transformResampleVolFromVol(volT1nobias, volT2, volT1)
        #aims.write(volT1nobiasT2, data_directory + 'volT1nobiasT2_%s.nii.gz' %(keyWord)) 
        
        ## read in the initial grey-white classification and transform it into T2
        #volT1_greyWhite = aims.read(pathToClassifFileOriginal)
        #volT1_greyWhiteT2 = highres_cortex.od_cutOutRois.transformResampleVolFromVol(volT1_greyWhite, volT2, volT1)
        #aims.write(volT1_greyWhiteT2, data_directory + 'volT1_GW_T2_%s.nii.gz' %(keyWord)) 

        ## mask the volT1_nobias_inT2 by the grey_white region
        #volT1nobiasT2 = aims.read(data_directory + 'volT1nobiasT2_%s.nii.gz' %(keyWord))
        ##volCortexT2 = aims.read(data_directory + 'volCortexT2_%s.nii.gz' %(keyWord))
        ##arrCortexT2 = np.array(volCortexT2, copy = False)
        #volT1_greyWhiteT2 = aims.read(data_directory + 'volT1_GW_T2_%s.nii.gz' %(keyWord)) 
        #arrT1nobiasT2 = np.array(volT1nobiasT2, copy = False)
        #arr_greyWhiteT2 = np.array(volT1_greyWhiteT2, copy = False)
        #arrT1nobiasT2[arr_greyWhiteT2 == 0] = 0
        #aims.write(volT1nobiasT2, data_directory + 'volT1nobiasT2_maskHemi_%s.nii.gz' %(keyWord))         
        ##volT1inT2 = highres_cortex.od_cutOutRois.transformResampleVolFromVol(volT1, volT2)
        ##aims.write(volT1inT2, data_directory + 'volT1inT2_%s.nii.gz' %(keyWord)) 
    
       ## subprocess.check_call(["VipSkeleton", "-i", data_directory + 'volCortexT2_%s.nii.gz' %(keyWord), "-so", data_directory + 'skeleton_%s.nii.gz' %(keyWord), "-vo", data_directory + 'roots_%s.nii.gz' %(keyWord), "-g", data_directory + 'volT1nobiasT2_maskHemi_%s.nii.gz' %(keyWord), "-w", "t"])     
        
##subprocess.check_call(["VipSkeleton", "-i", self.hemi_cortex, "-so", self.skeleton, "-vo", self.roots, "-g", data_directory + 'volT1inT2_%s.nii.gz' %(keyWord), "-w", "t"])     

print 'after all: keyWord', keyWord, ' and pathToClassifFile ', pathToClassifFile



## launch the distance map calculation # classif file must be here with 100, 200, 0 (no 50 and 150)
t0dist = timeit.default_timer()
subprocess.check_call(['time', 'od_distmapsMain.py', 
'-i', pathToClassifFile, '-d', data_directory, '-k', keyWord])
t1dist = timeit.default_timer()

### launch the heat map calculation # classif file must be here with 100, 200, 0 (no 50 and 150)
t0heat = timeit.default_timer()

# launch this process for images in initial T1 space:
if options.workOnT1inT2Space is not None:    
    print 'start heat calculation for images in T2 space resampled from T1'
    subprocess.check_call(['time', 'od_heatMain.py', '-i', pathToClassifFile, '-r', '-d', data_directory, '-k', keyWord])
else:
    print 'start heat calculation for images in initial T1 space'
    subprocess.check_call(['time', 'od_heatMain.py', '-i', pathToClassifFile, '-d', data_directory, '-k', keyWord])

t1heat = timeit.default_timer()

## launch the isovolumes calculation # classif file must be here with 100, 200, 0  (no 50 and 150)
t0iso = timeit.default_timer()
subprocess.check_call(['time', 'od_isovolumeMain.py', 
'-i', pathToClassifFile, '-d', data_directory, '-k', keyWord])
t1iso = timeit.default_timer()

## launch the calculation of the 'cortical columns' # classif file must be here with 100, 200, 0 , 50 and 150
t0col = timeit.default_timer()
pathToClassifWithBorders = data_directory + 'dist/classif_with_outer_boundaries_%s.nii.gz' %(keyWord)
subprocess.check_call(['time', 'od_column-regionsMain.py', 
'-i', pathToClassifWithBorders, '-d', data_directory, '-k', keyWord])
t1col = timeit.default_timer()

tDist = t1dist - t0dist
tHeat = t1heat - t0heat
tIso = t1iso - t0iso
tCol = t1col - t0col
header = ['sizeX', 'sizeY', 'sizeZ', 'voxelX', 'voxelY', 'voxelZ', 'tDist', 'tHeat', 'tIso', 'tCol', 'tTotal', 'keyWord', 'classifFile']
volClassif = aims.read(pathToClassifFile)
headerClassif = volClassif.header()
content = [str(headerClassif['sizeX']), str(headerClassif['sizeY']), str(headerClassif['sizeZ']), str(headerClassif['voxel_size'][0]), str(headerClassif['voxel_size'][1]), str(headerClassif['voxel_size'][2]), str(np.round(tDist)), str(np.round(tHeat)), str(np.round(tIso)), str(np.round(tCol)), str(np.round(tDist + tHeat + tIso + tCol)), keyWord, pathToClassifFile]
headerLine = ('\t').join(header)
contentLine = ('\t').join(content)
statFileName = data_directory + "statFile_%s" %(keyWord)

# note if the processing was performed on laptop:
if workOnLaptop:
    statFileName += '_laptop.txt'
else:
    statFileName += '.txt'

f = open(statFileName, "w")
f.write(headerLine + '\n')
f.write(contentLine + '\n')
f.close()
    
    
    
    
    
