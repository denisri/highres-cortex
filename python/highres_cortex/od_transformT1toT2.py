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
# this function will take the original T1 volume, the T2 volume, and 
# transform T1 into T2 space, respecting the bounding box
# update: it also takes the processing results of Morphologist (in T1 space) 
# and transforms them into T2 space (needed to compare and validate results
# of Morphologist on original images in the new T2 space)

# example how to run this file:
# python /volatile/od243208/brainvisa_sources/highres-cortex/python/highres_cortex/od_transformT1toT2.py -d /volatile/od243208/brainvisa_manual/testT1toT2_LinearWithHemi/ -r /neurospin/lnao/dysbrain/raw_niftis/

import random
from soma import aims, aimsalgo
import subprocess
from optparse import OptionParser
from scipy.stats import mode
import sys, glob, os, os.path, subprocess, sys, time, timeit
import numpy as np
import highres_cortex.od_cutOutRois
from soma.aims import volumetools


brainvisa_db_neurospin = '/neurospin/lnao/dysbrain/brainvisa_db_morphologist/dysbrain/'
brainvisa_raw_niftis = '/neurospin/lnao/dysbrain/raw_niftis/'
pathToTextures = '/neurospin/lnao/dysbrain/randomized_flipped_data/manual_work/'
# T2 images after bias correction with parameters (field rigidity = 5, sampling = 4) optimized on 16.02.2015
#pathToNobiasT2 = '/neurospin/lnao/dysbrain/nobiasT2_FR5_S4/'
pathToNobiasT2 = '/neurospin/lnao/dysbrain/nobiasT2_FR5_S16/'

#volT1 = aims.read('/volatile/od243208/brainvisa_db_morphologist/dysbrain/ac140155/t1mri/reversed_t1map_2/ac140155.nii.gz')
#volT2 = aims.read('/volatile/od243208/raw_niftis/ac140155/20140703_152737t2spctraiso05mmVOIs011a1001.nii.gz')
def transformT1toT2(volT1, volT2, transformOther):
    import numpy as np
    # get the bounding box of T1 (in pixels)
    size_x = volT1.getSizeX()
    size_y = volT1.getSizeY()
    size_z = volT1.getSizeZ()
    print 'size x, y, z of the volT1 is ', size_x, size_y, size_z, ' pixels'
    
    # + or - 0.5 voxel size to their coordinates
    voxSize1 = volT1.getVoxelSize()
    print 'voxSize1 : ', voxSize1[0], voxSize1[1], voxSize1[2]
    
    # we add a whole voxel to compensate both directions
    size_x = (size_x + 1) * voxSize1[0]         # (in mm)
    size_y = (size_y + 1) * voxSize1[1]
    size_z = (size_z + 1) * voxSize1[2]
    print 'size x, y, z of the volT1 is ', size_x, size_y, size_z, ' mm'

    # get coordinates of 8 points of the (extended) bounding box
    xCoords = (0, size_x, size_x, 0, 0, size_x, size_x, 0)
    yCoords = (0, 0, size_y, size_y, 0, 0, size_y, size_y)
    zCoords = (0, 0, 0, 0, size_z, size_z, size_z, size_z)

    # get transformation T1 to T2
    header1 = volT1.header()
    header2 = volT2.header()
    transform1 = header1['transformations'][0]
    affineTransform1 = aims.AffineTransformation3d(transform1)
    transform2 = header2['transformations'][0]
    dims2 = volT2.getSize()
    voxSize2 = header2['voxel_size']
    affineTransform2 = aims.AffineTransformation3d(transform2)
    t1_to_t2_original = affineTransform2.inverse() * affineTransform1
    t1_to_t2 = affineTransform2.inverse() * affineTransform1

    # transform all 8 points into T2 space
    newPoints = []    
    xCoords_new = []
    yCoords_new = []
    zCoords_new = []
    for i in range(0,8):
        p = aims.Point3df(xCoords[i], yCoords[i], zCoords[i])
        p1 = t1_to_t2.transform(p)
        newPoints.append(p1)
        xCoords_new.append(p1[0])
        yCoords_new.append(p1[1])
        zCoords_new.append(p1[2])
        #print 'old point coords ', xCoords[i], yCoords[i], zCoords[i]
        #print 'new point coords ', p1[0], p1[1], p1[2]
        #print '-------------------------------------------'

    # check if some of the new coordinates are negative, then translate the volume to make them all positive
    shiftX = 0
    shiftY = 0
    shiftZ = 0

    if min(xCoords_new) < 0 :
        shiftX = min(xCoords_new)
        #xCoords_new = [i - m for i in xCoords_new]    
        
    if min(yCoords_new) < 0 :
        shiftY = min(yCoords_new)
        #yCoords_new = [i - m for i in yCoords_new]
    
    if min(zCoords_new) < 0 :
        shiftZ = min(zCoords_new)
        #zCoords_new = [i - m for i in zCoords_new]
    
    print 'correction X, Y, Z: ', shiftX, shiftY, shiftZ
    print 'transformation t1_to_t2 : '
    print t1_to_t2
    
    # modify the transformation to compensate for these negative shifts
    t1_to_t2.translation()[0] -= shiftX
    t1_to_t2.translation()[1] -= shiftY
    t1_to_t2.translation()[2] -= shiftZ
    
    print 'modified transformation t1_to_t2 translation : '
    print t1_to_t2
    
    # todo! need to increase the size of the resulting volume!
    
    #for i in range(0,8):
        #print 'old point coords ', xCoords[i], yCoords[i], zCoords[i]
        #print 'new point coords ', xCoords_new[i], yCoords_new[i], zCoords_new[i]
        #print '-------------------------------------------'  
    
    # now we have coordinates of 8 points in T2 space. Create a "Volume" using them
    newMaxX = max(xCoords_new) - shiftX
    newMaxY = max(yCoords_new) - shiftY
    newMaxZ = max(zCoords_new) - shiftZ
    
    # check that these sizes are at least as large as the initial T2 images
    #print '    newMaxX = np.max(', newMaxX, ' or ', dims2[0], ' * ', voxSize2[0]
    #print newMaxX
    #print dims2[0]
    #print voxSize2[0]
    #z = dims2[0] * voxSize2[0]
    #print z
    #y = np.maximum(newMaxX, z)
    #print y
    #print (y - 1)
    newMaxX = np.maximum(newMaxX, dims2[0] * voxSize2[0])
    newMaxY = np.maximum(newMaxY, dims2[1] * voxSize2[1])
    newMaxZ = np.maximum(newMaxZ, dims2[2] * voxSize2[2])

    print 'volT1 max size in mm after transformation is ', newMaxX, newMaxY, newMaxZ
    #volNew = aims.Volume( newMaxX, newMaxY, newMaxZ, dtype='int32' )
    
    # transform the original T1 volumt into T2 space
    print 'voxSize2 : ', voxSize2[0], voxSize2[1], voxSize2[2]
    # modified! Resample volT1 with LINEAR interpolation!!! (1). and not with NN (0).
    resampler1 = aims.ResamplerFactory_S16().getResampler(1)
    resampler1.setRef(volT1)
    print '------------- RESAMPLER:  ---------  volT1.header()  ----------------------'
    #print volT1.header()
    vol_resamp = resampler1.doit(
        t1_to_t2, int(newMaxX/voxSize2[0]), int(newMaxY/voxSize2[1]),
        int(newMaxZ/voxSize2[2]), aims.Point3df(*voxSize2[:3])).volume()
    print 'given dims to resampler : ', newMaxX/voxSize2[0], newMaxY/voxSize2[1], newMaxZ/voxSize2[2]
    print 'voxSize2 : ', voxSize2[:]
    print 'type voxSize2 : ', type(voxSize2)
    print 'type: newMaxX/voxSize2[0] ', type(newMaxX/voxSize2[0])

    
    # now need to apply the same translation to the T2 volume
    t2Translation = aims.AffineTransformation3d()
    #t2Translation.translation()[:] = t1_to_t2.translation()[:]    
    #t2Translation.translation()[:] = [- t1_to_t2.translation()[0], - t1_to_t2.translation()[1], - t1_to_t2.translation()[2]]   
    t2Translation.translation()[:] = [- shiftX, - shiftY, - shiftZ]   
    
    # modified! Resample volT2 with LINEAR interpolation!!! (1). and not with NN (0).
    resampler2 = aims.ResamplerFactory_S16().getResampler(1)
    resampler2.setRef(volT2)
    print 'transformation to be applied to the T2 volume: '
    print t2Translation
    vol_resamp2 = resampler2.doit(
        t2Translation, int(newMaxX/voxSize2[0]), int(newMaxY/voxSize2[1]),
        int(newMaxZ/voxSize2[2]), aims.Point3df(*voxSize2[:3])).volume()
    
    # save the results
    resList = [vol_resamp, vol_resamp2, t1_to_t2_original, t1_to_t2, t2Translation]

    ##################################### transform skeletons, texture, GW, hemis if required #########################################
    vols = []
    names = []
 
    if transformOther:
        print 'start transforming other volumes : ' , str(transformOther)
         
        #for realSide in ['L', 'R']:            
            ## 1. GW classif
            #volsGWlist = glob.glob(brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/default_analysis/segmentation/%sgrey_white_%s.nii.gz' %(realSide, realPatientID))   # ok only if there is one file
            #if len(volsGWlist) == 1:
                #volGW = aims.read(volsGWlist[0])
                #print 'Took the volGW: ', volsGWlist[0]
                ## GW volume is resampled as the NN (0)
                #resampler3 = aims.ResamplerFactory_S16().getResampler(0)
                #resampler3.setRef(volGW)
                #vol_resamp3 = resampler3.doit(t1_to_t2, newMaxX/voxSize2[0], newMaxY/voxSize2[1], newMaxZ/voxSize2[2], voxSize2)
                #vols.append(vol_resamp3)
                #print 'resampled the initial GW %s ' %(realSide)
                #names.append('_GW_%s_T1inNewT2.nii.gz' %(realSide))
            
                
            ## 2. skeletons
            #sulcilist = glob.glob(brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/default_analysis/folds/3.1/default_session_auto/segmentation/%sSulci_%s_default_session_auto.nii.gz' %(realSide, realPatientID))   # ok only if there is one file
            #if len(sulcilist) == 1: 
                #sulci = aims.read(sulcilist[0])
                #print 'Took the sulci: ', sulcilist[0]         
                ## Skeletons are resampled as NN, too (0)
                #resampler4 = aims.ResamplerFactory_S16().getResampler(0)
                #resampler4.setRef(sulci)
                #vol_resamp4 = resampler4.doit(t1_to_t2, newMaxX/voxSize2[0], newMaxY/voxSize2[1], newMaxZ/voxSize2[2], voxSize2)    
                #vols.append(vol_resamp4)    
                #print 'resampled the initial T1 %s skeleton ' %(realSide)
                #names.append('_sulciSkel_%s_T1inNewT2.nii.gz' %(realSide))
                
            ## 3. textures        
            ## to transform texture: need to transform the corresponding mesh
            #fileHemi = glob.glob(brainvisa_db_neurospin + '%s/t1mri/reversed_t1map_2/default_analysis/segmentation/mesh/%s_%shemi.gii' %(realPatientID, realPatientID, realSide))
            #if len(fileHemi) == 1:
                #print 'found the hemisphere file : ', fileHemi 
                ## resample the hemisphere
                #volHemi = aims.read(fileHemi[0])
                #aims.SurfaceManip.meshTransform(volHemi, t1_to_t2)
                #print 'resampled the initial T1 %s hemisphere ' %(realSide)
                #vols.append(volHemi)    
                #names.append('_hemi_%s_T1inNewT2.gii' %(realSide))     
                
        ## 4. B1 map        
        #fileB1 = glob.glob(brainvisa_raw_niftis + '%s/*b1map*.nii.gz' %(realPatientID))           
        #if len(fileB1) == 1:
            #print 'found the hemisphere file : ', fileB1[0] 
            ## resample the hemisphere
            #volB1 = aims.read(fileB1[0])
            #resampler5 = aims.ResamplerFactory_S16().getResampler(1)
            #resampler5.setRef(volB1)
            #vol_resamp5 = resampler5.doit(t1_to_t2, newMaxX/voxSize2[0], newMaxY/voxSize2[1], newMaxZ/voxSize2[2], voxSize2)    
            #vols.append(vol_resamp5)    
            #print 'resampled the B1map '
            #names.append('_B1map_inNewT2.nii.gz')

 
         # 5. GRE...       
        fileGre = glob.glob(brainvisa_raw_niftis + '%s/*dpgreiso*.nii.gz' %(realPatientID))           
        if len(fileGre) == 1:
            print 'found the fileGre file : ', fileGre[0] 
            # resample the hemisphere
            volGre = aims.read(fileGre[0])
            resampler6 = aims.ResamplerFactory_S16().getResampler(1)
            resampler6.setRef(volGre)
            vol_resamp6 = resampler6.doit(t1_to_t2, newMaxX/voxSize2[0], newMaxY/voxSize2[1], newMaxZ/voxSize2[2], voxSize2)    
            vols.append(vol_resamp6)    
            print 'resampled the Gre '
            names.append('_Gre_inNewT2.nii.gz')

            
        ############################## problem: doing this outside of this function leads to errors!!! #########################
        #realSide = 'R'
        #pathToSulciFile = brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/default_analysis/folds/3.1/default_session_auto/segmentation/%sSulci_%s_default_session_auto.nii.gz' %(realSide, realPatientID)
        #print 'found the sulci skeletons file : ', pathToSulciFile
        #volSulci = aims.read(pathToSulciFile)
        
        #resampler3 = aims.ResamplerFactory_S16().getResampler(0)
        #resampler3.setRef(volSulci)
        #vol_resamp3 = resampler3.doit(t1_to_t2, newMaxX/voxSize2[0], newMaxY/voxSize2[1], newMaxZ/voxSize2[2], voxSize2)
        #print 'resampled the initial T1 %s skeleton ' %(realSide)
        #aims.write(vol_resamp3, resultDirectory + keyWord + '_%s_SulciSkelT1_inT2.nii.gz' %(realSide))
        ################################################################################################################
        
        resList.append(vols)
        resList.append(names)    
        
    return(resList)         

    
if __name__ == '__main__':
    
    volT1 = None
    volT2 = None
    outputFile = None
    realPatientID = None
    keyWord = None
    directory = None
    recursiveInDirectory = None
    resultDirectory = None
    transformOtherVolumes = False       # False, so transform and cut only the T2 volume. Else, work with skeletons, textures and hemispheres, too


    parser = OptionParser('Transform T1 volume into T2 space without cutting it')
    parser.add_option('-i', dest='volT1', help='volT1. Optional')   
    parser.add_option('-t', dest='volT2', help='volT2. Optional') 
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    parser.add_option('-k', dest='keyWord', help='keyWord')
    parser.add_option('-d', dest='directory', help='directory')
    parser.add_option('-r', dest='recursiveInDirectory', help='recursiveInDirectory')
    parser.add_option('-o', dest='transformOtherVolumes', action = 'store_true', help='Select if need to transform other volumes: skeletons, textures and hemispheres, too. False is default') 
    options, args = parser.parse_args(sys.argv)
    print options
    print args   
    
    subjectList=[]
    
    if options.directory is None:
        print >> sys.stderr, 'New: exit. no directory given'
        sys.exit(1)
    else:
        directory = options.directory           

    if options.realPatientID is not None:
        realPatientID = options.realPatientID     
        subjectList.append(realPatientID)
        
        if options.volT1 is None:
            # volume T1 can be optionally given. If nothing is give: take the standard volume
            volT1 = brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/%s.nii.gz' %(realPatientID)  
            print 'Took the standard T1 file: ', volT1
            #print >> sys.stderr, 'New: exit. no volT1 given'
            #sys.exit(1)
        else:
            volT1 = options.volT1
                
        if options.volT2 is None:
            # volume T21 can be optionally given. If nothing is give: take the standard volume
            finder2 = aims.Finder()
            
            #volT2list = glob.glob(brainvisa_raw_niftis + realPatientID + '/*t2*.nii.gz')   # ok only if there is one t2-file, or the first one of the t2-files is correct
          
            # work with T2 image after the bias correction!!
            volT2list = glob.glob(pathToNobiasT2 + 'nobias_%s.nii.gz' %(realPatientID))   # ok only if there is one t2-file, or the first one of the t2-files is correct
            
            if len(volT2list) == 1: #finder2.check(volT2)
                volT2 = volT2list[0]         
                print 'Took the standard T2 file: ', volT2
            else:
                print >> sys.stderr, 'found ', len(volT2list), ' suitable T2 files. Exit'
                sys.exit(1)                                
        else:
            volT2 = options.volT2     
        
        if options.keyWord is None:
            keyWord = realPatientID
        else:
            keyWord = options.keyWord             
    else:
        # check, maybe recursive processing is required
        if options.recursiveInDirectory is not None:
            # process all subjects in a given directory
            recursiveInDirectory = options.recursiveInDirectory
            subjectList = [o for o in os.listdir(recursiveInDirectory) if os.path.isdir(os.path.join(recursiveInDirectory,o))]            
        else:
            # exit. nothing was given        
            print >> sys.stderr, 'New: exit. neither realPatientID nor the folder for the processing was given'
            sys.exit(1)
    
    if options.transformOtherVolumes is not None:
        transformOtherVolumes = options.transformOtherVolumes
    
    # now process all the required subjects
    for realPatientID in subjectList:    
        print '####################### process subject  ', realPatientID, '  ##############################################'
        # check if the destination folder exists. If not - create it
        resultDirectory = directory + realPatientID + '/'

        if len(subjectList) > 1 or volT1 is None:
            volT1 = brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/%s.nii.gz' %(realPatientID)  
            print 'Took the standard T1 file: ', volT1
        
        if len(subjectList) > 1 or volT2 is None:
            finder2 = aims.Finder()
            #volT2list = glob.glob(brainvisa_raw_niftis + realPatientID + '/*t2*.nii.gz')   # ok only if there is one t2-file, or the first one of the t2-files is correct
    
            # work with T2 image after the bias correction!!
            volT2list = glob.glob(pathToNobiasT2 + 'nobias_%s.nii.gz' %(realPatientID))   # ok only if there is one t2-file, or the first one of the t2-files is correct
 
            if len(volT2list) == 1: # finder2.check(volT2):
                volT2 = volT2list[0]
                print 'Took the standard T2 file: ', volT2
            else:
                print 'found ', len(volT2list), ' of volT2 files. Continue with the next subject'
                continue
 
        if not os.path.exists(resultDirectory):
            os.makedirs(resultDirectory)

        keyWord = realPatientID
            
        vol1 = aims.read(volT1)
        vol2 = aims.read(volT2)    

        vols = transformT1toT2(vol1, vol2, transformOtherVolumes)
        vol1t = vols[0]
        pathToVol1 = resultDirectory + keyWord + '_T1inNewT2.nii.gz'
        aims.write(vol1t, pathToVol1)
        vol2t = vols[1]
        pathToVol2 = resultDirectory + keyWord + '_NewT2.nii.gz'
        aims.write(vol2t, pathToVol2)
        
        # save transformations
        t1_to_t2_original = vols[2]
        t1_to_t2 = vols[3]
        t2Translation = vols[4]
        pathTot1_to_t2_original = resultDirectory + keyWord + '_t1_to_t2_original.trm'
        aims.write(t1_to_t2_original, pathTot1_to_t2_original)
        
        pathTot1_to_t2 = resultDirectory + keyWord + '_t1_to_t2.trm'
        aims.write(t1_to_t2, pathTot1_to_t2)
        
        pathTot2 = resultDirectory + keyWord + '_t2.trm'
        aims.write(t2Translation, pathTot2)
        
        
        # Crop the transformed volumes using Denis's function
        volsToCrop = []         # read in the volumes that will be cropped in the next step
        volsToCrop.append(aims.read(pathToVol2))
        namesNewCropp = []
        namesNewCropp.append('_NewT2_cropped.nii.gz')

        # if other volumes were transformed : crop them, too
        print 'other volumes were transformed : crop them, too'
        if transformOtherVolumes:
            volsNew = vols[5]
            namesNew = vols[6]            
          
            print 'create new CroppedNames:'
            for i in range(len(volsNew)):
                # write out the volume transformed into ne new T2 space
                aims.write(volsNew[i], resultDirectory + keyWord + namesNew[i])
                # read it in for the cropping
                volsToCrop.append(aims.read(resultDirectory + keyWord + namesNew[i]))
                # get the name for the new cropped file
                nameCropped = namesNew[i].replace('_T1inNewT2', '_T1inNewT2_cropped')
                print 'initial file : ', resultDirectory + keyWord + namesNew[i], ' newCroppedName: ', nameCropped
                namesNewCropp.append(nameCropped)        
            
        # crop these volumes using Denis's function and save them
        print 'start cropping the volumes'
        vol1tVol = aims.read(pathToVol1)
        cropped = volumetools.crop_volumes(vol1tVol, volsToCrop, threshold=80, border=10)
        aims.write(cropped[0], resultDirectory + keyWord + '_T1inNewT2_cropped.nii.gz')
        
        for i in range(1, len(cropped)):
            print i, ' write file: ', resultDirectory + keyWord + namesNewCropp[i-1]
            aims.write(cropped[i], resultDirectory + keyWord + namesNewCropp[i-1])

        print 'done cropping the volumes'

 
    
    
    
    
    
    
    