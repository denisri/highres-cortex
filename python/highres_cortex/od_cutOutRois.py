#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright CEA (2014).
# Copyright Université Paris XI (2014).
#
# Contributor: Olga Domanova <olga.domanova@cea.fr>.

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

import numpy as np
from soma import aims, aimsalgo
import subprocess

def voronoiFromTexture(volGW_border, tex, hemi, stopLabel, directory, keyWord):
    """ This function takes a volume from which ROI will be cut out
    This volume should be read in with the border = 1 (for later dilation)
    It also takes a texture labelled (e.g. manually) from which seeds will be taken for the
    Voronoi classification.
   # DomainLabel - where to propagate the classification (is taken from the dilated volume)
    StopLabel - where to stop it
    directory - is a path where some intermediate file will be written to. Names will include the keyWord
    This function returns a volume with the selected region
    """
    vs = volGW_border.getVoxelSize()[:3]
    dims = volGW_border.getSize()
    c = aims.Converter(intype=volGW_border, outtype=aims.Volume('S16'))
    volShort = c(volGW_border)
    volGW_dilated = aimsalgo.AimsMorphoDilation(volShort, 1)

    # Had problem with some textures: instead of value '10', it was '9.9999999999'. -> use round
    for i, vertex in enumerate(hemi.vertex()):
        label = round(tex[0].item(i)) + 1
        posVox = (int(round(vertex[0] / vs[0])), int(round(vertex[1] / vs[1])), int(round(vertex[2] / vs[2])))
        # print posVox
        if posVox[0] >= 0 and posVox[0] < dims[0] and posVox[1] >= 0 and posVox[1] < dims[1] and posVox[2] >= 0 and posVox[2] < dims[2]:
            volGW_dilated.setValue(label, posVox[0], posVox[1], posVox[2])
                
    arrDilated = np.asarray(volGW_dilated.volume())  
    arrDilated[np.asarray(volGW_border) == 200] = 0    
    
    aims.write(volGW_dilated, directory + 'seedsNoWM_%s.nii.gz' %(keyWord))
    # Voronoi classification
    subprocess.call(['AimsVoronoi', '-i', directory + 'seedsNoWM_%s.nii.gz' %(keyWord), '-o', directory +  'voronoi_%s.nii.gz' %(keyWord), '-d', '32767', '-f', str(stopLabel)]) 
    volVoronoi = aims.read(directory +  'voronoi_%s.nii.gz' %(keyWord))
    
    # label '0' in texture was transformed into '1'. So we will set all voxels with value '1' to '0'
    arrVor = np.array(volVoronoi, copy = False)
    arrVor[arrVor == (stopLabel + 1)] = 0
    aims.write(volVoronoi, directory +  'voronoi_%s.nii.gz' %(keyWord))
    return(volVoronoi)

    


def excludeROI(volInitial, volLabel, fromLabel, excludeLabel):
    """this function that takes 2 volumes 
    from voxels of the first volume those are deleted that have label 'fromLabel' while in the 2nd volume they are labeled as 'excludeLabel'
    """
    arrInitial = np.array(volInitial, copy = False)
    arrLabel = np.array(volLabel, copy = False)
    arrInitial[(arrInitial == fromLabel) & (arrLabel == excludeLabel)] = 0
    return(volInitial)
    

def correctVoxelLabels(vol, pathToVol, directory, keyWord, minSize, connectivity):
    """ This function takes a volume, finds all unique labels, and a minimum size of a connected component required
    It takes all the too small connected components, and re-classifies their voxels depending on the neighborhood
    ignoreLabel - is the label for background (do not need to re-classify these pixels
    
    take connectivity 6 : it will better clean the classification than 26.
    for correction of classification: use the connectivity 6 
    """
    arr = np.asarray(vol)
    uniqueLabels = np.unique(arr)
    totalN = 0
    volCheckLabelMode = aims.read(pathToVol)
    
    for i in uniqueLabels:
        print '----------------------------------work with unique label %s. Re-classify small CCs. -----------------' %(i)
        
        # threshold the volume at the current label
        subprocess.call(['AimsThreshold', '-i', pathToVol, '-o', directory + 'voronoi_Thr%s_%s.nii.gz' %(str(i), keyWord), '-m', 'eq', '-t', str(i), '-b', 'true', '--fg', str(i)])

        # find connected components
        subprocess.call(['AimsConnectComp', '-i', directory + 'voronoi_Thr%s_%s.nii.gz' %(str(i), keyWord), '-o', directory + 'voronoi_Thr%s_CC%s_%s.nii.gz' %(str(i), str(connectivity), keyWord), '-c', str(connectivity)])

        # find big connected components
        subprocess.call(['AimsConnectComp', '-i', directory + 'voronoi_Thr%s_%s.nii.gz' %(str(i), keyWord), '-o', directory + 'voronoi_Thr%s_CC%sBigger%s_%s.nii.gz' %(str(i), str(connectivity), str(minSize), keyWord),  '-c', str(connectivity), '-s', str(minSize), '--verbose'])

        volAllCC = aims.read(directory + 'voronoi_Thr%s_CC%s_%s.nii.gz' %(str(i), str(connectivity), keyWord))
        arrAllCC = np.asarray(volAllCC)
        arrAllCC[arrAllCC != 0] = 1     # set all included voxels to "1"
        volBigCC = aims.read(directory + 'voronoi_Thr%s_CC%sBigger%s_%s.nii.gz' %(str(i), str(connectivity), str(minSize), keyWord))
        arrBigCC = np.asarray(volBigCC)
        arrBigCC[arrBigCC != 0] = 1     # set all included voxels to "1"
        verticesDel = np.where(arrAllCC != arrBigCC)   # vertices that are labeled initially but those that are in too small CCs
        numberChanged = 0

        print 'for label ', i, '  the number of voxels to check is  ', len(verticesDel[0])
        for j in range(len(verticesDel[0])):
            x = verticesDel[0][j]
            y = verticesDel[1][j]
            z = verticesDel[2][j]
            
            listOfNeighbors = []
            listOfNeighbors.append(vol.value(x, y, z + 1))
            listOfNeighbors.append(vol.value(x, y, z - 1))
            listOfNeighbors.append(vol.value(x + 1, y, z))
            listOfNeighbors.append(vol.value(x - 1, y, z))
            listOfNeighbors.append(vol.value(x, y + 1, z))
            listOfNeighbors.append(vol.value(x, y - 1, z))
            oldL = vol.value(x, y, z)
            newL = np.argmax(np.bincount(listOfNeighbors))
            
            # TODO later?: check how many neighbours have the most frequent label. If there is the equal number of labels e.g. '0' and '11', then '0' is assigned. In this case maybe prefer to set '11'?
            #if len(mode(listOfNeighbors)) > 1:
                #lab1 = int(mode(listOfNeighbors)[0][0])
                #lab2 = int(mode(listOfNeighbors)[1][0])
                #labNew = int(str(lab2) + str(lab1))
                #print 'several labels with equal frequency: ', lab1, ' and ', lab2
                #volCheckLabelMode.setValue
            
            if newL != oldL:
                print 'component number ', j, x, y, z, listOfNeighbors,  ' old label: ', oldL, ' newLabel: ', newL
                # substitute the value of this voxel with the mode value of the neighboring voxels
                vol.setValue(int(newL), x, y, z)
                numberChanged += 1
        print 'numberChanged ', numberChanged, 'from : ', len(verticesDel[0])
        totalN += numberChanged
    
    print 'totalN changed: ', totalN
    return(vol)



def labelDirect26Neighbors(vol, volDomain, label, labelWhereToExtend, domain, newLabel):
    """ this function takes a volume and a label to look for
    direct 26 neighbors of labeled voxels are found and are labeled with a newlabel
    it is performed in a given domain of the volDomain
    """            
    arr = np.array(vol, copy = False)   
    dims = vol.getSize()
    #print 'start labelDirect26Neighbors for label ', label
    verticesL = np.where(arr == label)
    #print 'len(verticesL[0]) ', len(verticesL[0])
    numberChanged = 0
    
    # iterate in vertices and look at their neighbors:
    for i in range(len(verticesL[0])):
        x = verticesL[0][i]
        y = verticesL[1][i]
        z = verticesL[2][i]
        for x1 in [x-1, x, x+1]:
            for y1 in [y-1, y, y+1]:    
                for z1 in [z-1, z, z+1]:
                    if x1 >= 0 and x1 < dims[0] and y1 >= 0 and y1 < dims[1] and z1 >= 0 and z1 < dims[2]:
                        if (vol.value(x1, y1, z1) == labelWhereToExtend) and (volDomain.value(x1, y1, z1) == domain):
                            #print 'this is a voxel to work with ?x1, y1, z1 ', x1, y1, z1, ' with value ', vol.value(x1, y1, z1), ' with GW value ', volDomain.value(x1, y1, z1)
                            vol.setValue(newLabel, x1, y1, z1)
                            numberChanged += 1
                                
    print ' the number of vertices set to ', newLabel, '  is : ', numberChanged
    return(vol)
                

                            
def dilateRoiConnectivity26(vol, volDomain, labelsToInclude, labelWhereToExtend, domain, iterNum):
    """ this function takes a volume and a list of labels. corresponding labeled regions are dilated 
    regions are dilated according to neighbors' connectivity. dilation is repeated a given number of times
    it is performed in a given domain using a previous function (labelDirect26Neighbors)
    """            
    arr = np.array(vol, copy = False)    
    for l in labelsToInclude:
        print '###################start dilateRoiConnectivity26 for label ', l, '  ###################### '
        # first step: find direct neighbors of the given label and set them to label + 1
        newLabel = l + 1
        vol = labelDirect26Neighbors(vol, volDomain, l, labelWhereToExtend, domain, newLabel)
        
        # now repeat this step a required number of times, with the only difference that the label assigned is the same as the target one
        for k in range(iterNum - 1):
            print 'iteration #  ', (k + 1)
            vol = labelDirect26Neighbors(vol, volDomain, newLabel, labelWhereToExtend, domain, newLabel)    
    return(vol)
        
        

                        
                        
def labelDirect6Neighbors(vol, volDomain, label, labelWhereToExtend, domain, newLabel):
    """ this function takes a volume and a label to look for direct 6 neighbors of labeled voxels 
    neighbors are found and are labeled with a newlabel. it is performed in a domain of the volDomain
    """            
    arr = np.array(vol, copy = False)   
    dims = vol.getSize()
    #print 'start labelDirect6Neighbors for label ', label
    verticesL = np.where(arr == label)
    #print 'len(verticesL[0]) ', len(verticesL[0])
    numberChanged = 0
    
    # iterate in vertices and look at their neighbors:
    for i in range(len(verticesL[0])):
        x = verticesL[0][i]
        y = verticesL[1][i]
        z = verticesL[2][i]
        
        # 3D 6 -connectivity: 
        #x1 = x, x-1, x+1, x, x, x 
        #y1 = y, y, y, y-1, y+1, y 
        #z1 = z-1, z, z, z, z, z +1
        
        for x1, y1, z1 in zip([x, x-1, x+1, x, x, x], [y, y, y, y-1, y+1, y], [z-1, z, z, z, z, z +1]):
            if x1 >= 0 and x1 < dims[0] and y1 >= 0 and y1 < dims[1] and z1 >= 0 and z1 < dims[2]:                
                if (vol.value(x1, y1, z1) == labelWhereToExtend) and (volDomain.value(x1, y1, z1) == domain):
                    vol.setValue(newLabel, x1, y1, z1)
                    numberChanged += 1                                
    print ' the number of vertices set to ', newLabel, '  is : ', numberChanged
    return(vol)
            
        
def dilateRoiConnectivity6(vol, volDomain, labelsToInclude, labelWhereToExtend, domain, iterNum):
    """ this function takes a volume and a list of labels. corresponding labeled regions are dilated 
    regions are dilated according to neighbors' connectivity. dilation is repeated a given number of times
    it is performed in a given domain using a previous function (labelDirect26Neighbors)
    """            
    arr = np.array(vol, copy = False)    
    for l in labelsToInclude:
        print 'start dilateRoiConnectivity6 for label ', l
        # first step: find direct neighbors of the given label and set them to label + 1
        newLabel = l + 1
        vol = labelDirect6Neighbors(vol, volDomain, l, labelWhereToExtend, domain, newLabel)
        
        # now repeat this step a required number of times, with the only difference that the label assigned is the same as the target one
        for k in range(iterNum - 1):
            print 'iteration #  ', (k + 1)
            vol = labelDirect6Neighbors(vol, volDomain, newLabel, labelWhereToExtend, domain, newLabel)    
    return(vol)
              
        
def transformResampleVolFromVol(vol, vol2, vol1 = None):
    """ this function takes an original volume, takes space coordinates of the second volume, transforms and resamples the original volume
    into the space of the third volume.
    If the second volume is not given, the original volume is directly transformed into the 2nd volume's space
    """            
    if vol1 is None:
        header1 = vol.header()
    else:
        header1 = vol1.header()

    header2 = vol2.header()
    transform1 = header1['transformations'][0]
    affineTransform1 = aims.AffineTransformation3d(transform1)

    transform2 = header2['transformations'][0]
    dims2 = header2['volume_dimension']    
    voxSize2 = header2['voxel_size']
    affineTransform2 = aims.AffineTransformation3d(transform2)
    
    t1_to_t2 = affineTransform2.inverse() * affineTransform1
    resampler1 = aims.ResamplerFactory_S16().getResampler(0)
    resampler1.setRef(vol)
    vol_resamp = resampler1.doit(t1_to_t2, dims2[0], dims2[1], dims2[2], voxSize2)
    return(vol_resamp)        
        
        
        

######################################


















