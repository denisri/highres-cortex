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
# this function will take the Morphologist's processing results of the original T1 volume, 
# the calculated transformation into the new T2 space.

import random
from soma import aims
import subprocess
from optparse import OptionParser
from soma import aims, aimsalgo
from scipy.stats import mode
import sys, glob, os, os.path, subprocess, sys, time, timeit
import numpy as np
import highres_cortex.od_cutOutRois

brainvisa_db_neurospin = '/neurospin/lnao/dysbrain/brainvisa_db_morphologist/dysbrain/'
brainvisa_raw_niftis = '/neurospin/lnao/dysbrain/raw_niftis/'
brainvisa_T1intoT2 = '/neurospin/lnao/dysbrain/testT1T2/'
pathToTextures = '/neurospin/lnao/dysbrain/randomized_flipped_data/manual_work/'

#volT1 = aims.read('/volatile/od243208/brainvisa_db_morphologist/dysbrain/ac140155/t1mri/reversed_t1map_2/ac140155.nii.gz')
#volT2 = aims.read('/volatile/od243208/raw_niftis/ac140155/20140703_152737t2spctraiso05mmVOIs011a1001.nii.gz')
def transformT1toT2(volT1, transfToT2, volT2):
    import numpy as np
    sizes = volT2.getSize()
    voxSize = volT2.getVoxelSize()
    print 'sizes : ', sizes[0], sizes[1], sizes[2]
    print 'voxSize : ', voxSize[0], voxSize[1], voxSize[2]
    resampler1 = aims.ResamplerFactory_S16().getResampler(0)    
    #resampler1 = aims.ResamplerFactory_FLOAT().getResampler(0)

    print 'type(sizes[0])', type(float(sizes[0]))
    print 'type(voxSize)', type(voxSize)
    print 'type(transfToT2)', type(transfToT2)
    resampler1.setRef(volT1)
    vol_resamp = resampler1.doit(transfToT2, float(sizes[0]), float(sizes[1]), float(sizes[2]), voxSize)    
    #vol_resamp = resampler1.doit(transfToT2, sizes[0], sizes[1], sizes[2], voxSize)

    return(vol_resamp)
        
if __name__ == '__main__':
    
    outputFile = None
    realPatientID = None
    keyWord = None
    directory = None    
    t2directory = None
    recursiveInDirectory = None
    resultDirectory = None

    parser = OptionParser('Transform Morphologist\'s results for T1 volume into the new T2 space')
    parser.add_option('-p', dest='realPatientID', help='realPatientID')
    parser.add_option('-k', dest='keyWord', help='keyWord')
    parser.add_option('-d', dest='directory', help='directory')    
    parser.add_option('-t', dest='t2directory', help='t2directory')
    parser.add_option('-r', dest='recursiveInDirectory', help='recursiveInDirectory')
    options, args = parser.parse_args(sys.argv)
    print options
    print args   
    
    subjectList=[]
    
    if options.directory is None:
        print >> sys.stderr, 'New: exit. no directory given'
        sys.exit(1)
    else:
        directory = options.directory     
        
    if options.t2directory is None:
        print >> sys.stderr, 'New: exit. no t2directory given'
        sys.exit(1)
    else:
        t2directory = options.t2directory      

    if options.realPatientID is not None:
        realPatientID = options.realPatientID     
        subjectList.append(realPatientID)
                
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
                 
    # now process all the required subjects

      
    for realPatientID in subjectList:    
        print '####################### process subject  ', realPatientID, '  ##############################################'
        # check if the destination folder exists. If not - create it
        resultDirectory = directory + realPatientID + '/'

        # check if the required transformation file is available
        finder2 = aims.Finder()
        transT2list = glob.glob(brainvisa_T1intoT2 + realPatientID + '/%s_t2.trm' %(realPatientID))   # ok only if there is one t2-file, or the first one of the t2-files is correct
        if len(transT2list) == 1: # finder2.check(volT2):
            transT2 = aims.read(transT2list[0])
            print 'Took the transformation into the new T2 space: ', transT2list[0]
        else:
            print 'found ', len(transT2list), ' of transT2 files. Continue with the next subject'
            continue        
        
        # find the volT2 into whose space we'll transform other volumes
        volT2 = aims.read(brainvisa_T1intoT2 + realPatientID + '/%s_NewT2.nii.gz' %(realPatientID))  
        
        # find all the volumes to be processed
        if not os.path.exists(resultDirectory):
            os.makedirs(resultDirectory)

        keyWord = realPatientID
        vols = []
        newVols = []
        
        # for the left and the right side
        for realSide in ['L', 'R']:            
            # 1. GW classif
            volsGWlist = glob.glob(brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/default_analysis/segmentation/%sgrey_white_%s.nii.gz' %(realSide, realPatientID))   # ok only if there is one file
            if len(volsGWlist) == 1:
                volGW = aims.read(volsGWlist[0])
                print 'Took the volGW: ', volsGWlist[0]
                vols.append(volGW)                
                volGWT2 = transformT1toT2(volGW, transT2, volT2)
                aims.write(volGWT2, resultDirectory + keyWord + '_GW_%s_T1inNewT2.nii.gz' %(realSide))
                
            # 2. skeletons
            sulcilist = glob.glob(brainvisa_db_neurospin + realPatientID + '/t1mri/reversed_t1map_2/default_analysis/folds/3.1/default_session_auto/segmentation/%sSulci_%s_default_session_auto.nii.gz' %(realSide, realPatientID))   # ok only if there is one file
            if len(sulcilist) == 1: 
                sulci = aims.read(sulcilist[0])
                print 'Took the sulci: ', sulcilist[0]
                vols.append(sulci)    
                sulciT2 = transformT1toT2(sulci, transT2, volT2)
                aims.write(sulciT2, resultDirectory + keyWord + '_sulciSkel_%s_T1inNewT2.nii.gz' %(realSide))
        ## transform these volumes into the given new T2 space
        #newVols = map(transformT1toT2, vols, transT2, volT2)
            
           
            
         
            
            
    
    
    
            
    