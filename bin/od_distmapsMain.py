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

import numpy as np
from soma import aims, aimsalgo
import highres_cortex.cortex_topo
import glob, os, subprocess, sys, time
from optparse import OptionParser

# for example: launch this script for patient ml140175 L hemisphere. For the whole cortex, or for the cut out regions
#od_distmapsMain.py -i /volatile/od243208/brainvisa_manual/ml140175/GWNoInnerSulciSkel_ml140175_L.nii.gz -d /volatile/od243208/brainvisa_manual/ml140175/ -k ml140175_L
#od_distmapsMain.py -i /volatile/od243208/brainvisa_manual/ml140175/VoronoiForDist_L_ml140175.nii.gz -d /volatile/od243208/brainvisa_manual/ml140175/ -k ml140175_L_Cut

pathToClassifFile = None
data_directory = None
result_directory = None
keyWord = None

parser = OptionParser('Calculate distance map in a cortical region')
parser.add_option('-i', dest='pathToClassifFile', help='Path to the volume with labeled cortex (100), and white matter (200)')   # if nothing is given: exit
parser.add_option('-d', dest='data_directory', help='directory for the results') # if nothing is given exit
parser.add_option('-k', dest='keyWord', help='KeyWord for the result files (including the patient ID and the hemisphere)') # if nothing is given exit

options, args = parser.parse_args(sys.argv)
print options
print args

if options.pathToClassifFile is None:
    print >> sys.stderr, 'New: exit. No classification volume was given'
    sys.exit(1)
else:
    pathToClassifFile = options.pathToClassifFile
    
if options.data_directory is None:
    print >> sys.stderr, 'New: exit. No directory for results was given'
    sys.exit(1)
else:
    data_directory = options.data_directory
    result_directory = data_directory + '/dist/'
    
if options.keyWord is None:
    print >> sys.stderr, 'New: exit. No keyword for results was given'
    sys.exit(1)
else:
    keyWord = options.keyWord
    
# in the given directory create the subdirectory for the results
if not os.path.exists(result_directory):
    os.makedirs(result_directory)

print '####################################### starting od_dystmaps.py ##############################################'
classif = aims.read(pathToClassifFile)
#dist_from_white = highres_cortex.cortex_topo.fastmarching_negative(classif, [100], [200], 150)
dist_from_white = highres_cortex.cortex_topo.fastmarching_negative(classif, [100], [200], 150, False)

aims.write(dist_from_white, result_directory + 'distwhite_%s.nii.gz' %(keyWord))
print '####################################### done : dist_from_white ###############################################'

# need to visualize the distance map. Therefore get all negative values and set them to some positive values
volDistFromWhiteCopy = aims.read(result_directory + 'distwhite_%s.nii.gz' %(keyWord))
arrDistFromWhiteCopy = np.array(volDistFromWhiteCopy, copy = False)
arrDistFromWhiteCopy[arrDistFromWhiteCopy < 0] = np.max(arrDistFromWhiteCopy) + 5
aims.write(volDistFromWhiteCopy, result_directory + 'distwhiteVisu_%s.nii.gz' %(keyWord))

# now calculate distance to CSF
#dist_from_CSF = highres_cortex.cortex_topo.fastmarching_negative(classif, [100], [0], 50)
dist_from_CSF = highres_cortex.cortex_topo.fastmarching_negative(classif, [100], [0], 50, False)
print '####################################### done : dist_from_CSF #################################################'
aims.write(dist_from_CSF, result_directory + 'distCSF_%s.nii.gz' %(keyWord))

volDistFromCSFCopy = aims.read(result_directory + 'distCSF_%s.nii.gz' %(keyWord))
arrDistFromCSFCopy = np.array(volDistFromCSFCopy, copy = False)

arrDistFromCSFCopy[arrDistFromCSFCopy < 0] = np.max(arrDistFromCSFCopy) + 5
aims.write(volDistFromCSFCopy, result_directory + 'distCSFVisu_%s.nii.gz' %(keyWord))
aims.write(classif, result_directory + 'classif_with_outer_boundaries_%s.nii.gz' %(keyWord))


#AimsThreshold -m di -t -1 \
    #-i /volatile/od243208/brainvisa_manual/ml140175/classif_with_outer_boundaries_ml140175_L.nii.gz \
    #-o classif_with_background_ml140175_L.nii
#subprocess.call(['AimsThreshold', '-m', 'di', '-t', '-1', '-i', result_directory + 'classif_with_outer_boundaries_%s.nii.gz' %(keyWord), '-o', result_directory + 'classif_with_background_%s.nii' %(keyWord)])
arrClassif = np.array(classif, copy = False)
arrClassif[arrClassif == -1] = 32767            # there are in theory no voxels with the value -1
aims.write(classif, result_directory + 'classif_with_background_%s.nii' %(keyWord))

#subprocess.call(['AimsMerge', '-m', 'sv', '-i', pathToClassifFile, '-M', result_directory + 'classif_with_background_%s.nii' %(keyWord), '-o', result_directory + 'classif_with_outer_boundaries_%s.nii.gz' %(keyWord)])  
#AimsMerge -m sv \
    #-i /volatile/od243208/brainvisa_manual/ml140175/GWNoInnerSulciSkel_ml140175_L.nii.gz \
    #-M /volatile/od243208/brainvisa_manual/ml140175/classif_with_background_ml140175_L.nii \
    #-o /volatile/od243208/brainvisa_manual/ml140175/classif_with_outer_boundaries_ml140175_L.nii.gz
classif = aims.read(pathToClassifFile)
arrClassifOriginal = np.array(classif, copy = False)
arrClassifOriginal[arrClassif != 0] = arrClassif[arrClassif != 0]
aims.write(classif, result_directory + 'classif_with_outer_boundaries_%s.nii' %(keyWord))




