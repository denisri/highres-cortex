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
#od_heatMain.py -i /volatile/od243208/brainvisa_manual/ml140175/VoronoiForDist_L_ml140175.nii.gz -d /volatile/od243208/brainvisa_manual/ml140175/ -k ml140175_L_Cut
#od_heatMain.py -i /volatile/od243208/brainvisa_manual/ml140175/GWNoInnerSulciSkel_ml140175_L.nii.gz -d /volatile/od243208/brainvisa_manual/ml140175/ -k ml140175_L


from soma import aims, aimsalgo
import sys, glob, os, subprocess, sys, time
import numpy as np
from optparse import OptionParser
import highres_cortex.cortex_topo, highres_cortex.div_gradn

#read in the path and the directory
pathToClassifFile = None 
data_directory = None
result_directory = None
keyWord = None
workOnT1inT2Space = False
n_iter = 500 # 500
time_step = 0.01
n_iter2 = 100 # 100
time_step2 = 0.001


parser = OptionParser('Calculate the heat map in a cortical region - OD')
parser.add_option('-i', dest='pathToClassifFile', help='Path to the volume with labeled cortex (100), and white matter (200)')   # if nothing is given: exit
parser.add_option('-d', dest='data_directory', help='directory for the results') # if nothing is given exit
parser.add_option('-k', dest='keyWord', help='KeyWord for the result files (including the patient ID and the hemisphere)') # if nothing is given exit
parser.add_option('-r', dest='workOnT1inT2Space', action = 'store_true', help='Select if working on with T1 images that were resampled into T2 space and then processed by Morphologist. False is default') 


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
    result_directory = data_directory + 'heat/'
    
if options.keyWord is None:
    print >> sys.stderr, 'New: exit. No keyword for results was given'
    sys.exit(1)
else:
    keyWord = options.keyWord
    
if options.workOnT1inT2Space is not None:
    # need to change parameters for the heat equation calculation
    n_iter = 200 # 500
    time_step = 0.04
    n_iter2 = 10 # 100
    time_step2 = 0.001
   
# in the given directory create the subdirectory for the results
if not os.path.exists(result_directory):
    os.makedirs(result_directory)

print '############################################# starting od_heatMain_NEW.py #####################################################'

############# new version by Yann (July 2015):
#ylLaplacian --classif ../classif.nii.gz --output heat.nii.gz


subprocess.check_call(['ylLaplacian', '--classif', pathToClassifFile, '--output', result_directory + 'heat_%s.nii.gz' %(keyWord)])
volHeat = aims.read(result_directory + 'heat_%s.nii.gz' %(keyWord))

# Normalized gradient's divergence
vol_divGrad = highres_cortex.div_gradn.divergence_gradient(volHeat)
aims.write(vol_divGrad, result_directory + "heat_div_gradn_%s.nii.gz" %(keyWord))





