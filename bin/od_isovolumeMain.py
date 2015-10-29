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
#od_isovolumeMain.py -i /volatile/od243208/brainvisa_manual/ml140175/GWNoInnerSulciSkel_ml140175_L.nii.gz -d /volatile/od243208/brainvisa_manual/ml140175/ -k ml140175_L

from soma import aims, aimsalgo
import sys, glob, os, subprocess, sys, time
import numpy as np
from optparse import OptionParser
import highres_cortex.cortex_topo, highres_cortex.div_gradn

#read in the path and the directory
pathToClassifFile = None #'/volatile/od243208/brainvisa_manual/%s/GWNoInnerSulciSkel_%s_%s.nii.gz' %(realPatientID, realPatientID, realSide)
data_directory = None
result_directory = None
heat_directory = None
keyWord = None

parser = OptionParser('Calculate iso-volumes in a cortical region')
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
    result_directory = data_directory + 'isovolume/'
    heat_directory = data_directory + 'heat/'

    
if options.keyWord is None:
    print >> sys.stderr, 'New: exit. No keyword for results was given'
    sys.exit(1)
else:
    keyWord = options.keyWord
    
# in the given directory create the subdirectory for the results
if not os.path.exists(result_directory):
    os.makedirs(result_directory)


#subprocess.check_call(['time', 'AimsThreshold', '-b', '--fg', '1', '-m', 'eq', '-t', '100', '-i', pathToClassifFile, '-o', result_directory + 'domain_%s.nii' %(keyWord)])
volClassif = aims.read(pathToClassifFile)
arrClassif = np.array(volClassif, copy = False) 
arrClassif[arrClassif != 100] = 0
arrClassif[arrClassif == 100] = 1
aims.write(volClassif, result_directory + 'domain_%s.nii' %(keyWord))
 
 
subprocess.check_call(['time', 'ylAdvectTubes', '--verbose', '--step', '0.05', '--domain', result_directory + 'domain_%s.nii' %(keyWord), '--grad-field', heat_directory + 'heat_%s.nii.gz' %(keyWord), '--divergence', heat_directory + 'heat_div_gradn_%s.nii.gz' %(keyWord), '--output-volumes', result_directory + 'white-tube-volumes_%s.nii.gz' % (keyWord), '--output-surfaces', result_directory + 'white-tube-surfaces_%s.nii.gz' % (keyWord)])
# time for the whole cortex : 6m48.759s

volWV = aims.read(result_directory + 'white-tube-volumes_%s.nii.gz' % (keyWord))
volWS = aims.read(result_directory + 'white-tube-surfaces_%s.nii.gz' % (keyWord))
volWVS = volWV / volWS
aims.write(volWVS, result_directory + 'white-tube-VoverS_%s.nii.gz' % (keyWord))
del volWS, volWVS

#time cartoLinearComb.py -f 'I1/I2' \
    #-i /volatile/od243208/brainvisa_manual/ml140175/white-tube-volumes_ml140175_L.nii.gz \
    #-i /volatile/od243208/brainvisa_manual/ml140175/white-tube-surfaces_ml140175_L.nii.gz \
    #-o /volatile/od243208/brainvisa_manual/ml140175/white-tube-VoverS_ml140175_L.nii.gz


subprocess.check_call(['time', 'ylAdvectTubes', '--verbose', '--step', '-0.05', '--domain', result_directory + 'domain_%s.nii' %(keyWord), '--grad-field', heat_directory + 'heat_%s.nii.gz' %(keyWord), '--divergence', heat_directory + 'heat_div_gradn_%s.nii.gz' %(keyWord), '--output-volumes', result_directory + 'pial-tube-volumes_%s.nii.gz' % (keyWord), '--output-surfaces', result_directory + 'pial-tube-surfaces_%s.nii.gz' % (keyWord)])
#time ylAdvectTubes --verbose \
    #--step -0.05 \
    #--domain /volatile/od243208/brainvisa_manual/ml140175/domain_ml140175_L.nii.gz \
    #--grad-field /volatile/od243208/brainvisa_manual/ml140175/heat_ml140175_L.nii.gz \
    #--divergence /volatile/od243208/brainvisa_manual/ml140175/heat_div_gradn_ml140175_L.nii.gz \
    #--output-volumes /volatile/od243208/brainvisa_manual/ml140175/pial-tube-volumes_ml140175_L.nii.gz \
    #--output-surfaces /volatile/od243208/brainvisa_manual/ml140175/pial-tube-surfaces_ml140175_L.nii.gz  
## time for the whole cortex : 5m48.343s

volPV = aims.read(result_directory + 'pial-tube-volumes_%s.nii.gz' % (keyWord))
volPS = aims.read(result_directory + 'pial-tube-surfaces_%s.nii.gz' % (keyWord))
volPVS = volPV / volPS
aims.write(volPVS, result_directory + 'pial-tube-VoverS_%s.nii.gz' % (keyWord))
del volPS, volPVS
#time cartoLinearComb.py -f 'I1/I2' \
    #-i /volatile/od243208/brainvisa_manual/ml140175/pial-tube-volumes_ml140175_L.nii.gz \
    #-i /volatile/od243208/brainvisa_manual/ml140175/pial-tube-surfaces_ml140175_L.nii.gz \
    #-o /volatile/od243208/brainvisa_manual/ml140175/pial-tube-VoverS_ml140175_L.nii.gz

    
volTV = volPV + volWV
aims.write(volTV, result_directory + 'total-tube-volumes_%s.nii.gz' % (keyWord))
#time cartoLinearComb.py -f 'I1+I2' \
    #-i /volatile/od243208/brainvisa_manual/ml140175/pial-tube-volumes_ml140175_L.nii.gz \
    #-i /volatile/od243208/brainvisa_manual/ml140175/white-tube-volumes_ml140175_L.nii.gz \
    #-o /volatile/od243208/brainvisa_manual/ml140175/total-tube-volumes_ml140175_L.nii.gz
    
volPVF = volPV / volTV
aims.write(volPVF, result_directory + 'pial-volume-fraction_%s.nii.gz' % (keyWord))
del volWV, volPV, volTV
#time cartoLinearComb.py -f 'I1/I2' \
    #-i /volatile/od243208/brainvisa_manual/ml140175/pial-tube-volumes_ml140175_L.nii.gz \
    #-i /volatile/od243208/brainvisa_manual/ml140175/total-tube-volumes_ml140175_L.nii.gzz \
    #-o /volatile/od243208/brainvisa_manual/ml140175/pial-volume-fraction_ml140175_L.nii.gz

    

    
    
