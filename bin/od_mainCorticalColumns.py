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
#time od_mainCorticalColumns.py -i /volatile/od243208/brainvisa_manual/ml140175_copyFrom18_11_2014/VoronoiForDist_L_ml140175.nii.gz -d /volatile/od243208/brainvisa_manual/ml140175_columns_Cut/ -k ml140175_L_Cut
# time od_mainCorticalColumns.py -i /volatile/od243208/brainvisa_db_morphologist/dysbrain/ad140157/t1mri/reversed_t1map_2/default_analysis/segmentation/Lgrey_white_ad140157.nii.gz -d /volatile/od243208/brainvisa_manual/ad140157_columns/ -k ad140157_L

# this is the main script to run on a classified GW volume
# it launches scripts by Yann Leprince: dist, heat, isovolume, column-regions to compute 'cortical columns'

from soma import aims, aimsalgo
import sys, glob, os, subprocess, sys, time, timeit
import numpy as np
from optparse import OptionParser

#read in the path and the directory
pathToClassifFile = None #'/volatile/od243208/brainvisa_manual/%s/GWNoInnerSulciSkel_%s_%s.nii.gz' %(realPatientID, realPatientID, realSide)
data_directory = None
keyWord = None

parser = OptionParser('Calculate iso-volumes and cortical columns')
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
    
if options.keyWord is None:
    print >> sys.stderr, 'New: exit. No keyword for results was given'
    sys.exit(1)
else:
    keyWord = options.keyWord
    
# in the given directory create the subdirectory for the results
if not os.path.exists(data_directory):
    os.makedirs(data_directory)

    
    
    
############# todo before this processing : ###############
###################### check: work in space T1 or T2 #######################################################################
###################### check: work with the whole cortex or with the 'cut' and 'dilated' region? ###########################
###################### according to this: launch the cutRoi script, adapt the keyWord and the destination data directory ###
############################################################################################################################
    
    

# launch the distance map calculation # classif file must be here with 100, 200, 0 (no 50 and 150)
t0dist = timeit.default_timer()
subprocess.check_call(['time', 'od_distmapsMain.py', 
'-i', pathToClassifFile, '-d', data_directory, '-k', keyWord])
t1dist = timeit.default_timer()

# launch the heat map calculation # classif file must be here with 100, 200, 0 (no 50 and 150)
t0heat = timeit.default_timer()
subprocess.check_call(['time', 'od_heatMain.py', 
'-i', pathToClassifFile, '-d', data_directory, '-k', keyWord])
t1heat = timeit.default_timer()

# launch the isovolumes calculation # classif file must be here with 100, 200, 0  (no 50 and 150)
t0iso = timeit.default_timer()
subprocess.check_call(['time', 'od_isovolumeMain.py', 
'-i', pathToClassifFile, '-d', data_directory, '-k', keyWord])
t1iso = timeit.default_timer()

# launch the calculation of the 'cortical columns' # classif file must be here with 100, 200, 0 , 50 and 150
t0col = timeit.default_timer()
pathToClassifWithBorders = data_directory + 'dist/classif_with_outer_boundaries_%s.nii.gz' %(keyWord)
subprocess.check_call(['time', 'od_column-regionsMain.py', 
'-i', pathToClassifWithBorders, '-d', data_directory, '-k', keyWord])
t1col = timeit.default_timer()

tDist = t1dist - t0dist
tHeat = t1heat - t0heat
tIso = t1iso - t0iso
tCol = t1col - t0col
header = ['sizeX', 'sizeY', 'sizeZ', 'voxelX', 'voxelY', 'voxelZ', 'tDist', 'tHeat', 'tIso', 'tCol', 'tTotal']
volClassif = aims.read(pathToClassifFile)
headerClassif = volClassif.header()
content = [str(headerClassif['sizeX']), str(headerClassif['sizeY']), str(headerClassif['sizeZ']), str(headerClassif['voxel_size'][0]), str(headerClassif['voxel_size'][1]), str(headerClassif['voxel_size'][2]), str(np.round(tDist)), str(np.round(tHeat)), str(np.round(tIso)), str(np.round(tCol)), str(np.round(tDist + tHeat + tIso + tCol))]
headerLine = ('\t').join(header)
contentLine = ('\t').join(content)
f = open(data_directory + "statFile_%s.txt" %(keyWord), "w")
f.write(headerLine + '\n')
f.write(contentLine + '\n')
f.close()
    
    
    
    
    
