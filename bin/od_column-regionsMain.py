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
# time od_column-regionsMain.py -i /volatile/od243208/brainvisa_manual/ml140175/dist/classif_with_outer_boundaries_ml140175_L.nii.gz -d /volatile/od243208/brainvisa_manual/ad140157_columns/ -k ad140157_L

# od_column-regionsMain.py -i /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/af140169/af140169_T1inT2_ColumnsCutNew20It/dist/classif_with_outer_boundaries_af140169_R_cut_noSulci_extended.nii.gz -d /neurospin/lnao/dysbrain/testBatchColumnsExtrProfiles/af140169/af140169_T1inT2_ColumnsCutNew20It/ -k af140169_R_cut_noSulci_extended

from soma import aims, aimsalgo
import sys, glob, os, subprocess, sys, time
import numpy as np
from optparse import OptionParser
import highres_cortex.cortex_topo, highres_cortex.div_gradn, highres_cortex.od_get_exchanged_propvol, highres_cortex.od_relabel_conjunction, highres_cortex.od_relabel, highres_cortex.od_randomize_labels

#read in the path and the directory
pathToClassifFile = None 
data_directory = None
result_directory = None
heat_directory = None
keyWord = None

parser = OptionParser('Calculate column-regions in a cortical region')
parser.add_option('-i', dest='pathToClassifFile', help='Path to the volume with labeled cortex (100), and white matter (200), as well as the borders (50 and 150)')   # if nothing is given: exit
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
    result_directory = data_directory + 'column_regions/'
    heat_directory = data_directory + 'heat/'
    iso_directory = data_directory + 'isovolume/'

    
if options.keyWord is None:
    print >> sys.stderr, 'New: exit. No keyword for results was given'
    sys.exit(1)
else:
    keyWord = options.keyWord
    
# in the given directory create the subdirectory for the results
if not os.path.exists(result_directory):
    os.makedirs(result_directory)


#AimsThreshold -b -m eq -t 50 \
    #-i /volatile/od243208/brainvisa_manual/ml140175/classif_with_outer_boundaries_ml140175_L.nii.gz \
    #-o /volatile/od243208/brainvisa_manual/ml140175/CSF_interface_ml140175_L.nii
volClassif = aims.read(pathToClassifFile)
arrSurfCSF = np.array(volClassif, copy = False)
arrSurfCSF[np.where(arrSurfCSF != 50)] = 0
arrSurfCSF[np.where(arrSurfCSF == 50)] = 32767
aims.write(volClassif, result_directory + 'CSF_interface_%s.nii' % (keyWord))   # OK
       
#AimsThreshold -b -m eq -t 150 \
    #-i ../classif_with_outer_boundaries.nii.gz \
    #-o white_interface.nii
volClassif = aims.read(pathToClassifFile)
arrSurfWhite = np.array(volClassif, copy = False)
arrSurfWhite[np.where(arrSurfWhite != 150)] = 0
arrSurfWhite[np.where(arrSurfWhite == 150)] = 32767
aims.write(volClassif, result_directory + 'white_interface_%s.nii' % (keyWord))   # OK
    

#ylLabelEachVoxel --verbose \
    #-i CSF_interface.nii.gz \
    #-o CSF_labelled_interface.nii \
    #--first-label 100000001
subprocess.check_call(['ylLabelEachVoxel', '--verbose', '-i', result_directory + 'CSF_interface_%s.nii' % (keyWord), '-o', result_directory + 'CSF_labelled_interface_%s.nii' % (keyWord), '--first-label', '100000001'])       # OK
    
#ylLabelEachVoxel --verbose \
    #-i white_interface.nii.gz \
    #-o white_labelled_interface.nii \
    #--first-label 200000001
subprocess.check_call(['ylLabelEachVoxel', '--verbose', '-i', result_directory + 'white_interface_%s.nii' % (keyWord), '-o', result_directory + 'white_labelled_interface_%s.nii' % (keyWord), '--first-label', '200000001'])      # OK

#AimsThreshold -b --fg -1 -m di -t 100 \
    #-i ../classif.nii.gz \                  # can take the classif with outer boundaries! as cortex is the same there
    #-o negative_outside_cortex.nii
volClassif = aims.read(pathToClassifFile)
arrNegOutCortex = np.array(volClassif, copy = False)
arrNegOutCortex[np.where(arrNegOutCortex != 100)] = -1
arrNegOutCortex[np.where(arrNegOutCortex == 100)] = 0
aims.write(volClassif, result_directory + 'negative_outside_cortex_%s.nii' % (keyWord))      # OK

#AimsFileConvert -t S32 \
    #-i negative_outside_cortex.nii \
    #-o negative_outside_cortex_S32.nii
c = aims.Converter(intype=volClassif, outtype=aims.Volume('S32'))
volNegOutCortex = c(volClassif)
aims.write(volNegOutCortex, result_directory + 'negative_outside_cortex_S32_%s.nii' % (keyWord))    # OK

#AimsMerge -m sv \
    #-i negative_outside_cortex_S32.nii \
    #-M CSF_labelled_interface.nii \
    #-o CSF_labelled_interface_negative_outside.nii
arrNegOutCortex = np.array(volNegOutCortex, copy = False)
volCSFLabelInt = aims.read(result_directory + 'CSF_labelled_interface_%s.nii' % (keyWord))
arrCSFLabelInt = np.array(volCSFLabelInt, copy = False)
arrNegOutCortex[arrCSFLabelInt != 0] = arrCSFLabelInt[arrCSFLabelInt != 0]
aims.write(volNegOutCortex, result_directory + 'CSF_labelled_interface_negative_outside_%s.nii' % (keyWord))   # OK
    
#AimsMerge -m ao -v 200000000 \
    #-i CSF_labelled_interface_negative_outside.nii \
    #-M white_labelled_interface.nii \
    #-o propvol_CSF_labels.nii.gz
volWhiteLabInt = aims.read(result_directory + 'white_labelled_interface_%s.nii' % (keyWord))
arrWhiteLabInt = np.array(volWhiteLabInt, copy = False)
arrNegOutCortex[arrWhiteLabInt != 0] = 200000000
aims.write(volNegOutCortex, result_directory + 'propvol_CSF_labels_%s.nii.gz' % (keyWord))    # OK

#AimsMerge -m sv \
    #-i negative_outside_cortex_S32.nii \
    #-M white_labelled_interface.nii \
    #-o white_labelled_interface_negative_outside.nii
volNegOutCortex = aims.read(result_directory + 'negative_outside_cortex_S32_%s.nii' % (keyWord))
arrNegOutCortex = np.array(volNegOutCortex, copy = False)
arrNegOutCortex[arrWhiteLabInt != 0] = arrWhiteLabInt[arrWhiteLabInt != 0]
aims.write(volNegOutCortex, result_directory + 'white_labelled_interface_negative_outside_%s.nii' % (keyWord))       # OK  
    
#AimsMerge -m ao -v 100000000 \
    #-i white_labelled_interface_negative_outside.nii \
    #-M CSF_labelled_interface.nii \
    #-o propvol_white_labels.nii.gz
arrNegOutCortex[np.where(arrCSFLabelInt != 0)] = 100000000
aims.write(volNegOutCortex, result_directory + 'propvol_white_labels_%s.nii.gz' % (keyWord))     # OK 

subprocess.check_call(['time', 'ylPropagateAlongField', '--verbose', '--grad-field', heat_directory + 'heat_%s.nii.gz' % (keyWord), '--seeds', result_directory + 'propvol_CSF_labels_%s.nii.gz' % (keyWord), '--step', '-0.05', '--target-label', '200000000', '--output', result_directory + 'heat_CSF_labels_on_white_%s.nii.gz' % (keyWord)])   # OK 
#ylPropagateAlongField --verbose \
    #--grad-field ../heat/heat.nii.gz \
    #--seeds propvol_CSF_labels.nii.gz \
    #--step -0.05 \
    #--target-label 200000000 \
    #--output heat_CSF_labels_on_white.nii.gz
#time for the whole cortex 1:27.7

subprocess.check_call(['time', 'ylPropagateAlongField', '--verbose', '--grad-field', heat_directory + 'heat_%s.nii.gz' % (keyWord), '--seeds', result_directory + 'propvol_white_labels_%s.nii.gz' % (keyWord), '--step', '0.05', '--target-label', '100000000', '--output', result_directory + 'heat_white_labels_on_CSF_%s.nii.gz' % (keyWord)])          # OK 
#ylPropagateAlongField --verbose \
    #--grad-field ../heat/heat.nii.gz \
    #--seeds propvol_white_labels.nii.gz \
    #--step 0.05 \
    #--target-label 100000000 \
    #--output heat_white_labels_on_CSF.nii.gz
#time for the whole cortex 1:43.87

volCSF_labels_on_white = aims.read(result_directory + 'heat_CSF_labels_on_white_%s.nii.gz' % (keyWord))
volwhite_labels_on_CSF = aims.read(result_directory + 'heat_white_labels_on_CSF_%s.nii.gz' % (keyWord))
volClassif = aims.read(pathToClassifFile)
volExchangedPropVol = highres_cortex.od_get_exchanged_propvol.getExchangedPropagationVolume(volCSF_labels_on_white, volwhite_labels_on_CSF, volClassif, result_directory, keyWord)
aims.write(volExchangedPropVol, result_directory + "exchanged_propvol_%s.nii.gz" %(keyWord))            
#python get_exchanged_propvol.py  # -> exchanged_propvol.nii.gz


    
# Why is the previous step necessary?
#
# The obvious alternative is to do exactly as described in the OHBM paper: do
# the projections on the original labels of each voxel.
#
# The previous case aggregates the adjacent voxels of one interface that point
# towards the same voxel on the other interface. This reduces
# over-segmentation.
#
# Another way of reducing over-segmentation would be to aggregate together
# voxels that have one projection in common, instead of both (see conjunction
# step later on). But this introduces the problem of transitivity. This was
# investigated previously on the ferret data (under the name Billiard), but was
# considered a dead-end and the above solution seems to solve this problem most
# efficiently.


# There is a problem with the propagation of labels: the step size is fixed,
# which means that sometimes the point can skip the corner of a voxel, and thus
# go directly from a bulk voxel to an outside voxel. In this case it is
# recorded as a "dead-end" advection path, no resulting label is recorded and
# it appears as zero in the result.
#
# This problem also appears in the previous "exchange" step, but is mitigated
# by the subsequent connex component detection (each failed propagation is
# assigned a different label).
#
# Quick fix: fix the conjunction step to not aggregate zeros.
#
# TODO: the proper way to fix this would be to force the advection path to
# respect the boundaries of voxels, so that the corner of voxels cannot be
# skipped over. This would also prevent the advection path from crossing the
# thin CSF surface within the sulcus (comes from skeleton).

# I could take into account the fake cortex–CSF interface that exists at the
# cut plane, by assigning it a special label (e.g. 500000000) in the
# exchanged_propvol label. It would then need to be treated specially: any
# voxel that projects onto this label would be excluded from the region list,
# and thus would not take part in the merging step. This would prevent the
# creation of regions that connect to this spurious surface, but this would not
# prevent the nearby regions from being deformed by the perturbation of the
# field. It would thus probably be overkill to implement this special case.
# Care is needed when dealing with regions close to the cut plane anyway.


#AimsMerge -m oo -l 150 -v 0 \
    #-i exchanged_propvol.nii.gz \
    #-M ../classif_with_outer_boundaries.nii.gz \
    #-o ./exchanged_labels_on_CSF.nii
arrExchangedPropVol = np.array(volExchangedPropVol, copy = False)   
arrClassif = np.array(volClassif, copy = False)
arrExchangedPropVol[arrClassif == 150] = 0
aims.write(volExchangedPropVol, result_directory + 'exchanged_labels_on_CSF_%s.nii' %(keyWord))         # OK

    
#AimsMerge -m oo -l 50 -v 0 \
    #-i ./exchanged_propvol.nii.gz \
    #-M ../classif_with_outer_boundaries.nii.gz \
    #-o ./exchanged_labels_on_white.nii
volExchangedPropVol = aims.read(result_directory + "exchanged_propvol_%s.nii.gz" %(keyWord))
arrExchangedPropVol = np.array(volExchangedPropVol, copy = False)   
arrExchangedPropVol[arrClassif == 50] = 0
aims.write(volExchangedPropVol, result_directory + 'exchanged_labels_on_white_%s.nii' %(keyWord))         # OK


#ylPropagateAlongField --verbose \
    #--grad-field ../heat/heat.nii.gz \
    #--seeds exchanged_labels_on_CSF.nii \
    #--step -0.05 \
    #--target-label 0 \
    #--output heat_CSF_on_bulk.nii.gz \
    #--dest-points heat_CSF_points_on_bulk.nii.gz
subprocess.check_call(['time', 'ylPropagateAlongField', '--verbose', '--grad-field', heat_directory + 'heat_%s.nii.gz' % (keyWord), '--seeds',result_directory + 'exchanged_labels_on_CSF_%s.nii' %(keyWord), '--step', '-0.05', '--target-label', '0', '--output', result_directory + 'heat_CSF_on_bulk_%s.nii.gz' % (keyWord), '--dest-points', result_directory + 'heat_CSF_points_on_bulk_%s.nii.gz' % (keyWord)])     
# time for the full cortex:      4:56.95
    
    
#ylPropagateAlongField --verbose \
    #--grad-field ../heat/heat.nii.gz \
    #--seeds exchanged_labels_on_white.nii \
    #--step 0.05 \
    #--target-label 0 \
    #--output heat_white_on_bulk.nii.gz \
    #--dest-points heat_white_points_on_bulk.nii.gz
subprocess.check_call(['time', 'ylPropagateAlongField', '--verbose', '--grad-field', heat_directory + 'heat_%s.nii.gz' % (keyWord), '--seeds',result_directory + 'exchanged_labels_on_white_%s.nii' %(keyWord), '--step', '0.05', '--target-label', '0', '--output', result_directory + 'heat_white_on_bulk_%s.nii.gz' % (keyWord), '--dest-points', result_directory + 'heat_white_points_on_bulk_%s.nii.gz' % (keyWord)])     
# time for the full cortex:    5:59.33

#python relabel_conjunction.py  # -> ./conjunction.nii.gz
vol1 = aims.read(result_directory + 'heat_CSF_on_bulk_%s.nii.gz' % (keyWord))
vol2 = aims.read(result_directory + 'heat_white_on_bulk_%s.nii.gz' % (keyWord))
volRelabeledConj = highres_cortex.od_relabel_conjunction.relabel_conjunctions(vol1, vol2)
aims.write(volRelabeledConj, result_directory + 'conjunction_%s.nii.gz' % (keyWord))


#ylMergeCortexColumnRegions --verbose 2 \
    #-i conjunction.nii.gz \
    #-o merged.nii \
    #--proj-csf heat_CSF_points_on_bulk.nii.gz \
    #--proj-white heat_white_points_on_bulk.nii.gz \
    #--goal-diameter 1
subprocess.check_call(['time', 'ylMergeCortexColumnRegions', '--verbose', '2', '-i', result_directory + 'conjunction_%s.nii.gz' % (keyWord), '-o',result_directory + 'merged_%s.nii' %(keyWord), '--proj-csf', result_directory + 'heat_CSF_points_on_bulk_%s.nii.gz' % (keyWord), '--proj-white', result_directory + 'heat_white_points_on_bulk_%s.nii.gz' % (keyWord), '--goal-diameter', '1'])     
# time for the full cortex : 0:58.83

#python relabel.py
vol1 = aims.read(result_directory + 'merged_%s.nii' %(keyWord))
vol2 = highres_cortex.od_relabel.relabel(vol1)
aims.write(vol2, result_directory + 'merged_relabelled_%s.nii.gz' % (keyWord))

#python randomize_labels.py
vol1 = highres_cortex.od_randomize_labels.relabel(vol2)
aims.write(vol1, result_directory + 'merged_randomized_%s.nii.gz' %(keyWord))

#print np.max(np.array(vol1)) # number of different columns 111067


## test for another diameter of cortical columns. E.g. of 3 mm, and 5 mm, and 9mm
diams = [3, 5, 7, 9]
for diam in diams:

    subprocess.check_call(['ylMergeCortexColumnRegions', '--verbose', '2', '-i', result_directory + 'conjunction_%s.nii.gz' % (keyWord), '-o',result_directory + 'merged_%s_diam%s.nii' %(keyWord, diam), '--proj-csf', result_directory + 'heat_CSF_points_on_bulk_%s.nii.gz' % (keyWord), '--proj-white', result_directory + 'heat_white_points_on_bulk_%s.nii.gz' % (keyWord), '--goal-diameter', str(diam)])     

    #python relabel.py
    vol1 = aims.read(result_directory + 'merged_%s_diam%s.nii' %(keyWord, diam))
    vol2 = highres_cortex.od_relabel.relabel(vol1)
    aims.write(vol2, result_directory + 'merged_relabelled_%s_diam%s.nii.gz' % (keyWord, diam))

    #python randomize_labels.py
    vol1 = highres_cortex.od_randomize_labels.relabel(vol2)
    aims.write(vol1, result_directory + 'merged_randomized_%s_diam%s.nii.gz' %(keyWord, diam))





    
    
    
    
