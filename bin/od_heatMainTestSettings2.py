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
import sys, glob, os, subprocess, sys, time, timeit
import numpy as np
from optparse import OptionParser
import highres_cortex.cortex_topo, highres_cortex.div_gradn, highres_cortex.od_cutOutRois
from scipy.stats import mode
 

#read in the path and the directory
pathToClassifFile = None 
data_directory = None
result_directory = None
keyWord = None

parser = OptionParser('Calculate the heat map in a cortical region - OD')
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
    result_directory = data_directory + 'heat/'
    
if options.keyWord is None:
    print >> sys.stderr, 'New: exit. No keyword for results was given'
    sys.exit(1)
else:
    keyWord = options.keyWord
    
# in the given directory create the subdirectory for the results
if not os.path.exists(result_directory):
    os.makedirs(result_directory)

print '############################################# starting od_heatMain.py #####################################################'

#subprocess.check_call(['AimsThreshold', '-b', '-m', 'di', '-t', '100', '-i', pathToClassifFile, '-o', result_directory + 'all_but_cortex_%s.nii' %(keyWord)])
#AimsThreshold -b -m di -t 100 \
    #-i ../classif.nii.gz \
    #-o ./all_but_cortex.nii
volClassif = aims.read(pathToClassifFile)
arrClassif = np.array(volClassif, copy = False)
arrClassif[arrClassif != 100] = 32767
arrClassif[arrClassif == 100] = 0
aims.write(volClassif, result_directory + 'all_but_cortex_%s.nii' %(keyWord))

# read in the classification file and convert it into float format    
#AimsFileConvert -t FLOAT \
    #-i ../classif.nii.gz \
    #-o heat.nii.gz
heatmap_before = aims.read(pathToClassifFile, 1)
c = aims.Converter(intype = heatmap_before, outtype = aims.Volume('FLOAT'))
heatmap_before = c(heatmap_before)
aims.write(heatmap_before, result_directory + 'heat_before_%s.nii.gz' %(keyWord))

statFileName = result_directory + "statFile_%s.txt" %(keyWord)
f = open(statFileName, "w")

# Each run refines the previous one
# python heat.py 500 0.01
#os.system("time python /volatile/od243208/brainvisa_sources/perso/domanova/trunk/to_sort/Test_lamination_randomized/heat/od_heat.py 500 0.01 result_directory + 'heat_%s.nii' %keyWord result_directory + 'all_but_cortex_%s.nii' %keyWord result_directory + 'heat_%s.nii' %keyWord")
# time python /volatile/od243208/python/Test_lamination_randomized/heat/od_heat.py 500 0.01 /volatile/od243208/brainvisa_manual/ml140175/heat_ml140175_L.nii.gz /volatile/od243208/brainvisa_manual/ml140175/all_but_cortex_ml140175_L.nii /volatile/od243208/brainvisa_manual/ml140175/heat_ml140175_L.nii.gz

# Yann has used the following settings
# I step:  n_iter = 500, time_step = 0.01
# II step: n_iter = 100, time_step = 0.001
# these settings were good for the resolution 1mm. For the new images with resolution = 0.5mm, they have to be adapted

# todo: test various combinations of n_iter and time step,s check for convergence
# parameters for the first step
#nIter = [200, 500, 800]
#tStep = [0.05, 0.02, 0.01]

#after the first try: adapted the parameters:
#nIter = [500, 800]
#tStep = [0.04, 0.03, 0.02, 0.01, 0.005]

## now the third try: 
#nIter = [200, 500, 800]
#tStep = [0.04, 0.03, 0.02]

# now the fourth try: 
# fifth: the same, as the values are good
#nIter = [300, 500]
#tStep = [0.04, 0.03]

#nIter = [300]   # just to test the diff between the 1st and 2nd step
#tStep = [0.04]  # just to test the diff between the 1st and 2nd step

# test for the initial images:
#nIter = [500]   
#tStep = [0.01]  

## test again: reduced number of iterations
#nIter = [50, 100, 150, 200, 250, 300]
#tStep = [0.04]

# test again: reduced number of iterations
#nIter = [150, 200, 250, 300]
#tStep = [0.04]

## test again: reduced number of iterations - CORRECTED convergence calculation!!
#nIter = [50, 100, 150, 200, 250, 300]
#tStep = [0.05, 0.04, 0.03, 0.02]

# test the SELECTED parameters
nIter = [200]
tStep = [0.04]

# parameters for the second step
#nIterII = [100, 200]
#tStepII = [0.001, 0.002]

#after the first try: adapted the parameters:
#nIterII = [50, 100, 200]
#tStepII = [0.0005, 0.001]

## and the third try for the 2nd step:
#nIterII = [10, 30, 50, 70]
#tStepII = [0.0005, 0.001, 0.002]

## and the fourth try for the 2nd step:
#nIterII = [10, 30, 50]
#tStepII = [0.0005, 0.001]

# and the fifth try for the 2nd step:
#nIterII = [10, 30, 50, 100, 150]
#tStepII = [0.0005, 0.001, 0.002, 0.004]

## test for the initial images
#nIterII = [100]
#tStepII = [0.001]

## test again: reduced number of iterations
#nIterII = [10, 50]
#tStepII = [0.04, 0.001, 0.002]

# test again: reduced number of iterations - CORRECTED convergence calculation!!
#nIterII = [10, 30, 50]
#tStepII = [0.001, 0.002, 0.004]

# test the SELECTED parameters
nIterII = [10]
tStepII = [0.001]

timesI = []


for n_iter in nIter:
    for time_step in tStep:
        print '################### Step I. n_iter = ', n_iter, ' time_step = ', time_step
        f.write('################### Step I. n_iter = ' + str(n_iter) + ' time_step = ' + str(time_step) + '\n')
        # the first step of heat calculation
        t0Step1 = timeit.default_timer()
        heatmap_before = aims.read(result_directory + 'heat_before_%s.nii.gz' %(keyWord), 1)
        mask = aims.read(result_directory + 'all_but_cortex_%s.nii' %(keyWord), 1)
        aimsmask = aims.AimsData(mask)  # Important for reference-counting!
        diff = aimsalgo.MaskedDiffusionSmoother_FLOAT(time_step)
        diff.setMask(aimsmask, 32767)
        heatmap = diff.doSmoothing(heatmap_before, n_iter, True)
        aims.write(heatmap, result_directory + 'heat_I_%s_%s_%s.nii.gz' %(n_iter, time_step, keyWord))
        t1Step1 = timeit.default_timer()
        f.write('Time for the Step I is ' + str(np.round(t1Step1 - t0Step1)) + '\n')
        vol1step = aims.read(result_directory + 'heat_I_%s_%s_%s.nii.gz' %(n_iter, time_step, keyWord))
        arr0 = np.array(vol1step, copy = False)
                        
        # the second step, finer step size
        for n_iter2 in nIterII:
            for time_step2 in tStepII:
                t0Step2 = timeit.default_timer()
                print '*************************  Step II. n_iter2 = ', n_iter2, ' time_step2 = ', time_step2
                f.write('*************************  Step II. n_iter2 = ' + str(n_iter2) + ' time_step2 = ' + str(time_step2) + '\n')
                diff = aimsalgo.MaskedDiffusionSmoother_FLOAT(time_step2)
                diff.setMask(aimsmask, 32767)
                heatmap = aims.read(result_directory + 'heat_I_%s_%s_%s.nii.gz' %(n_iter, time_step, keyWord))
                heatmap = diff.doSmoothing(heatmap, n_iter2, True)
                aims.write(heatmap, result_directory + 'heat_I_%s_%s_II_%s_%s_%s.nii.gz' %(n_iter, time_step, n_iter2, time_step2, keyWord))
                volHeat = aims.read(result_directory + 'heat_I_%s_%s_II_%s_%s_%s.nii.gz' %(n_iter, time_step, n_iter2, time_step2, keyWord))
                arr1 = np.array(volHeat, copy = False)
                t1Step2 = timeit.default_timer()
                f.write('Time for the Step II is ' + str(np.round(t1Step2 - t0Step2)) + '\n')
                # difference between the 1st step and the smoothing result:
                maxDiff1smoothing = np.max(np.abs(arr1 - arr0))
                print 'Max difference between the 1st and the 2nd step is  ', str(maxDiff1smoothing)
                f.write('Max difference between the 1st and the 2nd step is  ' + str(maxDiff1smoothing) + '\n')
 
                
                # Convergence test
                print '+++++++++++++++++++++++++  Convergence test after step II. n_iter2 = ', n_iter2, ' time_step2 = ', time_step2
                f.write('+++++++++++++++++++++++++  Convergence test after step II. n_iter2 = ' + str(n_iter2) + ' time_step2 = ' + str(time_step2) + '\n')
                diff = aimsalgo.MaskedDiffusionSmoother_FLOAT(time_step2)
                diff.setMask(aimsmask, 32767)
                heatmap2 = diff.doSmoothing(heatmap, 10, True)   # corrected!!! Convergence - test it with step = 10
                aims.write(heatmap2, result_directory + 'heat_I_%s_%s_II_%s_%s_Converg_%s.nii.gz' %(n_iter, time_step, n_iter2, time_step2, keyWord))
                volHeat2 = aims.read(result_directory + 'heat_I_%s_%s_II_%s_%s_Converg_%s.nii.gz' %(n_iter, time_step, n_iter2, time_step2, keyWord))
                arr2 = np.array(volHeat2, copy = False)
                maxDiff = np.max(np.abs(arr1 - arr2))
                print 'Max difference to the volume to test convergence is  ', str(maxDiff)
                f.write('Max difference to the volume to test convergence is  ' + str(maxDiff) + '\n')
                aims.write(volHeat - volHeat2, result_directory + 'heat_I_%s_%s_II_%s_%s_ConvergDiff_%s.nii.gz' %(n_iter, time_step, n_iter2, time_step2, keyWord))
                t2Step2 = timeit.default_timer()
                f.write('Time for the Step II-Convergence is ' + str(np.round(t2Step2 - t1Step2)) + '\n')


    
        
        
        
       

#n_iter = 500 # 500
#time_step = 0.01
#heatmap_before = aims.read(result_directory + 'heat_%s.nii.gz' %(keyWord), 1)
#mask = aims.read(result_directory + 'all_but_cortex_%s.nii' %(keyWord), 1)
#aimsmask = aims.AimsData(mask)  # Important for reference-counting!
#diff = aimsalgo.MaskedDiffusionSmoother_FLOAT(time_step)
#diff.setMask(aimsmask, 32767)
#heatmap = diff.doSmoothing(heatmap_before, n_iter, True)
#aims.write(heatmap, result_directory + 'heat_%s.nii.gz' %(keyWord))


#n_iter = 100 # 100
#time_step = 0.001
##mask = aims.read(result_directory + 'all_but_cortex_%s.nii' %(keyWord), 1)
##aimsmask = aims.AimsData(mask)  # Important for reference-counting!
#diff = aimsalgo.MaskedDiffusionSmoother_FLOAT(time_step)
#diff.setMask(aimsmask, 32767)
#heatmap = diff.doSmoothing(heatmap, n_iter, True)
#aims.write(heatmap, result_directory + 'heat_%s.nii.gz' %(keyWord))
#volHeat = aims.read(result_directory + 'heat_%s.nii.gz' %(keyWord))
#arr1 = np.array(volHeat, copy = False)

## does this calculation converge?
#n_iter = 100   # 100
#time_step = 0.001
##mask = aims.read(result_directory + 'all_but_cortex_%s.nii' %(keyWord), 1)
##aimsmask = aims.AimsData(mask)  # Important for reference-counting!
#diff = aimsalgo.MaskedDiffusionSmoother_FLOAT(time_step)
#diff.setMask(aimsmask, 32767)
#heatmap2 = diff.doSmoothing(heatmap, n_iter, True)
#aims.write(heatmap2, result_directory + 'heatTestConvergence_%s.nii.gz' %(keyWord))
#arr2 = np.array(aims.read(result_directory + 'heatTestConvergence_%s.nii.gz' %(keyWord)), copy = False)
#print 'Max difference to the volume to test convergence is  ', str(np.max(np.abs(arr1 - arr2)))

## Normalized gradient's divergence
##python div_gradn.py
#vol_divGrad = highres_cortex.div_gradn.divergence_gradient(volHeat)
#aims.write(vol_divGrad, result_directory + "heat_div_gradn_%s.nii.gz" %(keyWord))

f.close()


