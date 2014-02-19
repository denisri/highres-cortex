#!/bin/sh -e

AimsThreshold -b -m di -t 100 \
    -i ../classif.nii.gz \
    -o ./all_but_cortex.nii
AimsFileConvert -t FLOAT \
    -i ../classif.nii.gz \
    -o heat.nii.gz

# Each run refines the previous one
python heat.py 500 0.01
python heat.py 100 0.001

# Normalized gradient
ylNumpyComb.py -f 'sqrt(I1**2 + I2**2 + I3**2)' \
    -i heat_gradx.nii.gz \
    -i heat_grady.nii.gz \
    -i heat_gradz.nii.gz \
    -o heat_grad_norm.nii
cartoLinearComb.py -f 'I1/I2' \
    -i heat_gradx.nii.gz \
    -i heat_grad_norm.nii \
    -o heat_gradnx.nii
cartoLinearComb.py -f 'I1/I2' \
    -i heat_grady.nii.gz \
    -i heat_grad_norm.nii \
    -o heat_gradny.nii
cartoLinearComb.py -f 'I1/I2' \
    -i heat_gradz.nii.gz \
    -i heat_grad_norm.nii \
    -o heat_gradnz.nii

python div_gradn.py
