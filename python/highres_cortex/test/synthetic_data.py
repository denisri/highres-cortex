# -*- coding: utf-8 -*-
# Copyright Télécom ParisTech (2015).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
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

import os.path
import math

from soma import aims
import numpy

def _convert_to_float_triple(param):
    try:
        value = float(param)
        return (value, value, value)
    except TypeError:
        pass

    assert len(param) == 3
    return tuple(float(item) for item in param)

def make_cortex_sphere_classif(inner_radius, outer_radius,
                               voxel_size, margin=None):
    inner_radius = _convert_to_float_triple(inner_radius)
    outer_radius = _convert_to_float_triple(outer_radius)
    voxel_size = _convert_to_float_triple(voxel_size)
    if margin is None:
        margin = max(voxel_size)
    margin = _convert_to_float_triple(margin)

    for i in range(3):
        assert outer_radius[i] > inner_radius[i] > 0

    size = [int(math.ceil(2 * (ax_radius + ax_margin) / ax_voxel_size))
            for ax_radius, ax_margin, ax_voxel_size
            in zip(outer_radius, margin, voxel_size)]

    classif_volume = aims.Volume(size[0], size[1], size[2], dtype="S16")
    classif_volume.header()["voxel_size"] = list(voxel_size)
    emplace_cortex_sphere_classif(classif_volume,
                                  inner_radius, outer_radius, margin)
    return classif_volume

def make_centred_coord_grids(classif_volume):
    size = (classif_volume.getSizeX(),
            classif_volume.getSizeY(),
            classif_volume.getSizeZ())
    voxel_size = classif_volume.header()["voxel_size"][:3]

    s = [numpy.linspace(-ax_voxel_size * (ax_size // 2),
                        ax_voxel_size * (ax_size // 2),
                        ax_size)
         for ax_voxel_size, ax_size in zip(voxel_size, size)]
    return numpy.ix_(*s)


def emplace_cortex_sphere_classif(classif_volume,
                                  inner_radius, outer_radius, margin):
    inner_radius = _convert_to_float_triple(inner_radius)
    outer_radius = _convert_to_float_triple(outer_radius)
    margin = _convert_to_float_triple(margin)

    size = (classif_volume.getSizeX(),
            classif_volume.getSizeY(),
            classif_volume.getSizeZ())
    voxel_size = classif_volume.header()["voxel_size"][:3]
    for i in range(3):
        assert outer_radius[i] > inner_radius[i] > 0
        assert size[i] * voxel_size[i] >= 2 * outer_radius[i]

    grids = make_centred_coord_grids(classif_volume)

    np_classif = numpy.asarray(classif_volume)
    classif_volume.fill(0)
    np_classif[sum(grid ** 2 / radius ** 2
                   for radius, grid in zip(outer_radius, grids)) < 1] = 100
    np_classif[sum(grid ** 2 / radius ** 2
                   for radius, grid in zip(inner_radius, grids)) < 1] = 200

def _make_similar_volume(data_array, ref):
    volume = aims.Volume(data_array)
    volume.header().update(ref.header())
    return volume

def make_sphere_and_reference_result(inner_radius, outer_radius, voxel_size):
    inner_radius = float(inner_radius)
    outer_radius = float(outer_radius)
    assert outer_radius > inner_radius > 0
    voxel_size = _convert_to_float_triple(voxel_size)

    thickness = outer_radius - inner_radius

    classif_volume = make_cortex_sphere_classif(inner_radius, outer_radius,
                                                voxel_size)
    grids = make_centred_coord_grids(classif_volume)
    distance_to_centre = numpy.sqrt(sum(grid ** 2 for grid in grids))

    distance_to_white = distance_to_centre - inner_radius
    distance_to_CSF = outer_radius - distance_to_centre

    euclidean_metric = numpy.clip((outer_radius - distance_to_centre)
                                  / (outer_radius - inner_radius),
                                  0, 1)

    with numpy.errstate(divide="ignore"):
        laplacian_value = numpy.clip(
            inner_radius / (outer_radius - inner_radius) *
            (outer_radius / distance_to_centre - 1),
            0, 1)

    equivolumic_metric = numpy.clip(
        (outer_radius ** 3 - distance_to_centre ** 3) /
        (outer_radius ** 3 - inner_radius ** 3),
        0, 1)

    return (classif_volume,
            _make_similar_volume(distance_to_white, ref=classif_volume),
            _make_similar_volume(distance_to_CSF, ref=classif_volume),
            _make_similar_volume(euclidean_metric, ref=classif_volume),
            _make_similar_volume(laplacian_value, ref=classif_volume),
            _make_similar_volume(equivolumic_metric, ref=classif_volume))

def write_sphere_and_reference_result(inner_radius, outer_radius, voxel_size,
                                      dir="."):
    inner_radius = float(inner_radius)
    outer_radius = float(outer_radius)
    assert outer_radius > inner_radius > 0
    voxel_size = _convert_to_float_triple(voxel_size)

    (classif,
     distance_to_white, distance_to_CSF,
     euclidean_metric,
     laplacian_value,
     equivolumic_metric) = (make_sphere_and_reference_result(
                                inner_radius, outer_radius, voxel_size))

    aims.write(classif,
                os.path.join(dir, "classif.nii.gz"))
    aims.write(distance_to_white,
               os.path.join(dir, "reference_distwhite.nii.gz"))
    aims.write(distance_to_CSF,
               os.path.join(dir, "reference_distCSF.nii.gz"))
    aims.write(euclidean_metric,
               os.path.join(dir, "reference_euclidean.nii.gz"))
    aims.write(laplacian_value,
               os.path.join(dir, "reference_laplacian.nii.gz"))
    aims.write(equivolumic_metric,
               os.path.join(dir, "reference_equivolumic.nii.gz"))

if __name__ == "__main__":
    import os
    import shutil
    os.mkdir("sphere_3_6")
    write_sphere_and_reference_result(3, 6, 0.3, dir="sphere_3_6")
    os.mkdir("sphere_2_5")
    write_sphere_and_reference_result(2, 5, 0.3, dir="sphere_2_5")
    os.mkdir("sphere_1_4")
    write_sphere_and_reference_result(1, 4, 0.3, dir="sphere_1_4")
    os.mkdir("sphere_5_8")
    write_sphere_and_reference_result(5, 8, 0.3, dir="sphere_5_8")
    os.mkdir("sphere_10_13")
    write_sphere_and_reference_result(10, 13, 0.3, dir="sphere_10_13")
