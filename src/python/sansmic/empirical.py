# coding: utf-8
#
# Copyright (c) 2024 National Technology and Engineering Solutions of
# Sandia, LLC. Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.
#
# SPDX-License-Identifier: BSD-3-Clause.

"""
The values in this module are determined empirically and taken from
published research. Among the most significant data included here
is the table of water properties by weight-percent NaCl, temperature,
and density. See :cite:t:`kaufman60sodium` and :cite:t:`lide2000crc`.
"""

from numpy import array

rho_NaCl_s = float(2.16)  # density of rock salt
"""
Density of salt in the surrounding rock.
"""

a = array(
    [
        7.6091661297e-1,
        -3.8715516531,
        7.8254117638,
        -7.8395924067,
        3.8789476581,
        -0.7533873462,
    ],
    dtype=float,
)
"""
Coefficients used for wall recession rate equation (SAND2015-XXXX Ch 4).
"""

# Coefficients used for source term see (SAND2015-XXXX Ch2.3)
c = array(
    [
        1.0019795678,
        0.67459890538,
        0.32350531044,
    ],
    dtype=float,
)
"""
Coefficients used for converting between weight percent and specific gravity (SAND2015 Ch2.3).
"""

sg_matrix = array(
    [
        [
            0.99987,
            0.99973,
            0.99823,
            0.99707,
            0.99567,
            0.99224,
            0.98807,
            0.98324,
            0.97183,
            0.95838,
        ],
        [
            1.00747,
            1.00707,
            1.00534,
            1.00409,
            1.00261,
            0.99908,
            0.99482,
            0.99,
            0.9785,
            0.9651,
        ],
        [
            1.01509,
            1.01442,
            1.01246,
            1.01112,
            1.00957,
            1.00593,
            1.00161,
            0.9967,
            0.9852,
            0.9719,
        ],
        [
            1.03038,
            1.0292,
            1.0268,
            1.0253,
            1.02361,
            1.01977,
            1.01531,
            1.0103,
            0.9988,
            0.9855,
        ],
        [
            1.04575,
            1.04408,
            1.04127,
            1.03963,
            1.03781,
            1.03378,
            1.02919,
            1.0241,
            1.0125,
            0.9994,
        ],
        [
            1.06121,
            1.05907,
            1.05589,
            1.05412,
            1.05219,
            1.04798,
            1.04326,
            1.0381,
            1.0264,
            1.0134,
        ],
        [
            1.07677,
            1.07419,
            1.07068,
            1.06879,
            1.06676,
            1.06238,
            1.05753,
            1.0523,
            1.0405,
            1.0276,
        ],
        [
            1.09244,
            1.08946,
            1.08566,
            1.08365,
            1.08153,
            1.07699,
            1.07202,
            1.0667,
            1.0549,
            1.042,
        ],
        [
            1.10824,
            1.10491,
            1.10085,
            1.09872,
            1.09651,
            1.09182,
            1.08674,
            1.0813,
            1.0694,
            1.0565,
        ],
        [
            1.12419,
            1.12056,
            1.11621,
            1.11401,
            1.11171,
            1.10688,
            1.1017,
            1.0962,
            1.0842,
            1.0713,
        ],
        [
            1.14031,
            1.13643,
            1.1319,
            1.12954,
            1.12715,
            1.12218,
            1.11691,
            1.1113,
            1.0993,
            1.0864,
        ],
        [
            1.15663,
            1.15254,
            1.14779,
            1.14533,
            1.14285,
            1.13774,
            1.13238,
            1.1268,
            1.1146,
            1.1017,
        ],
        [
            1.17318,
            1.16891,
            1.16395,
            1.1614,
            1.15883,
            1.15358,
            1.14812,
            1.1425,
            1.1303,
            1.1172,
        ],
        [
            1.18999,
            1.18557,
            1.1804,
            1.17776,
            1.17511,
            1.16971,
            1.16414,
            1.1584,
            1.1463,
            1.1331,
        ],
        [
            1.20709,
            1.20254,
            1.19717,
            1.19443,
            1.1917,
            1.18614,
            1.18045,
            1.1747,
            1.1626,
            1.1492,
        ],
    ],
    dtype=float,
)
"""
Specific gravity data for sodium chloride :cite:`kaufman60sodium,lide2000crc`.
"""

wt_pct_vec = array(
    [
        [0.0],
        [1.0],
        [2.0],
        [4.0],
        [6.0],
        [8.0],
        [10.0],
        [12.0],
        [14.0],
        [16.0],
        [18.0],
        [20.0],
        [22.0],
        [24.0],
        [26.0],
    ],
    dtype=float,
)
"""
Weight-percent salt content data for sodium chloride :cite:`kaufman60sodium,lide2000crc`.
"""

T_degC_vec = array(
    [[0, 10, 20, 25, 30, 40, 50, 60, 80, 100]],
    dtype=float,
)
"""
Temperature (in degrees Celsius) data for sg vs wt-pct for sodium chloride :cite:`kaufman60sodium,lide2000crc`.
"""

wt_pct_max = array(
    [26.34, 26.35, 26.43, 26.48, 26.56, 26.71, 26.89, 27.09, 27.53, 28.12],
    dtype=float,
)
"""
Weight-percent salt content data for sodium chloride :cite:`kaufman60sodium,lide2000crc`.
"""

T_degC_max = array(
    [0, 10, 20, 25, 30, 40, 50, 60, 80, 100],
    dtype=float,
)
"""
Temperature (in degrees Celsius) data for sg vs wt pct for sodium chloride :cite:`kaufman60sodium,lide2000crc`.
"""

sg_max = array(
    [1.2093, 1.2044, 1.1999, 1.1978, 1.1957, 1.1914, 1.1872, 1.183, 1.1745, 1.166],
    dtype=float,
)
"""
Specific gravity data for sodium chloride :cite:`kaufman60sodium,lide2000crc`.
"""
