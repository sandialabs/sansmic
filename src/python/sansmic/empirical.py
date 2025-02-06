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
        7.609_166_129_7e-1,
        -3.871_551_653_1,
        7.825_411_763_8,
        -7.839_592_406_7,
        3.878_947_658_1,
        -0.753_387_346_2,
    ],
    dtype=float,
)
"""
Coefficients used for wall recession rate equation (SAND2015-XXXX Ch 4).
"""

# Coefficients used for source term see (SAND2015-XXXX Ch2.3)
c = array(
    [
        1.001_979_567_8,
        0.674_598_905_38,
        0.323_505_310_44,
    ],
    dtype=float,
)
"""
Coefficients used for converting between weight percent and specific gravity (SAND2015 Ch2.3).
"""

sg_matrix = array(
    [
        [
            0.999_87,
            0.999_73,
            0.998_23,
            0.997_07,
            0.995_67,
            0.992_24,
            0.988_07,
            0.983_24,
            0.971_83,
            0.958_38,
        ],
        [
            1.007_47,
            1.007_07,
            1.005_34,
            1.004_09,
            1.002_61,
            0.999_08,
            0.994_82,
            0.99,
            0.978_5,
            0.965_1,
        ],
        [
            1.015_09,
            1.014_42,
            1.012_46,
            1.011_12,
            1.009_57,
            1.005_93,
            1.001_61,
            0.996_7,
            0.985_2,
            0.971_9,
        ],
        [
            1.030_38,
            1.029_2,
            1.026_8,
            1.025_3,
            1.023_61,
            1.019_77,
            1.015_31,
            1.010_3,
            0.998_8,
            0.985_5,
        ],
        [
            1.045_75,
            1.044_08,
            1.041_27,
            1.039_63,
            1.037_81,
            1.033_78,
            1.029_19,
            1.024_1,
            1.012_5,
            0.999_4,
        ],
        [
            1.061_21,
            1.059_07,
            1.055_89,
            1.054_12,
            1.052_19,
            1.047_98,
            1.043_26,
            1.038_1,
            1.026_4,
            1.013_4,
        ],
        [
            1.076_77,
            1.074_19,
            1.070_68,
            1.068_79,
            1.066_76,
            1.062_38,
            1.057_53,
            1.052_3,
            1.040_5,
            1.027_6,
        ],
        [
            1.092_44,
            1.089_46,
            1.085_66,
            1.083_65,
            1.081_53,
            1.076_99,
            1.072_02,
            1.066_7,
            1.054_9,
            1.042,
        ],
        [
            1.108_24,
            1.104_91,
            1.100_85,
            1.098_72,
            1.096_51,
            1.091_82,
            1.086_74,
            1.081_3,
            1.069_4,
            1.056_5,
        ],
        [
            1.124_19,
            1.120_56,
            1.116_21,
            1.114_01,
            1.111_71,
            1.106_88,
            1.101_7,
            1.096_2,
            1.084_2,
            1.071_3,
        ],
        [
            1.140_31,
            1.136_43,
            1.131_9,
            1.129_54,
            1.127_15,
            1.122_18,
            1.116_91,
            1.111_3,
            1.099_3,
            1.086_4,
        ],
        [
            1.156_63,
            1.152_54,
            1.147_79,
            1.145_33,
            1.142_85,
            1.137_74,
            1.132_38,
            1.126_8,
            1.114_6,
            1.101_7,
        ],
        [
            1.173_18,
            1.168_91,
            1.163_95,
            1.161_4,
            1.158_83,
            1.153_58,
            1.148_12,
            1.142_5,
            1.130_3,
            1.117_2,
        ],
        [
            1.189_99,
            1.185_57,
            1.180_4,
            1.177_76,
            1.175_11,
            1.169_71,
            1.164_14,
            1.158_4,
            1.146_3,
            1.133_1,
        ],
        [
            1.207_09,
            1.202_54,
            1.197_17,
            1.194_43,
            1.191_7,
            1.186_14,
            1.180_45,
            1.174_7,
            1.162_6,
            1.149_2,
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
