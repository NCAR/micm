"""
Copyright (C) 2023
National Center for Atmospheric Research,
SPDX-License-Identifier: Apache-2.0

File:
    test_kpp_to_micm.py

Usage:
    pytest test_kpp_to_micm.py --log-cli-level=DEBUG
"""

import kpp_to_micm

def test_parse_kpp_arrhenius():
    """
    examples
    ARR_ab(1.0e-12, 2000.0)
    ARR_ac(1.0e-12, -3.0)
    ARR_abc(1.0e-12, 2000.0, -3.0)
    """

    kpp_A, kpp_B, kpp_C = 1.0e-12, 2000.0, -3.0

    arr_dict = kpp_to_micm.parse_kpp_arrhenius(
        'ARR_ab(%.2e, %.2f)' % (kpp_A, kpp_B))
    assert arr_dict['A'] == kpp_A
    assert arr_dict['B'] == kpp_B

    arr_dict = kpp_to_micm.parse_kpp_arrhenius(
        'ARR_ac(%.2e, %.2f)' % (kpp_A, kpp_C))
    assert arr_dict['A'] == kpp_A
    assert arr_dict['C'] == kpp_C

    arr_dict = kpp_to_micm.parse_kpp_arrhenius(
        'ARR_abc(%.2e, %.2f, %.2f)' % (kpp_A, kpp_B, kpp_C))
    assert arr_dict['A'] == kpp_A
    assert arr_dict['B'] == kpp_B
    assert arr_dict['C'] == kpp_C
