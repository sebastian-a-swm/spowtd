"""Test code for specific yield classes

"""

import numpy as np

import pytest

import yaml

import spowtd.plot_specific_yield as plot_mod
import spowtd.specific_yield as specific_yield_mod
from spowtd.test import conftest


@pytest.mark.parametrize(
    'sy_type, expected_sy',
    [
        ('peatclsm', None),  # Expected values retrieved from R script
        (
            'spline',
            [
                0.1358,
                0.1358,
                0.1358,
                0.1358,
                0.1358,
                0.13794492,
                0.16751247,
                0.21664407,
                0.26417862,
                0.6857,
            ],
        ),
    ],
)
def test_specific_yield(sy_type, expected_sy):
    """Test specific yield functions"""
    with open(conftest.get_parameter_file_path(sy_type), 'rt') as sy_file:
        sy_parameters = yaml.safe_load(sy_file)['specific_yield']
    specific_yield = specific_yield_mod.create_specific_yield_function(
        sy_parameters
    )
    if sy_type == 'spline':
        zeta_mm = np.linspace(-0.9, 0.2, 10) * 1000
        assert np.allclose(specific_yield(zeta_mm), expected_sy)
    else:
        assert sy_type == 'peatclsm'
        # In this case, expected_sy is provided by an R script, and is
        # large (201 values), so we load it from a file produced by
        # the R script instead of including it in the pytest
        # decorator.
        del expected_sy
        sy_table = conftest.peatclsm_specific_yield_table()
        zeta_m = np.linspace(-0.995, 1.005, 201)
        zeta_m_ref, expected_sy = zip(*sy_table)
        assert np.allclose(zeta_m, zeta_m_ref)
        assert np.allclose(specific_yield(zeta_m * 1000), expected_sy)
