"""Test code for generating PEST files

"""

import io

import pytest

import numpy as np

import spowtd.pestfiles as pestfiles_mod
import spowtd.recession as recession_mod
import spowtd.rise as rise_mod
from spowtd.test import conftest


reference_text = {
    (
        'rise',
        'tpl',
        'peatclsm',
    ): """ptf @
specific_yield:
  type: peatclsm
  sd: @sd                      @
  theta_s: @theta_s                 @
  b: @b                       @
  psi_s: @psi_s                   @
transmissivity:
  type: peatclsm
  Ksmacz0: 7.3  # m/s
  alpha: 3  # dimensionless
  zeta_max_cm: 1.0""",
    (
        'rise',
        'tpl',
        'spline',
    ): """ptf @
specific_yield:
  type: spline
  zeta_knots_mm:
    - -291.7
    - -183.1
    - -15.74
    - 10.65
    - 38.78
    - 168.3
  sy_knots:  # Specific yield, dimensionless
    - @sy_knot_1               @
    - @sy_knot_2               @
    - @sy_knot_3               @
    - @sy_knot_4               @
    - @sy_knot_5               @
    - @sy_knot_6               @
transmissivity:
  type: spline
  zeta_knots_mm:
    - -291.7
    - -5.167
    - 168.3
    - 1000
  K_knots_km_d:  # Conductivity, km /d
    - 0.005356
    - 1.002
    - 6577.0
    - 8430.0
  minimum_transmissivity_m2_d: 7.442  # Minimum transmissivity, m2 /d""",
    (
        'curves',
        'tpl',
        'peatclsm',
    ): """ptf @
specific_yield:
  type: peatclsm
  sd: @sd                      @
  theta_s: @theta_s                 @
  b: @b                       @
  psi_s: @psi_s                   @
transmissivity:
  type: peatclsm
  Ksmacz0: @Ksmacz0                 @  # m/s
  alpha: @alpha                   @  # dimensionless
  zeta_max_cm: 1.0""",
    (
        'curves',
        'tpl',
        'spline',
    ): """ptf @
specific_yield:
  type: spline
  zeta_knots_mm:
    - -291.7
    - -183.1
    - -15.74
    - 10.65
    - 38.78
    - 168.3
  sy_knots:  # Specific yield, dimensionless
    - @sy_knot_1               @
    - @sy_knot_2               @
    - @sy_knot_3               @
    - @sy_knot_4               @
    - @sy_knot_5               @
    - @sy_knot_6               @
transmissivity:
  type: spline
  zeta_knots_mm:
    - -291.7
    - -5.167
    - 168.3
    - 1000
  K_knots_km_d:  # Conductivity, km /d
    - @K_knot_1                @
    - @K_knot_2                @
    - @K_knot_3                @
    - @K_knot_4                @
  minimum_transmissivity_m2_d: @T_min                   @  # Minimum transmissivity, m2 /d""",
}


# A way of identifying the sample:
#   mapping from number of distinct zetas to sample number
SAMPLE_NUMBER = {141: 1, 201: 2}


@pytest.mark.parametrize(
    ('parameterization', 'reference_file_contents'),
    [
        ('peatclsm', reference_text[('rise', 'tpl', 'peatclsm')]),
        ('spline', reference_text[('rise', 'tpl', 'spline')]),
    ],
)
def test_generate_rise_tpl(
    classified_connection, parameterization, reference_file_contents
):
    """Generation of template files for rise curve calibration"""
    outfile = io.StringIO()
    with open(
        conftest.get_parameter_file_path(parameterization), 'rt'
    ) as parameter_file:
        pestfiles_mod.generate_rise_pestfiles(
            classified_connection,
            parameter_file=parameter_file,
            outfile_type='tpl',
            configuration_file=None,
            outfile=outfile,
        )
    assert (
        outfile.getvalue().splitlines() == reference_file_contents.splitlines()
    )


def test_generate_rise_ins(classified_connection):
    """Generation of instruction files for rise curve calibration"""
    # XXX Use a fixture
    rise_mod.find_rise_offsets(classified_connection)
    cursor = classified_connection.cursor()
    cursor.execute(
        """
    SELECT count(distinct zeta_number)
    FROM rising_interval_zeta"""
    )
    n_zeta = cursor.fetchone()[0]
    cursor.close()
    assert n_zeta in SAMPLE_NUMBER
    sample = SAMPLE_NUMBER[n_zeta]
    outfile = io.StringIO()
    with open(
        # Parameterization doesn't matter for ins file
        conftest.get_parameter_file_path('peatclsm'),
        'rt',
    ) as parameter_file:
        pestfiles_mod.generate_rise_pestfiles(
            classified_connection,
            parameter_file=parameter_file,
            outfile_type='ins',
            configuration_file=None,
            outfile=outfile,
        )
    with open(
        conftest.get_sample_file_path('rise_calibration', sample, 'ins'), 'rt'
    ) as ref_file:
        assert outfile.getvalue().splitlines() == ref_file.read().splitlines()


@pytest.mark.parametrize('parameterization', ['peatclsm', 'spline'])
def test_generate_rise_pst(classified_connection, parameterization):
    """Generation of control files for rise curve calibration"""
    # XXX Use a fixture
    rise_mod.find_rise_offsets(classified_connection)
    cursor = classified_connection.cursor()
    cursor.execute(
        """
    SELECT count(distinct zeta_number)
    FROM rising_interval_zeta"""
    )
    n_zeta = cursor.fetchone()[0]
    cursor.close()
    sample = SAMPLE_NUMBER[n_zeta]
    outfile = io.StringIO()
    with open(
        conftest.get_parameter_file_path(parameterization), 'rt'
    ) as parameter_file:
        pestfiles_mod.generate_rise_pestfiles(
            classified_connection,
            parameter_file=parameter_file,
            outfile_type='pst',
            configuration_file=None,
            outfile=outfile,
        )
    with open(
        conftest.get_sample_file_path(
            '{}_rise_calibration'.format(parameterization), sample, 'pst'
        ),
        'rt',
    ) as ref_file:
        out_lines = outfile.getvalue().splitlines()
        ref_lines = ref_file.read().splitlines()
        obs_start = ref_lines.index('* observation data')
        assert ref_lines[obs_start] == '* observation data'
        assert out_lines[: obs_start + 1] == ref_lines[: obs_start + 1]
        obs_end = ref_lines.index('* model command line')
        assert ref_lines[obs_end] == '* model command line'
        assert out_lines[obs_end:] == ref_lines[obs_end:]
        (ref_obs_values, ref_obs_lines) = extract_observations(ref_lines)
        (out_obs_values, out_obs_lines) = extract_observations(out_lines)
        assert ref_obs_lines == out_obs_lines
        assert np.allclose(ref_obs_values, out_obs_values)


@pytest.mark.parametrize(
    ('parameterization', 'reference_file_contents'),
    [
        ('peatclsm', reference_text[('curves', 'tpl', 'peatclsm')]),
        ('spline', reference_text[('curves', 'tpl', 'spline')]),
    ],
)
def test_generate_curves_tpl(
    classified_connection, parameterization, reference_file_contents
):
    """Generation of template files for master curves calibration"""
    outfile = io.StringIO()
    with open(
        conftest.get_parameter_file_path(parameterization), 'rt'
    ) as parameter_file:
        pestfiles_mod.generate_curves_pestfiles(
            classified_connection,
            parameter_file=parameter_file,
            outfile_type='tpl',
            configuration_file=None,
            outfile=outfile,
        )
    assert (
        outfile.getvalue().splitlines() == reference_file_contents.splitlines()
    )


def test_generate_curves_ins(classified_connection):
    """Generation of instruction files for master curves calibration"""
    # XXX Use a fixture
    rise_mod.find_rise_offsets(classified_connection)
    recession_mod.find_recession_offsets(classified_connection)
    cursor = classified_connection.cursor()
    cursor.execute(
        """
    SELECT count(distinct zeta_number)
    FROM rising_interval_zeta"""
    )
    n_zeta = cursor.fetchone()[0]
    cursor.close()
    assert n_zeta in SAMPLE_NUMBER
    sample = SAMPLE_NUMBER[n_zeta]
    outfile = io.StringIO()
    with open(
        # Parameterization doesn't matter for ins file
        conftest.get_parameter_file_path('peatclsm'),
        'rt',
    ) as parameter_file:
        pestfiles_mod.generate_curves_pestfiles(
            classified_connection,
            parameter_file=parameter_file,
            outfile_type='ins',
            configuration_file=None,
            outfile=outfile,
        )
    with open(
        conftest.get_sample_file_path('curves_calibration', sample, 'ins'),
        'rt',
    ) as ref_file:
        assert outfile.getvalue().splitlines() == ref_file.read().splitlines()


@pytest.mark.parametrize('parameterization', ['peatclsm', 'spline'])
def test_generate_curves_pst(classified_connection, parameterization):
    """Generation of control files for curves curve calibration"""
    # XXX Use a fixture
    rise_mod.find_rise_offsets(classified_connection)
    recession_mod.find_recession_offsets(classified_connection)
    cursor = classified_connection.cursor()
    cursor.execute(
        """
    SELECT count(distinct zeta_number)
    FROM rising_interval_zeta"""
    )
    n_rise_zeta = cursor.fetchone()[0]
    cursor.execute(
        """
    SELECT count(distinct zeta_number)
    FROM recession_interval_zeta"""
    )
    n_recession_zeta = cursor.fetchone()[0]
    cursor.close()
    sample = SAMPLE_NUMBER[n_rise_zeta]
    outfile = io.StringIO()
    with open(
        conftest.get_parameter_file_path(parameterization), 'rt'
    ) as parameter_file:
        pestfiles_mod.generate_curves_pestfiles(
            classified_connection,
            parameter_file=parameter_file,
            outfile_type='pst',
            configuration_file=None,
            outfile=outfile,
        )
    with open(
        conftest.get_sample_file_path(
            '{}_curves_calibration'.format(parameterization), sample, 'pst'
        ),
        'rt',
    ) as ref_file:
        out_lines = outfile.getvalue().splitlines()
        ref_lines = ref_file.read().splitlines()
        obs_start = ref_lines.index('* observation data')
        assert ref_lines[obs_start] == '* observation data'
        assert out_lines[: obs_start + 1] == ref_lines[: obs_start + 1]
        obs_end = ref_lines.index('* model command line')
        assert ref_lines[obs_end] == '* model command line'
        assert out_lines[obs_end:] == ref_lines[obs_end:]
        (ref_obs_values, ref_obs_lines) = extract_observations(ref_lines)
        (out_obs_values, out_obs_lines) = extract_observations(out_lines)
        assert ref_obs_lines == out_obs_lines
        assert np.allclose(ref_obs_values, out_obs_values)


def extract_observations(pestfile_lines):
    """Extract observations from a PEST control file

    Splits the observation section of a PEST control file and returns
    a Numpy array of floating-point observation values and a list of
    lists of whitespace-separated items from each row in the
    observation section (excluding the observation values), in that
    order.

    """
    obs_start = pestfile_lines.index('* observation data')
    obs_end = pestfile_lines.index('* model command line')
    observation_rows = [
        line.strip().split()
        for line in pestfile_lines[obs_start + 1 : obs_end]
    ]
    return (
        np.array([float(row[1]) for row in observation_rows], dtype='float64'),
        [row[:1] + row[2:] for row in observation_rows],
    )
