"""Determine time offsets that align pieces of a hydrograph

"""

import logging

import numpy.linalg as linalg_mod
import numpy as np

import spowtd.regrid as regrid_mod


LOG = logging.getLogger('spowtd.fit_offsets')


def get_series_time_offsets(series_list, head_step):
    """Find a time offset that minimizes difference in head crossing times

    This function is used in assembly of recession curves.  See further
    documentation under get_series_offsets.

    """
    head_mapping, index_mapping = build_connected_head_mapping(
        series_list, head_step
    )
    del series_list, head_step
    return get_series_offsets(head_mapping, index_mapping)


def get_series_offsets(head_mapping, index_mapping, covariance=None, weight=None):
    """Find offsets that minimizes difference in head crossing times

    Given a sequence of (x, head) data series, rediscretize to get an
    average crossing x of each integer multiple of head_step for each
    series.  Then find an x offset for each data series that minimizes the
    sum of squared differences in times when all data series cross those head
    values.  This offset *replaces* the existing offset of the series,
    x_new = x - min(x) + time_offset

    The series in the list with the largest initial head is treated as the
    reference series, and assigned an offset of zero.

    Returns (indices, offsets, head_mapping) connecting indices of series in
    series_list to their corresponding offset.  In general, these will be a
    subset of series_list if there were some heads at which there was no
    overlap.  The head mapping is a mapping between discrete head ids and
    sequences of (series_id, x_mean) pairs, representing the mean x at which
    that series crossed that head.  Head ids are integers which, when
    multiplied by head_step, will give once more a head with units.

    This function is used in assembly of both recession and rise curves.

    """
    series_ids, offsets = find_offsets(head_mapping, covariance)
    assert len(series_ids) == len(offsets)
    original_indices = [index_mapping[series_id] for series_id in series_ids]
    # We need to map series ids in head_mapping back to their original
    #   indices
    output_mapping = {}
    for head_id, crossings in list(head_mapping.items()):
        output_mapping[head_id] = [
            (index_mapping[series_id], t_mean)
            for series_id, t_mean in crossings
        ]
    return (original_indices, offsets, output_mapping, weight)


def build_connected_head_mapping(series_list, head_step):
    """Construct a mapping between head ids and series crossing times

    Calls build_exhaustive_head_mapping to create an initial mapping between
    head ids and series crossing times, and then identifies and retains just
    the largest connected component of overlapping events.  Heads with only one
    series are also removed as these are uninformative.  Series are identified
    by an id that sorts the series by initial head.

    The output is two dicts:

    head_mapping: maps head_id to a sequence of (series_id, time) pairs
                  indicating when each series crossed that head.
    index_mapping: maps series_id to the index of the series in series_list.

    """
    if not series_list:
        raise ValueError('empty series list')
    # We need to retain the indices in series_list so that the caller knows
    # which offsets go with which series, but we also need to sort by initial
    # head; so, retain a mapping from the series_id we use for finding offsets
    # to index in sorted_list
    dec = sorted(
        (
            (t - t.min(), H, index)
            for (index, (t, H)) in enumerate(series_list)
        ),
        key=lambda t_H_index: t_H_index[1][0],
    )
    del series_list
    sorted_list = []
    index_mapping = {}
    for new_index, (t, H, original_index) in enumerate(dec):
        sorted_list.append((t, H))
        index_mapping[new_index] = original_index

    del dec
    head_mapping = build_exhaustive_head_mapping(sorted_list, head_step)
    del sorted_list
    series_at_head = dict(
        (head, set(series_id for (series_id, t_mean) in value))
        for (head, value) in list(head_mapping.items())
    )
    connected_components = get_connected_components(series_at_head)
    del series_at_head
    if len(connected_components) > 1:
        LOG.info(
            '%s connected sets of head of sizes '
            '%s; will keep only largest component.',
            len(connected_components),
            tuple(len(cc) for cc in connected_components),
        )
    head_mappings = split_mapping_by_keys(
        head_mapping, connected_components[:1]
    )
    del connected_components
    assert len(head_mappings) == 1
    head_mapping = head_mappings[0]
    del head_mappings
    # Eliminate all heads with only one series, these are uninformative
    for head_id, seq in list(head_mapping.items()):
        # Don't use "assert seq" here, this is an ndarray
        assert len(seq) > 0  # pylint: disable=len-as-condition
        if len(seq) == 1:
            del head_mapping[head_id]
    return head_mapping, index_mapping


def build_exhaustive_head_mapping(series, head_step=1):
    """Construct a mapping between head ids and series crossing times

    series is a sequence of (time, head) data series.  Each series is
    regridded via interpolation to instead give times at which the head time
    series crosses an integer multiple of head_step.  That integer multiple
    (head id) is the key to the returned mapping, which maps between head_id
    and a sequence of (series_id, time) pairs indicating when each series
    crossed that head.  The series id is just the index of the series passed
    in.

    """
    head_mapping = {}
    for series_id, (t, H) in enumerate(series):
        # take averages of t at distinct H
        all_times = {}
        for head_id, time in regrid_mod.regrid(t, H, head_step):
            all_times.setdefault(head_id, []).append(time)
        for head_id, time in list(all_times.items()):
            t_mean = np.mean(time)
            head_mapping.setdefault(head_id, []).append((series_id, t_mean))
    return head_mapping


def find_offsets(head_mapping, covariance=None):
    """Find the time offsets that align the series in head_mapping

    Finds the set of time offsets that minimize the sum of squared differences
    in times at which each series crosses a particular head.  Input is a
    mapping of head id (a hashable value corresponding to a head, normally an
    integer) to a sequence of (series_id, time) pairs wherein series_id is an
    identifier for a sequence and time is the time at which the series crossed
    the corresponding head value.

    The series with the series_id that is largest (last in sort order) is
    treated as the reference and given an offset of zero; all other offsets
    are relative to that one.

    Returns series_ids, offsets where series_ids are the identifiers

    """
    # Assemble mapping of series ids to row numbers for the offset-finding
    #   problem
    series_ids = sorted(
        set().union(
            *(
                (series_id for series_id, _ in seq)
                for seq in list(head_mapping.values())
            )
        )
    )
    series_indices = dict(zip(series_ids, range(len(series_ids))))
    A, b = assemble_linear_system(
        head_mapping,
        series_indices,
    )
    if covariance is not None:
        assert covariance.shape == (len(b), len(b))
        offsets = gls_solve(A, b, covariance)
    else:
        offsets = ols_solve(A, b)
    assert offsets.shape == (A.shape[1],)
    # This was the boundary condition, zero offset for reference (last) id
    offsets = np.concatenate((offsets, [0]))
    # Offsets are by index, but reverse mapping is trivial because series ids
    #   are sorted
    assert len(series_ids) == len(offsets), '{} != {}'.format(
        len(series_ids), len(offsets)
    )
    return (series_ids, offsets)


def assemble_linear_system(
    head_mapping,
    series_indices,
):
    """Assemble linear system for fitting offsets

    Creates and populates matrix A and vector b representing the overdetermined
    system Ax = b that will be solved.

    """
    number_of_equations = sum(
        len(series_at_head) for series_at_head in list(head_mapping.values())
    )
    number_of_unknowns = len(series_indices) - 1
    LOG.info(
        '%s equations, %s unknowns', number_of_equations, number_of_unknowns
    )
    # Reference series corresponds to the highest series id; it has the
    #   largest initial head, because we sorted them
    reference_index = max(series_indices)
    LOG.info('Reference index: %s', reference_index)
    A = np.zeros((number_of_equations, number_of_unknowns))
    b = np.zeros((number_of_equations,))
    row_template = np.zeros((number_of_unknowns,))
    row_index = 0
    for head_id, series_at_head in sorted(head_mapping.items()):
        row_template[:] = 0
        sids, times = list(zip(*series_at_head))
        number_of_series_at_head = len(sids)
        indices = [
            series_indices[index] for index in sids if index != reference_index
        ]
        row_template[indices] = 1.0 / number_of_series_at_head
        mean_time = np.mean(times)
        for series_id, t in series_at_head:
            A[row_index] = row_template
            # !!! some redundancy here
            if series_id != reference_index:
                series_index = series_indices[series_id]
                A[row_index, series_index] -= 1
            b[row_index] = t - mean_time
            row_index += 1

    assert row_index == number_of_equations, row_index
    return A, b


def ols_solve(A, b):
    """Solve Ax = b by ordinary least squares

    Returns x, the solution to the linear system

    """
    number_of_unknowns = A.shape[1]
    ATA = np.dot(A.transpose(), A)
    assert ATA.shape == (number_of_unknowns, number_of_unknowns), ATA.shape
    ATb = np.dot(A.transpose(), b)
    return linalg_mod.solve(ATA, ATb)  # pylint: disable=E1101


def gls_solve(A, b, O):
    """Solve Ax = b with covariance matrix O by generalized least squares

    Returns x, the solution to the linear system

    """
    number_of_unknowns = A.shape[1]
    assert O.shape == (A.shape[0], A.shape[0])
    Oinv = np.linalg.inv(O)
    assert Oinv.shape == (A.shape[0], A.shape[0])
    ATOinvA = A.T @ Oinv @ A
    assert ATOinvA.shape == (
        number_of_unknowns,
        number_of_unknowns,
    ), ATOinvA.shape
    ATOinvb = A.T @ Oinv @ b
    assert ATOinvb.shape == (number_of_unknowns,)
    return linalg_mod.solve(ATOinvA, ATOinvb)  # pylint: disable=E1101


def split_mapping_by_keys(mapping, key_lists):
    """Split up a mapping (dict) according to connected components

    Each connected component is a sequence of keys from mapping; returns a
    corresponding sequence of head mappings, each with only the keys from that
    connected component.

    """
    mappings = []
    for seq in key_lists:
        mappings.append(
            dict(
                (head_id, value)
                for head_id, value in list(mapping.items())
                if head_id in seq
            )
        )
    return mappings


def get_connected_components(head_mapping):
    """Find all overlapping sequences of heads in head_mapping

    Head_mapping is a mapping between head_ids and sets of ids for series that
    have data at that head_id.  Connected components are tuples of series ids
    that overlap in head.

    Series id sequences are returned sorted from longest to shortest.

    Head ids rather than heads are used because they must be hashable.

    """
    groups = {}
    for head_id, series_at_head in list(head_mapping.items()):
        matches = [
            keys
            for keys, group in list(groups.items())
            if not series_at_head.isdisjoint(group)
        ]
        new_keys = (head_id,) + sum(matches, ())
        others = [groups.pop(keys) for keys in matches]
        new_group = series_at_head.union(*others)
        groups[new_keys] = new_group
    connected_components = sorted(list(groups.keys()), key=len, reverse=True)
    # Sanity check: union should include all head_id ids
    assert sum(len(cc) for cc in connected_components) == len(head_mapping)
    assert set().union(*[set(cc) for cc in connected_components]) == set(
        head_mapping
    )
    return connected_components
