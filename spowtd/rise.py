"""Determine rising curves

"""

import logging

import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt

from spowtd.fit_offsets import (
    get_series_offsets,
    build_connected_head_mapping,
)


LOG = logging.getLogger('spowtd.rise')


def find_rise_offsets(
    connection, reference_zeta_mm=None, recharge_error_weight=0
):
    """Determine rising curves

    If recharge_error_weight is provided, an estimated variance-covariance
    matrix is assembled and used when solving the rise curve assembly problem.
    """
    cursor = connection.cursor()
    observed_recharge = True
    Sy_previous = 0
    master_rise_crossing_depth_mm_previous = 0
    compute_rise_offsets(
        cursor, reference_zeta_mm, recharge_error_weight=recharge_error_weight, Sy_previous=Sy_previous, master_rise_crossing_depth_mm_previous=master_rise_crossing_depth_mm_previous, observed_recharge=observed_recharge
    )
    cursor.close()
    connection.commit()


def get_series_storage_offsets(
    series_list, head_step, recharge_error_weight=0, Sy_previous=0, master_rise_crossing_depth_mm_previous=0, observed_recharge=0
):
    """Find a storage offset that minimizes difference in head crossing times

    This function is used in assembly of rise curves.  See further
    documentation under get_series_offsets.

    If recharge_error_weight is given, an estimated variance-covariance matrix
    is assembled and used when solving the rise curve assembly problem.

    """
    head_mapping, index_mapping = build_connected_head_mapping(
        series_list, head_step
    )
    if not recharge_error_weight:
        return get_series_offsets(head_mapping, index_mapping)
    covariance,weight = assemble_rise_covariance(
        head_mapping,
        index_mapping,
        series_list,
        head_step,
        recharge_error_weight,
        observed_recharge=observed_recharge,
        master_rise_crossing_depth_mm_previous=master_rise_crossing_depth_mm_previous,
    )
    return get_series_offsets(
        head_mapping, index_mapping, covariance=covariance, weight=weight
    )


def get_rise_covariance(connection, recharge_error_weight):
    """Use database connection to build covariance of rise event errors"""
    cursor = connection.cursor()
    series, _, _, _ = assemble_rise_series(cursor)
    head_step = get_head_step(cursor)
    cursor.close()
    head_mapping, index_mapping = build_connected_head_mapping(
        series, head_step
    )
    observed_recharge=True
    return assemble_rise_covariance(
        head_mapping, index_mapping, series, head_step, recharge_error_weight, observed_recharge=observed_recharge
    )


def assemble_rise_covariance(
    head_mapping, index_mapping, series, head_step, recharge_error_weight=1e3, master_rise_crossing_depth_mm_previous=0, observed_recharge=0,
):
    """Assemble covariance of rise event errors

    Creates and populates matrix omega representing the errors in the
    overdetermined system Ax = b for the rise curve assembly problem.

    head_mapping: maps head_id to a sequence of (series_id, time) pairs
                  indicating when each series crossed that head.
    index_mapping: maps series_id to the index of the series in series_list.
    series: list of (depth, rise) pairs where depth is a (0, storm_depth)
            duple and rise is a (zeta_initial, zeta_final) duple.
    head_step: factor by which to multiply integer head_id to get water
               level zeta.
    zeta_error_factor: ratio of error in zeta measurement to error induced
                       by recharge depth mismeasurement.  If zero or None, the
                       K x K identity matrix is returned where K is the total
                       number of equations.

    The covariance of rise event errors is a symmetric, positive semidefinite
    matrix with shape (n_equations, n_equations) characterizing the covariance
    among errors in each equation of the rise curve assembly problem.

    Recharge_error_factor is the relative weight for errors arising from
    mismeasurement of recharge depth.  The identity matrix, divided by
    recharge_error_weight, is added to the rise event error covariance to make
    it positive definite.

    This function does not verify that the computed covariance matrix is
    positive definite.  However, it does verify that the matrix is symmetric,
    and therefore successful Cholesky decomposition of the matrix
    (numpy.linalg.cholesky) will confirm that it is positive definite.  See
    Press et al. (2021) Numerical Recipes, 3rd ed.

    """
    # XXX Partly duplicates code in find_offsets

    # Assemble mapping of series ids to row numbers (unknowns) for the offset-
    #   finding problem.  Note that these are not the same as indices in the
    #   series list, which must be obtained via index_mapping:
    #   a_series = series[index_mapping[series_id]]
    series_ids = sorted(
        set().union(
            *(
                (series_id for series_id, _ in seq)
                for seq in list(head_mapping.values())
            )
        )
    )

    # Check input
    for head_id, series_at_head in head_mapping.items():
        zeta = head_id * head_step
        for series_id, depth_at_head in series_at_head:
            depth, rise = series[index_mapping[series_id]]
            assert len(depth) == 2
            assert len(rise) == 2
            assert depth[0] == 0, f'depth[0] = {depth[0]}'
            assert depth_at_head >= 0, depth_at_head
            assert depth_at_head <= depth[1], f'{depth_at_head} > {depth[1]}'
            assert (
                rise[0] <= zeta and rise[1] >= zeta
            ), f'Rise {rise} does not span {zeta}'
            del depth, rise, series_id, depth_at_head
        del head_id, zeta, series_at_head

    number_of_equations = sum(
        len(series_at_head) for series_at_head in list(head_mapping.values())
    )
    if not recharge_error_weight:
        return np.identity(number_of_equations)
    number_of_unknowns = len(series_ids) - 1
    # XXX /Duplicate code
    # (i, j)(k)
    ij = [
        (head_id, series_id)
        for head_id, series_at_head in sorted(head_mapping.items())
        for series_id, _ in series_at_head
    ]
    assert len(ij) == number_of_equations
    # k[i, j]
    eqn_no = {
        (head_id, series_id): eqn
        for eqn, (head_id, series_id) in enumerate(ij)
    }
    assert max(eqn_no.values()) == number_of_equations - 1

    omega = np.zeros(
        (number_of_equations, number_of_equations), dtype=np.float64
    )

    # R[j]
    rain_depth = {
        # Retrieving rain depth: obtain index, retrieve series;
        #   total_depth is the 2nd element of the first tuple
        series_id: series[index_mapping[series_id]][0][1]
        for series_id in series_ids
    }

    mean_rain_depth = sum(value for value in rain_depth.values()) / len(
        rain_depth
    )
    normalized_rain_depth = {
        key: value / mean_rain_depth for key, value in rain_depth.items()
    }

    # SA! outcommented this because it is still needed
    # del rain_depth

    # zeta_j^* - zeta_j^o
    rise = {
        # Retrieving rise: obtain index, retrieve series;
        #   initial_zeta, final_zeta is the second tuple
        series_id: (
            series[index_mapping[series_id]][1][1]
            - series[index_mapping[series_id]][1][0]
        )
        for series_id in series_ids
    }

    # zeta_bar
    mean_zeta = {
        series_id: sum(series[index_mapping[series_id]][1]) / 2
        for series_id in series_ids
    }

    # f[i, j]m
    coef = {
        (head_id, series_id): (head_id * head_step - mean_zeta[series_id])
        / rise[series_id]
        for head_id, series_at_head in head_mapping.items()
        for series_id, _ in series_at_head
    }
    #SA! added tolerance to cope with rounding errors of np in python and allow assertion to pass
    tolerance = 1e-10
    assert min(coef.values()) >= -0.5 - tolerance, min(coef.values())
    assert max(coef.values()) <= 0.5, max(coef.values())

    # To compute the terms in the covariance matrix, we need sets of heads
    # crossed by each series
    # I[j]
    series_heads = {}
    for head_id, series_at_head in head_mapping.items():
        for sid, _ in series_at_head:
            series_heads.setdefault(sid, []).append(head_id)
    # We also need sets of series crossing each head
    # J[i]
    Ji = {}
    for head_id, series_at_head in head_mapping.items():
        for sid, _ in series_at_head:
            Ji.setdefault(head_id, []).append(sid)

    #SA! added to calculate the weight variable for the weighted plotting of the master rise curve
    weight = {}
    for series_id, head_ids in series_heads.items():
        for head_id in head_ids:
            weight[series_id, head_id] = normalized_rain_depth[series_id]*coef[(head_id,series_id)]

    #SA! to calculate the average Sy of the master rise curve for the range of water levels that an event j spans
    #SA! only applied if it is NOT the first weight calculation, i.e., observed recharge is False
    if observed_recharge == False:

        # Average Sy of the master rise curve for all individual events
        Sy_avg = {}
        for series_id, head_ids in series_heads.items():
            min_WL = min(head_ids)
            max_WL = max(head_ids)
            Sy_avg[series_id] = abs(
                next((item[1] for item in master_rise_crossing_depth_mm_previous if item[0] == min_WL), None)-next((item[1] for item in master_rise_crossing_depth_mm_previous if item[0] == max_WL), None)
            )/(abs(min_WL-max_WL))

        #SA! calculate Rff using the 'adjusted recharge' based on the average Sy of the curve at those WL depths
        # divide by mean_rain_depth to get normalized rain
        Rff_app = {}
        for series_id, head_ids in series_heads.items():
            for head_id_1 in head_ids:
                for head_id_2 in head_ids:
                    Rff_app[(series_id, head_id_1, head_id_2)] = (
                        ((Sy_avg[series_id]*rise[series_id])/mean_rain_depth)** 2
                        * coef[(head_id_1, series_id)]
                        * coef[(head_id_2, series_id)]
                )

    # Terms
    # Rff[(j, i1, i2)] = R[j]^2 f[i_1, j] f[i_2, j]
    Rff = {}
    for series_id, head_ids in series_heads.items():
        for head_id_1 in head_ids:
            for head_id_2 in head_ids:
                Rff[(series_id, head_id_1, head_id_2)] = (
                    normalized_rain_depth[series_id] ** 2
                    * coef[(head_id_1, series_id)]
                    * coef[(head_id_2, series_id)]
                )

    # SA! switch approximated Rff ON or OFF
    if observed_recharge == False:
        Rff = Rff_app
    else:
        print('First iteration, use observed recharge')


    # Populate covariance matrix
    for k1 in range(number_of_equations):
        i1, j1 = ij[k1]
        Ji1 = set(Ji[i1])
        for k2 in range(number_of_equations):
            i2, j2 = ij[k2]
            Ji2 = set(Ji[i2])
            if j1 == j2:
                omega[k1, k2] += Rff[(j1, i1, i2)]
            if j1 in Ji2:
                omega[k1, k2] -= Rff[(j1, i1, i2)] / len(Ji2)
            if j2 in Ji1:
                omega[k1, k2] -= Rff[(j2, i2, i1)] / len(Ji1)
            omega[k1, k2] += sum(
                Rff[(j, i1, i2)] for j in Ji1.intersection(Ji2)
            ) / (len(Ji1) * len(Ji2))

    # Basic checks
    assert omega.shape == (
        number_of_equations,
        number_of_equations,
    ), f'Covariance matrix shape {omega.shape} != (# eqs, #eqs)'
    assert np.isfinite(omega).all(), (
        'Covariance matrix contains non-finite values: '
        f'{omega[~np.isfinite(omega)]}'
    )
    # Matrix may not be exactly symmetric due to roundoff
    assert np.allclose(
        omega, omega.T
    ), 'Covariance matrix not nearly symmetric'
    # Make exactly symmetric
    i_lower = np.tril_indices(number_of_equations, -1)
    omega[i_lower] = omega.T[i_lower]
    assert (omega == omega.T).all(), 'Covariance matrix not symmetric'
    # Normalize
    omega /= np.linalg.norm(omega)
    omega += np.identity(number_of_equations) / recharge_error_weight
    # Check that omega is positive definite
    L = np.linalg.cholesky(omega)
    assert np.allclose(np.dot(L, L.T), omega)
    return omega, weight


def compute_rise_offsets(cursor, reference_zeta_mm, recharge_error_weight=0, Sy_previous=0, master_rise_crossing_depth_mm_previous=0, observed_recharge=0):
    """Compute time offsets to populate rising_interval_zeta

    If a reference zeta is not None, the crossing-depth of this water level is
    used as the origin of the axis.  Otherwise, the highest water level in the
    longest assembled rise is used.

    """
    series, series_outlier, epoch, zeta_intervals = assemble_rise_series(cursor)
    delta_z_mm = get_head_step(cursor)

    # Solve for offsets
    indices, offsets, zeta_mapping, weight, = get_series_storage_offsets(
        series, delta_z_mm, recharge_error_weight=recharge_error_weight, Sy_previous=Sy_previous,
        master_rise_crossing_depth_mm_previous=master_rise_crossing_depth_mm_previous, observed_recharge=observed_recharge
    )

    indices_old=indices
    zeta_mapping_old = zeta_mapping

# START ------------------------------ #

    #SA! The general approach of the changes included below is the following:
    #SA! First an additional series is created including information on the specific yield
    #SA! Second a dictionary is created from the series to hold the original index as a key throughout the elimination process
    #SA! Third is the removal process: 1) events that are excluded during the offset calculation are removed from the series,
    #SA! 2) unrealistic events (Sy>1) are removed, 3) outlier detection and removal, and 4) events excluded during recalculation of offsets are removed

    # create a series_dictionary from the list to maintain the indices as a key value throughout the elimination process from the series
    indices_outlier = list(range(len(series)))
    series_dictionary = {}
    series_copy = copy.deepcopy(series)
    for key in indices_outlier:
        for value in series_copy:
            series_dictionary[key] = value
            series_copy.remove(value)
            #zeta_intervals.remove(value)
            break
    series_dictionary

    if len(indices) != len(series):
        drop_1 = list(sorted(set(range(len(series)))-set(indices)))
        for element in sorted(drop_1, reverse=True):
            del series_dictionary[element]
            #del zeta_intervals[element]

    # SA! remove unrealistic and outlier rise events
    outliers_removal = 0
    if outliers_removal == 1 and observed_recharge == True:
        #SA! This part removes the unrealistic rise events by filtering for series_id with Sy > 1 and removing them.
        #SA! Create an ndarray with the index of all series_ids to remove
        drop_2 = np.empty((0), dtype=int)
        for index, (t, H, Sy) in enumerate(series_outlier):
            if Sy > 1:
                drop_2 = np.append(drop_2, index)

        #SA! loop over ndarray in reverse order to maintain index position while removing
        #SA! remove from variables: series, indices and zeta_mapping (via crossing)
        for k in sorted(drop_2, reverse=True):
            try:
                #del zeta_intervals[indices.index(k)]
                indices.remove(k)
                del series_dictionary[k]
            except:
                print('This index is already removed from series')
            for discrete_zeta, crossings in zeta_mapping.items():
                for series_id, mean_crossing_depth_mm in crossings:
                    if k == series_id:
                        for i in range(len(crossings)):
                            if len(crossings) <= i:
                                continue
                            elif k in crossings[i]:
                                del crossings[i]

        #SA! This part identifies outliers at each zeta_mm increment based on the Sy values and adds the series_id of the outlier to a counting array
        #SA! At each zeta_mm increment dataframe with series_id and Sy value are created
        outlier_count = np.empty((0, 2), dtype=float)
        zeta_step_sy = []
        occurence_count = []
        for discrete_zeta, crossings in zeta_mapping.items():
            for series_id, mean_crossing_depth_mm in crossings:
                occurence_count.append(series_id)
                for index, (t, H, Sy) in enumerate(series_outlier):
                    if index == series_id:
                        for i in range(len(crossings)):
                            if index in crossings[i]:
                                specific_yield = np.array([series_id,Sy])
                zeta_step_sy.append(specific_yield.copy())
            df= pd.DataFrame(list(map(np.ravel, zeta_step_sy)), columns=['series_id','Sy'])
            zeta_step_sy=[]

            #SA! At each zeta_mm increment determine Sy outliers and add the series_id to the counting array
            #SA! If there are less than 5 events on the zeta_m apply boxplot outlier detection
            #SA! if there are 5 or more apply the 95th and 5th percentile outlier detection
            if len(df) >= 5:
                highest_allowed = df['Sy'].quantile(0.95)
                lowest_allowed = df['Sy'].quantile(0.05)
            else:
                IQR = df['Sy'].quantile(0.75)- df['Sy'].quantile(0.25)
                highest_allowed = df['Sy'].quantile(0.75) + 1.5*(IQR)
                lowest_allowed = df['Sy'].quantile(0.25) - 1.5*(IQR)
            df_outliers = df[(df['Sy'] > highest_allowed) | (df['Sy'] < lowest_allowed)]
            outlier_count = np.append(outlier_count, df_outliers.to_numpy(), axis =0)

        #SA! This part creates an ndarray with the outlier series_id based on the count
        drop_3 = np.empty((0), dtype=int)
        for j in range(len(series_dictionary)):
            how_many = np.count_nonzero(outlier_count[:,0] == j)
            total = np.array(occurence_count)
            if (how_many/(np.count_nonzero(total[:] == j)+0.000000001) > 1/2):
                drop_3 = np.append(drop_3, j)

        # SA! loop over ndarray in reverse order to maintain index position while removing
        # SA! remove from variables: series, indices and zeta_mapping (via crossing)
        for z in sorted(drop_3, reverse=True):
            #del zeta_intervals[indices.index(z)]
            del series_dictionary[z]
            indices.remove(z)
            for discrete_zeta, crossings in zeta_mapping.items():
                for series_id, mean_crossing_depth_mm in crossings:
                    if z == series_id:
                        for i in range(len(crossings)):
                            if len(crossings) <= i:
                                continue
                            elif z in crossings[i]:
                                del crossings[i]

    series_copy = list(series_dictionary.values())
    series=series_copy

    #SA! This solves for the offsets again after removing unrealistic and outlier events from series.
    #SA! Creates a copy of zeta_mapping and indices because the original with removed events are used.
    indices_copy, offsets, zeta_mapping_copy, weight = get_series_storage_offsets(
        series_copy, delta_z_mm, recharge_error_weight=recharge_error_weight, Sy_previous=Sy_previous, master_rise_crossing_depth_mm_previous=master_rise_crossing_depth_mm_previous, observed_recharge=observed_recharge)

    #SA! If additional series_ids were removed from offsets (necessary overlap of events)
    #SA! Check this series_id and drop all zeta_mm from zeta_mappings with this series_id and from the indices.
    for key in (zeta_mapping.keys() - zeta_mapping_copy.keys()):
        zeta_mapping.pop(key,None)
    unique_ids = []
    if len(indices) != len(offsets):
        for discrete_zeta, crossings in zeta_mapping.items():
            for series_id, mean_crossing_depth_mm in crossings:
                unique_ids.append(series_id)
        missing = list(sorted((set(indices))-(set(unique_ids))))
        for ele in missing:
            indices.remove(ele)


# STOP ------------------------------ #


    reference_zeta_off_grid = (
        reference_zeta_mm is not None
        and not np.allclose(reference_zeta_mm % delta_z_mm, 0)
    )
    if reference_zeta_off_grid:
        raise ValueError(
            'Reference zeta {} mm not evenly divisible by '
            'zeta step {} mm'.format(reference_zeta_mm, delta_z_mm)
        )
    if reference_zeta_mm is not None:
        reference_index = int(reference_zeta_mm / delta_z_mm)
    else:
        reference_index = max(zeta_mapping.keys())

    mean_zero_crossing_depth_mm = np.array(
        [
            offsets[indices.index(series_id)] + depth_mean_mm
            for series_id, depth_mean_mm in zeta_mapping[reference_index]
        ]
    ).mean()


    if recharge_error_weight == 0:
        print('no weighting applied')
    else:
        Sy, master_rise_crossing_depth_mm = get_master_rise_curve(
            indices, offsets, mean_zero_crossing_depth_mm, zeta_mapping, weight
        )

        all_master_rise = []
        all_master_rise.append(master_rise_crossing_depth_mm)

        Sy_previous=Sy
        master_rise_crossing_depth_mm_previous=master_rise_crossing_depth_mm


        #SA! SWITCH TO TURN OF ITERATION --> TRUE = OFF     FALSE = ON
        observed_recharge = True


        count = 1
        while observed_recharge == False:
            # Solve for offsets
            indices, offsets, zeta_mapping, weight = get_series_storage_offsets(
                series, delta_z_mm, recharge_error_weight=recharge_error_weight, Sy_previous=Sy_previous,
                master_rise_crossing_depth_mm_previous=master_rise_crossing_depth_mm_previous,
                observed_recharge=observed_recharge
            )

            reference_zeta_off_grid = (
                    reference_zeta_mm is not None
                    and not np.allclose(reference_zeta_mm % delta_z_mm, 0)
            )
            if reference_zeta_off_grid:
                raise ValueError(
                    'Reference zeta {} mm not evenly divisible by '
                    'zeta step {} mm'.format(reference_zeta_mm, delta_z_mm)
                )
            if reference_zeta_mm is not None:
                reference_index = int(reference_zeta_mm / delta_z_mm)
            else:
                reference_index = max(zeta_mapping.keys())

            mean_zero_crossing_depth_mm = np.array(
                [
                    offsets[indices.index(series_id)] + depth_mean_mm
                    for series_id, depth_mean_mm in zeta_mapping[reference_index]
                ]
            ).mean()

            Sy, master_rise_crossing_depth_mm = get_master_rise_curve(
                indices, offsets, mean_zero_crossing_depth_mm, zeta_mapping, weight
            )

            all_master_rise.append(master_rise_crossing_depth_mm)

            count = count + 1

            Sy_diff = sum((Sy_previous[i][1]-Sy[i][1])**2 for i in range(len(Sy)-1))
            print('Sy_diff = ' + str(Sy_diff))
            if Sy_diff >= 0.5:
                Sy_previous = Sy
                master_rise_crossing_depth_mm_previous = master_rise_crossing_depth_mm
                print('This is iteration number ' + str(count) + ', not converged yet.')
            else:
                observed_recharge = None
                print('An equilibrium (=converged) was reached after ' + str(count) + ' iterations.')

        for i, data_list in enumerate(all_master_rise):
            column1_data = [point[0] for point in data_list]
            column2_data = [point[1] for point in data_list]

            plt.plot(column2_data, column1_data, label=f'Iteration {i + 1} ')

        # Add labels and legend to the plot
        plt.xlabel('Dynamic storage (mm)')
        plt.ylabel('WL (mm)')
        plt.title('Master rise curve iterations')
        plt.legend()
        #plt.show()

    for i, series_id in enumerate(indices_old):
        interval = zeta_intervals[series_id]
        del series_id
        cursor.execute(
            """
        INSERT INTO rising_interval (
            start_epoch, rain_depth_offset_mm)
        SELECT :start_epoch,
                :rain_depth_offset_mm""",
            {
                'start_epoch': epoch[interval[0]],
                'rain_depth_offset_mm': (
                    offsets[i] - mean_zero_crossing_depth_mm
                ),
            },
        )
        del interval

    #SA! added to write out the weighted master rise curve for plotting, only if recharge error weight is not 0
    if recharge_error_weight != 0:
        for discrete_zeta, crossings in zeta_mapping_old.items():
            master_rise_crossing_depth = next((arr[0] for t, arr in master_rise_crossing_depth_mm if t == discrete_zeta), None)
            cursor.execute(
                """ 
            INSERT INTO rising_interval_zeta_weighted (
                zeta_number,
                master_rise_crossing_depth_mm)
            SELECT  :discrete_zeta,
                    :master_rise_crossing_depth_mm""",
                {
                    'discrete_zeta': discrete_zeta,
                    'master_rise_crossing_depth_mm': master_rise_crossing_depth,
                },
            )
            del discrete_zeta, crossings

    for discrete_zeta, crossings in zeta_mapping_old.items():
        for series_id, mean_crossing_depth_mm in crossings:
            interval = zeta_intervals[series_id]
            del series_id
            cursor.execute(
                """
            INSERT INTO rising_interval_zeta (
                start_epoch, zeta_number,
                mean_crossing_depth_mm)
            SELECT :start_epoch,
                    :discrete_zeta,
                    :mean_crossing_depth_mm""",
                {
                    'start_epoch': epoch[interval[0]],
                    'discrete_zeta': discrete_zeta,
                    'mean_crossing_depth_mm': mean_crossing_depth_mm,
                },
            )
            del mean_crossing_depth_mm, interval
        del discrete_zeta, crossings
    cursor.close()


def get_head_step(cursor):
    """Get grid interval in mm from database cursor

    If grid_interval_mm is not yet populated in zeta_grid, ValueError is raised

    """
    try:
        delta_z_mm = cursor.execute(
            """
        SELECT (grid_interval_mm) FROM zeta_grid
        """
        ).fetchone()[0]
    except TypeError:
        raise ValueError(  # pylint: disable=raise-missing-from
            "Discrete water level interval not yet set"
        )
    return delta_z_mm

def get_master_rise_curve(indices, offsets, mean_zero_crossing_depth_mm, zeta_mapping, weight):

    #SA! calculate the master rise - from that the specific yield - and loop over this until converges
    indices_sorted = np.argsort(indices)
    offsets_sorted = offsets[indices_sorted]
    rain_depth_offset_mm = offsets - mean_zero_crossing_depth_mm
    rain_depth_offset_mm = np.stack((np.array(indices).T,rain_depth_offset_mm.T), axis = 0)

    #unique_values = set(item[0] for sublist in zeta_mapping.values() for item in sublist)

    master_rise_crossing_depth_mm = []
    for discrete_zeta, crossings in zeta_mapping.items():
        value_sum = 0
        weight_sum = 0
        for series_id, mean_crossing_depth_mm in crossings:
            # find the weight linked a (discrete_zeta, series_id) combination
            if (indices.index(series_id), discrete_zeta) in weight:
                weight_app = 1/abs(weight[(indices.index(series_id),discrete_zeta)])
                if weight_app == float('inf'):
                    weight_app = 1
            value = rain_depth_offset_mm[1, np.where(rain_depth_offset_mm[0] == series_id)[0]] + mean_crossing_depth_mm
            weight_sum += weight_app
            value_sum += value*weight_app
        master_zeta_value=value_sum/weight_sum
        master_rise_crossing_depth_mm.append([discrete_zeta,master_zeta_value])
    #print(master_rise_crossing_depth_mm)

    #master_rise_crossing_depth_mm = np.array(master_rise_crossing_depth_mm)
    Sy = []
    for i in range(len(master_rise_crossing_depth_mm)-1):
        i=i+1
        Sy_value= (abs(master_rise_crossing_depth_mm[i][1])-abs(master_rise_crossing_depth_mm[i-1][1]))/(abs(master_rise_crossing_depth_mm[i][0])-abs(master_rise_crossing_depth_mm[i-1][0]))
        Sy.append((master_rise_crossing_depth_mm[i-1][0],Sy_value))

    master_rise_crossing_depth_mm = [tuple(x) for x in master_rise_crossing_depth_mm]

    return Sy, master_rise_crossing_depth_mm


def assemble_rise_series(cursor):
    """Assemble rise series

    Returns:
      series:  A list of pairs of arrays
               ((0, total_depth), (initial_zeta, final_zeta)),
               each representing a distinct rise event.
      epoch:  Time vector as seconds since the UNIX epoch
      zeta_intervals:  A list of (start, thru) indices in the epoch vector,
                       one for each series

    """
    epoch, zeta_mm = [
        np.array(v, dtype='float64')
        for v in zip(
            *cursor.execute(
                """
     SELECT epoch, zeta_mm FROM water_level
     ORDER BY epoch"""
            )
        )
    ]
    assert np.isfinite(zeta_mm).all()
    cursor.execute(
        """
SELECT s.start_epoch,
       s.thru_epoch,
       zi.start_epoch,
       zi.thru_epoch
FROM storm AS s
JOIN zeta_interval_storm AS zis
  ON s.start_epoch = zis.storm_start_epoch
JOIN zeta_interval AS zi
  ON zi.start_epoch = zis.interval_start_epoch
ORDER BY s.start_epoch"""
    )
    series = []
    # SA! added series_outlier for removal
    series_outlier = []
    rain_intervals = []
    zeta_intervals = []
    for (
        storm_start_epoch,
        storm_thru_epoch,
        zeta_start_epoch,
        zeta_thru_epoch,
    ) in cursor.fetchall():
        rain_start = np.argwhere(epoch == storm_start_epoch)[0, 0]
        # Epoch associated with rainfall intensities are *start*
        # epoch for the interval, so the time slice that *starts*
        # at the storm thru_epoch is not included.
        rain_stop = np.argwhere(epoch == storm_thru_epoch)[0, 0]
        zeta_start = np.argwhere(epoch == zeta_start_epoch)[0, 0]
        zeta_thru = np.argwhere(epoch == zeta_thru_epoch)[0, 0]
        cursor.execute(
            """
    SELECT total_depth_mm
    FROM storm_total_rain_depth
    WHERE storm_start_epoch = :storm_start_epoch""",
            {'storm_start_epoch': storm_start_epoch},
        )
        total_depth = cursor.fetchone()[0]
        assert zeta_thru > zeta_start
        zeta_seq = zeta_mm[zeta_start : zeta_thru + 1]
        assert (
            len(zeta_seq) >= 2
        ), 'A jump is defined by at least two zeta values'
        assert (
            np.diff(zeta_seq) > 0
        ).all(), '{} is not strictly increasing'.format(zeta_seq)
        initial_zeta = zeta_seq[0]
        final_zeta = zeta_seq[-1]
        assert (
            len(zeta_seq) > 0
        ), 'empty sequence'  # pylint:disable=len-as-condition
        rain_intervals.append((rain_start, rain_stop))
        zeta_intervals.append((zeta_start, zeta_thru + 1))
        series.append(
            (np.array((0, total_depth)), np.array((initial_zeta, final_zeta)))
        )
        # SA! Create an ndarray that includes Sy to filter out rise events based on Sy
        series_outlier.append(
            (np.array((0, total_depth)), np.array((initial_zeta, final_zeta)), total_depth/(final_zeta-initial_zeta))
        )


    del zeta_mm
    return series, series_outlier, epoch, zeta_intervals
