import datetime
import sqlite3
from spowtd.load import *
from spowtd.classify import *
from spowtd.zeta_grid import *
from spowtd.rise import *
from spowtd.set_curvature import *
from spowtd.recession import *
from spowtd.plot_rise import *

from collections import defaultdict

import numpy as np
import pandas as pd

load = 1
if load == 1:
    connection = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_4C_OG.sqlite3")  # "/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3"
    precipitation_data_file = open('/data/leuven/324/vsc32460/AC/spowtd/src/precipitation_Brunei_IMERG_CDF.txt', mode='r') #'/data/leuven/324/vsc32460/AC/spowtd/src/throughfall_Brunei_30min.txt' precipitation_Brunei_IMERG_CDF.txt
    evapotranspiration_data_file = open('/data/leuven/324/vsc32460/AC/spowtd/src/evapotranspiration_no_intercep_Brunei_30min.txt') #'/data/leuven/324/vsc32460/AC/spowtd/src/evapotranspiration_no_intercep_Brunei_30min.txt'
    water_level_data_file = open('/data/leuven/324/vsc32460/AC/spowtd/src/waterlevel_Brunei_20min_3times_windowaverage_GS.txt') #'/data/leuven/324/vsc32460/AC/spowtd/src/waterlevel_Brunei_20min_3times_windowaverage_GS.txt'
    time_zone_name = 'Asia/Brunei' #'Africa/Lagos' #'Asia/Brunei'
    cursor = connection.cursor()
    load_data(connection,precipitation_data_file,evapotranspiration_data_file,water_level_data_file,time_zone_name,)
    storm_rain_threshold_mm_h = 4.0
    rising_jump_threshold_mm_h = 4.0
    classify_intervals(connection, storm_rain_threshold_mm_h, rising_jump_threshold_mm_h)
    #df_zeta_grid = pd.read_sql_query("SELECT * from zeta_grid", connection)
    #grid_interval = df_zeta_grid.grid_intervall_mm
    populate_zeta_grid(connection, grid_interval_mm=1.0)
    find_rise_offsets(connection,reference_zeta_mm=None)
    set_curvature(connection, curvature_m_km2=1.0)
    find_recession_offsets(connection, reference_zeta_mm=None)
    #plot_rise(connection, parameters=None)

# APPLY ensemble of IMERG data over SPOWTD to
elif load == 2:
    for i in np.arange(0, 64):
        print(i)
        connection = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_throughfall_4C"+ "_mem_" + str(i) + "SD-1.sqlite3")  # "/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3"
        precipitation_data_file = open('/data/leuven/324/vsc32460/AC/spowtd/src/EM/throughfall_mem_' + str(i) + '/throughfall_Brunei_30min_mem_' + str(i) + 'SD-1.txt', mode='r')  # '/data/leuven/324/vsc32460/AC/spowtd/src/precipitation_Brunei_IMERG_CDF.txt'
        evapotranspiration_data_file = open('/data/leuven/324/vsc32460/AC/spowtd/src/evapotranspiration_no_intercep_Brunei_30min.txt')  # '/data/leuven/324/vsc32460/AC/spowtd/src/evapotranspiration_no_intercep_Brunei_30min.txt'
        water_level_data_file = open('/data/leuven/324/vsc32460/AC/spowtd/src/waterlevel_Brunei_20min_3times_windowaverage_GS.txt')  # '/data/leuven/324/vsc32460/AC/spowtd/src/waterlevel_Brunei_20min_3times_windowaverage_GS.txt'
        time_zone_name = 'Asia/Brunei' #'Africa/Lagos'
        cursor = connection.cursor()
        load_data(connection, precipitation_data_file, evapotranspiration_data_file, water_level_data_file,
                  time_zone_name, )
        storm_rain_threshold_mm_h = 4.0
        rising_jump_threshold_mm_h = 4.0
        classify_intervals(connection, storm_rain_threshold_mm_h, rising_jump_threshold_mm_h)
        # df_zeta_grid = pd.read_sql_query("SELECT * from zeta_grid", connection)
        # grid_interval = df_zeta_grid.grid_intervall_mm
        populate_zeta_grid(connection, grid_interval_mm=1.0)
        find_rise_offsets(connection, reference_zeta_mm=None)
        set_curvature(connection, curvature_m_km2=1.0)
        find_recession_offsets(connection, reference_zeta_mm=None)
        # plot_rise(connection, parameters=None)

elif load == 0:
    connection = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_classify_remove.sqlite3")
    cursor = connection.cursor()
    storm_rain_threshold_mm_h = 4.0
    rising_jump_threshold_mm_h = 4.0
    classify_intervals(connection, storm_rain_threshold_mm_h, rising_jump_threshold_mm_h)



print('1+1=2')
