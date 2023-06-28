import numpy as np
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt

# Read sqlite query results into a pandas DataFrame
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3")
df1_average_rising_depth = pd.read_sql_query("SELECT * from average_rising_depth", con)
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3")
df1_rising_interval_zeta = pd.read_sql_query("SELECT * from rising_interval_zeta", con)
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3")
df1_rising_interval = pd.read_sql_query("SELECT * from rising_interval", con)
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3")
df1_zeta_grid = pd.read_sql_query("SELECT * from zeta_grid", con)
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3")
df1_storm_total_rain_depth = pd.read_sql_query("SELECT * from storm_total_rain_depth", con)
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3")
df1_storm_total_rise = pd.read_sql_query("SELECT * from storm_total_rise", con)
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3")
df1_storm = pd.read_sql_query("SELECT * from storm", con)
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3")
df1_rising_curve_line_segment = pd.read_sql_query("SELECT * from rising_curve_line_segment", con)
#con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_load_remove.sqlite3")
#df1_median_rising_depth = pd.read_sql_query("SELECT * from median_rising_depth", con)

df1_median_rising_depth = pd.merge(df1_rising_interval,df1_rising_interval_zeta, on="start_epoch")
df1_median_rising_depth["median_crossing_depth_mm"]=df1_median_rising_depth["mean_crossing_depth_mm"]+df1_median_rising_depth["rain_depth_offset_mm"]
df_manipulation = df1_median_rising_depth[["zeta_number","median_crossing_depth_mm"]].copy().groupby('zeta_number').median('median_crossing_depth_mm').reset_index().rename(columns = {'zeta_number' : 'zeta_mm'})

con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG.sqlite3")
df2_average_rising_depth = pd.read_sql_query("SELECT * from average_rising_depth", con)
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG.sqlite3")
df2_rising_interval_zeta = pd.read_sql_query("SELECT * from rising_interval_zeta", con)
con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG.sqlite3")
df2_discrete_zeta = pd.read_sql_query("SELECT * from discrete_zeta", con)

con = sqlite3.connect("/data/leuven/324/vsc32460/AC/spowtd/Brunei100_44_GS_smooth_IMERG_stretch_remove.sqlite3")
df_test = pd.read_sql_query("SELECT * from rising_curve_line_segment", con)
#df_rec_time = pd.read_sql_query("SELECT * from recession_interval", con)
#df_rec_et = pd.read_sql_query("SELECT * from evapotranspiration_staging", con)
#df_rec_wl = pd.read_sql_query("SELECT * from water_level_staging", con)
#df_rec=df_rec_time.copy()
#df_rec["end_epoch"] = df_rec_time["start_epoch"]+df_rec_time["time_offset_s"]
#df_rec =df_rec.drop(df_rec.columns[[1,2]],axis = 1)

#df_rec["intial_zeta_mm"]=np.where((df_rec["start_epoch"] == df_rec_wl["epoch"]),df_rec_wl["zeta_mm"])
_df = pd.read_sql_query("SELECT * from rising_curve_line_segment", con)
df2 = df.drop(df.columns[[0, 1]],axis = 1)
df2["Sy"] = df2["rain_total_depth_mm"]/(df2["final_zeta_mm"]-df2["initial_zeta_mm"])
df2["mean_zeta_mm"] = df2["final_zeta_mm"]-((df2["final_zeta_mm"]-df2["initial_zeta_mm"])/2)

df2["depth"] = df2["mean_zeta_mm"].apply(lambda x: -25.0 if x>=-100.0 else (-125 if -100.0>x>=-150.0 else (-175 if -150.0>x>=-200.0 else (-225 if -200.0>x>=-250.0 else (-275 if -250.0>x>=-300.0 else (-325 if -300.0>x>=-350.0 else (-375 if -350.0>x>=-400.0 else (-450 if -400.0>x>=-500.0 else (-550 if -500.0>x>=-600.0 else -650)))))))))
df_median=df2.groupby(["depth"])["Sy"].median()
df_median=df_median.to_frame()
df_median.reset_index(inplace=True)


# Verify that result of SQL query is stored in the dataframe
#pd.display(df)
df2['length'] = df2[["final_zeta_mm", "initial_zeta_mm"]].apply(tuple, axis=1)
x = list(df2["Sy"])
y = tuple(df2[["final_zeta_mm", "initial_zeta_mm"]].apply(tuple, axis=1))

plt.plot(x, [i for (i,j) in y], 'rs', markersize = 4)
plt.plot(x, [j for (i,j) in y], 'bo', markersize = 4)
plt.plot((x,x),([i for (i,j) in y], [j for (i,j) in y]),c='black')
plt.plot(df_median["Sy"], df_median["depth"], 'go-', markersize = 7)
plt.axhline(y=0.0, color='grey', linestyle='-')
plt.axhline(y=-100, color='grey', linestyle='-')
plt.axhline(y=-200, color='grey', linestyle='-')
plt.axhline(y=-300, color='grey', linestyle='-')
plt.axhline(y=-400, color='grey', linestyle='-')
plt.axhline(y=-500, color='grey', linestyle='-')
plt.axhline(y=-600, color='grey', linestyle='-')
plt.axhline(y=-700, color='grey', linestyle='-')


plt.xlim(xmin=0, xmax=1.5)
plt.ylim(ymin=-1000, ymax=100)

plt.ylabel("zeta interval [mm]")
plt.xlabel("Sy [-]")
plt.title("Sy figure ")
plt.show()
plt.tight_layout()
plt.close()


#assign abcde --> apply
#groupby abcde kolom en bereken gemiddelde per group
# smoothing window 1H is 3 timesteps 2H is 6 timestepsdf2["depth"] = np.where(df2.final_zeta_mm > 0.0, '0', 'NA')
print(df2)
Sy_median=df2["final_zeta_mm"].rolling(10, center=True).median()



con.close()