import os
import earthaccess
import cdsapi
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import xarray as xr
import glob
import nctoolkit as nc
import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns
from haversine import haversine
from scipy.spatial import KDTree
import geopy.distance
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import colors as c
from mpl_toolkits.basemap import Basemap, shiftgrid

#
# #copernicus dataset fail
# c = cdsapi.Client()
#
# c.retrieve(
#     'satellite-surface-radiation-budget',
#     {
#         'format': 'zip',
#         'product_family': 'clara_a3',
#         'month': [
#             '01', '02', '03',
#             '04', '05', '06',
#             '07', '08', '09',
#             '10', '11', '12',
#         ],
#         'year': [
#             '1979', '1980', '1981',
#             '1982', '1983', '1984',
#             '1985', '1986', '1987',
#             '1988', '1989', '1990',
#             '1991', '1992', '1993',
#             '1994', '1995', '1996',
#             '1997', '1998', '1999',
#             '2000', '2001', '2002',
#             '2003', '2004', '2005',
#             '2006', '2007', '2008',
#             '2009', '2010', '2011',
#             '2012', '2013', '2014',
#             '2015', '2016', '2017',
#             '2018', '2019', '2020',
#         ],
#         'time_aggregation': 'monthly_mean',
#         'climate_data_record_type': 'thematic_climate_data_record',
#         'variable': [
#             'surface_net_downward_radiative_flux', 'surface_net_downward_shortwave_flux',
#         ],
#         'origin': 'eumetsat',
#     },
#     '/Users/leonardobertini/Desktop/NetSurfaceRadiation/download.zip')
#

# DOWNLOAD files and check for updates
# earthaccess.login()
# timechunks = [("1980-01-01T00:00", "2019-12-01T00:00:00")]
# files_long = []
# download_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/DSR_Flux'
# for pair in timechunks:
#     results = earthaccess.search_data(
#         short_name='M2TMNXOCN',
#         cloud_hosted=True,
#         temporal=pair,
#         count=2000)
#
#     files = earthaccess.download(results, download_path)
#     files_long.append(files)
#



download_path='/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/DSR_Flux/'
files=[os.path.join(download_path,item) for item in os.listdir(download_path)]
files.sort()

example = xr.open_mfdataset(paths=files, engine='netcdf4', )  # combine multiple files into single file
varname='SWGNTWTR'
DSR=example[varname]


# Annual mean (i.e., the averaged value over a year each grid)
DSR_annual_mean = DSR.resample(time="1YE").mean(dim='time', skipna=True)
DSR_annual_sd = DSR.resample(time="1YE").std(dim='time', skipna=True)


########## exploratory plots
pmap = DSR_annual_sd.isel(time=[0, 1]).plot(transform=ccrs.PlateCarree(),  # the data's projection
                                                 col='time', col_wrap=2, robust=True,  # multiplot settings
                                                 cmap=plt.cm.Spectral_r,
                                                 cbar_kwargs={
                                                     "orientation": "horizontal",
                                                     "shrink": 0.9,
                                                     "aspect": 40,
                                                     "pad": 0.1,
                                                 },
                                                 subplot_kws={'projection': ccrs.PlateCarree(
                                                     central_longitude=180)})  # the plot's projection # shift the original central longitude from 0 to 180

# We have to set the map's options on all axes
for ax in pmap.axes.flat:
    ax.coastlines(resolution="110m", linewidth=0.5)

# Plot main title
main_title = "{}{} weekly_max".format("oi", "oi")
plt.suptitle(main_title, fontweight='bold')
plt.show()


fig, ax = plt.subplots(figsize=(16,5.5))
annual_mean_global = DSR_annual_mean.groupby('time').mean(dim=['lat','lon'],skipna=True)

# Convert to dataframe
plotdata = annual_mean_global.to_dataframe()
# List the statistics
stat = plotdata.describe()
print("stat:")
print(stat)

# Plot time series
ax.plot(plotdata,label='annual_global_mean', marker="o", linewidth=2)
ax.legend(shadow=True, fancybox=True)


########## extracting closest DSR from my sites

def extract_DSR(ncfile, Dataframe, TimeRangeCol):
    # initiating new columns to append data to
    Dataframe["SST_HadSST_ann"] = np.nan
    Dataframe["SST_HadSST_min"] = np.nan
    Dataframe["SST_HadSST_max"] = np.nan
    Dataframe["SST_HadSST_sd"] = np.nan
    Dataframe["SST_HadSST_LatGridCell"] = np.nan
    Dataframe["SST_HadSST_LonGridCell"] = np.nan
    Dataframe["SST_HadSST_TimeSpan"] = np.nan
    Dataframe['SST_HadSST_TimeSpan'] = Dataframe['SST_HadSST_TimeSpan'].astype(object)

    for i in range(0, len(Dataframe)):

        ##todo change to function
        print(f"processing data entry row no. : {Dataframe['index'][i]}")

        ilat = Dataframe['Lat'][i]
        ilon = Dataframe['Lon'][i]

        nclats = np.asarray(ncfile.variables['latitude'][:])
        nclons = np.asarray(ncfile.variables['longitude'][:])

        time_start = '1980-01-01T00:00'
        time_end = '2019-12-31-01T00:00'

        # calculating closest grid point to target location
        grid_point = []
        d = []
        grid_point_idx = []
        for lat_index in range(0, len(nclats)):
            for lon_index in range(0, len(nclons)):
                grid_point.append((nclats[lat_index], nclons[lon_index]))
                grid_point_idx.append((lat_index, lon_index))
                d.append(geopy.distance.geodesic((ilat, ilon), (nclats[lat_index], nclons[lon_index])).km)

        out = list(zip(d, grid_point,
                       grid_point_idx))  # zipped tuple of ((distance), (latgrid, longrid), (latgrid_idx, longrid_idx) )
        out.sort(key=lambda x: x[0])  # sorting so first is the closest point relative to my target point

        # now find valid grid cell with value closest to target location iterating over tuple of point distances
        for close_point_idx in range(0, len(out)):
            # slice data
            check_series = ncfile.sel(latitude=out[close_point_idx][1][0],
                                      longitude=out[close_point_idx][1][1],
                                      time=slice(time_start, time_end))

            if np.isnan(np.nanmean(check_series.variables['sst'][
                                   :])):  # get average ignoring if there's nan in the series and so return true if there's at least one valid value
                print("Empty grid cell: Finding next valid point...")

            else:
                print(
                    f"found closest grid point with valid value {(out[close_point_idx][1][0], out[close_point_idx][1][1])}")
                print(f"original site has coords {ilat, ilon}")
                print(f"time considered is FROM {time_start} TO {time_end} ")

                values = np.asarray(check_series.variables['sst'][:])

                # update values on dataframe
                Dataframe.loc[i, "SST_HadSST_sd"] = np.std(values)
                Dataframe.loc[i, "SST_HadSST_ann"] = np.mean(values)
                Dataframe.loc[i, "SST_HadSST_min"] = values.min()
                Dataframe.loc[i, "SST_HadSST_max"] = values.max()
                Dataframe.loc[i, "SST_HadSST_LatGridCell"] = out[close_point_idx][1][0]
                Dataframe.loc[i, "SST_HadSST_LonGridCell"] = out[close_point_idx][1][1]
                Dataframe.at[i, "SST_HadSST_TimeSpan"] = list((time_start, time_end))

                # todo get column for closes valid grid point coordinates
                break

    return Dataframe

lit_rev_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 3 - Results Chapter - Leo - General/Coral_QGIS_Spatial_Analyses/Revised_Literature_Points_Summarised.xlsx'
Dataframe_Literature = pd.read_excel(lit_rev_path, sheet_name="Aggregated_MapLayer")
Dataframe_Literature = Dataframe_Literature.reset_index()
TEST = extract_DSR(ncfile=DSR_annual_mean, Dataframe=Dataframe_Literature, TimeRangeCol='Year_range')
