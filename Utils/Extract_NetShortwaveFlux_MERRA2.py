import os
import earthaccess
import cdsapi
import numpy as np
import xarray as xr
import pandas as pd
import geopy.distance
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


def extract_DSR(ncfile, Dataframe, TimeRangeCol, varname):
    # initiating new columns to append data to
    Dataframe["DSR_MERRA2_mean"] = np.nan
    Dataframe["DSR_MERRA2_sd"] = np.nan
    Dataframe["DSR_MERRA2_LatGridCell"] = np.nan
    Dataframe["DSR_MERRA2_LonGridCell"] = np.nan
    Dataframe['DSR_MERRA2_TimeSpan'] = np.nan
    Dataframe['DSR_MERRA2_TimeSpan'] = Dataframe['DSR_MERRA2_TimeSpan'].astype(object)
    Dataframe['DSR_MERRA2_Series'] = np.nan
    Dataframe['DSR_MERRA2_Series'] = Dataframe['DSR_MERRA2_Series'].astype(object)

    for i in range(0, len(Dataframe)):

        ##todo change to function
        print(f"processing data entry row no. : {Dataframe['index'][i]}")

        ilat = Dataframe['Lat'][i]
        ilon = Dataframe['Lon'][i]

        nclats = np.asarray(ncfile.variables['lat'][:])
        nclons = np.asarray(ncfile.variables['lon'][:])

        time_start = '1980-01-01'
        time_end = '2019-12-31-01'

        # calculating closest grid point to target location
        grid_point = []
        d = []
        grid_point_idx = []
        for lat_index in range(0, len(nclats)):
            for lon_index in range(0, len(nclons)):
                # add 3 degree buffer for lat long to speed things up
                if ilat - 3 < nclats[lat_index] < ilat + 3 and ilon - 3 < nclons[lon_index] < ilon + 3:
                    grid_point.append((nclats[lat_index], nclons[lon_index]))
                    grid_point_idx.append((lat_index, lon_index))
                    d.append(geopy.distance.geodesic((ilat, ilon), (nclats[lat_index], nclons[lon_index])).km)

        out = list(zip(d, grid_point,
                       grid_point_idx))  # zipped tuple of ((distance), (latgrid, longrid), (latgrid_idx, longrid_idx) )
        out.sort(key=lambda x: x[0])  # sorting so first is the closest point relative to my target point

        # now find valid grid cell with value closest to target location iterating over tuple of point distances
        for close_point_idx in range(0, len(out)):
            # slice data
            check_series = ncfile.sel(lat=out[close_point_idx][1][0],
                                      lon=out[close_point_idx][1][1],
                                      time=slice(time_start, time_end))

            if np.isnan(np.nanmean(check_series.variables[varname][:])):  # get average ignoring if there's nan in the series and so return true if there's at least one valid value
                print("Empty grid cell: Finding next valid point...")

            else:
                print(
                    f"found closest grid point with valid value {(out[close_point_idx][1][0], out[close_point_idx][1][1])}")
                print(f"original site has coords {ilat, ilon}")
                print(f"time considered is FROM {time_start} TO {time_end} ")

                values = np.asarray(check_series.variables[varname][:])

                # update values on dataframe
                Dataframe.loc[i, "DSR_MERRA2_sd"] = np.std(values)
                Dataframe.loc[i, "DSR_MERRA2_mean"] = np.mean(values)
                Dataframe.loc[i, "DSR_MERRA2_LatGridCell"] = out[close_point_idx][1][0]
                Dataframe.loc[i, "DSR_MERRA2_LonGridCell"] = out[close_point_idx][1][1]
                Dataframe.at[i, "DSR_MERRA2_TimeSpan"] = list((time_start, time_end))
                Dataframe.at[i, 'DSR_MERRA2_Series']=list(values)

                break

    return Dataframe


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

download_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/DSR_Flux_MERRA2/'
files = [os.path.join(download_path, item) for item in os.listdir(download_path)]
files.sort()

MERRA_dataset = xr.open_mfdataset(paths=files, engine='netcdf4', )  # combine multiple files into single file
varname = 'SWGNTWTR'
DSR_dataset = MERRA_dataset[[varname]]
DSR=DSR_dataset.variables[varname]


########## extracting closest DSR from my sites
lit_rev_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 3 - Results Chapter - Leo - General/Coral_QGIS_Spatial_Analyses/Revised_Literature_Points_Summarised.xlsx'
Dataframe_Literature = pd.read_excel(lit_rev_path, sheet_name="Aggregated_MapLayer")
Dataframe_Literature = Dataframe_Literature.reset_index()
TEST1 = extract_DSR(ncfile=DSR_dataset,varname='SWGNTWTR',
                   Dataframe=Dataframe_Literature, TimeRangeCol='Year_range')
outfile_name = os.path.join(os.path.dirname(lit_rev_path), 'DSR_MERRA2_Extracted_Results_LitReview.xlsx')
TEST1.to_excel(outfile_name)

#My Data
excel_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Tables_and_Regional_Sectors_Averages.xlsx'
Dataframe_MyData = pd.read_excel(excel_path, sheet_name="Table1")
# remove first row
Dataframe_MyData = Dataframe_MyData.drop([0]).reset_index()
TEST2 = extract_DSR(ncfile=DSR_dataset,varname='SWGNTWTR',
                   Dataframe=Dataframe_MyData, TimeRangeCol='Year_range')
outfile_name = os.path.join(os.path.dirname(excel_path), 'DSR_MERRA2_Extracted_Results_MyData.xlsx')
TEST2.to_excel(outfile_name)

#
# #doing some basic geospatial averaging for exploratory purposes and gif animations
# # Annual mean (i.e., the averaged value over a year each grid)
# cropped_ds = DSR_dataset.sel(lat=slice(-30,30), lon=slice(80,160))
#
# DSR_annual_mean = cropped_ds.resample(time="1YE").mean(dim='time', skipna=True)
# DSR_annual_sd = cropped_ds.resample(time="1YE").std(dim='time', skipna=True)
#
# ########## exploratory plots and animation
# t = 1
# p = DSR_annual_mean[varname].isel(time=[t]).plot(col='time',
#                                                  cmap=plt.cm.Spectral_r,
#                                                  cbar_kwargs={
#                                                      "orientation": "horizontal",
#                                                      "shrink": 0.9,
#                                                      "aspect": 40,
#                                                      "pad": 0.1, }
#                                                  )
# for ax in p.axes.flat:
#     ax.coastlines(resolution="110m", linewidth=0.5)
#
# pmap = DSR_annual_mean[varname].isel(time=t).plot(transform=ccrs.PlateCarree()),  # the data's projection
#                                             col='time', robust=True, cmap=plt.cm.Spectral_r,
#                                             cbar_kwargs={
#                                                 "orientation": "horizontal",
#                                                 "shrink": 0.9,
#                                                 "aspect": 40,
#                                                 "pad": 0.1,
#                                             },
#                                             subplot_kws={'projection': ccrs.PlateCarree(central_longitude=180)})  # the plot's projection # shift the original central longitude from 0 to 180
#
# # We have to set the map's options on all axes
# for ax in pmap.axes.flat:
#     ax.coastlines(resolution="110m", linewidth=0.5)
#
# # Plot main title
# main_title = "{}{} weekly_max".format("oi", "oi")
# plt.suptitle(main_title, fontweight='bold')
# plt.show()
#
# fig, ax = plt.subplots(figsize=(16, 5.5))
# annual_mean_global = DSR_annual_mean.groupby('time').mean(dim=['lat', 'lon'], skipna=True)
#
# # Convert to dataframe
# plotdata = annual_mean_global.to_dataframe()
# # List the statistics
# stat = plotdata.describe()
# print("stat:")
# print(stat)
#
# # Plot time series
# ax.plot(plotdata, label='annual_mean lat (-10:10) lon (106:136)', marker="o", linewidth=2)
# ax.legend(shadow=True, fancybox=True)
# plt.show()
#

