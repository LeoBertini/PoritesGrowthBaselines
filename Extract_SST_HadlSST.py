import os
import earthaccess
import netCDF4
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

def fix_calendar(ds, timevar='T'):
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds


def extract_HadSST(ncfile, Dataframe, TimeRangeCol):
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

        time_start = (Dataframe[TimeRangeCol][i]).split('-')[0] + '-01-01T00:00'
        time_end = (Dataframe[TimeRangeCol][i]).split('-')[1] + '-12-01T00:00'

        # if data point is too old for series. adopt earliest 1961-1990 mean and then subtract based on expected anomaly
        if int(Dataframe[TimeRangeCol][i].split('-')[0]) <= 1870:
            print(
                f"Time range not in HadSST series. Taking 1961-1990 average instead to then apply proxy-based anomalies from D'Arrigo")
            time_start = '1961-01-01T00:00'
            time_end = '1990-12-01T00:00'

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
                                   :])):  # get average temperature ignoring if there's nan in the series and so return true if there's at least one valid value
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


def extract_ERSSTv5(ncfile, Dataframe, TimeRangeCol):
    Dataframe["SST_ERSSTv5_sd"] = np.nan
    Dataframe["SST_ERSSTv5_ann"] = np.nan
    Dataframe["SST_ERSSTv5_min"] = np.nan
    Dataframe["SST_ERSSTv5_max"] = np.nan
    Dataframe["SST_ERSSTv5_LatGridCell"] = np.nan
    Dataframe["SST_ERSSTv5_LonGridCell"] = np.nan
    Dataframe["SST_ERSSTv5_TimeSpan"] = np.nan
    Dataframe['SST_ERSSTv5_TimeSpan'] = Dataframe['SST_ERSSTv5_TimeSpan'].astype(object)

    for i in range(0, len(Dataframe)):

        print(f"processing data entry row no. : {Dataframe['index'][i]}")

        ilat = Dataframe['Lat'][i]
        ilon = Dataframe['Lon'][i] + 180  # because longitude is gridded 0-360 in ERSST,

        nclats = np.asarray(ncfile.variables['lat'][:])
        nclons = np.asarray(ncfile.variables['lon'][:])

        time_start = (Dataframe[TimeRangeCol][i]).split('-')[0] + '-01-01T00:00'
        time_end = (Dataframe[TimeRangeCol][i]).split('-')[1] + '-12-01T00:00'

        # if data point is too old for series. adopt earliest 1961-1990 mean and then subtract based on expected anomaly
        if int(Dataframe[TimeRangeCol][i].split('-')[0]) <= 1870:
            print(
                f"Time range not in HadSST series. Taking 1961-1990 average instead to then apply proxy-based anomalies from D'Arrigo")
            time_start = '1961-01-01T00:00'
            time_end = '1990-12-01T00:00'

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
            check_series = ncfile.sel(lat=out[close_point_idx][1][0],
                                      lon=out[close_point_idx][1][1],
                                      time=slice(time_start, time_end))

            if np.isnan(np.nanmean(check_series.variables['sst'][:,
                                   :])):  # get average temperature ignoring if there's nan in the series and so return true if there's at least one valid value
                print("Empty grid cell: Finding next valid point...")

            else:
                print(
                    f"found closest grid point with valid value {(out[close_point_idx][1][0], out[close_point_idx][1][1] - 180)}")
                print(f"original site has coords {ilat, ilon - 180}")
                print(f"time considered is FROM {time_start} TO {time_end} ")

                values = np.asarray(check_series.variables['sst'][:])

                # update values on dataframe
                Dataframe.loc[i, "SST_ERSSTv5_sd"] = np.std(values)
                Dataframe.loc[i, "SST_ERSSTv5_ann"] = np.mean(values)
                Dataframe.loc[i, "SST_ERSSTv5_min"] = values.min()
                Dataframe.loc[i, "SST_ERSSTv5_max"] = values.max()
                Dataframe.loc[i, "SST_ERSSTv5_LatGridCell"] = out[close_point_idx][1][0]
                Dataframe.loc[i, "SST_ERSSTv5_LonGridCell"] = out[close_point_idx][1][
                                                                  1] - 180  # revert longitude back for right comparisson
                Dataframe.at[i, "SST_ERSSTv5_TimeSpan"] = list((time_start, time_end))
                break

    return Dataframe


def extract_COADS(ncfille, Dataframe, TimeRangeCol, LatVarName, LonVarName, SlideLatVal_Up, SlideLonVal_Up,
                  SlideLatVal_Btm, SlideLonVal_Btm):
    for i in range(0, len(Dataframe)):

        latrange = (round(Dataframe[LatVarName][i] + SlideLatVal_Up),
                    round(Dataframe[LatVarName][i] + SlideLatVal_Btm))

        longrange = (
            round(Dataframe[LonVarName][i] + SlideLonVal_Up) + 180,  # because longitude in the model is gridded 0-360,
            round(Dataframe[LonVarName][i] + SlideLonVal_Btm) + 180)

        time_start = (Dataframe[TimeRangeCol][i]).split('-')[0] + '-01-01T00:00'
        time_end = (Dataframe[TimeRangeCol][i]).split('-')[1] + '-12-01T00:00'

        dummy = ncfille.sel(X=slice(longrange[0], longrange[1]), Y=slice(latrange[0], latrange[1]),
                            T=slice(time_start, time_end))

        # for regions with NaN SST and between 1800-1850, get the 1854-1854 value to then subtract from anomalies based on proxy studies
        if np.isnan(float(dummy.variables['sst'].mean())):
            dummy = ncfille.sel(X=slice(longrange[0], longrange[1]), Y=slice(latrange[0], latrange[1]),
                                T=slice('1854-01-01T00:00', '1860-12-01T00:00'))

        year_jump = 12
        values = []
        valmax = []
        valmin = []
        month = 0
        while month < int(len(dummy.variables['sst'])):
            values.append(float(dummy.variables['sst'][month:month + year_jump, :, :].mean()))
            valmin.append(float(dummy.variables['sst'][month:month + year_jump, :, :].min()))
            valmax.append(float(dummy.variables['sst'][month:month + year_jump, :, :].max()))
            month = month + year_jump

        # then a global mean between those years
        Dataframe.loc[i, "SST_COADS_sd"] = np.std(values)
        Dataframe.loc[i, "SST_COADS_ann"] = np.mean(values)
        Dataframe.loc[i, "SST_COADS_range"] = np.mean((np.asarray(valmax) - np.asarray(valmin)))

    return Dataframe


# DOWNLOAD ERSSTv5 files and check for updates
earthaccess.login()
BBOX = {'SW_lat': -20,
        'SW_lon': 100,
        'NE_lat': 10,
        'NE_lon': 150}
timechunks = [("1854-01-01T00:00", "1870-01-12T00:00:00"),
              ("1880-01-01T00:00", "1910-01-12T00:00:00"),
              ("1910-01-01T00:00", "1930-01-12T00:00:00"),
              ("1940-01-01T00:00", "1990-01-12T00:00:00"),
              ("1843-01-01T00:00", "1855-01-12T00:00:00"),
              ("1990-01-01T00:00", "1999-01-12T00:00:00"),
              ("1950-01-01T00:00", "2018-01-12T00:00:00"),
              ]

files_long = []
download_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/ERSSTv5'
# for pair in timechunks:
#     results = earthaccess.search_data(
#         short_name='REYNOLDS_NCDC_L4_MONTHLY_V5',
#         cloud_hosted=True,
#         bounding_box=(BBOX['SW_lon'],
#                       BBOX['SW_lat'],
#                       BBOX['NE_lon'],
#                       BBOX['NE_lat']),
#         temporal=pair,
#         count=2000)
#
#     files = earthaccess.download(results, download_path)
#     files_long.append(files)


##TODO reading ERSSTv5 data until 2007
ERSSTv5_until2007 = xr.open_mfdataset(paths=os.path.join(download_path, '*.nc'),
                                      engine='netcdf4', )  # combine multiple files into single file

# TODO reading ERSSTv5 data from 2008 because calendar is defined differently
download_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/ERSSTv5_2008_onwards'
ERSSTv5_from2008 = xr.open_mfdataset(paths=os.path.join(download_path, '*.nc'),
                                     engine='netcdf4', )  # combine multiple files into single file
ERSSTv5_from2008 = ERSSTv5_from2008.convert_calendar('360_day',
                                                     align_on='year')  # calendar needs conversion to be able to concat array of nc files
ERSSTv5 = xr.merge([ERSSTv5_until2007, ERSSTv5_from2008])

##TODO reading COADS data
coads_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/COADS_SST/COADS_data.nc'
coads = xr.open_dataset(filename_or_obj=coads_path, engine='netcdf4', decode_times=False)
COADS_fixed = fix_calendar(coads)
COADS_fixed = xr.decode_cf(COADS_fixed)

# TODO reading HadlSST data
hadlsst_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/HadlSST/HadISST_sst.nc'
HadlSST = xr.open_dataset(filename_or_obj=hadlsst_path, engine='netcdf4', decode_times=True)
HadlSST = HadlSST.convert_calendar('360_day', align_on='year')

#My Data
# ead spreadsheet and feed data into function to get SSTs
excel_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Tables_and_Regional_Sectors_Averages.xlsx'
Dataframe_MyData = pd.read_excel(excel_path, sheet_name="Table1")
# remove first row
Dataframe_MyData = Dataframe_MyData.drop([0]).reset_index()
TEST_1 = extract_HadSST(ncfile=HadlSST, Dataframe=Dataframe_MyData, TimeRangeCol='MGA year range')
outfile_name = os.path.join(os.path.dirname(excel_path), 'SST_Extracted_Results_MyData.xlsx')
TEST_1.to_excel(outfile_name)


# Literature review + Australian Corals
lit_rev_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 2 - Leo - General/Coral_QGIS_Spatial_Analyses/Revised_Literature_Points_Summarised.xlsx'
Dataframe_Literature = pd.read_excel(lit_rev_path, sheet_name="Aggregated_MapLayer")
Dataframe_Literature = Dataframe_Literature.reset_index()
TEST_3 = extract_HadSST(ncfile=HadlSST, Dataframe=Dataframe_Literature, TimeRangeCol='Year_range')
outfile_name = os.path.join(os.path.dirname(lit_rev_path), 'SST_Extracted_Results_LitReview.xlsx')
TEST_3.to_excel(outfile_name)

