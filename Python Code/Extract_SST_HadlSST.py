import os
import numpy as np
import xarray as xr
import pandas as pd
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


# TODO reading HadlSST data
hadlsst_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/HadlSST/HadISST_sst.nc'
HadlSST = xr.open_dataset(filename_or_obj=hadlsst_path, engine='netcdf4', decode_times=True)
HadlSST = HadlSST.convert_calendar('360_day', align_on='year')

#My Data
# read spreadsheet and feed data into function to get SSTs
excel_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Tables_and_Regional_Sectors_Averages.xlsx'
Dataframe_MyData = pd.read_excel(excel_path, sheet_name="Table 2")
# remove first row
Dataframe_MyData = Dataframe_MyData.drop([0]).reset_index()
TEST_1 = extract_HadSST(ncfile=HadlSST, Dataframe=Dataframe_MyData, TimeRangeCol='MGA year range')
outfile_name = os.path.join(os.path.dirname(excel_path), 'SST_Extracted_Results_MuseumSpecimens.xlsx')
TEST_1.to_excel(outfile_name)


# Literature review
lit_rev_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 3 - Results Chapter - Leo - General/Coral_QGIS_Spatial_Analyses/Revised_Literature_Points_Summarised.xlsx'
Dataframe_Literature = pd.read_excel(lit_rev_path, sheet_name="To_Extract")
Dataframe_Literature = Dataframe_Literature.reset_index()
TEST_3 = extract_HadSST(ncfile=HadlSST, Dataframe=Dataframe_Literature, TimeRangeCol='Year_range')
outfile_name = os.path.join(os.path.dirname(lit_rev_path), 'SST_Extracted_Results_LitReview_SouthChinaSea.xlsx')
TEST_3.to_excel(outfile_name)

