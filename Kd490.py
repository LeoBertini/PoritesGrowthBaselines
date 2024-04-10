import os
import time

import earthaccess
import cdsapi
import numpy as np
import xarray as xr
import pandas as pd
import geopy.distance
import requests
import matplotlib.pyplot as plt


def extract_Kd490(ncfile, Dataframe, TimeRange, varname):
    # initiating new columns to append data to
    Dataframe["Kd490_mean"] = np.nan
    Dataframe["Kd490_sd"] = np.nan
    Dataframe["Kd490_LatGridCell"] = np.nan
    Dataframe["Kd490_LonGridCell"] = np.nan
    Dataframe['Kd490_TimeSpan'] = np.nan
    Dataframe['Kd490_TimeSpan'] = Dataframe['Kd490_TimeSpan'].astype(object)
    Dataframe['Kd490_Series'] = np.nan
    Dataframe['Kd490_Series'] = Dataframe['Kd490_Series'].astype(object)

    for i in range(0, len(Dataframe)):

        ##todo change to function
        print(f"processing data entry row no. : {Dataframe['index'][i]}")

        ilat = Dataframe['Lat'][i]
        ilon = Dataframe['Lon'][i]

        nclats = np.asarray(ncfile.variables['lat'][:])
        nclons = np.asarray(ncfile.variables['lon'][:])

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
                                      lon=out[close_point_idx][1][1])

            if np.isnan(np.nanmean(check_series.variables[varname][
                                   :])):  # get average ignoring if there's nan in the series and so return true if there's at least one valid value
                print("Empty grid cell: Finding next valid point...")

            else:
                print(
                    f"found closest grid point with valid value {(out[close_point_idx][1][0], out[close_point_idx][1][1])}")
                print(f"original site has coords {ilat, ilon}")

                values = np.asarray(check_series.variables[varname][:])

                # update values on dataframe
                Dataframe.loc[i, "Kd490_mean"] = np.nanstd(values)
                Dataframe.loc[i, "Kd490_sd"] = np.nanmean(values)
                Dataframe.loc[i, "Kd490_LatGridCell"] = out[close_point_idx][1][0]
                Dataframe.loc[i, "Kd490_LonGridCell"] = out[close_point_idx][1][1]
                Dataframe.at[i, 'Kd490_TimeSpan'] = list(TimeRange)
                Dataframe.at[i, 'Kd490_Series'] = list(values)

                break

    return Dataframe


text_file = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Kd490/Kd490_download_links.txt'
with open(text_file, 'r') as file:
    # Create an empty list to store the lines
    links = []

    # Iterate over the lines of the file
    for line in file:
        # Remove the newline character at the end of the line
        line = line.strip()

        # Append the line to the list
        links.append(line)

# download files
download_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Kd490/'
for link in links:
    file_name = link.split('/')[-1]
    file_path = os.path.join(download_path, file_name)
    if not os.path.isfile(file_path):
        r = requests.get(link, allow_redirects=True)
        open(file_path, 'wb').write(r.content)
    else:
        print('The following file was already downloaded')
        print(file_path)

files = [os.path.join(download_path, item) for item in os.listdir(download_path) if item.startswith('AQUA_MODIS')]
files.sort()

Kd490 = xr.open_mfdataset(paths=files, engine='netcdf4', concat_dim="Time",
                          combine="nested")  # combine multiple files into single file
varname = 'Kd_490'
KD_dataset = Kd490[[varname]]

time_start = files[0].split("/")[-1].split('.')[1].split('_')[0]
time_end = files[-1].split("/")[-1].split('.')[1].split('_')[0]

########## extracting closest Kd490 from my sites
lit_rev_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 3 - Results Chapter - Leo - General/Coral_QGIS_Spatial_Analyses/Revised_Literature_Points_Summarised.xlsx'
Dataframe_Literature = pd.read_excel(lit_rev_path, sheet_name="Aggregated_MapLayer")
Dataframe_Literature = Dataframe_Literature.reset_index()
TEST1 = extract_Kd490(ncfile=KD_dataset, varname='Kd_490',
                      Dataframe=Dataframe_Literature, TimeRange=(time_start,time_end))
outfile_name = os.path.join(os.path.dirname(lit_rev_path), 'Kd490_Extracted_Results_LitReview.xlsx')
TEST1.to_excel(outfile_name)


#My Data
excel_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Tables_and_Regional_Sectors_Averages.xlsx'
Dataframe_MyData = pd.read_excel(excel_path, sheet_name="Table1")
# remove first row
Dataframe_MyData = Dataframe_MyData.drop([0]).reset_index()
TEST2 = extract_Kd490(ncfile=KD_dataset,varname='Kd_490',
                   Dataframe=Dataframe_MyData, TimeRange=(time_start,time_end))
outfile_name = os.path.join(os.path.dirname(excel_path), 'Kd490_Extracted_Results_MyData.xlsx')
TEST2.to_excel(outfile_name)

# #vizualising my data
# example_plot=KD_dataset[varname].isel(Time=0).plot()

