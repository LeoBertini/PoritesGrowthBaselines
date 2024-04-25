import os
import cdsapi
import numpy as np
import xarray as xr
import pandas as pd
import geopy.distance



def extract_DSR(ncfile, Dataframe, TimeRangeCol, varname):
    # initiating new columns to append data to
    Dataframe["DSR_ERB_mean"] = np.nan
    Dataframe["DSR_ERB_sd"] = np.nan
    Dataframe["DSR_ERB_LatGridCell"] = np.nan
    Dataframe["DSR_ERB_LonGridCell"] = np.nan
    Dataframe['DSR_ERB_TimeSpan'] = np.nan
    Dataframe['DSR_ERB_TimeSpan'] = Dataframe['DSR_ERB_TimeSpan'].astype(object)
    Dataframe['DSR_ERB_Series'] = np.nan
    Dataframe['DSR_ERB_Series'] = Dataframe['DSR_ERB_Series'].astype(object)

    for i in range(0, len(Dataframe)):

        ##todo change to function
        print(f"processing data entry row no. : {Dataframe['index'][i]}")

        ilat = Dataframe['Lat'][i]
        ilon = Dataframe['Lon'][i]

        nclats = np.asarray(ncfile.variables['lat'][:])
        nclons = np.asarray(ncfile.variables['lon'][:])

        time_start = '1979-01-01T00:00:00.000000000'
        time_end = '2020-12-01T00:00:00.000000000'

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
                Dataframe.loc[i, "DSR_ERB_sd"] = np.nanstd(values)
                Dataframe.loc[i, "DSR_ERB_mean"] = np.nanmean(values)
                Dataframe.loc[i, "DSR_ERB_LatGridCell"] = out[close_point_idx][1][0]
                Dataframe.loc[i, "DSR_ERB_LonGridCell"] = out[close_point_idx][1][1]
                Dataframe.at[i, "DSR_ERB_TimeSpan"] = list((time_start, time_end))
                Dataframe.at[i, 'DSR_ERB_Series']=list(values)
                break

    return Dataframe

#download data
#
# #copernicus dataset download
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


download_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/DSR_Flux_ERB/'
files = [os.path.join(download_path, item) for item in os.listdir(download_path) if item.startswith('SNS')]
files.sort()

ERB = xr.open_mfdataset(paths=files, engine='netcdf4', )  # combine multiple files into single file
varname = 'SNS'
DSR_dataset = ERB[[varname]]
DSR = DSR_dataset.variables[varname]


########## extracting closest DSR from my sites

lit_rev_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 3 - Results Chapter - Leo - General/Coral_QGIS_Spatial_Analyses/Revised_Literature_Points_Summarised.xlsx'
Dataframe_Literature = pd.read_excel(lit_rev_path, sheet_name="Aggregated_MapLayer")
Dataframe_Literature = Dataframe_Literature.reset_index()
TEST1 = extract_DSR(ncfile=DSR_dataset,varname='SNS',
                   Dataframe=Dataframe_Literature, TimeRangeCol='Year_range')
outfile_name = os.path.join(os.path.dirname(lit_rev_path), 'DSR_ERB_Extracted_Results_LitReview.xlsx')
TEST1.to_excel(outfile_name)

#My Data
excel_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Tables_and_Regional_Sectors_Averages.xlsx'
Dataframe_MyData = pd.read_excel(excel_path, sheet_name="Table 2")
# remove first row
Dataframe_MyData = Dataframe_MyData.drop([0]).reset_index()
TEST2 = extract_DSR(ncfile=DSR_dataset,varname='SNS',
                   Dataframe=Dataframe_MyData, TimeRangeCol='Year_range')
outfile_name = os.path.join(os.path.dirname(excel_path), 'DSR_ERB_Extracted_Results_MyData.xlsx')
TEST2.to_excel(outfile_name)


