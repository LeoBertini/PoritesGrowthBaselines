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


def extract_ERSSTv5(ncfile, Dataframe, TimeRangeCol, LatVarName, LonVarName, SlideLatVal_Up, SlideLonVal_Up,
                    SlideLatVal_Btm, SlideLonVal_Btm):

    for i in range(0, len(Dataframe)):

       ##todo change to function
       print(f"processing data entry row no. : {Dataframe['index'][i]}")

       ilat= Dataframe['Lat'][i]
       ilon = Dataframe['Lon'][i]+180 #because longitude is gridded 0-360,

       nclats=np.asarray(ncfile.variables['lat'][:])
       nclons=np.asarray(ncfile.variables['lon'][:])

       time_start = (Dataframe[TimeRangeCol][i]).split('-')[0] + '-01-01T00:00'
       time_end = (Dataframe[TimeRangeCol][i]).split('-')[1] + '-12-01T00:00'

       # if data point is too old for series. adopt earliest 5 years of ERSST record
       if int(Dataframe[TimeRangeCol][i].split('-')[0]) <= 1854:
           time_start = '1854-01-01T00:00'
           time_end = '1860-12-01T00:00'

        #calculating closest grid point to target location
        grid_point=[]
        d = []
        grid_point_idx = []
        for lat_index in range(0,len(nclats)):
            for lon_index in range(0, len(nclons)):
                grid_point.append((nclats[lat_index],nclons[lon_index]))
                grid_point_idx.append((lat_index,lon_index))
                d.append(geopy.distance.geodesic((ilat,ilon), (nclats[lat_index], nclons[lon_index])).km)

        out=list(zip(d,grid_point, grid_point_idx)) #zipped tuple of ((distance), (latgrid, longrid), (latgrid_idx, longrid_idx) )
        out.sort(key=lambda x: x[0]) #sorting so first is the closest point relative to my target point

        #now find valid grid cell with value closest to target location iterating over tuple of point distances
        for close_point_idx in range(0, len(out)):
            #slice data
            check_series = ncfile.sel(lat=out[close_point_idx][1][0],
                                      lon=out[close_point_idx][1][1],
                                      time=slice(time_start, time_end))

            if np.isnan(np.nanmean(check_series.variables['sst'][:,:])): #get average temperature ignoring if there's nan in the series and so return true if there's at least one valid value
                print("Empty grid cell: Finding next valid point")

            else:
                print(f"found closest grid point with valid value {(out[close_point_idx][1][0],out[close_point_idx][1][1]-180) }")
                print(f"original grid point is {ilat,ilon-180}")

                values = np.asarray(check_series.variables['sst'][:,:])

                # update values on dataframe
                Dataframe.loc[i, "SST_ERSSTv5_sd"] = np.std(values)
                Dataframe.loc[i, "SST_ERSSTv5_ann"] = np.mean(values)
                Dataframe.loc[i, "SST_ERSSTv5_range"] = values.max()-values.min()
                #todo get column for min
                #todo get column for max
                #todo get column for closes valid grid point coordinates
                break
       #######







        latrange = (round(Dataframe[LatVarName][i] + SlideLatVal_Up),
                    round(Dataframe[LatVarName][i] + SlideLatVal_Btm))

        longrange = (round(Dataframe[LonVarName][i] + SlideLonVal_Up)+180,
                     round(Dataframe[LonVarName][i] + SlideLonVal_Btm)+180)

        time_start = (Dataframe[TimeRangeCol][i]).split('-')[0] + '-01-01T00:00'
        time_end = (Dataframe[TimeRangeCol][i]).split('-')[1] + '-12-01T00:00'

        dummy = ncfille.sel(lon=slice(longrange[0], longrange[1]), lat=slice(latrange[0], latrange[1]),
                            time=slice(time_start, time_end), method='nearest')



        # for regions with NaN SST and between 1800-1850, get the 1854-1854 value to then subtract from anomalies based on proxy studies
        if np.isnan(float(dummy.variables['sst'].mean())):
            dummy = ncfille.sel(lon=slice(longrange[0], longrange[1]), lat=slice(latrange[0], latrange[1]),
                                time=slice('1854-01-01T00:00', '1860-12-01T00:00'))

        year_jump = 12
        values = []
        valmax = []
        valmin = []
        month = 0
        while month < int(len(dummy.variables['sst'])):
            values.append(float(dummy.variables['sst'][month:month + year_jump, :, :, :].mean()))
            valmin.append(float(dummy.variables['sst'][month:month + year_jump, :, :, :].min()))
            valmax.append(float(dummy.variables['sst'][month:month + year_jump, :, :, :].max()))
            month = month + year_jump

        # then a global mean between those years
        Dataframe.loc[i, "SST_ERSSTv5_sd"] = np.std(values)
        Dataframe.loc[i, "SST_ERSSTv5_ann"] = np.mean(values)
        Dataframe.loc[i, "SST_ERSSTv5_range"] = np.mean((np.asarray(valmax) - np.asarray(valmin)))

    return Dataframe


def extract_COADS(ncfille, Dataframe, TimeRangeCol, LatVarName, LonVarName, SlideLatVal_Up, SlideLonVal_Up,
                  SlideLatVal_Btm, SlideLonVal_Btm):
    for i in range(0, len(Dataframe)):

        latrange = (round(Dataframe[LatVarName][i] + SlideLatVal_Up),
                    round(Dataframe[LatVarName][i] + SlideLatVal_Btm))

        longrange = (round(Dataframe[LonVarName][i] + SlideLonVal_Up)+180, #because longitude in the model is gridded 0-360,
                     round(Dataframe[LonVarName][i] + SlideLonVal_Btm)+180)

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
ERSSTv5_from2008 = ERSSTv5_from2008.convert_calendar('360_day', align_on='year') #calendar needs conversion to be able to concat array of nc files

ERSSTv5 = xr.merge([ERSSTv5_until2007,ERSSTv5_from2008])


##TODO reading COADS data
coads_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/COADS_SST/COADS_data.nc'
coads = xr.open_dataset(filename_or_obj=coads_path, engine='netcdf4', decode_times=False)
COADS_fixed = fix_calendar(coads)
COADS_fixed = xr.decode_cf(COADS_fixed)


#TODO reading HadlSST data



##Indonesian samples
# TODO read spreadsheet and feed data into function to get averages
excel_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Tables_and_Regional_Sectors_Averages.xlsx'
Dataframe_MyData = pd.read_excel(excel_path, sheet_name="Table1")
# remove first row
Dataframe_MyData = Dataframe_MyData.drop([0]).reset_index()


#TODO find closest available gridcell where value is not NaN

Dataframe_MyData = extract_ERSSTv5(ncfille=ERSSTv5, Dataframe=Dataframe_MyData, TimeRangeCol='MGA year range',
                                   LatVarName='Lat', LonVarName='Lon',
                                   SlideLatVal_Up=-1, SlideLonVal_Up=-1,
                                   SlideLatVal_Btm=0, SlideLonVal_Btm=0)

Dataframe_MyData = extract_COADS(ncfille=COADS_fixed, Dataframe=Dataframe_MyData, TimeRangeCol='MGA year range',
                                 LatVarName='Lat', LonVarName='Lon',
                                 SlideLatVal_Up=-2, SlideLonVal_Up=-2,
                                 SlideLatVal_Btm=0, SlideLonVal_Btm=0)

# get nan rows to try again expanding the grid search
Dataframe_MyData_1 = Dataframe_MyData.loc[Dataframe_MyData['SST_COADS_ann'].notnull()]
NotfoundDF_1 = Dataframe_MyData.loc[Dataframe_MyData['SST_COADS_ann'].isnull()]
if not NotfoundDF_1.empty:
    NotfoundDF_1 = NotfoundDF_1.reset_index(drop=True)
    print('expanding grid search for COADS')
    NotfoundDF_1 = extract_COADS(ncfille=COADS_fixed, Dataframe=NotfoundDF_1, TimeRangeCol='MGA year range',
                                 LatVarName='Lat', LonVarName='Lon',
                                 SlideLatVal_Up=-1, SlideLonVal_Up=-1,
                                 SlideLatVal_Btm=2, SlideLonVal_Btm=2)

NotfoundDF_2 = NotfoundDF_1.loc[NotfoundDF_1['SST_COADS_ann'].isnull()]
Dataframe_MyData_2 = NotfoundDF_1.loc[NotfoundDF_1['SST_COADS_ann'].notnull()]
if not NotfoundDF_2.empty:
    NotfoundDF_2 = NotfoundDF_2.reset_index(drop=True)
    print('expanding grid search for COADS')
    NotfoundDF_2 = extract_COADS(ncfille=COADS_fixed, Dataframe=NotfoundDF_2, TimeRangeCol='MGA year range',
                                 LatVarName='Lat', LonVarName='Lon',
                                 SlideLatVal_Up=-3, SlideLonVal_Up=-3,
                                 SlideLatVal_Btm=3, SlideLonVal_Btm=3)

Dataframe_MyData_3 = NotfoundDF_2.loc[NotfoundDF_2['SST_COADS_ann'].notnull()]
NotfoundDF_3 = NotfoundDF_2.loc[NotfoundDF_2['SST_COADS_ann'].isnull()]

Final_DF_MyData = pd.concat([Dataframe_MyData_1, Dataframe_MyData_2, Dataframe_MyData_3, NotfoundDF_3])

# saving results to file
file_name_results = os.path.join(os.path.dirname(excel_path), 'MyDATA_SST_extracted_results.xlsx')
Final_DF_MyData.to_excel(file_name_results, index=False)

# Australian Corals
excel_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 2 - Leo - General/Coral_QGIS_Spatial_Analyses/Australian_Coral_Growth_Datasets.xlsx'
Dataframe_Australia = pd.read_excel(excel_path)

Dataframe_Australia = extract_ERSSTv5(ncfille=ERSSTv5, Dataframe=Dataframe_Australia, TimeRangeCol='Year_range',
                                      LatVarName='Lat', LonVarName='Lon',
                                      SlideLatVal_Up=-1, SlideLonVal_Up=-1,
                                      SlideLatVal_Btm=0, SlideLonVal_Btm=0)

Dataframe_Australia_1 = Dataframe_Australia.loc[Dataframe_Australia['SST_ERSSTv5_ann'].notnull()]
Australia_notfound1 = Dataframe_Australia.loc[Dataframe_Australia['SST_ERSSTv5_ann'].isnull()].reset_index(drop=True)
# expand grid search
Australia_notfound1 = extract_ERSSTv5(ncfille=ERSSTv5, Dataframe=Australia_notfound1, TimeRangeCol='Year_range',
                                      LatVarName='Lat', LonVarName='Lon',
                                      SlideLatVal_Up=-1, SlideLonVal_Up=-1,
                                      SlideLatVal_Btm=1, SlideLonVal_Btm=1)

Dataframe_Australia_2 = Australia_notfound1.loc[Australia_notfound1['SST_ERSSTv5_ann'].notnull()]
Australia_notfound2 = Australia_notfound1.loc[Australia_notfound1['SST_ERSSTv5_ann'].isnull()].reset_index(drop=True)

Australia_notfound2 = extract_ERSSTv5(ncfille=ERSSTv5, Dataframe=Australia_notfound2, TimeRangeCol='Year_range',
                                      LatVarName='Lat', LonVarName='Lon',
                                      SlideLatVal_Up=-1, SlideLonVal_Up=-1,
                                      SlideLatVal_Btm=2, SlideLonVal_Btm=2)

Dataframe_Australia = pd.concat([Dataframe_Australia_1, Dataframe_Australia_2, Australia_notfound2]).reset_index(
    drop=True)  # all ERSSTv5 extracted

Dataframe_Australia = extract_COADS(ncfille=COADS_fixed, Dataframe=Dataframe_Australia, TimeRangeCol='Year_range',
                                    LatVarName='Lat', LonVarName='Lon',
                                    SlideLatVal_Up=-2, SlideLonVal_Up=-2,
                                    SlideLatVal_Btm=0, SlideLonVal_Btm=0)

AustraliaCOADS_1 = Dataframe_Australia.loc[Dataframe_Australia['SST_COADS_ann'].notnull()]
AustraliaCOADS_notfound_1 = Dataframe_Australia.loc[Dataframe_Australia['SST_COADS_ann'].isnull()].reset_index(
    drop=True)

AustraliaCOADS_notfound_1 = extract_COADS(ncfille=COADS_fixed, Dataframe=AustraliaCOADS_notfound_1,
                                          TimeRangeCol='Year_range',
                                          LatVarName='Lat', LonVarName='Lon',
                                          SlideLatVal_Up=-2, SlideLonVal_Up=-2,
                                          SlideLatVal_Btm=1, SlideLonVal_Btm=1)

Dataframe_Australia = pd.concat([AustraliaCOADS_1, AustraliaCOADS_notfound_1])

# saving results to file
file_name_results = os.path.join(os.path.dirname(excel_path), 'Australia_SST_extracted_results.xlsx')
Dataframe_Australia.to_excel(file_name_results, index=False)

# literature review SST extraction
lit_rev_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 2 - Leo - General/Coral_QGIS_Spatial_Analyses/Revised_Literature_Points_Summarised.xlsx'
Dataframe_Literature = pd.read_excel(lit_rev_path, sheet_name="Revised_Literature_Points")

try:
    Dataframe_Literature = extract_ERSSTv5(ncfille=ERSSTv5, Dataframe=Dataframe_Literature, TimeRangeCol='Year_range',
                                          LatVarName='Lat', LonVarName='Lon',
                                          SlideLatVal_Up=-1, SlideLonVal_Up=-1,
                                          SlideLatVal_Btm=0, SlideLonVal_Btm=0)
finally:
    file_name_results = os.path.join(os.path.dirname(lit_rev_path), 'LitReview_SST_extracted_results.xlsx')

Dataframe_Literature_1 = Dataframe_Literature.loc[Dataframe_Literature['SST_ERSSTv5_ann'].notnull()]

Dataframe_Literature_notfound1 = Dataframe_Literature.loc[Dataframe_Literature['SST_ERSSTv5_ann'].isnull()].reset_index(drop=True)

try:
    Dataframe_Literature_notfound1 = extract_ERSSTv5(ncfille=ERSSTv5, Dataframe=Dataframe_Literature_notfound1, TimeRangeCol='Year_range',
                                          LatVarName='Lat', LonVarName='Lon',
                                          SlideLatVal_Up=-2, SlideLonVal_Up=-2,
                                          SlideLatVal_Btm=2, SlideLonVal_Btm=2)
finally:
    print('Check file')


Dataframe_Literature.to_excel(file_name_results, index=False)
