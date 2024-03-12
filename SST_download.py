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

def fix_calendar(ds, timevar='T'):
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds



#DOWNLOAD ERSSTv5 files and check for updates
earthaccess.login()
BBOX = {'SW_lat': -20,
        'SW_lon': 100,
        'NE_lat' : 10,
        'NE_lon' : 150}
timechunks = [("1854-01-01T00:00", "1870-01-12T00:00:00"),
              ("1880-01-01T00:00", "1910-01-12T00:00:00"),
              ("1910-01-01T00:00", "1930-01-12T00:00:00"),
              ("1940-01-01T00:00", "1990-01-12T00:00:00"),
              ("1843-01-01T00:00", "1855-01-12T00:00:00"),
              ("1990-01-01T00:00", "1999-01-12T00:00:00"),
              ("1950-01-01T00:00", "2008-01-12T00:00:00"),
              ]

files_long = []
download_path='/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/SST_downloads'
for pair in timechunks:
    results = earthaccess.search_data(
        short_name='REYNOLDS_NCDC_L4_MONTHLY_V5',
        cloud_hosted=True,
        bounding_box = (BBOX['SW_lon'],
                        BBOX['SW_lat'],
                        BBOX['NE_lon'],
                        BBOX['NE_lat']),
        temporal=pair,
        count=2000)

    files = earthaccess.download(results, download_path)
    files_long.append(files)

ds2 = xr.open_mfdataset(paths=os.path.join(download_path,'*.nc'), engine='netcdf4') #combine multiple files into single file

## getting sst for my Indonesian samples
#TODO read spreadsheet and feed data into function to get averages
excel_path='/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Tables_and_Regional_Sectors_Averages.xlsx'
Dataframe=pd.read_excel(excel_path, sheet_name="Table1")
#remove first row
Dataframe=Dataframe.drop([0]).reset_index()

for i in range(0, len(Dataframe)):
    latrange =(round(Dataframe['Lat'][i]-1),
               round(Dataframe['Lat'][i]))

    longrange =(round(Dataframe['Lon'][i]-1),
                round(Dataframe['Lon'][i]))

    time_start = (Dataframe['MGA year range'][i]).split('-')[0]+'-01-01T00:00'
    time_end = (Dataframe['MGA year range'][i]).split('-')[1]+'-01-01T00:00'

    dummy = ds2.sel(lon = slice(longrange[0],longrange[1]), lat = slice(latrange[0], latrange[1]), time=slice(time_start,time_end))

    # for regions with NaN SST and between 1800-1850, get the 1854-1854 value to then subtract from anomalies based on proxy studies
    if np.isnan(float(dummy.variables['sst'].mean())):
        dummy = ds2.sel(lon = slice(longrange[0],longrange[1]), lat = slice(latrange[0], latrange[1]), time=slice('1854-01-01T00:00','1854-12-01T00:00'))

    Dataframe.loc[i, "SST_stdv"] = float(dummy.variables['sst'].std())
    Dataframe.loc[i, "SST_recon"] = float(dummy.variables['sst'].mean())
    print(float(dummy.variables['sst'].std()))
    print(float(dummy.variables['sst'].mean()))

#saving results to file
file_name_results=os.path.join(os.path.dirname(excel_path), 'SST_extracted_results.xlsx')
Dataframe.to_excel(file_name_results, index=False )


#getting sst for Australian sites
excel_path='/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 2 - Leo - General/Coral_QGIS_Spatial_Analyses/Australian_Coral_Growth_Datasets.xlsx'
Dataframe_Australia=pd.read_excel(excel_path)

#FIRST ERSSTv5
for i in range(0, len(Dataframe_Australia)):
    latrange =(round(Dataframe_Australia['Lat'][i]-1),
               round(Dataframe_Australia['Lat'][i]+1))

    longrange =(round(Dataframe_Australia['Lon'][i]-1), #further from coast to avoid falling into blank grid point
                round(Dataframe_Australia['Lon'][i]+1))

    time_start = (Dataframe_Australia['Year_range'][i]).split('-')[0]+'-01-01T00:00'
    time_end = (Dataframe_Australia['Year_range'][i]).split('-')[1]+'-12-01T00:00'

    dummy = ds2.sel(lon = slice(longrange[0],longrange[1]), lat = slice(latrange[0], latrange[1]), time=slice(time_start,time_end))

    #TODO-COMPLETE get annual mean per year

    year_jump=12
    values=[]
    valmax=[]
    valmin=[]
    for year in range(0, int(len(dummy.variables['sst'])/year_jump)):
        year=year+year+year_jump+1
        values.append(float(dummy.variables['sst'][year:year+year_jump,:,:,:].mean()))
        valmin.append(float(dummy.variables['sst'][year:year+year_jump,:,:,:].min()))
        valmax.append(float(dummy.variables['sst'][year:year + year_jump, :, :, :].max()))

    # then a global mean between those years
    Dataframe_Australia.loc[i, "SST_ERSSTv5_sd"] = np.std(values)
    Dataframe_Australia.loc[i, "SST_ERSSTv5_ann"] = np.mean(values)
    Dataframe_Australia.loc[i, "SST_ERSSTv5_range"] = np.mean((np.asarray(valmax)-np.asarray(valmin)))


#TODO reading COADS data
coads_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/COADS_SST/COADS_data.nc'
coads= xr.open_dataset(filename_or_obj=coads_path, engine='netcdf4', decode_times=False)

coads_fxed = fix_calendar(coads)
coads_fxed = xr.decode_cf(coads_fxed)

Dataframe_Australia["SST_COADS_sd"] = np.zeros(len(Dataframe_Australia))
Dataframe_Australia["SST_COADS_ann"] = np.zeros(len(Dataframe_Australia))
Dataframe_Australia["SST_COADS_range"] = np.zeros(len(Dataframe_Australia))


for i in range(0, len(Dataframe_Australia)):
    latrange =(round(Dataframe_Australia['Lat'][i]-2),
               round(Dataframe_Australia['Lat'][i]))

    longrange =(round(Dataframe_Australia['Lon'][i]-2), #further from coast to avoid falling into blank grid point
                round(Dataframe_Australia['Lon'][i]))

    time_start = (Dataframe_Australia['Year_range'][i]).split('-')[0]+'-01-01T00:00'
    time_end = (Dataframe_Australia['Year_range'][i]).split('-')[1]+'-12-01T00:00'

    dummy = coads_fxed.sel(X = slice(longrange[0],longrange[1]), Y = slice(latrange[0], latrange[1]), T=slice(time_start,time_end))

    year_jump=12
    values=[]
    valmax=[]
    valmin=[]
    for year in range(0, int(len(dummy.variables['sst'])/year_jump)):
        year=year+year+year_jump+1
        values.append(float(dummy.variables['sst'][year:year+year_jump,:,:].mean()))
        valmin.append(float(dummy.variables['sst'][year:year+year_jump,:,:].min()))
        valmax.append(float(dummy.variables['sst'][year:year + year_jump, :, :].max()))

    # then a global mean between those years
    Dataframe_Australia.loc[i, "SST_COADS_sd_sd"] = np.std(values)
    Dataframe_Australia.loc[i, "SST_COADS_sd_ann"] = np.mean(values)
    Dataframe_Australia.loc[i, "SST_COADS_sd_range"] = np.mean((np.asarray(valmax)-np.asarray(valmin)))


#saving extracted results
file_name_results=os.path.join(os.path.dirname(excel_path), 'SST_extracted_results_Australia.xlsx')
Dataframe_Australia.to_excel(file_name_results, index=False )


