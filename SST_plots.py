import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns
import os
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns
from scipy import stats


def reg_coef(x,y,label=None,color=None,**kwargs):
    ax = plt.gca()
    r,p = pearsonr(x,y)
    slope, intercept, r_value, pv, se = stats.linregress(x,y)
    ax.annotate(f" r = {(round(r_value,2))} | Y = {round(slope,3)} * X + {round(intercept,2)}", xy=(0.5,0.5), xycoords='axes fraction', ha='center')
    ax.set_axis_off()


#scatter plot matrix
# importing

excel_path='/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 2 - Leo - General/Coral_QGIS_Spatial_Analyses/Australian_Coral_Growth_Datasets.xlsx'
Dataframe_Australia=pd.read_excel(excel_path)

Dataframe_Australia = Dataframe_Australia.filter(items=['ReefProvince','Calcification','SST_HadlSST_ann','SST_COADS_ann','SST_ERSSTv5_ann'])
#Removing NaNs from one column where GBR lacks HAdlSST
Dataframe_Australia = Dataframe_Australia[Dataframe_Australia['SST_HadlSST_ann'].notna()]
pd.plotting.scatter_matrix(Dataframe_Australia, figsize=(10,10))
plt.show()



#filter only WA corals
WA_DF = Dataframe_Australia.loc[Dataframe_Australia['ReefProvince'] == 'WA']
pd.plotting.scatter_matrix(WA_DF, figsize=(10,10))
plt.show()

#doing regressions
g = sns.PairGrid(Dataframe_Australia)
g.map_diag(sns.distplot)
g.map_lower(sns.regplot)
g.map_upper(reg_coef)
fig_name=os.path.join(os.path.dirname(excel_path),'scatter_plot_fig.png')
g.savefig(fname=fig_name, dpi=300)

g = sns.PairGrid(WA_DF)
g.map_diag(sns.distplot)
g.map_lower(sns.regplot)
g.map_upper(reg_coef)
fig_name=os.path.join(os.path.dirname(excel_path),'scatter_plot_fig_WA.png')
g.savefig(fname=fig_name, dpi=300)

#plot map
map = Basemap(projection='merc',
              llcrnrlon=100+10,
              llcrnrlat=-20-10,
              urcrnrlon=150+10,
              urcrnrlat=10+10,
              resolution='i') # projection, lat/lon extents and resolution of polygons to draw
map.drawcoastlines()
#parallels = np.arange(-6,0, 1) # make latitude lines ever 5 degrees from 30N-50N
#meridians = np.arange(100, 110, 1) # make longitude lines every 5 degrees from 95W to 70W
#map.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
#map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
map.drawlsmask(land_color='Linen', ocean_color='#ffffff') # can use HTML names or codes for colors
lons,lats= np.meshgrid(lon,lat) #for this dataset, longitude is 0 through 360, so you need to subtract 180 to properly display on map
x,y = map(lons,lats)
temp = map.contourf(x,y,sst[0,0,:,:])
cb = map.colorbar(temp,"bottom", size="5%", pad="2%")
plt.title('SST')
cb.set_label('SST (degC)')
plt.show()
