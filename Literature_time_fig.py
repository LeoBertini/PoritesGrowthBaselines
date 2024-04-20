import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Tables_and_Regional_Sectors_Averages.xlsx'
litedata = pd.read_excel(path, sheet_name="Table1_Literature")

dict = {'Year': list(range(1850, 2019)),
        'Count': np.zeros(len(range(1850, 2019))),
        'Ext': np.zeros(len(range(1850, 2019))),
        'Den': np.zeros(len(range(1850, 2019))),
        'Cal': np.zeros(len(range(1850, 2019)))}

COUNT_TABLE = pd.DataFrame(dict)
COUNT_TABLE['Ext'] = COUNT_TABLE['Ext'].replace(0, np.nan)
COUNT_TABLE['Den'] = COUNT_TABLE['Den'].replace(0, np.nan)
COUNT_TABLE['Cal'] = COUNT_TABLE['Cal'].replace(0, np.nan)
COUNT_TABLE['Ext'] = COUNT_TABLE['Ext'].astype(object)
COUNT_TABLE['Den'] = COUNT_TABLE['Den'].astype(object)
COUNT_TABLE['Cal'] = COUNT_TABLE['Cal'].astype(object)

for i in range(0, len(litedata)):
    text = litedata['Record year range'][i]
    y1 = pd.eval(text.split("-")[0])
    y2 = pd.eval(text.split("-")[1])
    ext = litedata['Extension'][i]
    cal = litedata['Calcification'][i]
    den = litedata['Density'][i]

    for j in range(0, len(COUNT_TABLE)):
        if y1 <= COUNT_TABLE['Year'][j] <= y2:
            COUNT_TABLE['Count'][j] = COUNT_TABLE['Count'][j] + 1

            Ext_list = np.asarray(COUNT_TABLE['Ext'][j]).flatten().tolist()
            Ext_list.append(ext)
            COUNT_TABLE.at[j, 'Ext'] = Ext_list

            Den_list = np.asarray(COUNT_TABLE['Den'][j]).flatten().tolist()
            Den_list.append(den)
            COUNT_TABLE.at[j, 'Den'] = Den_list

            Cal_list = np.asarray(COUNT_TABLE['Cal'][j]).flatten().tolist()
            Cal_list.append(cal)
            COUNT_TABLE.at[j, 'Cal'] = Cal_list

Ext_avg = []
Den_avg = []
Cal_avg = []
for j in range(0, len(COUNT_TABLE)):
    e = np.nanmean(COUNT_TABLE['Ext'][j])
    Ext_avg.append(e)
    d = np.nanmean(COUNT_TABLE['Den'][j])
    Den_avg.append(d)
    c = np.nanmean(COUNT_TABLE['Cal'][j])
    Cal_avg.append(c)

COUNT_TABLE['Ext_avg'] = Ext_avg
COUNT_TABLE['Den_avg'] = Den_avg
COUNT_TABLE['Cal_avg'] = Cal_avg

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.70)

twin1 = ax.twinx()
twin2 = ax.twinx()
twin3 = ax.twinx()

# Offset the right spine of twin2.  The ticks and label have already been
# placed on the right by twinx above.
twin1.spines.right.set_position(("axes", 1.01))
twin2.spines.right.set_position(("axes", 1.15))
twin3.spines.right.set_position(("axes", 1.35))

p1, = ax.plot(COUNT_TABLE['Year'].tolist(), COUNT_TABLE['Count'].tolist(), "C0", label="Number of Observations")
p2, = twin1.plot(COUNT_TABLE['Year'].tolist(), COUNT_TABLE['Ext_avg'].tolist(), "C1", label="Extension")
p3, = twin2.plot(COUNT_TABLE['Year'].tolist(), COUNT_TABLE['Den_avg'].tolist(), "C2", label="Density")
p4, = twin3.plot(COUNT_TABLE['Year'].tolist(), COUNT_TABLE['Cal_avg'].tolist(), "C4", label="Calcification")

ax.set(xlim=(1850, 2018), ylim=(0, 150), xlabel="Record Year", ylabel="Number of observations")
twin1.set(ylim=(0, 30), ylabel="Extension")
twin2.set(ylim=(0.5, 2.5), ylabel="Density")
twin3.set(ylim=(1, 2), ylabel="Calcification")


ax.yaxis.label.set_color(p1.get_color())
twin1.yaxis.label.set_color(p2.get_color())
twin2.yaxis.label.set_color(p3.get_color())
twin3.yaxis.label.set_color(p4.get_color())

ax.tick_params(axis='y', colors=p1.get_color())
twin1.tick_params(axis='y', colors=p2.get_color())
twin2.tick_params(axis='y', colors=p3.get_color())
twin3.tick_params(axis='y', colors=p4.get_color())

ax.legend(handles=[p1, p2, p3, p4])

plt.show()
