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
Ext_sd_upper = []
Ext_sd_lower = []

Den_avg = []
Den_sd_upper = []
Den_sd_lower = []

Cal_avg = []
Cal_sd_upper = []
Cal_sd_lower = []

for j in range(0, len(COUNT_TABLE)):
    e = np.nanmean(COUNT_TABLE['Ext'][j])
    e_sd = np.nanstd(COUNT_TABLE['Ext'][j])
    Ext_avg.append(e)
    Ext_sd_upper.append(e + e_sd)
    Ext_sd_lower.append(e - e_sd)

    d = np.nanmean(COUNT_TABLE['Den'][j])
    d_sd = np.nanstd(COUNT_TABLE['Den'][j])
    Den_avg.append(d)
    Den_sd_upper.append(d + d_sd)
    Den_sd_lower.append(d - d_sd)

    c = np.nanmean(COUNT_TABLE['Cal'][j])
    c_sd = np.nanstd(COUNT_TABLE['Cal'][j])
    Cal_avg.append(c)
    Cal_sd_upper.append(c + c_sd)
    Cal_sd_lower.append(c - c_sd)

COUNT_TABLE['Ext_avg'] = Ext_avg
COUNT_TABLE['Den_avg'] = Den_avg
COUNT_TABLE['Cal_avg'] = Cal_avg
COUNT_TABLE['Ext_upper'] = Ext_sd_upper
COUNT_TABLE['Ext_lower'] = Ext_sd_lower
COUNT_TABLE['Den_upper'] = Den_sd_upper
COUNT_TABLE['Den_lower'] = Den_sd_lower
COUNT_TABLE['Cal_upper'] = Cal_sd_upper
COUNT_TABLE['Cal_lower'] = Cal_sd_lower

cutoffpoint=88

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.65)

twin1 = ax.twinx()
twin2 = ax.twinx()
twin3 = ax.twinx()

# Offset the right spine of twin2.  The ticks and label have already been
# placed on the right by twinx above.
twin1.spines.right.set_position(("axes", 1.00))
twin1.spines.right.set_color('C1')
twin2.spines.right.set_position(("axes", 1.20))
twin2.spines.right.set_color('C2')
twin3.spines.right.set_position(("axes", 1.45))
twin3.spines.right.set_color('C4')

p1, = ax.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Count'][cutoffpoint:].tolist(), color='black', label="observations")
p11, = ax.plot(COUNT_TABLE['Year'][:cutoffpoint].tolist(), COUNT_TABLE['Count'][:cutoffpoint].tolist(), label="Number of Observations", linestyle='dashed', color='black')

p2, = twin1.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Ext_avg'][cutoffpoint:].tolist(), "C1", label="Extension")
p21, = twin1.plot(COUNT_TABLE['Year'][:cutoffpoint].tolist(), COUNT_TABLE['Ext_avg'][:cutoffpoint].tolist(), "C1", label="Extension", linestyle='dashed')
p22, = twin1.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Ext_upper'][cutoffpoint:].tolist(), "C1" ,linestyle='dotted')
p23, = twin1.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Ext_lower'][cutoffpoint:].tolist(), "C1", linestyle='dotted')

p3, = twin2.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Den_avg'][cutoffpoint:].tolist(), "C2", label="Density")
p31, = twin2.plot(COUNT_TABLE['Year'][:cutoffpoint].tolist(), COUNT_TABLE['Den_avg'][:cutoffpoint].tolist(), "C2", label="Density", linestyle='dashed')
p32, = twin2.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Den_upper'][cutoffpoint:].tolist(), "C2", linestyle='dotted')
p33, = twin2.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Den_lower'][cutoffpoint:].tolist(), "C2", linestyle='dotted')

p4, = twin3.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Cal_avg'][cutoffpoint:].tolist(), "C4", label="Calcification")
p41, = twin3.plot(COUNT_TABLE['Year'][:cutoffpoint].tolist(), COUNT_TABLE['Cal_avg'][:cutoffpoint].tolist(), "C4", label="Calcification", linestyle='dashed')
p42, = twin3.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Cal_upper'][cutoffpoint:].tolist(), "C4", linestyle='dotted')
p43, = twin3.plot(COUNT_TABLE['Year'][cutoffpoint:].tolist(), COUNT_TABLE['Cal_lower'][cutoffpoint:].tolist(), "C4", linestyle='dotted')

ax.set(xlim=(1850, 2018), ylim=(0, 150), xlabel="Record Year", ylabel="Number of observations")
ax.yaxis.set_ticks(np.arange(0, 150, 15))
ax.grid(True, linestyle='--', color='grey', linewidth=0.5)

twin1.set(ylim=(-10, 50), ylabel="Extension [$\mathregular{mm.yr^{-1}}$]")
twin1.yaxis.set_ticks(np.arange(8, 22, 2))

twin2.set(ylim=(0.8, 5), ylabel="Density [$\mathregular{g.cm^{-3}}$]")
twin2.yaxis.set_ticks(np.arange(0.8, 1.8, 1/3.75))

twin3.set(ylim=(-1.4, 2.4), ylabel="Calcification [$\mathregular{g.cm^{-2}.yr^{-1}}$]")
twin3.yaxis.set_ticks(np.arange(1, 2.5, 1.5/7.5))

ax.yaxis.label.set_color(p1.get_color())
twin1.yaxis.label.set_color(p2.get_color())
twin2.yaxis.label.set_color(p3.get_color())
twin3.yaxis.label.set_color(p4.get_color())

ax.tick_params(axis='y', colors=p1.get_color())
twin1.tick_params(axis='y', colors=p2.get_color())
twin2.tick_params(axis='y', colors=p3.get_color())
twin3.tick_params(axis='y', colors=p4.get_color())

ax.legend(handles=[p1, p2, p3, p4])
plt.rcParams["font.family"] = "Arial"
#plt.show()
savepath='/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Figures/'
plt.savefig(savepath+'test_literature_time.svg')