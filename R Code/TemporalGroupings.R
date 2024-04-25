library(readxl)
library(ggplot2)
library(lme4)
library(lmerTest)
library(reshape2)
library(ggrepel)
library(randomcoloR)
library(ggnewscale)
library(cowplot)
library(RColorBrewer)
library(svglite)
library(ggnewscale)
library(sjPlot)

# Getting the path of current R file.. this is where figures will be saved by default
setwd('/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/RScripts')

#My data
Mydata_path="/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Extracted_Results_MuseumSpecimens_Final.xlsx"
Mydata = read_excel(Mydata_path, sheet = 'Grouped')
Mydata$Calcification_sd=as.numeric(Mydata$Calcification_sd)
Mydata$Extension_sd=as.numeric(Mydata$Extension_sd)
Mydata$Density_sd=as.numeric(Mydata$Density_sd)
Mydata['MidYear'] = Mydata$YearMin + (Mydata$YearMax-Mydata$YearMin)/2

#Literature data
Lit_path="/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Extracted_Results_LitReview_Final.xlsx"
Litdata = read_excel(Lit_path, sheet = 'Grouped')
Litdata$YearMax=as.numeric(Litdata$YearMax)
Litdata$YearMin=as.numeric(Litdata$YearMin)
Litdata$Calcification_sd=as.numeric(Litdata$Calcification_sd)
Litdata$Extension_sd=as.numeric(Litdata$Extension_sd)
Litdata$Density_sd=as.numeric(Litdata$Density_sd)
Litdata['MidYear'] = Litdata$YearMin + (Litdata$YearMax-Litdata$YearMin)/2

#plotting 
color_scheme=c( '#f4cccc','#42d4f4', '#469990', 'red','blue', 'green', '#c294cf','#c21c45' ,
                '#f58231', '#911eb4','#A3AABE','#32CD32ff', '#c5e513ff', '#ff00ccfa', 
                '#5A5A5A', '#c0ccc0','#808000', '#c97f7f', '#11afccff', '#FFC000','#817567','#6c8dc5', '#FFD59A',
                '#a97947','#81b781','#ffe599', '#97ebdb' )

shape_codes = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)

stroke_colors=c( 'red','black', 'blue','#f58231','#f58231','black', 'black','black' ,
                        '#911eb4','red','red','red', 'black','black', 
                        'blue','blue','green','green', '#469990', '#469990', 'red','#a2d8f2','#AFCCFF')

fill_colors=c('red','black','blue','#f58231','#f58231','black', 'black','black' ,
                 '#911eb4','red','red','red', 'black','black', 
                 'blue','blue','black','green', 'green','#469990', '#469990', 'red','#a2d8f2','#AFCCFF')

# calcification rates
# only keep data from Indonesian studies
Litdata_filtered = subset(Litdata, 
                          Study == "Cahyarini 2008 - Jakarta" | 
                          Study == "Razak et al 2019 - Togean Is." |
                          Study ==  "Razak et al 2019 - West Papua"|
                          Study == "Edinger 2000 - Central Java")



#calcification
Fig1 = ggplot()+
  geom_pointrange(data = Mydata ,aes(x = MidYear,
                                     y = Calcification_mean,
                                     ymin=Calcification_mean-Calcification_sd, 
                                     ymax=Calcification_mean+Calcification_sd), 
                  alpha=0.3, size=0.01 ) +
  
  geom_linerange(data = Mydata,aes(x = MidYear, 
                                   y = Calcification_mean,
                                   xmin= YearMin,
                                   xmax = YearMax), alpha=0.3)+
  
  geom_point(data = Mydata ,aes(x = MidYear,
                                y = Calcification_mean, 
                                shape=Regional_Sector, 
                                fill=Regional_Sector, 
                                color=Regional_Sector),size=2)+
  
  scale_shape_manual(values = shape_codes)+
  scale_fill_manual(values = fill_colors)+
  scale_color_manual(values = stroke_colors)+

  new_scale_color() +
  ####literature
  geom_pointrange(data = Litdata_filtered,aes(x = MidYear,
                                     y = Calcification_mean,
                                     ymin=Calcification_mean-Calcification_sd, 
                                     ymax=Calcification_mean+Calcification_sd),
                  alpha=.3, size=0.01)+

  geom_linerange(data = Litdata_filtered, aes(x = MidYear, 
                                  xmin= YearMin,
                                  xmax = YearMax,
                                  y = Calcification_mean), alpha=.3)+
  
  geom_point(data = Litdata_filtered ,aes(x = MidYear,
                                 y = Calcification_mean, 
                                 color=Study), size=1.5)+

  scale_color_manual(values = c('red','black','#469990','#911eb4'))+
  
  scale_x_continuous(limits=c(1800,2015), breaks=seq(1800,2015,20))+
  scale_y_continuous(limits=c(.5,3), breaks=seq(.5,3,.5))+
  ylab(bquote(atop('Calcification ± 1sd', '(g.'~cm^-2~yr^-1~')')))+
  xlab('Year range (min-max)')+
  theme_bw() + 
  theme(axis.text = element_text(size = 5.5, color = 'black'), 
        axis.title = element_text(size = 7),
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
        legend.position ='')


#density
Fig2 = ggplot()+
  geom_pointrange(data = Mydata ,aes(x = MidYear,
                                     y = Density_mean,
                                     ymin=Density_mean-Density_sd, 
                                     ymax=Density_mean+Density_sd), 
                  alpha=0.3, size=0.01 ) +
  
  geom_linerange(data = Mydata,aes(x = MidYear, 
                                   y = Density_mean,
                                   xmin= YearMin,
                                   xmax = YearMax), alpha=0.3)+
  
  geom_point(data = Mydata ,aes(x = MidYear,
                                y = Density_mean, 
                                shape=Regional_Sector, 
                                fill=Regional_Sector, 
                                color=Regional_Sector),size=2)+
  
  scale_shape_manual(values = shape_codes)+
  scale_fill_manual(values = fill_colors)+
  scale_color_manual(values = stroke_colors)+
  
  new_scale_color() +
  ####literature
  geom_pointrange(data = Litdata_filtered,aes(x = MidYear,
                                              y = Density_mean,
                                              ymin=Density_mean-Density_sd, 
                                              ymax=Density_mean+Density_sd),
                  alpha=.3, size=0.01)+
  
  geom_linerange(data = Litdata_filtered, aes(x = MidYear, 
                                              xmin= YearMin,
                                              xmax = YearMax,
                                              y = Density_mean), alpha=.3)+
  
  
  geom_point(data = Litdata_filtered ,aes(x = MidYear,
                                          y = Density_mean, 
                                          color=Study), size=1.5)+
  scale_color_manual(values = c('red','black','#469990','#911eb4'))+
  
  scale_x_continuous(limits=c(1800,2015), breaks=seq(1800,2015,20))+
  scale_y_continuous(limits=c(1 ,1.7), breaks=seq(1,1.7,.1))+
  ylab(bquote(atop('Density ± 1sd', '(g.'~cm^-3~')')))+
  xlab('Year range (min-max)')+
  theme_bw() + 
  theme(axis.text = element_text(size = 5.5, color = 'black'), 
        axis.title = element_text(size = 7),
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
        legend.position ='')


#extension
Fig3 = ggplot()+
  geom_pointrange(data = Mydata ,aes(x = MidYear,
                                     y = Extension_mean,
                                     ymin=Extension_mean-Extension_sd, 
                                     ymax=Extension_mean+Extension_sd), 
                  alpha=0.3, size=0.01 ) +
  
  geom_linerange(data = Mydata,aes(x = MidYear, 
                                   y = Extension_mean,
                                   xmin= YearMin,
                                   xmax = YearMax), alpha=0.3)+
  
  geom_point(data = Mydata ,aes(x = MidYear,
                                y = Extension_mean, 
                                shape=Regional_Sector, 
                                fill=Regional_Sector, 
                                color=Regional_Sector),size=2)+
  
  scale_shape_manual(values = shape_codes)+
  scale_fill_manual(values = fill_colors)+
  scale_color_manual(values = stroke_colors)+
  
  new_scale_color() +
  ####literature
  geom_pointrange(data = Litdata_filtered,aes(x = MidYear,
                                              y = Extension_mean,
                                              ymin=Extension_sd-Extension_sd, 
                                              ymax=Extension_sd+Extension_sd),
                  alpha=.3, size=0.01)+
  
  geom_linerange(data = Litdata_filtered, aes(x = MidYear, 
                                              xmin= YearMin,
                                              xmax = YearMax,
                                              y = Extension_mean), alpha=.3)+
  
  
  geom_point(data = Litdata_filtered ,aes(x = MidYear,
                                          y = Extension_mean, 
                                          color=Study), size=1.5)+
  scale_color_manual(values = c('red','black','#469990','#911eb4'))+
  
  scale_x_continuous(limits=c(1800,2015), breaks=seq(1800,2015,20))+
  scale_y_continuous(limits=c(4 ,20), breaks=seq(4,20,2))+
  ylab(bquote(atop('Extension ± 1sd', '(mm.'~yr^-1~')')))+
  xlab('Year range (min-max)')+
  theme_bw() + 
  theme(axis.text = element_text(size = 5.5, color = 'black'), 
        axis.title = element_text(size = 7),
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
        legend.position ='')


lgd =cowplot::get_legend(Fig1 + guides(fill = guide_legend(byrow = TRUE)) +
                                theme(legend.spacing.y = unit(0.0, "cm"),
                                legend.key.size = unit(0.4, 'cm'),
                                legend.position = "right",
                                legend.title = element_text(size = 8),
                                legend.text =element_text(size = 6)))
ggpubr::as_ggplot(lgd)
ggsave("legend_test.png", width = 20, height = 20, units = "cm")


PLT = plot_grid(Fig1, Fig3, Fig2, nrow=3, labels=c('a)', 'b)', 'c)'), label_size = 10)
PLT2 = plot_grid(PLT, lgd, nrow=1, ncol=2)
PLT2
ggsave("test.png", width = 16, height = 16, units = "cm", dpi=300, bg = "white")
sjPlot::save_plot("test.svg", fig = PLT2, width = 16, height = 16)








############################ clean code after this
# importing datasets
datapath="/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Results/Coral_Growth_Data_Chapter4.xlsx" 
Lough_DF = read_excel(datapath, sheet = 'DataArrangedForPlot')
Lough_DF$Ext_cmyr = Lough_DF$Ext_mmyr*0.1

AMR_DF = read_excel(datapath, sheet = 'CalciRates')

Complete_DF = merge(Lough_DF,AMR_DF, by='CoralColony')
Complete_DF = within(Complete_DF, rm('ExtType_MinOrMax'))
Complete_DF = distinct(Complete_DF)
Complete_DF$TotalAgeMax = as.numeric(Complete_DF$TotalAgeMax)
Complete_DF$TotalAgeMax = as.numeric(Complete_DF$TotalAgeMin)
Complete_DF$AMR_Ext_min = Complete_DF$AMR_MeanDistance_cm*10 / Complete_DF$TotalAgeMax
Complete_DF$AMR_Ext_max = Complete_DF$AMR_MeanDistance_cm*10 / Complete_DF$TotalAgeMin
Complete_DF$AMR_Ext_avg = (Complete_DF$AMR_Ext_min + Complete_DF$AMR_Ext_max)/2
Complete_DF$Calci_AMR_avg = (Complete_DF$Calci_AMR_max + Complete_DF$Calci_AMR_min)/2


# Averaging Lough Track data for GIS and other Plots -------------------------------------------------------------------------
DF_GIS = Lough_DF %>%
  group_by(CoralColony) %>%
  summarize(MeanExt = mean(Ext_mmyr, na.rm=TRUE),
            MeanDensity = mean(TrackDensity, na.rm=TRUE),
            MeanCalcification = mean(TrackCalcification, na.rm=TRUE))


DF_GIS_merged = merge(DF_GIS, Lough_DF, by="CoralColony")
DF_GIS_merged = DF_GIS_merged[c("CoralColony",
                                "Location",
                                "MeanExt",
                                "MeanDensity",
                                "MeanCalcification",
                                "Collected_in")]

DF_GIS_merged = distinct(DF_GIS_merged)
write.csv(DF_GIS_merged, "Lough_averaged_values.csv")

#MERGING AMR dataset with LoughDataset ---------------------------------
Complete_DF = merge(Complete_DF,DF_GIS_merged, by='CoralColony')
Complete_DF = within(Complete_DF, rm("Ext_mmyr", "Ext_cmyr","TrackLength","TrackDensity","TrackCalcification", "Track_index",
                                     "TrackDurationMin","TrackDurationMax" ,"YearRangeMin","YearRangeMax", "Location.y", "Collected_in","Collected_in.y", "Location"))
Complete_DF = distinct(Complete_DF)

Complete_DF  = Complete_DF %>%
  rename(
    Collected_in = Collected_in.x,
    Location = Location.x,
    Lough_MeanExt = MeanExt,
    Lough_MeanDensity = MeanDensity,
    Lough_MeanCalcification = MeanCalcification
  )


# get colorpallete based on Locations ------------------------------------
colourcount =length(unique(Complete_DF$Location))
getPallete = colorRampPalette(brewer.pal(8,"Set1"))
color_scheme=c( '#e6194B', '#bd7dbd', '#f58231','#0000ffff',
                '#911eb4', '#42d4f4','#006400',
                '#32CD32ff', '#c5e513ff', '#469990','#ff00ccfa',
                '#000000', '#5A5A5A', '#C78752','#FFC000',
                '#800020', '#11afccff', '#808000','#4363d8')




# MY DATA  
Complete_DF_V = Complete_DF[Complete_DF$Slab_Orientation=='Vertical',]
Complete_DF_H = Complete_DF[Complete_DF$Slab_Orientation=='Horizontal',]

# Wrangling My Data -------------------------------------------------------

#Assigning Group categories on dataframe to calculate group averages over time
#Group 1 : 1800-1850 (Java - unknown sites)
#Group 2 : 1850-1900 ()
#Group 3 : 1910-1930 (Jakarta)
#Group 4 : 1910-1930 (Elsewhere)
#Group 5 : 1950-1980 
#Group 6 : Jakarta 1980-2000
#Group 7 : 2000's onwards


#VERTICAL
AMR_DF_Vertical=Complete_DF_V %>% separate_wider_delim(Collected_in, delim = "-", names = c("YearMin", "YearMax"), too_few = "align_start")
for(i in 1:nrow(AMR_DF_Vertical)) {
  
  if(is.na(AMR_DF_Vertical$YearMax[i])){
    AMR_DF_Vertical$YearMax[i] = AMR_DF_Vertical$YearMin[i]
  }
}


#Add group column for temporal analysis to my data
AMR_DF_Vertical$Temporal_Group = NaN
AMR_DF_Vertical$SampleType = 'This study'

for(i in 1:nrow(AMR_DF_Vertical)) {
  
  if(AMR_DF_Vertical$YearMax[i]>=1800 & AMR_DF_Vertical$YearMax[i]<=1825){
    AMR_DF_Vertical$Temporal_Group[i] = '1800-1825 (Java, metadata deficient)'
  }
  
  if(AMR_DF_Vertical$YearMax[i]>=1830 & AMR_DF_Vertical$YearMax[i]<=1860){
    AMR_DF_Vertical$Temporal_Group[i] = '1830-1860 (Java & Maluku, metadata deficient)'
  }
  
  if(AMR_DF_Vertical$YearMax[i]>=1880 & AMR_DF_Vertical$YearMax[i]<=1900){
    AMR_DF_Vertical$Temporal_Group[i] = '1880-1910 (Romang Is. & Selayar Is. & Binongko Is.)'
  }
  
  if(AMR_DF_Vertical$YearMax[i]>1910 & AMR_DF_Vertical$YearMax[i]<=1931 &
     grepl('Jakarta',AMR_DF_Vertical$Location[i]) ){
    AMR_DF_Vertical$Temporal_Group[i] = '1910-1930 (Jarkarta)'
  }
  
  if(AMR_DF_Vertical$YearMax[i]>1915 & AMR_DF_Vertical$YearMax[i]<=1935 &
     !grepl('Jakarta',AMR_DF_Vertical$Location[i]) ){
    AMR_DF_Vertical$Temporal_Group[i] = '1915-1935 (Leksoela Is. & SW Timor & Tanimbar Is. & Karakelong Is.)'
  }
  
  if(AMR_DF_Vertical$YearMax[i]>=1940 & AMR_DF_Vertical$YearMax[i]<=1955){
    AMR_DF_Vertical$Temporal_Group[i] = '1940-1960 (Biak, Papua)'
  }
  
#  if(AMR_DF_Vertical$YearMax[i]>=1975 & AMR_DF_Vertical$YearMax[i]<=1985){
#    AMR_DF_Vertical$Temporal_Group[i] = '1975-1985 (SW Sulawesi & Tukang Besi Is. & Taka Bonerate)'
#  }
  
  if(AMR_DF_Vertical$YearMax[i]>=1975 & AMR_DF_Vertical$YearMax[i]<=1985 &
     grepl('Sulawesi', AMR_DF_Vertical$Location[i])) {
    AMR_DF_Vertical$Temporal_Group[i] = '1976-1984 (SW Sulawesi)'
  }
  
  if(AMR_DF_Vertical$YearMax[i]>=1975 & AMR_DF_Vertical$YearMax[i]<=1985 &
     grepl(paste(c('Taka Bonerate','Tukang Besi'), collapse='|'), AMR_DF_Vertical$Location[i])) {
    AMR_DF_Vertical$Temporal_Group[i] = '1976-1984 (Taka Bonerate & Tukang Besi Is.)'
  }
  
  
}

#HORIZONTAL
AMR_DF_Horizontal=Complete_DF_H %>% separate_wider_delim(Collected_in, delim = "-", names = c("YearMin", "YearMax"), too_few = "align_start")
for(i in 1:nrow(AMR_DF_Horizontal)) {
  
  if(is.na(AMR_DF_Horizontal$YearMax[i])){
    AMR_DF_Horizontal$YearMax[i] = AMR_DF_Horizontal$YearMin[i]
  }
}

#Add group column for temporal analysis to my data
AMR_DF_Horizontal$Temporal_Group = NaN
AMR_DF_Horizontal$SampleType = 'This study'

for(i in 1:nrow(AMR_DF_Horizontal)) {
  
  if(AMR_DF_Horizontal$YearMax[i]>=1800 & AMR_DF_Horizontal$YearMax[i]<=1825){
    AMR_DF_Horizontal$Temporal_Group[i] = '1800-1825 (Java, metadata deficient)'
  }
  
  if(AMR_DF_Horizontal$YearMax[i]>=1830 & AMR_DF_Horizontal$YearMax[i]<=1860){
    AMR_DF_Horizontal$Temporal_Group[i] = '1830-1860 (Java & Maluku, metadata deficient)'
  }
  
  if(AMR_DF_Horizontal$YearMax[i]>=1880 & AMR_DF_Horizontal$YearMax[i]<=1900){
    AMR_DF_Horizontal$Temporal_Group[i] = '1880-1910 (Romang Is. & Selayar Is. & Binongko Is.)'
  }
  
  if(AMR_DF_Horizontal$YearMax[i]>1910 & AMR_DF_Horizontal$YearMax[i]<=1931 &
     grepl('Jakarta',AMR_DF_Horizontal$Location[i]) ){
    AMR_DF_Horizontal$Temporal_Group[i] = '1910-1930 (Jarkarta)'
  }
  
  if(AMR_DF_Horizontal$YearMax[i]>1915 & AMR_DF_Horizontal$YearMax[i]<=1935 &
     !grepl('Jakarta',AMR_DF_Horizontal$Location[i]) ){
    AMR_DF_Horizontal$Temporal_Group[i] = '1915-1935 (Leksoela Is. & SW Timor & Tanimbar Is. & Karakelong Is.)'
  }
  
  if(AMR_DF_Horizontal$YearMax[i]>=1940 & AMR_DF_Horizontal$YearMax[i]<=1955){
    AMR_DF_Horizontal$Temporal_Group[i] = '1940-1960 (Biak, Papua)'
  }
  
  if(AMR_DF_Horizontal$YearMax[i]>=1975 & AMR_DF_Horizontal$YearMax[i]<=1985 &
     grepl('Sulawesi', AMR_DF_Horizontal$Location[i])) {
    AMR_DF_Horizontal$Temporal_Group[i] = '1976-1984 (SW Sulawesi)'
  }
  
  if(AMR_DF_Horizontal$YearMax[i]>=1975 & AMR_DF_Horizontal$YearMax[i]<=1985 &
     grepl(paste(c('Taka Bonerate','Tukang Besi'), collapse='|'), AMR_DF_Horizontal$Location[i])) {
    AMR_DF_Horizontal$Temporal_Group[i] = '1976-1984 (Taka Bonerate & Tukang Besi Is.)'
  }
  

#if(AMR_DF_Horizontal$YearMax[i]>=1975 & AMR_DF_Horizontal$YearMax[i]<=1985){
#    AMR_DF_Horizontal$Temporal_Group[i] = '1975-1985 (SW Sulawesi & Tukang Besi Is. & Taka Bonerate)'
#  }
  
}

#### CREATING DF WITH GROUP AVERAGES #### #### 
#Removing unnecessary columns
AMR_DF_Vertical=within(AMR_DF_Vertical, rm('YearRangeMin.x',
                                           'YearRangeMax.x', 
                                           'YearRangeMin.y',
                                           'YearRangeMax.y',
                                           'PercentageChange_EXT',
                                           'PercentageChange_Density',
                                           'PercentageChange_Calci', 
                                           'SIO'))
#Removing duplicate rows
AMR_DF_Vertical=AMR_DF_Vertical[!duplicated(AMR_DF_Vertical), ]

#Removing unnecessary columns
AMR_DF_Horizontal=within(AMR_DF_Horizontal, rm('YearRangeMin.x',
                                           'YearRangeMax.x', 
                                           'YearRangeMin.y',
                                           'YearRangeMax.y',
                                           'PercentageChange_EXT',
                                           'PercentageChange_Density',
                                           'PercentageChange_Calci', 
                                           'SIO'))
#Removing duplicate rows
AMR_DF_Horizontal=AMR_DF_Horizontal[!duplicated(AMR_DF_Horizontal), ]



group_avg_df_AMR = function(Dataframe){
  
  Temporal_Averages_MyData = Dataframe %>%
    group_by(Temporal_Group) %>%
    summarise(Ext_rate_mean = mean(AMR_Ext_avg, na.rm=TRUE),
              Ext_rate_sd= sd(AMR_Ext_avg, na.rm=TRUE),
              Density_mean = mean(AMR_MeanDensity_gcm3, na.rm=TRUE),
              Density_sd = sd(AMR_MeanDensity_gcm3, na.rm=TRUE),
              Calci_mean = mean(Calci_AMR_avg, na.rm=TRUE),
              Calci_sd = sd(Calci_AMR_avg, na.rm=TRUE))
  
  Temporal_Averages_MyData=Temporal_Averages_MyData %>% separate_wider_delim(Temporal_Group, 
                                                                             delim = "-", 
                                                                             names = c("YearMin", "YearMaxx"), 
                                                                             too_few = "align_start",   
                                                                             cols_remove = FALSE)
  
  Temporal_Averages_MyData=Temporal_Averages_MyData %>% separate_wider_delim(YearMaxx, 
                                                                             delim = "(", 
                                                                             names = c("YearMax", "YearMaxxx"), 
                                                                             too_few = "align_start",   
                                                                             cols_remove = FALSE)
  
  Temporal_Averages_MyData = Temporal_Averages_MyData[ , !names(Temporal_Averages_MyData) %in% 
                                                         c("YearMaxx","YearMaxxx")]
  
  Temporal_Averages_MyData$YearMin = as.numeric(Temporal_Averages_MyData$YearMin)
  Temporal_Averages_MyData$YearMax = as.numeric(Temporal_Averages_MyData$YearMax)
  Temporal_Averages_MyData$MidYear = (Temporal_Averages_MyData$YearMin+ Temporal_Averages_MyData$YearMax)/2
  Temporal_Averages_MyData$SampleType='Coral head'
  Temporal_Averages_MyData$DataSource='This study'
  #TODO 
  # double check groupings to see if make sense and add number of colonies in each group
  Temporal_Averages_MyData$N_colonies_analysed = NaN
  #Group 1 : 1800-1850
  #Group 2 : 1850-1900
  #Group 3 : 1910-1930 (Jakarta)
  #Group 4 : 1910-1930 (Elsewhere)
  #Group 5 : 1950-1980 
  #Group 6 : Jakarta 1980-2000
  return(Temporal_Averages_MyData)
}

group_avg_df_Lough = function(Dataframe){
  
  Temporal_Averages_MyData = Dataframe %>%
    group_by(Temporal_Group) %>%
    summarise(Ext_rate_mean = mean(Lough_MeanExt, na.rm=TRUE),
              Ext_rate_sd= sd(Lough_MeanExt, na.rm=TRUE),
              Density_mean = mean(Lough_MeanDensity, na.rm=TRUE),
              Density_sd = sd(Lough_MeanDensity, na.rm=TRUE),
              Calci_mean = mean(Lough_MeanCalcification, na.rm=TRUE),
              Calci_sd = sd(Lough_MeanCalcification, na.rm=TRUE))
  
  Temporal_Averages_MyData=Temporal_Averages_MyData %>% separate_wider_delim(Temporal_Group, 
                                                                             delim = "-", 
                                                                             names = c("YearMin", "YearMaxx"), 
                                                                             too_few = "align_start",   
                                                                             cols_remove = FALSE)
  
  Temporal_Averages_MyData=Temporal_Averages_MyData %>% separate_wider_delim(YearMaxx, 
                                                                             delim = "(", 
                                                                             names = c("YearMax", "YearMaxxx"), 
                                                                             too_few = "align_start",   
                                                                             cols_remove = FALSE)
  
  Temporal_Averages_MyData = Temporal_Averages_MyData[ , !names(Temporal_Averages_MyData) %in% 
                                                         c("YearMaxx","YearMaxxx")]
  
  Temporal_Averages_MyData$YearMin = as.numeric(Temporal_Averages_MyData$YearMin)
  Temporal_Averages_MyData$YearMax = as.numeric(Temporal_Averages_MyData$YearMax)
  Temporal_Averages_MyData$MidYear = (Temporal_Averages_MyData$YearMin+ Temporal_Averages_MyData$YearMax)/2
  Temporal_Averages_MyData$SampleType='Coral head'
  Temporal_Averages_MyData$DataSource='This study'
  #TODO 
  # double check groupings to see if make sense and add number of colonies in each group
  Temporal_Averages_MyData$N_colonies_analysed = NaN
  #Group 1 : 1800-1850
  #Group 2 : 1850-1900
  #Group 3 : 1910-1930 (Jakarta)
  #Group 4 : 1910-1930 (Elsewhere)
  #Group 5 : 1950-1980 
  #Group 6 : Jakarta 1980-2000
  return(Temporal_Averages_MyData)
}


Temporal_Averages_MyData_V = group_avg_df_AMR(AMR_DF_Vertical)
Temporal_Averages_MyData_H = group_avg_df_AMR(AMR_DF_Horizontal)
Temporal_Averages_MyData_Lough = group_avg_df_Lough(AMR_DF_Vertical)



# Wrangling Literature Data -----------------------------------------------
lit_path = '/Users/leonardobertini/Library/CloudStorage/OneDrive-UniversityofBristol/Coral_QGIS_Spatial_Analyses/Revised_Literature_Points_Summarised.xlsx'
Literature_data = read_excel(lit_path, sheet = 'Revised_Literature_Points')

#FilteringData
#Only comparing museum coral data with measurements from other coral heads from the following regions that can be considered coral triangle
# Coral Triangle - Indonesia
# Coral Triangle - PNG
# Indonesia - Natuna Islands
# Indonesia - Jakarta

Lit_Filtered = Literature_data %>%
  filter(Reef_Province %in% c('Coral Triangle - Indonesia','Indonesia - Natuna Islands', 'Coral Triangle - PNG', 'Indonesia - Jakarta'),
         Core_Head ==  'Head')
Lit_Filtered$Ext_rate = as.numeric(Lit_Filtered$Ext_rate)
Lit_Filtered$Density = as.numeric(Lit_Filtered$Density)
Lit_Filtered$Calcification = Lit_Filtered$Ext_rate*Lit_Filtered$Density

#Adding min and max years to Literature dataset
Lit_Filtered=Lit_Filtered %>% separate_wider_delim(Year, delim = "-", names = c("YearMin", "YearMax"), too_few = "align_start")

for(i in 1:nrow(Lit_Filtered)) {
  
  if(is.na(Lit_Filtered$YearMax[i])){
    Lit_Filtered$YearMax[i] = Lit_Filtered$YearMin[i]
  }
  # do stuff with row
}

#Add group column for temporal analysis to lit review data
#Group 01 : 1970-1985 (South Sulawesi, Taka Bonerate)
#Group 02 : 1970-1985 (Jakarta)
#Group 03 : 1970-1990 (Central Java)
#Group 04 : 1990-1995 (Central Java) = Karimunjawa + Panjang Jepara + Bondo
#Group 05 : 1990-1995 (Sulawesi) = Edinger 2000
#Group 06 : 1990-1995 (Maluku) = Tanjung Setan + Ambon, Hila + Ambon, Wayame + Ambon, Wailiha
#Group 07 : 1996-2000 (PNG) 
#Group 08 : 1999-2005 (East Kalimantan) 
#Group 09 : 2005-2010 (Central Java)
#Group 10 :2001-2002 Wakatopbi (Kaledupa and Hoga sites)

Lit_Filtered$Temporal_Group = NaN
Lit_Filtered$SampleType = 'Literature data'

for(i in 1:nrow(Lit_Filtered)) {
  
  if(Lit_Filtered$YearMin[i]>1970 & Lit_Filtered$YearMax[i]<=1985 &
     grepl('Taka Bonerate', Lit_Filtered$Site[i] )){
    Lit_Filtered$Temporal_Group[i] = '1970-1985: Taka Bonerate South Sulawesi; Maier et al. (2004)'
  }
  
  if(Lit_Filtered$YearMin[i]>1970 & Lit_Filtered$YearMax[i]<=1985 &
     grepl('Thousand Islands', Lit_Filtered$Site[i] )){
    Lit_Filtered$Temporal_Group[i] = '1980-1985: Jakarta; Scoffin (1986)'
  }
  
  if(Lit_Filtered$YearMin[i]>1970 & Lit_Filtered$YearMax[i]<=1990 &
     grepl(c('Central'), Lit_Filtered$Site[i] )){
    Lit_Filtered$Temporal_Group[i] = '1970-1990: Central Java; Supriharyono (1998)'
  }
  
  if(Lit_Filtered$YearMin[i]>=1990 & Lit_Filtered$YearMax[i]<=1995 &
     grepl(paste(c('Karimunjawa','Panjang Jepara', 'Bondo'), collapse='|'), Lit_Filtered$Site[i])){
    Lit_Filtered$Temporal_Group[i] = '1990-1995: Central Java; Edinger et al. (2000)'
  }
  
  if(Lit_Filtered$YearMin[i]>=1990 & Lit_Filtered$YearMax[i]<=1995 &
     grepl(c('Sulawesi'), Lit_Filtered$Site[i] )){
    Lit_Filtered$Temporal_Group[i] = '1990-1995: SW Sulawesi; Edinger et al. (2000)'
  }
  
  if(Lit_Filtered$YearMin[i]>=1990 & Lit_Filtered$YearMax[i]<=1995 &
     grepl(paste(c('Ambon', 'Tanjung Setan'), collapse='|'), Lit_Filtered$Site[i])){
    Lit_Filtered$Temporal_Group[i] = '1990-1995: Ambon Maluku; Edinger et al. (2000)'
  }
  
  if(Lit_Filtered$YearMin[i]>=1999 & Lit_Filtered$YearMax[i]<=2005 &
     grepl('East Kalimantan', Lit_Filtered$Site[i] )){
    Lit_Filtered$Temporal_Group[i] = '1999-2005: East Kalimantan; Supriharyono (2004)'
  }
  
  if(Lit_Filtered$YearMin[i]>=2001 & Lit_Filtered$YearMax[i]<=2007 &
     grepl(paste(c('Karimunjawa', 'Bangkalan'), collapse='|'), Lit_Filtered$Site[i])){
    Lit_Filtered$Temporal_Group[i] = '2001-2007: Central and East Java; Nurgraha (2008)'
  }
  
  if(Lit_Filtered$YearMin[i]>=2001 & Lit_Filtered$YearMax[i]<=2002 &
     grepl('Wakatobi', Lit_Filtered$Site[i]) & 
     !grepl('Sampela', Lit_Filtered$Sample_ID[i]))  { #Sampela site is outlier, very low growth rate due to sedimentation impacts
    Lit_Filtered$Temporal_Group[i] = '2001-2002: Wakatobi SE Sulawesi; Crabbe et al. (2006)'
  }
  
}

Lit_Filtered = Lit_Filtered[Lit_Filtered$Temporal_Group != NaN,]

Lit_Filtered$Temporal_Group = as.factor(Lit_Filtered$Temporal_Group)
Temporal_Averages_Lit = Lit_Filtered %>%
  group_by(Temporal_Group) %>%
  summarise(Ext_rate_mean = mean(Ext_rate, na.rm=TRUE),
            Ext_rate_sd= sd(Ext_rate, na.rm=TRUE),
            Density_mean = mean(Density, na.rm=TRUE),
            Density_sd = sd(Density, na.rm=TRUE),
            Calci_mean = mean(Ext_rate*Density/10, na.rm=TRUE),
            Calci_sd = sd(Ext_rate*Density/10, na.rm=TRUE))

#adding number of samples per group n=
# 1 1970-1985 (Jakarta)             61  
# 2 1970-1985 (Taka Bonerate)      12
# 3 1970-1990 (Central Java)        15 
# 4 1990-1995 (Central Java)        26
# 5 1990-1995 (Maluku)              21
# 6 1990-1995 (Sulawesi)            27
# 8 1999-2005 (East Kalimantan)     8
# 9 2001-2002 (Wakatobi)            6        
# 10 2001-2007 (Central Java)       20

#Number_of_colonies_analysed
Temporal_Averages_Lit$N_colonies_analysed = c(61, 12, 15, 26, 21, 27, 8, 6, 20)
Temporal_Averages_Lit$SampleType='Coral head'
Temporal_Averages_Lit$DataSource='Literature data'
Temporal_Averages_Lit=Temporal_Averages_Lit %>% separate_wider_delim(Temporal_Group, 
                                                                     delim = "-", 
                                                                     names = c("YearMin", "YearMaxx"), 
                                                                     too_few = "align_start",   
                                                                     cols_remove = FALSE)

Temporal_Averages_Lit=Temporal_Averages_Lit %>% separate_wider_delim(YearMaxx, 
                                                                     delim = ":", 
                                                                     names = c("YearMax", "YearMaxxx"), 
                                                                     too_few = "align_start",   
                                                                     cols_remove = FALSE)


Temporal_Averages_Lit = Temporal_Averages_Lit[ , !names(Temporal_Averages_Lit) %in% 
                                                 c("YearMaxx","YearMaxxx")]


Temporal_Averages_Lit$YearMin = as.numeric(Temporal_Averages_Lit$YearMin)
Temporal_Averages_Lit$YearMax = as.numeric(Temporal_Averages_Lit$YearMax)
Temporal_Averages_Lit$MidYear = (Temporal_Averages_Lit$YearMin+ Temporal_Averages_Lit$YearMax)/2



# Now wrangling Lit data for Coral CORES and merge to figure --------------------------------------------------------------------
Lit_Filtered_CORES = Literature_data %>%
  filter(Reef_Province %in% c('Coral Triangle - Indonesia',
                              'Indonesia - Natuna Islands', 
                              'Coral Triangle - PNG', 
                              'Indonesia - Jakarta'), 
         Core_Head ==  'Core')
Lit_Filtered_CORES = Lit_Filtered_CORES[Lit_Filtered_CORES$Ext_rate!='NULL',]

#Groupings for cores
# 1970-2012 Indonesia - Nusa Penida 
# 1940-1981 Philippines
# 1994-1996 PNG (Hydrotermal)
# 1993-2011 PNG
# 1981-2012 Jakarta
# 1980-2010 Malaysia
# 1980-2010 Singapore
# 1992-2011 Indonesia - Nusa Tengara
# 1957-2011 Wakatobi

Lit_Filtered_CORES=Lit_Filtered_CORES %>% separate_wider_delim(Year, delim = "-", names = c("YearMin", "YearMax"), too_few = "align_start")
Lit_Filtered_CORES$Temporal_Group=NaN

for(i in 1:nrow(Lit_Filtered_CORES)) {
  
  if(Lit_Filtered_CORES$YearMin[i]>=1949 & Lit_Filtered_CORES$YearMax[i]<=2012 &
     grepl('Nusa Penida', Lit_Filtered_CORES$Site[i] )){
    Lit_Filtered_CORES$Temporal_Group[i] = '1949-2012: Indonesia Nusa Penida; Tito et al. (2016)'
  }
  
  if(Lit_Filtered_CORES$YearMin[i]>=1993 & Lit_Filtered_CORES$YearMax[i]<=2011 &
     grepl('Biak', Lit_Filtered_CORES$Site[i])){
    Lit_Filtered_CORES$Temporal_Group[i] = '1993-2011: Biak Papua Indonesia; Suharsono & Cahyarini (2012)'
  }
  
  if(Lit_Filtered_CORES$YearMin[i]>=1992 & Lit_Filtered_CORES$YearMax[i]<=2011 &
     grepl('Maumere', Lit_Filtered_CORES$Site[i])){
    Lit_Filtered_CORES$Temporal_Group[i] = '1992-2011: Maumere Nusa Tengara Indonesia; Suharsono & Cahyarini (2012)'
  }
  
  if(Lit_Filtered_CORES$YearMin[i]>=1957 & Lit_Filtered_CORES$YearMax[i]<=2011 &
     grepl('Wakatobi', Lit_Filtered_CORES$Site[i] )){
    Lit_Filtered_CORES$Temporal_Group[i] = '1957-2011: SE Sulawesi Wakatobi; Suharsono & Cahyarini (2012)'
  }
  
  if(Lit_Filtered_CORES$YearMin[i]>=1981 & Lit_Filtered_CORES$YearMax[i]<=2012 &
     grepl('Tunda', Lit_Filtered_CORES$Site[i])){
    Lit_Filtered_CORES$Temporal_Group[i] ='1981-2012: Banten Province Jakarta; Zamani & Arman (2016)'
  }
  
  if(Lit_Filtered_CORES$YearMin[i]>=1992 & Lit_Filtered_CORES$YearMax[i]<=2011 & Lit_Filtered_CORES$YearMin[i]!=1990 & 
     grepl('Thousand', Lit_Filtered_CORES$Site[i])){
    Lit_Filtered_CORES$Temporal_Group[i] ='1992-2011: Jakarta Thousand Islands; Arman et al. (2013)'
  }
  
  
  if(Lit_Filtered_CORES$YearMin[i]>=1992 & Lit_Filtered_CORES$YearMax[i]<=2005 & 
     grepl(paste(c('Jukung','Air', 'Bidadari'), collapse='|'), Lit_Filtered_CORES$Site[i])){
    Lit_Filtered_CORES$Temporal_Group[i] ='1992-2005: Jakarta Thousand Islands; Cahyarini (2008)'

  }
  
  }
  
Lit_Filtered_CORES = Lit_Filtered_CORES[Lit_Filtered_CORES$Temporal_Group != 'NaN',]
Lit_Filtered_CORES$Ext_rate = as.numeric(Lit_Filtered_CORES$Ext_rate)
Lit_Filtered_CORES$Density = as.numeric(Lit_Filtered_CORES$Density)

Temporal_Averages_Lit_CORES = Lit_Filtered_CORES %>%
  group_by(Temporal_Group) %>%
  summarise(Ext_rate_mean = mean(Ext_rate, na.rm=TRUE),
            Ext_rate_sd= sd(Ext_rate, na.rm=TRUE),
            Density_mean = mean(Density, na.rm=TRUE),
            Density_sd = sd(Density, na.rm=TRUE),
            Calci_mean = mean(Ext_rate*Density/10, na.rm=TRUE),
            Calci_sd = sd(Ext_rate*Density/10, na.rm=TRUE))

Temporal_Averages_Lit_CORES=Temporal_Averages_Lit_CORES %>% separate_wider_delim(Temporal_Group, 
                                                                     delim = "-", 
                                                                     names = c("YearMin", "YearMaxx"), 
                                                                     too_few = "align_start", 
                                                                     too_many = "merge",
                                                                     cols_remove = FALSE)

Temporal_Averages_Lit_CORES=Temporal_Averages_Lit_CORES %>% separate_wider_delim(YearMaxx, 
                                                                                 delim = ":", 
                                                                                 names = c("YearMax", "YearMaxxxx"), 
                                                                                 too_few = "align_start", 
                                                                                 too_many = "merge",
                                                                                 cols_remove = FALSE)
Temporal_Averages_Lit_CORES = Temporal_Averages_Lit_CORES[ , !names(Temporal_Averages_Lit_CORES) %in% 
                                                 c("YearMaxx","YearMaxxxx")]

Temporal_Averages_Lit_CORES$YearMin = as.numeric(Temporal_Averages_Lit_CORES$YearMin)
Temporal_Averages_Lit_CORES$YearMax = as.numeric(Temporal_Averages_Lit_CORES$YearMax)
Temporal_Averages_Lit_CORES$MidYear=(Temporal_Averages_Lit_CORES$YearMin + Temporal_Averages_Lit_CORES$YearMax)/2
Temporal_Averages_Lit_CORES$N_colonies_analysed= NaN
Temporal_Averages_Lit_CORES$SampleType='Coral core'
Temporal_Averages_Lit_CORES$DataSource='Literature data'

  
  
#Spatio/temporal plots for groupings  -----------------------------------

#filtering core and literature data sets to include only relevant areas
unique(Temporal_Averages_Lit$Temporal_Group)
remove=c('1994-1996 (Hydrothermal _ Ambitle Is., PNG)')
Temporal_Averages_Lit = Temporal_Averages_Lit[!(Temporal_Averages_Lit$Temporal_Group %in% remove),]

unique(Temporal_Averages_Lit_CORES$Temporal_Group)
remove=c('1994-1996 PNG (Hydrotermal)', '1940-1981 Philippines')
Temporal_Averages_Lit_CORES=Temporal_Averages_Lit_CORES[!(Temporal_Averages_Lit_CORES$Temporal_Group %in% remove),]


#TODO Fix legend label names 'Spatio-temporal groupings' / 'Data source' :Literature, this study
#### MERGING CORAL HEADS AND CORES 
Temporal_Averages_MyData_V$Measurement_Style='AMR-V'
Temporal_Averages_MyData_Lough$Measurement_Style='MGA'
Temporal_Averages_Lit$Measurement_Style='Literature'
Temporal_GROUPED_final = rbind(Temporal_Averages_MyData_V, Temporal_Averages_MyData_Lough, Temporal_Averages_Lit)
Temporal_GROUPED_final$Measurement_Style=as.factor(Temporal_GROUPED_final$Measurement_Style)

Fig1 = ggplot()+
  geom_pointrange(data = Temporal_GROUPED_final,aes(x = MidYear,
                                                    y = Ext_rate_mean,
                                                    ymin=Ext_rate_mean-Ext_rate_sd, 
                                                    ymax=Ext_rate_mean+Ext_rate_sd,
                                                    color = Temporal_Group,
                                                    fill = Temporal_Group,
                                                    shape = DataSource,
                                                    alpha= Measurement_Style))+
  
  scale_color_manual(values = color_scheme)+
  scale_fill_manual(values = color_scheme)+
  scale_shape_manual(values = c(21,24))+
  scale_alpha_manual(values=c(.4,1,1))+
  geom_linerange(data = Temporal_GROUPED_final,aes(x = MidYear, 
                                                   xmin= YearMin,
                                                   xmax = YearMax,
                                                   y = Ext_rate_mean,
                                                   color=Temporal_Group,
                                                   alpha=Measurement_Style))+
  scale_alpha_manual(values=c(.4,1,1))+
  
  scale_x_continuous(limits=c(1800,2015), breaks=seq(1800,2015,10))+
  scale_y_continuous(limits=c(0,20), breaks=seq(0,20,2))+
  ylab(bquote(atop('Mean skeletal extension ± 1sd', '(mm'~y^-1~')')))+
  xlab('Year range (min-max)')+
  theme_bw() + 
  theme(axis.text = element_text(size = 9, color = 'black'), 
        axis.title = element_text(size = 9), 
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .1),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .1),
        legend.position ='none')+
  
  ####cores
  geom_pointrange(data = Temporal_Averages_Lit_CORES,aes(x = MidYear,
                                                         y = Ext_rate_mean,
                                                         ymin=Ext_rate_mean-Ext_rate_sd, 
                                                         ymax=Ext_rate_mean+Ext_rate_sd,
                                                         shape=Temporal_Group), alpha=.5)+
  scale_shape_manual(values = c(0, 3, 4, 5, 6, 13, 14, 21, 24))+
  
  geom_linerange(data = Temporal_Averages_Lit_CORES,aes(x = MidYear, 
                                                        xmin= YearMin,
                                                        xmax = YearMax,
                                                        y = Ext_rate_mean), alpha=.5)

  #Plotting density 
  Fig2 = ggplot()+
    
    geom_pointrange(data = Temporal_GROUPED_final,aes(x = MidYear,
                                                      y = Density_mean,
                                                      ymin=Density_mean-Density_sd, 
                                                      ymax=Density_mean+Density_sd,
                                                      color = Temporal_Group,
                                                      fill = Temporal_Group,
                                                      shape = DataSource, 
                                                      alpha=Measurement_Style))+
    scale_alpha_manual(values=c(.4,1,1))+
    
    scale_color_manual(values = color_scheme)+
    scale_fill_manual(values = color_scheme)+
    scale_shape_manual(values = c(21,24))+
    geom_linerange(data = Temporal_GROUPED_final,aes(x = MidYear, 
                                                     xmin= YearMin,
                                                     xmax = YearMax,
                                                     y = Density_mean,
                                                     color=Temporal_Group,
                                                     alpha=Measurement_Style))+
    scale_alpha_manual(values=c(.4,1,1))+
    
    scale_x_continuous(limits=c(1800,2015), breaks=seq(1800,2015,10))+
    scale_y_continuous(limits=c(.8,1.8), breaks=seq(.8,1.8,.20))+
    ylab(bquote(atop('Mean density ± 1sd', '(g.'~cm^-3~')')))+
    xlab('Year range (min-max)')+
    theme_bw() + 
    theme(axis.text = element_text(size = 9, color = 'black'), 
          axis.title = element_text(size = 9), 
          panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .1),
          panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .1),
          legend.position ='none')+
    
    
    ####cores
    geom_pointrange(data = Temporal_Averages_Lit_CORES,aes(x = MidYear,
                                                           y = Density_mean,
                                                           ymin=Density_mean-Density_sd, 
                                                           ymax=Density_mean+Density_sd,
                                                           shape=Temporal_Group), alpha=.5)+
    scale_shape_manual(values = c(0, 3, 4, 5, 6, 13, 14, 21, 24))+
    
    
    geom_linerange(data = Temporal_Averages_Lit_CORES,aes(x = MidYear, 
                                                          xmin= YearMin,
                                                          xmax = YearMax,
                                                          y = Density_mean), alpha=.5)
  
  # calcification rates
  Fig3 = ggplot()+
    
    geom_pointrange(data = Temporal_GROUPED_final,aes(x = MidYear,
                                                      y = Calci_mean,
                                                      ymin=Calci_mean-Calci_sd, 
                                                      ymax=Calci_mean+Calci_sd,
                                                      color = Temporal_Group,
                                                      fill = Temporal_Group,
                                                      shape = DataSource,
                                                      alpha=Measurement_Style))+
    scale_alpha_manual(values=c(.4,1,1))+
    
    scale_color_manual(values = color_scheme)+
    scale_fill_manual(values = color_scheme)+
    scale_shape_manual(values = c(21,24))+
    geom_linerange(data = Temporal_GROUPED_final,aes(x = MidYear, 
                                                     xmin= YearMin,
                                                     xmax = YearMax,
                                                     y = Calci_mean,
                                                     color=Temporal_Group,
                                                     alpha=Measurement_Style))+
    scale_alpha_manual(values=c(.4,1,1))+
    
    scale_x_continuous(limits=c(1800,2015), breaks=seq(1800,2015,10))+
    scale_y_continuous(limits=c(.5,3), breaks=seq(.5,3,.5))+
    ylab(bquote(atop('Mean Calcification ± 1sd', '(g.'~cm^-2~yr^-1~')')))+
    xlab('Year range (min-max)')+
    theme_bw() + 
    theme(axis.text = element_text(size = 9, color = 'black'), 
          axis.title = element_text(size = 9), 
          panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .1),
          panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .1),
          legend.position ='none')+
    

    ####cores
    geom_pointrange(data = Temporal_Averages_Lit_CORES,aes(x = MidYear,
                                                           y = Calci_mean,
                                                           ymin=Calci_mean-Calci_sd, 
                                                           ymax=Calci_mean+Calci_sd,
                                                           shape=Temporal_Group), alpha=.5)+
    scale_shape_manual(values = c(0, 3, 4, 5, 6, 13, 14, 21, 24))+
    
    geom_linerange(data = Temporal_Averages_Lit_CORES,aes(x = MidYear, 
                                                          xmin= YearMin,
                                                          xmax = YearMax,
                                                          y = Calci_mean), alpha=.5)
  
  #Single Figure
  PLT=plot_grid(Fig3,Fig1, Fig2, nrow=3, labels=c('a)', 'b)', 'c)'))

  
  

# 
# #save fig
# ggsave(filename='/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Figures/Temporal_Groupings_draft.svg',
#        plot = MERGED_Fig,
#        device = svglite, 
#        width = 13 ,
#        height = 13 ,
#        units = "cm")
# 
# #save lgd
# ggsave(filename='/Users/leonardobertini/Library/CloudStorage/OneDrive-SharedLibraries-UniversityofBristol/grp-Chapter 4 - Leo - General/Figures/Temporal_Groupings_legend.svg',
#        plot = legend_b,
#        device = svglite, 
#        width = 30 ,
#        height = 30 ,
#        units = "cm")

