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
Mydata_path="/Users/leonardobertini/Library/CloudStorage/Extracted_Results_SST_DSR_Kd490_MuseumSpecimens.xlsx"
Mydata = read_excel(Mydata_path, sheet = 'Grouped')
Mydata$Calcification_sd=as.numeric(Mydata$Calcification_sd)
Mydata$Extension_sd=as.numeric(Mydata$Extension_sd)
Mydata$Density_sd=as.numeric(Mydata$Density_sd)
Mydata['MidYear'] = Mydata$YearMin + (Mydata$YearMax-Mydata$YearMin)/2

#Literature data
Lit_path="/Users/leonardobertini/Library/CloudStorage/Extracted_Results_SST_DSR_Kd490_Literature.xlsx"
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