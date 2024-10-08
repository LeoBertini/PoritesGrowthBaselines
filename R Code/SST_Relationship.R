library(readxl)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(reshape2)
library(ggrepel)
library(randomcoloR)
library(ggnewscale)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(svglite)
library(ggpubr)
library(ggpmisc)
library(scales)
library(sjPlot)


# SST_relationship --------------------------------------------------------

# importing datasets
Mydata_path="/Users/leonardobertini/Library/CloudStorage/Extracted_Results_SST_DSR_Kd490_MuseumSpecimens_Final.xlsx"
Mydata = read_excel(Mydata_path, sheet = 'Sheet1')
Mydata$Regional_Sector=as.factor(Mydata$Regional_Sector)
Mydata$Endobionts=as.factor(Mydata$Endobionts)
Mydata$Corallith=as.factor(Mydata$Corallith)

Mydata['Source'] ='This Study'

#basic lm to check relationships
summary(lm(Mydata$MGA_Calcification~Mydata$MGA_Extension))
plot(Mydata$MGA_Calcification,Mydata$MGA_Extension)

#basic stats of min, max, average 
tapply(Mydata$MGA_Calcification, Mydata$Source, summary) 
tapply(Mydata$MGA_Extension, Mydata$Source, summary) 
tapply(Mydata$MGA_Density, Mydata$Source, summary) 
tapply(Mydata$Lat, Mydata$Source, summary) 
tapply(Mydata$SST_HadlSST_ann, Mydata$Source, summary)

Mydata['Calcification_Predicted'] = Mydata$SST_HadlSST_ann*0.33 -6.98

sd((Mydata$MGA_Calcification- Mydata$Calcification_Predicted)/Mydata$Calcification_Predicted)

#basic
summary(lm(Mydata$MGA_Calcification~Mydata$SST_HadlSST_ann))
plot(Mydata$SST_HadlSST_ann,Mydata$MGA_Calcification)

#basic2
summary(lm(Mydata$MGA_Calcification~Mydata$MGA_Extension))

#basic3
fig0 = ggplot()+
        geom_smooth(inherit.aes = FALSE, data=Mydata, method='lm', 
                    aes(x=SST_HadlSST_ann, 
                    y=MGA_Calcification), 
                    color='black', linetype='dashed', alpha=0.2, se=FALSE)+
        geom_point(data=Mydata, aes(x=SST_HadlSST_ann,
                              y=MGA_Calcification,
                              shape=Regional_Sector), 
                   size=3)+
       
  annotate(geom="text", x = 28, y = 2.5, label = 'Cal = 0.70*SST -18.35 | R2=0.16; p < 0.002')+
  scale_shape_manual(values= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))+
  xlab(bquote(atop('HadlSSTv1.1 ('~degree*C~')')))+
  ylab(bquote(atop('Calcification [g'~cm^-2~~y^-1~']')))+
  theme_bw()+
  theme(legend.position = 'right')
  

#basic averages of growth variables
mean(Mydata$MGA_Extension)
sd(Mydata$MGA_Extension)

mean(Mydata$MGA_Density)
sd(Mydata$MGA_Density)

mean(Mydata$MGA_Calcification)
sd(Mydata$MGA_Calcification)

mean(Mydata$SST_HadlSST_ann)
sd(Mydata$SST_HadlSST_ann)

mean(Mydata$DSR_ERB_mean)
sd(Mydata$DSR_ERB_mean)

mean(Mydata$Kd490_mean)
sd(Mydata$Kd490_mean)


#aov
res.aov = aov(MGA_Calcification ~ Corallith + Endobionts, data = Mydata)
summary(res.aov)

#basic stats for corallith groups
bxp1 <- ggplot(data=Mydata, aes(x = Corallith, y = MGA_Calcification))+geom_boxplot()+xlab('')
bxp2 <- ggplot(data=Mydata, aes(x = Corallith, y = MGA_Extension))+geom_boxplot()+xlab('')
bxp3 <- ggplot(data=Mydata, aes(x = Corallith, y = MGA_Density))+geom_boxplot()+xlab('Corallith? (Y/N)')

#basic stats for endobiont groups
bxp4 <- ggplot(data=Mydata, aes(x = factor(Endobionts, level=c('L','M','H')), y = MGA_Calcification))+geom_boxplot()+xlab('')
bxp5 <- ggplot(data=Mydata, aes(x = factor(Endobionts, level=c('L','M','H')), y = MGA_Extension))+geom_boxplot()+xlab('')
bxp6 <- ggplot(data=Mydata, aes(x = factor(Endobionts, level=c('L','M','H')), y = MGA_Density))+geom_boxplot()+xlab('Macro-endobiont abundance (Low-Medium-High)')

B1 = cowplot::plot_grid(bxp1, bxp2, bxp3, bxp4, bxp5, bxp6,  nrow=3, ncol=1,
                  labels = c('a)', 'b)', 'c)'))

B2 = cowplot::plot_grid(bxp4, bxp5, bxp6,  nrow=3, ncol=1,
               labels = c('d)', 'e)', 'f)'))
BOXES= cowplot::plot_grid(B1, B2, ncol=2)
BOXES


# filtering for unpaired t_test
LowEndo=Mydata[Mydata$Endobionts=='L',]
mean(LowEndo$MGA_Calcification)
sd(LowEndo$MGA_Calcification)

MediumEndo=Mydata[Mydata$Endobionts=='M',]
mean(MediumEndo$MGA_Calcification, na.rm=TRUE)
sd(MediumEndo$MGA_Calcification,na.rm=TRUE)

HighEndo=Mydata[Mydata$Endobionts=='H',]
mean(HighEndo$MGA_Calcification, na.rm=TRUE)
sd(HighEndo$MGA_Calcification,na.rm=TRUE)


# filtering for unpaired t_test
Coralliths=Mydata[Mydata$Corallith=='Y',]
mean(Coralliths$MGA_Calcification, na.rm=TRUE)
sd(Coralliths$MGA_Calcification, na.rm=TRUE)

NoCoralliths=Mydata[Mydata$Corallith!='Y',]
mean(NoCoralliths$MGA_Calcification, na.rm=TRUE)
sd(NoCoralliths$MGA_Calcification, na.rm=TRUE)

#checking for homoscedasticity and if p>0.05 in this test, then variances of both samples are homogenous
var.test(Coralliths$MGA_Calcification, NoCoralliths$MGA_Calcification)
t.test(Coralliths$MGA_Calcification,NoCoralliths$MGA_Calcification, var.equal = TRUE)

#checking for homoscedasticity and if p>0.05 in this test, then variances of both samples are homogenous
var.test(Coralliths$MGA_Extension, NoCoralliths$MGA_Extension)
t.test(Coralliths$MGA_Extension,NoCoralliths$MGA_Extension, var.equal = TRUE)

#checking for homoscedasticity and if p>0.05 in this test, then variances of both samples are homogenous
var.test(Coralliths$MGA_Density, NoCoralliths$MGA_Density)
t.test(Coralliths$MGA_Density,NoCoralliths$MGA_Density, var.equal = TRUE)
#conclusion: None of growth variables is significantly different between corallith and non-corallith groups


# Group by mean using dplyr
require("dyplr")     
MyData_grouped <- Mydata %>% group_by(Regional_Sector) %>% 
  summarise(Extension_mean=mean(MGA_Extension, na.rm=TRUE),
            Extension_sd=sd(MGA_Extension, na.rm=TRUE),
            Density_mean=mean(MGA_Density, na.rm=TRUE),
            Density_sd=sd(MGA_Density, na.rm=TRUE),
            Calcification_mean=mean(MGA_Calcification, na.rm=TRUE),
            Calcification_sd=sd(MGA_Calcification, na.rm=TRUE),
            Calcification_pred_mean=mean(Calcification_Predicted, na.rm=TRUE),
            Calcification_pred_sd=sd(Calcification_Predicted, na.rm=TRUE),
            SST_ann_mean=mean(SST_HadlSST_ann, na.rm=TRUE),
            SST_ann_sd=sd(SST_HadlSST_ann, na.rm=TRUE),
            DSR_ERB_mean_g=mean(DSR_ERB_mean, na.rm=TRUE),
            DSR_ERB_sd_g=sd(DSR_ERB_mean, na.rm=TRUE),
            Kd490_mean_g=mean(Kd490_mean, na.rm=TRUE),
            Kd490_sd_g=sd(Kd490_mean, na.rm=TRUE),
            Age_mean_mean =mean(Age_mean, na.rm=TRUE ))

#overall supression across all samples
mean(Mydata$MGA_Calcification/Mydata$Calcification_Predicted)

subset_in = subset(MyData_grouped, 
                Regional_Sector == 'Jakata Bay'|
                Regional_Sector == 'Jakata Bay - Kapal Is.'|
                Regional_Sector == 'Spermonde - SW Sulawesi early 20th century'|
                Regional_Sector == 'Spermonde - SW Sulawesi late 20th century'|
                Regional_Sector == 'Togian Is. - Gulf of Tomini, North Sulawesi'
                  )

subset_out = subset(MyData_grouped, 
                   Regional_Sector != 'Jakata Bay'|
                     Regional_Sector != 'Jakata Bay - Kapal Is.'|
                     Regional_Sector != 'Spermonde - SW Sulawesi early 20th century'|
                     Regional_Sector != 'Spermonde - SW Sulawesi late 20th century'|
                     Regional_Sector != 'Togian Is. - Gulf of Tomini, North Sulawesi'
)



#getting year min and year max
Regional_Sectors=unique(Mydata$Regional_Sector)
YearMin=vector()
YearMax=vector()
library(dplyr)
library(tidyr)

  for (j in sort(Regional_Sectors)) {
    print(j)
    sub= Mydata[Mydata$Regional_Sector == j,]
    sub2= sub %>% separate(`MGA year range`, c('YearMin', 'YearMax'))
    YearMin[j]= min(sub2$YearMin)
    YearMax[j]= max(sub2$YearMax)
}
MyData_grouped['Source'] = 'This Study'
MyData_grouped['YearMin'] = YearMin
MyData_grouped['YearMax'] = YearMax

#some basic comparisons
aa = MyData_grouped[MyData_grouped$Regional_Sector=='Spermonde - SW Sulawesi early 20th century',]
bb = MyData_grouped[MyData_grouped$Regional_Sector=='Spermonde - SW Sulawesi late 20th century',]

aa$Extension_mean/bb$Extension_mean
aa$Density_mean/bb$Density_mean


#Literature data
Lit_path="/Users/leonardobertini/Library/CloudStorage/Extracted_Results_SST_DSR_Kd490_Literature.xlsx"
Litdata = read_excel(Lit_path, sheet = 'Sheet1')
Litdata$Reef_Province=as.factor(Litdata$Reef_Province)
Litdata$Extension=as.numeric(Litdata$Extension)
Litdata$Density=as.numeric(Litdata$Density)
Litdata$Calcification=as.numeric(Litdata$Calcification)
Litdata$Study=factor(Litdata$Study)
Litdata = Litdata[!is.na(Litdata$Calcification),] #removing points from which there's no calcification
Litdata['Reference'] = 'Literature'
Litdata['Calcification_Predicted'] = Litdata$SST_HadlSST_ann*0.33 -6.98


#Grouped Lit by Study and then get year ranges
LitGrouped =  Litdata %>% group_by(Study) %>% 
  summarise(Extension_mean=mean(Extension, na.rm=TRUE),
            Extension_sd=sd(Extension, na.rm=TRUE),
            Density_mean=mean(Density, na.rm=TRUE),
            Density_sd=sd(Density, na.rm=TRUE),
            Calcification_mean=mean(Calcification, na.rm=TRUE),
            Calcification_sd=sd(Calcification, na.rm=TRUE),
            SST_ann_mean=mean(SST_HadlSST_ann, na.rm=TRUE),
            SST_ann_sd=sd(SST_HadlSST_ann, na.rm=TRUE),
            DSR_ERB_mean_g=mean(DSR_ERB_mean, na.rm=TRUE),
            DSR_ERB_sd_g=sd(DSR_ERB_mean, na.rm=TRUE),
            Kd490_mean_g=mean(Kd490_mean, na.rm=TRUE),
            Kd490_sd_g=sd(Kd490_mean, na.rm=TRUE),
            )

#getting year min and year max
Study=unique(Litdata$Study)
YearMin=vector()
YearMax=vector()
library(dplyr)
library(tidyr)

for (j in sort(Study)) {
  print(j)
  sub= Litdata[Litdata$Study == j,]
  sub2= sub %>% separate(Year_range, c('YearMin', 'YearMax'))
  YearMin[j]= min(sub2$YearMin)
  YearMax[j]= max(sub2$YearMax)
}

LitGrouped['YearMin'] = YearMin
LitGrouped['YearMax'] = YearMax
LitGrouped['Reference'] = 'Literature'



color_scheme=c( '#f4cccc','#42d4f4', '#469990', 'red','blue', 'green', '#c294cf','#c21c45' ,
                '#f58231', '#911eb4','#A3AABE','#32CD32ff', '#c5e513ff', '#ff00ccfa', 
                '#5A5A5A', '#c0ccc0','#808000', '#c97f7f', '#11afccff','#179FE0','#a2d8f2','#AFCCFF',
                '#FFC000','#817567','#6c8dc5', '#f785c9',
                '#a97947','#81b781','#ffe599', '#97ebdb' )


####figures   
  fig1 = ggplot()+
    geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                   y = Calcification,
                                   color = Study),
               size=1)+
    
    geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                  y = Calcification_mean,
                                  shape = Regional_Sector),
               size=1.5, color = 'black')+
    geom_abline(slope = 0.33, intercept = -6.98, linetype='solid', size=0.25, color='black')+ #coefficients from Lough(2008)

    scale_shape_manual(values= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))+
    scale_color_manual(values=color_scheme)+
    
    xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Calcification [g'~cm^-2~~y^-1~']')))+
    theme_bw()+
    theme_bw()+
    theme(axis.text = element_text(size = 5.5, color = 'black'), 
          axis.title = element_text(size = 7),
          panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
          panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
          legend.position ='')
  

  
#sjPlot::save_plot("Fig4aLarge.svg", fig = fig1, width = 18, height = 10)

  fig2 = ggplot()+
    geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                  y = Extension,
                                  color = Study),
               size=1)+
    
    geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                         y = Extension_mean,
                                         shape = Regional_Sector),
               size=1.5, color = 'black')+
    geom_abline(slope = 2.975, intercept = -65.459, linetype='solid', size=0.25, color='black')+
    scale_shape_manual(values= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))+
    scale_color_manual(values=color_scheme)+
    
    xlab(bquote(atop('SST ('~degree*C~')')))+
    ylab(bquote(atop('Extension','(mm'~y^-1~')')))+
    theme_bw()+
    theme(axis.text = element_text(size = 5.5, color = 'black'), 
          axis.title = element_text(size = 7),
          panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
          panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
          legend.position ='')
  
  
  fig3 = ggplot()+
    geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                  y = Density,
                                  color = Study),
               size=1)+
    
    geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                         y = Density_mean,
                                         shape = Regional_Sector),
               size=1.5, color = 'black')+
    geom_abline(slope = -0.089, intercept = 3.653, linetype='solid', size=0.25, color='black')+
    scale_shape_manual(values= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))+
    scale_color_manual(values=color_scheme)+
    
    xlab(bquote(atop('SST ('~degree*C~')')))+
    ylab(bquote(atop('Density','(g'~cm^-3~')')))+
    theme_bw()+
    theme(axis.text = element_text(size = 5.5, color = 'black'), 
          axis.title = element_text(size = 7),
          panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
          panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .025),
          legend.position ='')

lgd =cowplot::get_legend(fig1 + guides(fill = guide_legend(byrow = TRUE)) +
                           theme(legend.spacing.y = unit(0.1, "cm"),
                                 legend.key.size = unit(10, 'pt'),
                                 legend.position = "right",
                                 legend.title = element_text(size = 8),
                                 legend.text =element_text(size = 6)))
ggpubr::as_ggplot(lgd)
ggsave("legend_test.png", width = 20, height = 20, units = "cm")

FIG = cowplot::plot_grid(fig1, fig2, fig3, nrow=3, ncol=1, labels = c("a)", "b)", "c)"), label_size = 10)
ggsave("test.png", width = 10, height = 16, units = "cm", dpi=300, bg = "white")

FIGG = cowplot::plot_grid(FIG, lgd, nrow=1, ncol=2)
sjPlot::save_plot("test.svg", fig = FIGG, width = 16, height = 16)


########################
#### now plotting calci:SST with point fill as colour gradient based on remove-sensing variable

#normalizing Kd490 --> as proxy of percent of DSR absorbed within the first 1m
Litdata$KdNorm=(1-(1/exp(Litdata$Kd490_mean))) #proportion of initial irradiance that is absorbed after penetrating 1 m. 
MyData_grouped$KdNorm=(1-(1/exp(MyData_grouped$Kd490_mean_g)))

fig4 = ggplot()+
  geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                y = Calcification,
                                fill = DSR_ERB_mean),
             size=1.5, shape=21, color='grey',stroke=0.1)+
  
  geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                       y = Calcification_mean,
                                       shape = Regional_Sector,
                                       fill = DSR_ERB_mean_g),
             size=1.5, shape=22 ,stroke=0.2)+  
  geom_abline(slope = 0.33, intercept = -6.98, linetype='solid', size=0.5, color='black')+
  #scale_fill_gradient2(low = "blue", high = "red", midpoint=midval, mid='white')+
  scale_fill_distiller(palette = "Spectral", limits = c(140,240), oob=squish)+
  xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Calcification [g'~cm^-2~~y^-1~']')))+
  theme_bw()+
  theme(legend.position = c(0.15, 0.75), 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8),
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))+
        labs(fill=bquote('DSR [W'~m^-2~']'))

fig5 = ggplot()+
  geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                y = Calcification,
                                fill = Kd490_mean),
             size=1.5, shape=21, color='grey',stroke=0.1)+
  
  geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                       y = Calcification_mean,
                                       shape = Regional_Sector,
                                       fill = Kd490_mean_g),
             size=1.5, shape=22 ,stroke=0.2)+
  geom_abline(slope = 0.33, intercept = -6.98, linetype='solid', size=0.5, color='black')+
  new_scale_color() +
  scale_fill_distiller(palette = "RdBu", limits = c(0,0.2), oob=squish)+
  xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Calcification [g'~cm^-2~~y^-1~']')))+
  theme_bw()+
  theme(legend.position = c(0.15, 0.75), 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8),
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))+
  labs(fill=bquote(K[d]~'490 ['~m^-1~']'))


#low turbidity and high radiation --> DSR/Kd490
fig6 = ggplot()+
  geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                y = Calcification,
                                fill = DSR_ERB_mean*(1-KdNorm)),
             size=1.5, shape=21, color='grey',stroke=0.1)+
  
  geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                       y = Calcification_mean,
                                       shape = Regional_Sector,
                                       fill = MyData_grouped$DSR_ERB_mean_g*(1-MyData_grouped$KdNorm)),
             size=1.5, shape=22, stroke=0.2)+
  geom_abline(slope = 0.33, intercept = -6.98, linetype='solid', size=0.5, color='black')+
  scale_fill_distiller(palette = "Spectral", limits=c(125, 225), oob=squish)+
  xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Calcification [g'~cm^-2~~y^-1~']')))+
  theme_bw()+
  theme(legend.position = c(0.20, 0.75), 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8), 
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))+
  labs(fill=bquote('DSR/e'^(K[d]~'490')))


FIG_GRADIENTS1 = cowplot::plot_grid(fig4, fig5, fig6, nrow=3, ncol=1, labels = c("a)", "b)", "c)"), label_size=10)
ggsave(filename = '/Users/leonardobertini/Desktop/Light_Turbidity_Calci.png', plot=FIG_GRADIENTS1, dpi=300, height=18, width=9, units='cm')


#############################

fig7 = ggplot()+
  geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                y = Extension,
                                fill = DSR_ERB_mean),
             size=1.5, shape=21, color='grey',stroke=0.1)+
  
  geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                       y = Extension_mean,
                                       shape = Regional_Sector,
                                       fill = DSR_ERB_mean_g),
             size=1.5, shape=22 ,stroke=0.2)+  
  geom_abline(slope = 2.975, intercept = -65.459, linetype='solid', size=0.5, color='black')+
  scale_fill_distiller(palette = "Spectral", limits = c(140,240), oob=squish)+
  xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Extension', '[mm.'~yr^-1~']')))+
  theme_bw()+
  theme(legend.position = c(0.15, 0.75), 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8),
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))+
  labs(fill=bquote('DSR [W'~m^-2~']'))



fig8 = ggplot()+
  geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                y = Extension,
                                fill = Kd490_mean),
             size=1.5, shape=21, color='grey',stroke=0.1)+
  
  geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                       y = Extension_mean,
                                       shape = Regional_Sector,
                                       fill = Kd490_mean_g),
             size=1.5, shape=22 ,stroke=0.2)+
  geom_abline(slope = 2.975, intercept = -65.459, linetype='solid', size=0.5, color='black')+
  new_scale_color() +
  scale_fill_distiller(palette = "RdBu", limits = c(0,0.2), oob=squish)+
  xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Extension', '[mm.'~yr^-1~']')))+
  theme_bw()+
  theme(legend.position = c(0.15, 0.75), 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8),
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))+
  labs(fill=bquote(K[d]~'490 ['~m^-1~']'))

fig9= ggplot()+
    geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                  y = Extension,
                                  fill = DSR_ERB_mean/exp(Kd490_mean)),
               size=1.5, shape=21, color='grey',stroke=0.1)+
    
    geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                         y = Extension_mean,
                                         shape = Regional_Sector,
                                         fill = DSR_ERB_mean_g/exp(Kd490_mean_g)),
               size=1.5, shape=22, stroke=0.2)+
    geom_abline(slope = 2.975, intercept = -65.459, linetype='solid', size=0.5, color='black')+
    scale_fill_distiller(palette = "Spectral", limits=c(130, 230), oob=squish)+
    xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Extension', '[mm.'~yr^-1~']')))+
  theme_bw()+
    theme(legend.position = c(0.20, 0.75), 
          legend.text=element_text(size=7), 
          legend.title = element_text(size=8),
          legend.key.height= unit(0.22, 'cm'), 
          legend.key.width = unit(0.4, 'cm'), 
          text = element_text(size=8), 
          panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
          panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))+
  labs(fill=bquote('Ratio DSR/'~e^(K[d]~'490')))
  
FIG_GRADIENTS2 = cowplot::plot_grid(fig7, fig8, fig9, nrow=3, ncol=1, labels = c("d)", "e)", "f)"), label_size=10)
ggsave(filename = '/Users/leonardobertini/Desktop/Light_Turbidity_Extension.png', plot=FIG_GRADIENTS2, dpi=300, height=18, width=9, units='cm')


#############################

fig10 = ggplot()+
  geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                y = Density,
                                fill = DSR_ERB_mean),
             size=1.5, shape=21, color='grey',stroke=0.1)+
  
  geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                       y = Density_mean,
                                       shape = Regional_Sector,
                                       fill = DSR_ERB_mean_g),
             size=1.5, shape=22 ,stroke=0.2)+  
  geom_abline(slope = -0.089, intercept = 3.653, linetype='solid', size=0.5, color='black')+
  scale_fill_distiller(palette = "Spectral", limits = c(140,240), oob=squish)+
  xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Density','(g.'~cm^-3~')')))+
  theme_bw()+
  theme(legend.position = c(0.15, 0.25), 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8),
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))+
  labs(fill=bquote('DSR [W'~m^-2~']'))



fig11 = ggplot()+
  geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                y = Density,
                                fill = Kd490_mean),
             size=1.5, shape=21, color='grey',stroke=0.1)+
  
  geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                       y = Density_mean,
                                       shape = Regional_Sector,
                                       fill = Kd490_mean_g),
             size=1.5, shape=22 ,stroke=0.2)+
  geom_abline(slope = -0.089, intercept = 3.653, linetype='solid', size=0.5, color='black')+
  new_scale_color() +
  scale_fill_distiller(palette = "RdBu", limits = c(0,0.2), oob=squish)+
  xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Density','(g.'~cm^-3~')')))+
  theme_bw()+
  theme(legend.position = c(0.15, 0.25),
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8),
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))+
  labs(fill=bquote(K[d]~'490 ['~m^-1~']'))

fig12= ggplot()+
  geom_point(data = Litdata,aes(x = SST_HadlSST_ann,
                                y = Density,
                                fill = DSR_ERB_mean/exp(Kd490_mean)),
             size=1.5, shape=21, color='grey',stroke=0.1)+
  
  geom_point(data = MyData_grouped,aes(x = SST_ann_mean,
                                       y = Density_mean,
                                       shape = Regional_Sector,
                                       fill = DSR_ERB_mean_g/exp(Kd490_mean_g)),
             size=1.5, shape=22, stroke=0.2)+
  geom_abline(slope = -0.089, intercept = 3.653, linetype='solid', size=0.5, color='black')+
  scale_fill_distiller(palette = "Spectral", limits=c(130, 230), oob=squish)+
  xlab(bquote(atop('SST ('~degree*C~')')))+
  ylab(bquote(atop('Density','(g.'~cm^-3~')')))+
  theme_bw()+
  theme(legend.position = c(0.20, 0.25), 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8), 
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))+
  labs(fill=bquote('Ratio DSR/'~e^(K[d]~'490')))

FIG_GRADIENTS3 = cowplot::plot_grid(fig10, fig11, fig12, nrow=3, ncol=1, labels = c("g)", "h)", "i)"), label_size=10)
ggsave(filename = '/Users/leonardobertini/Desktop/Light_Turbidity_Density.png', plot=FIG_GRADIENTS3, dpi=300, height=18, width=9, units='cm')


FIG_GRADIENTS34=cowplot::plot_grid(fig4, fig5, fig7, fig8, fig10, fig11, nrow=3, ncol=2, labels=c("a)", "b)", "c)", "d)", "e)", "f)"), label_size=10)

########################

# correlations between density and Kd490 -----------------------------------------------------------------------
calcification_binded = c(Litdata$Calcification,MyData_grouped$Calcification_mean)
source_binded = c(Litdata$Reference,MyData_grouped$Source)
extension_binded = c(Litdata$Extension,MyData_grouped$Extension_mean)
density_binded = c(Litdata$Density,MyData_grouped$Density_mean)
sst_binded = c(Litdata$SST_HadlSST_ann,MyData_grouped$SST_ann_mean)
kd490_binded = c(Litdata$Kd490_mean,MyData_grouped$Kd490_mean_g)
dsr_binded = c(Litdata$DSR_ERB_mean,MyData_grouped$DSR_ERB_mean_g)
ratio_binded = c(dsr_binded/kd490_binded)
calc_offset_binded =c(abs((Litdata$Calcification/Litdata$Calcification_Predicted)-1), 
                     abs((MyData_grouped$Calcification_mean/MyData_grouped$Calcification_pred_mean)-1))

binded_dataset = data.frame(calcification_binded,extension_binded,density_binded,sst_binded,kd490_binded, dsr_binded, source_binded, ratio_binded, calc_offset_binded)

#some outliers?
Q = quantile(kd490_binded, probs=c(.25, .75), na.rm = FALSE)
iqr = IQR(kd490_binded)
up =  Q[2]+1.5*iqr # Upper Range  
low = Q[1]-1.5*iqr # Lower Range

eliminated = subset(binded_dataset, binded_dataset$kd490_binded > low &  binded_dataset$kd490_binded < up)

#linear models
m1=lm(binded_dataset$density_binded~binded_dataset$kd490_binded)
summary(m1)

m2=lm(binded_dataset$density_binded~binded_dataset$dsr_binded)
summary(m2)

m3=lm(eliminated$density_binded~eliminated$ratio_binded)
summary(m3)

m4=lm(binded_dataset$calc_offset_binded~binded_dataset$dsr_binded)
summary(m4)
plot(binded_dataset$dsr_binded,binded_dataset$calc_offset_binded, ylim=c(0,1))

#figures
figx1= ggplot(data=binded_dataset, aes(x=kd490_binded, 
                                   y=density_binded))+
  geom_point(aes(shape=source_binded,
                 fill=source_binded))+
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=c('white','black'))+
  geom_smooth(method='lm', color='black',linetype='dashed', fullrange=TRUE )+
  xlab(bquote(~K[d]~'490'~'['~m^-1~']'))+
  ylab(bquote(atop('Density [g.'~cm^-3~']')))+
  theme_bw()+
  theme(legend.position = '', 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8), 
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))
  

figx2= ggplot(data=binded_dataset, aes(x=dsr_binded, 
                                   y=density_binded))+
  geom_point(aes(shape=source_binded,
                 fill=source_binded))+
  scale_fill_manual(values=c('white','black'))+
  
  scale_shape_manual(values=c(21,22))+
  geom_smooth(method='lm', color='black',linetype='dashed', fullrange=TRUE )+
  xlab(bquote('DSR [W'~m^-2~']'))+
  ylab(bquote(atop('Density [g.'~cm^-3~']')))+
  theme_bw()+
  theme(legend.position = '', 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8), 
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))

figx3= ggplot(data=binded_dataset, aes(x=ratio_binded, 
                                       y=density_binded))+
  geom_point(aes(shape=source_binded, 
                 fill=source_binded))+
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=c('white','black'))+
  geom_smooth(method='lm', color='black',linetype='dashed', fullrange=TRUE )+
  xlab(bquote('Ratio DSR/'~K[d]~'490'))+
  ylab(bquote(atop('Density [g.'~cm^-3~']')))+
  theme_bw()+
  theme(legend.position = '', 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8), 
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))


PLTx=cowplot::plot_grid(figx2, figx1, figx3, nrow=1, ncol=3, labels=c('a)', 'b)', 'c)'), label_size = 10)

figx4= ggplot(data=binded_dataset, aes(x=dsr_binded, 
                                       y=calc_offset_binded*100))+
  geom_point(aes(shape=source_binded, 
                 fill=source_binded))+
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=c('white','black'))+
  geom_smooth(method='lm', color='black',linetype='dashed', fullrange=TRUE )+
  ylim(0,80)+
  xlab(bquote('DSR [W'~m^-2~']'))+
  ylab(bquote(atop('Absolute Calcification Offset [%]','w.r.t  Lough(2008)'~SST[ann]~':Calc predictions')))+
  theme_bw()+
  theme(legend.position = '', 
        legend.text=element_text(size=7), 
        legend.title = element_text(size=8),
        legend.key.height= unit(0.22, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        text = element_text(size=8), 
        panel.grid.major = element_line(linetype = 'dotted', colour = "black", linewidth = .05),
        panel.grid.minor = element_line(linetype = 'dotted', colour = "black", linewidth = .05))

PLTx2=cowplot::plot_grid(fig4, figx4, fig5,figx1, nrow=2, ncol=2, labels=c('a)', 'b)', 'c)', 'd)'), label_size = 10)
PLTx2
ggsave(filename = '/Users/leonardobertini/Desktop/FigExtra.png', plot=PLTx2, dpi=300, height=10, width=15, units='cm')
sjPlot::save_plot("FigExtra.svg", fig = PLTx2, width = 18, height = 10)

