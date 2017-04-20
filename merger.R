# Kota Miura April, 2017

library("ggplot2")
library("dplyr")

gppath <- '/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/'
folders = c(
  "161122 ctrl croped and 16 frames",
  "161206 ctrl croped and 16 frames",  
  "161129 ROCKinh croped and 16 frames",
  "161205 10umblebb croped and 16 frames",
  "161222 claD0.06 croped and 16 frames",
  "170130 CK666 croped and 16 frames",
  "170206 KOCLCAB croped and 50 frames",
  "170215 SMIFH2 croped and 16 frames",
  "170217 Jasplakinolide croped and 16 frames",
  "170221 Jasplakinolide croped and 16 frames",
  "170222 Genistein croped and 16 frames",
  "170302 bcyclo croped and 16 frames",
  "170307 a5b1peptidomimetic croped and 16 frames",
  "170328 beta1ABp5d2 croped and 14 frames",
  "170405 Hela croped and 16 frames"
)

datalist = list(
  "ctrl1" = seq(5),
  "ctrl2"= seq(8),
  "ROCKinh" = c(0, 1, 2, 5, 6) + 1,
  "Blebb10um" = c(0, 2, 3, 6, 7) + 1,
  "ClaD06" = c(0, 1, 2, 3) + 1,
  "CK666"= c(0, 1, 2, 3, 7, 8) + 1, 
  "KOCLCAB" = c(1, 2, 4, 5, 6, 7, 8, 9) + 1,
  "SMIFH2" = c(1, 2, 3, 4, 5, 6, 7) + 1,
  "Jaspla" = c(4) + 1,
  "Jaspla2" = c(1, 2, 3, 4) + 1,
  "Genistein" = seq(11),
  "bcyclo" = seq(9),
  "a5b1" = c(0, 1, 2, 4, 5, 6, 7, 8, 9, 10) + 1,
  "beta1ABp5d2" = c(1, 2, 3, 4, 5, 6, 9) + 1,
  "Hela" = c(0, 1, 2, 3, 4, 5, 6, 8, 9) + 1
)

datanames = names(datalist)

## box plot order
xordering = c("ctrl1", "ctrl2", "ROCKinh", "Blebb10um", "ClaD06", "CK666", "KOCLCAB", "SMIFH2", "Jaspla", "Jaspla2", "Genistein", "bcyclo", "a5b1", "beta1ABp5d2", "Hela")
  
#   c(
#   "ctrl1",
#   "ctrl2",
#   "ROCKinh",
#   "Blebb10um",
#   "ClaD06",
#   "CK666", 
#   "KOCLCAB",
#   "SMIFH2",
#   "Jaspla",
#   "Genistein",
#   "bcyclo",
#   "a5b1",
#   "beta1ABp5d2",
#   "Hela"
# )

# prepare an empty list, which will hold all imported data as elements.
# key = experiment type
# value = dataframe, imported from ImageJ results
dflist <- vector("list", length(datanames))
# fill names first
names(dflist) <- datanames
# adding values (data) to each keys
for (ind in 1:length(datanames)){
  dflist[[datanames[ind]]] <- read.csv2( file.path(gppath, folders[[ind]], "results.csv"), sep=",", dec=".")
  dflist[[datanames[ind]]]$type = rep(datanames[[ind]], nrow(dflist[[datanames[ind]]]))
}

# merge all dataframes in the above dictionary in a single dataframe
for (i in 1:length(datanames)){
  if (i == 1){
    alldata = merge(dflist[[1]], dflist[[2]], all = TRUE)
  } else if (i > 2){
    alldata = merge(alldata, dflist[[i]], all = TRUE)
  }
}

# specifically create a data frame for all controls. 
controlAll = merge(dflist[[1]], dflist[[2]], all = TRUE)

# compute pre / post drug treatment RipOff density. 
alldata$RipOff_Density <- alldata$RipOff.counts5_10 / alldata$Area.um2. * 100
alldata$RipOffPreDensity <- alldata$RipOff.counts1_4 / alldata$Area.um2. * 100

#assessing changes of individuals (ratio of post/pre on per frame bases)
#alldata$RipOffChange <- alldata$RipOff_Density/alldata$RipOffPreDensity/6*4
alldata$RipOffChange <- alldata$RipOff_Density/alldata$RipOffPreDensity

# Correction of Rip-Off density, as post-drug rip-offs are based on less virus. 
alldata$RipOffChangeCorrected <- alldata$RipOffChange * ( alldata$Dots.Total / (alldata$Dots.Total - alldata$RipOff.counts1_4))

### time series

timecourseAll = data.frame()
for (i in 1:length(folders)){
  curfolder = file.path(gppath, folders[i])
  exptype = names(datalist)[i]
  items = length(datalist[[exptype]])
  pathvector = character(length(items))
  for (j in 1:items){
    expnum =  datalist[[exptype]][j]
    datapath = file.path(curfolder, paste('cell', expnum, '_virus_median.tif_counts.csv', sep=''))
    curdata = read.csv(datapath)
    curdata$Type = rep(names(datalist)[i], nrow(curdata))
    curdata$CellID = rep(as.character(j), nrow(curdata))
    cellarea = alldata$Area.um2.[alldata$CellID==expnum & alldata$type==exptype]
    
    curdata$RipOffDensity = curdata$Counts / cellarea * 100
    curdata$VirusDensity = curdata$TotalCounts / cellarea * 100
    timecourseAll = rbind(timecourseAll, curdata)
    #pathvector[j] = datapath
  }
  #datafiles = lapply(pathvector, read.csv)
}

##### PLOTTINGS #####

# checking the dependency of rip-off density to overall virus density, only control
ggplot(subset(alldata, type %in% c("ctrl1", "ctrl2")), aes(x=Density, y = RipOff_Density)) + geom_point()
ggsave(paste('Control_RipOffDensityVSdensity.png'), width=25, height=20, units="cm")

# checking the dependency of rip-off density to overall virus density, all data
ggplot(alldata, aes(x=Density, y = RipOff_Density)) + geom_point()
ggsave(paste('All_RipOffDensityVSdensity.png'), width=25, height=20, units="cm")

# checking the dependency of rip-off (pre-drug) density to overall virus density, all data
ggplot(alldata, aes(x=Density, y = RipOffPreDensity)) + geom_point()
ggsave(paste('All_RipOffPreDensityVSdensity.png'), width=25, height=20, units="cm")


ylabeltext = "Rip-Off-Density [counts/100um^2]"
ggplot(alldata, aes(x=type, y=RipOff_Density)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext)
ggsave(paste('RipOff_Density_allExperiments.png'), width=25, height=15, units="cm")

ylabeltext2 = "Rip-Off-Density (Before Adding Drug)\n[counts/100um^2]"
ggplot(alldata, aes(x=type, y=RipOffPreDensity)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext2)
ggsave(paste('RipOff_Density_allExperiments_BeforeDrug.png'), width=25, height=15, units="cm")

ylabeltext3 = "RipOff Change \npost / pre treatment per frame"
ggplot(alldata, aes(x=type, y=RipOffChange)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext3)
ggsave(paste('RipOffDensityChanges.png'), width=25, height=20, units="cm")

# only controls
ggplot(subset(alldata, type %in% c("ctrl1", "ctrl2")), aes(x=type, y=RipOffChange)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments") + scale_y_continuous(name=ylabeltext3)
ggsave(paste('RipOffDensityChanges_OnlyControls.png'), width=25, height=20, units="cm")

# checking the density for each type (VirusDensity_in_each_Experiment)
ylabeltext4 = "Virus Density in 2nd frame\n [counts/um^2]"
ggplot(alldata, aes(x=type, y = Density)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext4)
ggsave(paste('virusDensity_allExperiments.png'), width=25, height=15, units="cm")

ylabeltext5 = "RipOff Change (Density Corrected) \npost / pre treatment per frame"
ggplot(alldata, aes(x=type, y=RipOffChangeCorrected)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext5)
ggsave(paste('RipOffDensityChanges_DensityCorrected.png'), width=25, height=20, units="cm")

### plotting time series

ggplot(data = filter(timecourseAll, Frame > 2))  + 
  geom_path(aes(x=Frame, y=Counts, color=CellID)) + 
  scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + 
  scale_y_continuous(name='Rip-Off [Counts]') + 
  facet_wrap(~ Type, nrow = 4)
ggsave('TimeCourse_Counts.png', width=25, height=20, units="cm")

for (thename in names(datalist)){
  typemax = max(timecourseAll$Counts[timecourseAll$Type==thename])
  annoY = typemax/2
  ggplot(data=filter(timecourseAll, Type==thename, Frame > 2)) + 
    geom_path(aes(x=Frame, y=Counts, color=CellID)) + 
    scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + 
    scale_y_continuous(name='Rip-Off [Counts]') +
    annotate("text", x=12, y=annoY, label= thename, size=30, color='gray')
  ggsave(paste("TimeCourse_Counts_",thename,".png", sep=''), width=25, height=20, units="cm")
}

ggplot(data=filter(timecourseAll, Frame > 2)) + geom_path(aes(x=Frame, y=RipOffDensity, color=CellID)) + scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + scale_y_continuous(name='Rip-Off Density \n[counts/100um^2]') + facet_wrap(~ Type, nrow = 4)
ggsave(paste('TimeCourse_RipOffDensity.png'), width=25, height=20, units="cm")

for (thename in names(datalist)){
  typemax = max(timecourseAll$RipOffDensity[timecourseAll$Type==thename])
  annoY = typemax/2
  ggplot(data=filter(timecourseAll, Type==thename, Frame > 2)) + 
    geom_path(aes(x=Frame, y=RipOffDensity, color=CellID)) + 
    scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + 
    scale_y_continuous(name='Rip-Off Density \n[counts/100um^2]') +
    annotate("text", x=12, y=annoY, label= thename, size=30, color='gray')
  ggsave(paste("TimeCourse_RipOffDensity_",thename,".png", sep=''), width=25, height=20, units="cm")
}

ggplot(data=filter(timecourseAll, Frame > 2)) + geom_path(aes(x=Frame, y=VirusDensity, color=CellID)) + scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + scale_y_continuous(name='Virus Density \n[counts/100um^2]') + facet_wrap(~ Type, nrow = 4)
ggsave(paste('TimeCourse_VirusDensity.png'), width=25, height=20, units="cm")

for (thename in names(datalist)){
  typemax = max(timecourseAll$VirusDensity[timecourseAll$Type==thename])
  annoY = typemax/2
  ggplot(data=filter(timecourseAll, Type==thename, Frame > 2)) + 
    geom_path(aes(x=Frame, y=VirusDensity, color=CellID)) + 
    scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + 
    scale_y_continuous(name='Virus Density \n[counts/100um^2]', limits=c(0, NA)) +
    annotate("text", x=12, y=annoY, label= thename, size=30, color='gray')
  ggsave(paste("TimeCourse_VirusDensity_",thename,".png", sep=''), width=25, height=20, units="cm")
}

                                                          
