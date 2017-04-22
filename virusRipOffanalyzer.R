# Kota Miura April, 2017

library("ggplot2")
library("dplyr")

# inputcsv = read.csv('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/data_lists/test.csv')
inputcsvpath = file.choose()
inputcsv = read.csv( inputcsvpath )

### load "results.csv" 
gppath = as.character(inputcsv$root[1])
folders = as.character(inputcsv$datafolder)
datanames = as.character(inputcsv$type)

# prepare an empty list, which will hold all imported data as elements.
# key = experiment type
# value = dataframe, imported from ImageJ results
dflist <- vector("list", length(datanames))
# fill names first
names(dflist) <- datanames
datalist <- vector("list", length(datanames))
names(datalist) <- datanames
# adding values (data) to each keys
for (ind in 1:length(datanames)){
  print(datanames[ind])
  print(file.path(gppath, folders[[ind]], "results.csv"))
  dflist[[datanames[ind]]] <- read.csv( file.path(gppath, folders[[ind]], "results.csv"))
  dflist[[datanames[ind]]]$type = rep(datanames[[ind]], nrow(dflist[[datanames[ind]]]))
  ids = inputcsv[ind, 4:ncol(inputcsv)]
  ids = ids[!is.na(ids)]
  datalist[[datanames[ind]]] <- ids
}

#datanames = names(datalist)

## box plot order
xordering = c("ctrl1", "ctrl2", "ROCKinh", "Blebb10um", "ClaD06", "CK666", "KOCLCAB", "SMIFH2", "Jaspla", "Jaspla2", "Genistein", "bcyclo", "a5b1", "beta1ABp5d2", "Hela")
  
# dflist <- vector("list", length(datanames))
# names(dflist) <- datanames
# 
# for (ind in 1:length(datanames)){
#   dflist[[datanames[ind]]] <- read.csv2( file.path(gppath, folders[[ind]], "results.csv"), sep=",", dec=".")
#   dflist[[datanames[ind]]]$type = rep(datanames[[ind]], nrow(dflist[[datanames[ind]]]))
# }

# merge all dataframes in the above dictionary in a single dataframe
if ( length( datanames) > 1){
  for (i in 1:length(datanames)){
    if (i == 2){
      alldata = merge(dflist[[1]], dflist[[2]], all = TRUE)
    } else if (i > 2){
      alldata = merge(alldata, dflist[[i]], all = TRUE)
    }
  }
} else {
  alldata = dflist[[1]]
}

# specifically create a data frame for all controls. 
# controlAll = merge(dflist[[1]], dflist[[2]], all = TRUE)

### time series

timecourseAll = data.frame()
for (i in 1:length(folders)){
  curfolder = file.path(gppath, folders[i])
  exptype = names(datalist)[i]
  items = length(datalist[[exptype]])
  pathvector = character(length(items))
  print(exptype)
  for (j in 1:items){
    expnum =  datalist[[exptype]][j]
    print(paste('Cell', expnum))
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

### compute paramters
# compute pre / post drug treatment RipOff density. 
alldata$RipOffPreDensity <- alldata$RipOff.counts2_4 / alldata$Area.um2. * 100
alldata$RipOff_Density <- alldata$RipOff.counts5_7 / alldata$Area.um2. * 100
alldata$RipOff2_Density <- alldata$RipOff.counts8_10 / alldata$Area.um2. * 100

# alternative
# sum(timecourseAll$Counts[((timecourseAll$Frame >= 2) & (timecourseAll$Frame <= 4) & (timecourseAll$CellID==2))])

#assessing changes of individuals (ratio of post/pre on per frame bases)
#alldata$RipOffChange <- alldata$RipOff_Density/alldata$RipOffPreDensity/6*4
alldata$RipOffChange <- alldata$RipOff_Density/alldata$RipOffPreDensity
alldata$RipOff2Change <- alldata$RipOff2_Density/alldata$RipOffPreDensity

# Correction of Rip-Off density, as post-drug rip-offs are based on less virus. 
alldata$RipOffChangeCorrected <- alldata$RipOffChange * ( alldata$Dots.Total / (alldata$Dots.Total - alldata$RipOff.counts2_4))
alldata$RipOff2ChangeCorrected <- alldata$RipOffChange * ( alldata$Dots.Total / (alldata$Dots.Total - alldata$RipOff.counts2_4 - alldata$RipOff.counts5_7))

##### PLOTTINGS #####

if ( !dir.exists( file.path( gppath, 'plots') )) {
  dir.create( file.path( gppath, 'plots') )
}

outfoldername = strftime(Sys.time(), "%Y%m%d-%H%M")
outfolderpath = file.path(gppath, "plots", outfoldername)
if ( !dir.exists( outfolderpath )) {
  dir.create( outfolderpath )
}

outsubfoldername = file.path(outfolderpath, 'time_courses')
if ( !dir.exists(outsubfoldername)) {
  dir.create( outsubfoldername )
}

# checking the dependency of rip-off density to overall virus density, only control
if (('ctrl1' %in% datanames) & ('ctrl2' %in% datanames)) {
  ggplot(subset(alldata, type %in% c("ctrl1", "ctrl2")), aes(x=Density, y = RipOff_Density)) + geom_point()
  ggsave(file.path(outfolderpath , 'Control_RipOffDensityVSdensity.png'), width=25, height=20, units="cm")
}

# checking the dependency of rip-off density to overall virus density, all data
ggplot(alldata, aes(x=Density, y = RipOff_Density)) + geom_point()
ggsave(file.path(outfolderpath , 'All_RipOffDensityVSdensity.png'), width=25, height=20, units="cm")

# checking the dependency of rip-off (pre-drug) density to overall virus density, all data
ggplot(alldata, aes(x=Density, y = RipOffPreDensity)) + geom_point()
ggsave(file.path(outfolderpath , 'All_RipOffPreDensityVSdensity.png'), width=25, height=20, units="cm")


ylabeltext = "Rip-Off-Density\nF5-F7 [counts/100um^2]"
ggplot(alldata, aes(x=type, y=RipOff_Density)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext)
ggsave(file.path(outfolderpath , 'RipOff_Density_allExperiments.png'), width=25, height=15, units="cm")

ylabeltext = "Rip-Off-Density\nF8-F10 [counts/100um^2]"
ggplot(alldata, aes(x=type, y=RipOff_Density)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext)
ggsave(file.path(outfolderpath , 'RipOff_Density2_allExperiments.png'), width=25, height=15, units="cm")

ylabeltext2 = "Rip-Off-Density (Before Adding Drug)\n[counts/100um^2]"
ggplot(alldata, aes(x=type, y=RipOffPreDensity)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext2)
ggsave(file.path(outfolderpath , 'RipOff_Density_allExperiments_BeforeDrug.png'), width=25, height=15, units="cm")

ylabeltext3 = "RipOff Change \npost[5-7] / pre treatment per frame"
ggplot(alldata, aes(x=type, y=RipOffChange)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext3)
ggsave(file.path(outfolderpath , 'RipOffDensityChanges.png'), width=25, height=20, units="cm")

ylabeltext3 = "RipOff Change \npost[8-10] / pre treatment per frame"
ggplot(alldata, aes(x=type, y=RipOff2Change)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext3)
ggsave(file.path(outfolderpath , 'RipOff2DensityChanges.png'), width=25, height=20, units="cm")

# only controls
ggplot(subset(alldata, type %in% c("ctrl1", "ctrl2")), aes(x=type, y=RipOffChange)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments") + scale_y_continuous(name=ylabeltext3)
ggsave(file.path(outfolderpath , 'RipOffDensityChanges_OnlyControls.png'), width=25, height=20, units="cm")

ggplot(subset(alldata, type %in% c("ctrl1", "ctrl2")), aes(x=type, y=RipOff2Change)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments") + scale_y_continuous(name=ylabeltext3)
ggsave(file.path(outfolderpath , 'RipOff2DensityChanges_OnlyControls.png'), width=25, height=20, units="cm")

# RipOff changes, with corrections
ylabeltext5 = "RipOff Change (Density Corrected) \npost[5-7] / pre treatment per frame"
ggplot(alldata, aes(x=type, y=RipOffChangeCorrected)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext5)
ggsave(file.path(outfolderpath , 'RipOffDensityChanges_DensityCorrected.png'), width=25, height=20, units="cm")

ylabeltext5 = "RipOff Change (Density Corrected) \npost[8-10] / pre treatment per frame"
ggplot(alldata, aes(x=type, y=RipOff2ChangeCorrected)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext5)
ggsave(file.path(outfolderpath , 'RipOff2DensityChanges_DensityCorrected.png'), width=25, height=20, units="cm")

# checking the density for each type (VirusDensity_in_each_Experiment)
ylabeltext4 = "Virus Density in 2nd frame\n [counts/um^2]"
ggplot(alldata, aes(x=type, y = Density)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext4)
ggsave(file.path(outfolderpath , 'virusDensity_allExperiments.png'), width=25, height=15, units="cm")


### plotting time series


timecourseAllF2 = filter(timecourseAll, Frame >= 2, Frame <= 16)

ggplot(data = timecourseAllF2)  + 
  geom_path(aes(x=Frame, y=Counts, color=CellID)) + 
  scale_x_continuous(name = 'Time [Frame]', limits=c(2, 16)) + 
  scale_y_continuous(name='Rip-Off [Counts]', limits = NULL) + 
  facet_wrap(~ Type, nrow = 4)
ggsave(file.path(outsubfoldername , 'TimeCourse_Counts.png'), width=25, height=20, units="cm")

for (thename in names(datalist)){
  typemax = max(timecourseAllF2$Counts[timecourseAllF2$Type==thename])
  annoY = typemax/2
  ggplot(data= filter(timecourseAllF2, Type==thename) ) + 
    geom_path(aes(x=Frame, y=Counts, color=CellID)) + 
    scale_x_continuous(name = 'Time [Frame]', limits=c(2, 16)) + 
    scale_y_continuous(name='Rip-Off [Counts]', limits = NULL) +
    annotate("text", x=12, y=annoY, label= thename, size=30, color='gray')
  ggsave(file.path(outsubfoldername , paste("TimeCourse_Counts_",thename,".png", sep='')), width=25, height=20, units="cm")
}

ggplot(data= timecourseAllF2) + 
  geom_path(aes(x=Frame, y=RipOffDensity, color=CellID)) + 
  scale_x_continuous(name = 'Time [Frame]', limits=c(2, 16)) + 
  scale_y_continuous(name='Rip-Off Density \n[counts/100um^2]', limits = NULL) + 
  facet_wrap(~ Type, nrow = 4)
ggsave(file.path(outsubfoldername , paste('TimeCourse_RipOffDensity.png')), width=25, height=20, units="cm")

for (thename in names(datalist)){
  typemax = max(timecourseAllF2$RipOffDensity[timecourseAllF2$Type==thename])
  annoY = typemax/2
  ggplot(data= filter(timecourseAllF2, Type==thename)) + 
    geom_path(aes(x=Frame, y=RipOffDensity, color=CellID)) + 
    scale_x_continuous(name = 'Time [Frame]', limits=c(2, 16)) + 
    scale_y_continuous(name='Rip-Off Density \n[counts/100um^2]', limits = NULL) +
    annotate("text", x=12, y=annoY, label= thename, size=30, color='gray')
  ggsave(file.path(outsubfoldername , paste("TimeCourse_RipOffDensity_",thename,".png", sep='')), width=25, height=20, units="cm")
}

ggplot(data= timecourseAllF2) + 
  geom_path(aes(x=Frame, y=VirusDensity, color=CellID)) + 
  scale_x_continuous(name = 'Time [Frame]', limits=c(2, 16)) + 
  scale_y_continuous(name='Virus Density \n[counts/100um^2]', limits = NULL) + 
  facet_wrap(~ Type, nrow = 4)
ggsave(file.path(outsubfoldername , paste('TimeCourse_VirusDensity.png')), width=25, height=20, units="cm")

for (thename in names(datalist)){
  typemax = max(timecourseAllF2$VirusDensity[timecourseAllF2$Type==thename])
  annoY = typemax/2
  print(thename)
  ggplot(data= filter(timecourseAllF2, Type==thename)) + 
    geom_path(aes(x=Frame, y=VirusDensity, color=CellID)) + 
    scale_x_continuous(name = 'Time [Frame]', limits=c(2, 16)) + 
    scale_y_continuous(name='Virus Density \n[counts/100um^2]', limits = NULL ) +
    annotate("text", x=12, y=annoY, label= thename, size=30, color='gray')
  ggsave(file.path(outsubfoldername , paste("TimeCourse_VirusDensity_",thename,".png", sep='')), width=25, height=20, units="cm")
}

                                                          
