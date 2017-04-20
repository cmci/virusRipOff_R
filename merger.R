# Kota Miura April, 2017

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
  "170222 Genistein croped and 16 frames",
  "170302 bcyclo croped and 16 frames",
  "170307 a5b1peptidomimetic croped and 16 frames",
  "170328 beta1ABp5d2 croped and 14 frames",
  "170405 Hela croped and 16 frames"
)

datanames = c(
  "ctrl1",
  "ctrl2",
  "ROCKinh",
  "Blebb10um",
  "ClaD06",
  "CK666", 
  "KOCLCAB",
  "SMIFH2",
  "Jaspla",
  "Genistein",
  "bcyclo",
  "a5b1",
  "beta1ABp5d2",
  "Hela"
)

# prepare an empty list, which will hold all imported data as elements. 
dflist <- vector("list", length(datanames))
# fill names first
names(dflist) <- datanames
# adding elements
for (ind in 1:length(datanames)){
  dflist[[datanames[ind]]] <- read.csv2( file.path(gppath, folders[[ind]], "results.csv"), sep=",", dec=".")
  dflist[[datanames[ind]]]$type = rep(datanames[[ind]], nrow(dflist[[datanames[ind]]]))
}

# merge all data in a single dataframe
for (i in 1:length(datanames)){
  if (i == 1){
    alldata = merge(dflist[[1]], dflist[[2]], all = TRUE)
  } else if (i > 2){
    alldata = merge(alldata, dflist[[i]], all = TRUE)
  }
}
controlAll = merge(dflist[[1]], dflist[[2]], all = TRUE)

alldata$RipOff_Density <- alldata$RipOff.counts5_10 / alldata$Area.um2. * 100
alldata$RipOffPreDensity <- alldata$RipOff.counts1_4 / alldata$Area.um2. * 100

        
library("ggplot2")

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
xordering = c("ctrl1", "ctrl2", "ROCKinh", "Blebb10um", "ClaD06", "CK666", "KOCLCAB", "SMIFH2", "Jaspla", "Genistein", "bcyclo", "a5b1", "beta1ABp5d2", "Hela")
ggplot(alldata, aes(x=type, y=RipOff_Density)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext)
ggsave(paste('RipOff_Density_allExperiments.png'), width=25, height=15, units="cm")

ylabeltext2 = "Rip-Off-Density (Before Adding Drug)\n[counts/100um^2]"
ggplot(alldata, aes(x=type, y=RipOffPreDensity)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext2)
ggsave(paste('RipOff_Density_allExperiments_BeforeDrug.png'), width=25, height=15, units="cm")

#assessing changes of individuals (ratio of post/pre on per frame bases)
#alldata$RipOffChange <- alldata$RipOff_Density/alldata$RipOffPreDensity/6*4
alldata$RipOffChange <- alldata$RipOff_Density/alldata$RipOffPreDensity

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

alldata$RipOffChangeCorrected <- alldata$RipOffChange * ( alldata$Dots.Total / (alldata$Dots.Total - alldata$RipOff.counts1_4))

ylabeltext5 = "RipOff Change (Density Corrected) \npost / pre treatment per frame"
ggplot(alldata, aes(x=type, y=RipOffChangeCorrected)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext5)
ggsave(paste('RipOffDensityChanges_DensityCorrected.png'), width=25, height=20, units="cm")

                                                          
