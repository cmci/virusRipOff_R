# import time course of rip-off counts. 

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


datalist = list(
  "ctrl1" = seq(5),
  "ctrl2"= seq(8),
  "ROCKinh" = c(0, 1, 2, 5, 6) + 1,
  "Blebb10um" = c(0, 2, 3, 6, 7) + 1,
  "ClaD06" = c(0, 1, 2) + 1,
  "CK666"= c(0, 1, 2, 3, 7, 8) + 1, 
  "KOCLCAB" = c(1, 2, 4, 5, 6, 7, 8, 9) + 1,
  "SMIFH2" = c(1, 2, 3, 4, 5, 6, 7) + 1,
  "Jaspla" = c(4) + 1,
  "Genistein" = seq(11),
  "bcyclo" = seq(9),
  "a5b1" = c(0, 1, 9) + 1,
  "beta1ABp5d2" = c(1, 2, 3, 4, 5, 6, 9) + 1,
  "Hela" = c(0, 1, 2, 3, 4, 5, 6, 8, 9) + 1
)

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
# 
# typemax = max(timecourseAll$RipOffDensity[timecourseAll$Type=='ctrl2'])
# annoY = typemax/2
# ggplot(data=filter(timecourseAll, Type=="ctrl2")) + 
#   geom_path(aes(x=Frame, y=RipOffDensity, color=CellID)) + 
#   scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + 
#   scale_y_continuous(name='Rip-Off Density \n[counts/100um^2]') + 
#   annotate("text", x=12, y=annoY, label= "ctrl2", size=30, color='gray')
# ggsave("TimeCourse_RipOffDensity_ctrl2.png", width=25, Height=20, units="cm")

ggplot(data=timecourseAll) + 
  geom_path(aes(x=Frame, y=Counts, color=CellID)) + 
  scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + 
  scale_y_continuous(name='Rip-Off [Counts]') + 
  facet_wrap(~ Type, nrow = 4)
ggsave('TimeCourse_Counts.png', width=25, height=20, units="cm")

for (thename in names(datalist)){
  typemax = max(timecourseAll$Counts[timecourseAll$Type==thename])
  annoY = typemax/2
  ggplot(data=filter(timecourseAll, Type==thename)) + 
    geom_path(aes(x=Frame, y=Counts, color=CellID)) + 
    scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + 
    scale_y_continuous(name='Rip-Off [Counts]') +
    annotate("text", x=12, y=annoY, label= thename, size=30, color='gray')
  ggsave(paste("TimeCourse_Counts_",thename,".png", sep=''), width=25, height=20, units="cm")
}

ggplot(data=timecourseAll) + geom_path(aes(x=Frame, y=RipOffDensity, color=CellID)) + scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + scale_y_continuous(name='Rip-Off Density \n[counts/100um^2]') + facet_wrap(~ Type, nrow = 4)
ggsave(paste('TimeCourse_RipOffDensity.png'), width=25, height=20, units="cm")

for (thename in names(datalist)){
  typemax = max(timecourseAll$RipOffDensity[timecourseAll$Type==thename])
  annoY = typemax/2
  ggplot(data=filter(timecourseAll, Type==thename)) + 
    geom_path(aes(x=Frame, y=RipOffDensity, color=CellID)) + 
    scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + 
    scale_y_continuous(name='Rip-Off Density \n[counts/100um^2]') +
    annotate("text", x=12, y=annoY, label= thename, size=30, color='gray')
  ggsave(paste("TimeCourse_RipOffDensity_",thename,".png", sep=''), width=25, height=20, units="cm")
}

ggplot(data=timecourseAll) + geom_path(aes(x=Frame, y=VirusDensity, color=CellID)) + scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + scale_y_continuous(name='Virus Density \n[counts/100um^2]') + facet_wrap(~ Type, nrow = 4)
ggsave(paste('TimeCourse_VirusDensity.png'), width=25, height=20, units="cm")

for (thename in names(datalist)){
  typemax = max(timecourseAll$VirusDensity[timecourseAll$Type==thename])
  annoY = typemax/2
  ggplot(data=filter(timecourseAll, Type==thename)) + 
    geom_path(aes(x=Frame, y=VirusDensity, color=CellID)) + 
    scale_x_continuous(name = 'Time [Frame]', limits=c(1, 16)) + 
    scale_y_continuous(name='Virus Density \n[counts/100um^2]', limits=c(0, NA)) +
    annotate("text", x=12, y=annoY, label= thename, size=30, color='gray')
  ggsave(paste("TimeCourse_VirusDensity_",thename,".png", sep=''), width=25, height=20, units="cm")
}


#ggplot(data=filter(timecourseAll, Type=="ctrl2")) + geom_path(aes(x=Frame, y=RipOffDensity, color=CellID))


