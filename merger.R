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
  "170328 beta1ABp5d2 croped and 14 frames"
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
  "beta1ABp5d2"
)

dflist <- vector("list", length(datanames))
names(dflist) <- datanames
for (ind in 1:length(datanames)){
  dflist[[datanames[ind]]] <- read.csv2( file.path(gppath, folders[[ind]], "results.csv"), sep=",", dec=".")
  dflist[[datanames[ind]]]$type = rep(datanames[[ind]], nrow(dflist[[datanames[ind]]]))
}

ctrl1 = read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/161122 ctrl croped and 16 frames/results.csv', sep=",", dec=".")
ctrl2 <- read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/161206 ctrl croped and 16 frames/results.csv', sep=',', dec=".")
ctrl1$type = rep("Control1", 5)
ctrl2$type = rep("Control2", 8)
ROCKinh <- read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/161129 ROCKinh croped and 16 frames/results.csv', sep=',', dec='.')
ROCKinh$type <- rep("ROCKinhibitor", 5)
Blebb10um <- read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/161205 10umblebb croped and 16 frames/results.csv', sep=',', dec='.')
Blebb10um$type = rep("Blebb10um", 5)
ClaD06 = read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/161222 claD0.06 croped and 16 frames/results.csv', sep=',', dec='.')
ClaD06$type = rep('ClaD06', 3)
CK666=read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/170130 CK666 croped and 16 frames/results.csv', sep=',', dec='.')
CK666$type=rep('CK666', 6)
KOCLCAB = read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/170206 KOCLCAB croped and 50 frames/results.csv', sep=',', dec='.')
KOCLCAB$type = rep('KOCLCAB', 8)
SMIFH2 = read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/170215 SMIFH2 croped and 16 frames/results.csv', sep=',', dec = '.')
SMIFH2$type=rep('SMIFH2', 7)
Jaspla = read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/170217 Jasplakinolide croped and 16 frames/results.csv', sep=',', dec='.')
Jaspla$type = rep('Jaspla', 1)
Genistein = read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/170222 Genistein croped and 16 frames/results.csv', sep=',', dec='.')
Genistein$type = rep('Genistein', 11)
bcyclo=read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/170302 bcyclo croped and 16 frames/results.csv', sep = ',', dec = '.')
bcyclo$type = rep('bcyclo', 9)
a5b1 = read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/170307 a5b1peptidomimetic croped and 16 frames/results.csv', sep=',', dec = '.')
a5b1$type = rep('a5b1', 3)
beta1ABp5d2 = read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/170328 beta1ABp5d2 croped and 14 frames/results.csv', sep = ',', dec='.')
beta1ABp5d2$type = rep('beta1ABp5d2', 7)

alldata <- merge(ctrl1, ctrl2, all=TRUE)
alldata <- merge(alldata, ROCKinh, all = TRUE)
alldata <- merge(alldata, Blebb10um, all = TRUE)
alldata <- merge(alldata, ClaD06, all = TRUE)
alldata <- merge(alldata, CK666, all = TRUE)
alldata <- merge(alldata, KOCLCAB, all = TRUE)
alldata <- merge(alldata, SMIFH2, all = TRUE)
alldata <- merge(alldata, Jaspla, all = TRUE)
alldata <- merge(alldata, Genistein, all = TRUE)
alldata <- merge(alldata, bcyclo, all = TRUE)
alldata <- merge(alldata, a5b1, all = TRUE)
alldata <- merge(alldata, beta1ABp5d2, all = TRUE)
