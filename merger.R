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