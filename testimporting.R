inputcsv = file.choose()
tt = read.csv(inputcsv)
gppath = as.character(tt$root[1])
folders = as.character(tt$datafolder)
datanames = as.character(tt$type)
dflist <- vector("list", length(datanames))
names(dflist) <- datanames
datalist <- vector("list", length(datanames))
names(datalist) <- datanames
for (ind in 1:length(datanames)){
  print(datanames[ind])
  print(file.path(gppath, folders[[ind]], "results.csv"))
  dflist[[datanames[ind]]] <- read.csv( file.path(gppath, folders[[ind]], "results.csv"))
  dflist[[datanames[ind]]]$type = rep(datanames[[ind]], nrow(dflist[[datanames[ind]]]))
  ids = tt[ind, 4:ncol(tt)]
  ids = ids[!is.na(ids)]
  datalist[[datanames[ind]]] <- ids
}
