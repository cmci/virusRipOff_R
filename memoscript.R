ctrl1 = read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/161122 ctrl croped and 16 frames/results.csv', sep=",", dec='.')
ctrl1$type <- rep('Control1', 5)
ctrl2 <- read.csv2('/Users/miura/Dropbox/people/Tina/shared_Tina_Kota/Data Tina/161206 ctrl croped and 16 frames/results.csv', sep=',', dec='.')
ctrl2$type <- rep('Control2', 8)
controlAll = merge(ctrl1, ctrl2, all=TRUE)
p <- ggplot(controlAll, aes(type, RipOff.Ratio))
pp <- ggplot(controlAll, aes(x=Density,y=RipOff.Ratio))
pp + geom_point() + scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25))

