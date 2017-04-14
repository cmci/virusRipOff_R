ylabeltext = "Rip-Off-Density [counts/100um^2]"
xordering = c("Control1", "Control2", "ROCKinhibitor", "Blebb10um", "ClaD06", "CK666", "KOCLCAB", "SMIFH2", "Jaspla", "Genistein", "bcyclo", "a5b1", "beta1ABp5d2")
ggplot(alldata, aes(x=type, y=RipOff_Density)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext)

ylabeltext2 = "Rip-Off-Density (Before Adding Drug)\n[counts/100um^2]"
ggplot(alldata, aes(x=type, y=RipOffPreDensity)) + geom_boxplot() + geom_jitter(width=0.2) + scale_x_discrete(name="Experiments", limits=xordering) + scale_y_continuous(name=ylabeltext2, limits = c(0, 7))
