# Load packages
# install.packages("multcomp")
library("multcomp")

# load files and turn into factors
logFCs<-read.csv('species_logFC.csv', header=T)
logFCs$Species<-factor(logFCs$Species, 
                       levels=c("Enterobacter ludwigii",
                                "Gardnerella vaginalis",
                                "Rothia mucilaginosa",
                                "Enterobacter roggenkampii"))

logFCs$Temperature<-factor(logFCs$Temperature,levels=c("18","20","25","28","32","35","36"))

sumlogFCs<-summary(logFCs)

# 
plot(logFC~Species+Temperature,data=logFCs)

# Run anova model for temperature. 
logFCsaov<-aov(logFC~Temperature, data=logFCs)

# TUKEY HSD
tuktab<-TukeyHSD(logFCsaov)

# Plot HSD 
with(par(mai=c(1,1,1,1)),{plot(tuktab, las=1, cex.axis=1, srt=45)})
