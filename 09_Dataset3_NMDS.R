library(vegan)
library(ggplot2)

x<-read.csv('Murrieta_nrnt_notEuk.csv', sep=",")
# still some NAs
x[, 6:186][is.na(x[, 6:186])] <- 0

namesD <- as.data.frame(names(x[6:186]))
names(namesD)<-'sample'

rownames(x) <- x[,1] 
scounts <- x[,c(6:186)]

meta<- read.csv('Murrieta.metadata.csv', sep= ',', header=T)
meta <- meta[,c(1,12)]
names(meta)<-c('sample', 'temperature')
table(meta$temperature)
tempTreat <- factor(meta$temperature)

#remove 25-35 variable treatment
meta <- meta[,c(1,2)]
names(meta)<-c('sample', 'temperature')
meta2=meta[meta$temperature != '25-35',]
names(meta2)
table(meta2$temperature)
keep = scounts[, names(scounts) %in% meta2$sample] # here we drop all the columns that aren't in our new metadata
metakeep <- meta2[meta2$sample %in% names(keep),] # here we drop all the rows that aren't in our new metadata

#transpose so samples are rows
t_scounts <- as.data.frame(t(keep))
#rarefy data
min_depth = min(colSums(keep))
t_scounts_rarefied <- as.data.frame(round(rrarefy(t_scounts, min_depth)))
#sqr root transform rarefied table and determine best method calculating distance matrix
sqrt_t_scounts_rarefied = sqrt(t_scounts_rarefied)
rank.tscounts <- rankindex(as.matrix(sqrt_t_scounts_rarefied), t_scounts_rarefied, indices = c("bray", "euclid", "manhattan", "horn"), method = "spearman")
print(paste("The highest rank was given by the", names(sort(rank.tscounts, decreasing = TRUE)[1]), "method."))
scounts_dist = as.matrix((vegdist(t_scounts_rarefied, "horn")))

#NMDS
NMDS = metaMDS(scounts_dist)

#build a data frame with NMDS coordinates and metadata
#but first factor time and temp
meta2$temperature=factor(meta2$temperature,labels=c("25","28","32","35"))
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, temperature=meta2$temperature)
head(NMDS)
#how do samples fall in ordination space?
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=temperature)) +
  geom_point() +
  stat_ellipse()+
  theme_bw() +
  labs(title = "NMDS Plot")
ggsave("plots/Murrieta_NMDS_infection.png", width = 5, height = 5, units = "in", dpi = 300)


#ANOSIM
#temperature
anosim_temperature = anosim(scounts_dist, meta2$temperature)
anosim_temperature # take a look at results
summary(anosim_temperature)
plot(anosim_temperature)
#all of these showed significance which I think means taht our communities were different frome achother
pdf("plots/Murrieta_ANOSIM_temperature.pdf")
dev.off()

#permANOVA
#temperature
adonis_temperature = adonis2(scounts_dist ~ temperature, metakeep)
adonis_temperature 
summary(adonis_temperature)
plot(adonis_temperature)
#significant for temperature

tempTab <- knitr::kable(adonis_temperature, format = "pipe", digits = 3, caption = "permANOVA table for Murrieta temperature")
kableExtra::save_kable(tempTab, file = "tables/Murrieta_temp_peranova_table.txt")
