library(vegan)
library(ggplot2)

x<-read.csv('Wimalasiri_totalCounts_notEuk_taxLevel1.csv', sep=",")
# replace NAs with 0
x[, 7:75][is.na(x[, 7:75])] <- 0

namesD <- as.data.frame(names(x[7:75]))
names(namesD)<-'sample'

rownames(x) <- x[,1] 
scounts <- x[,c(7:75)]

meta <- read.csv("Wimalasiri-Yapa.metadata.csv", sep = ",", header = T)
meta <- meta[, c(1, 11, 12, 13, 15)]
names(meta) <- c("sample", "time", "temperature", "infection", "treatment")
table(meta$time)
table(meta$temperature)
infection <- factor(meta$infection, levels = c("Control", "Chikungunya"))
table(infection)
tempTreat <- factor(meta$temperature)

#transpose so samples are rows
t_scounts <- as.data.frame(t(scounts))
#rarefy data
min_depth = min(colSums(scounts))
t_scounts_rarefied <- as.data.frame(round(rrarefy(t_scounts, min_depth)))
#sqr root transform rarefied table and determine best method calculating distance matrix
sqrt_t_scounts_rarefied = sqrt(t_scounts_rarefied)
rank.tscounts <- rankindex(as.matrix(sqrt_t_scounts_rarefied), t_scounts_rarefied, indices = c("bray", "euclid", "manhattan", "horn"), method = "spearman")
print(paste("The highest rank was given by the", names(sort(rank.tscounts, decreasing = TRUE)[1]), "method."))
scounts_dist = as.matrix((vegdist(t_scounts_rarefied, "bray")))

#NMDS
NMDS = metaMDS(scounts_dist)

#build a data frame with NMDS coordinates and metadata
#but first factor time and temp
meta$time=factor(meta$time,labels=c("3","7"))
meta$temperature=factor(meta$temperature,labels=c("18","28","32"))
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, infection=infection, time=meta$time, temperature=meta$temperature)
head(NMDS)
#how do samples fall in ordination space?
ggplot(NMDS, aes(x=MDS1, y=MDS2, col=infection)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "NMDS Plot")
ggsave("plots/Wimalasiri_NMDS_infection.png", width = 5, height = 5, units = "in", dpi = 300)

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=time)) +
  geom_point() +
  stat_ellipse()+
  theme_bw() +
  labs(title = "NMDS Plot")
ggsave("plots/Wimalasiri_NMDS_time.png", width = 5, height = 5, units = "in", dpi = 300)

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=temperature)) +
  geom_point() +
  stat_ellipse()+
  theme_bw() +
  labs(title = "NMDS Plot")
ggsave("plots/Wimalasiri_NMDS_temperature.png", width = 5, height = 5, units = "in", dpi = 300)

#ANOSIM
#infection
anosim_infection = anosim(scounts_dist, infection)
anosim_infection # take a look at results
summary(anosim_infection)
plot(anosim_infection, main = "ANOSIM for infection", xlab = "Infection", ylab = "ANOSIM Rank")
pdf("plots/Wimalasiri_ANOSIM_infection.pdf")
dev.off()

#time
anosim_time = anosim(scounts_dist, meta$time)
anosim_time # take a look at results
summary(anosim_time)
plot(anosim_time, main = "ANOSIM for time", xlab = "Time", ylab = "ANOSIM Rank")
pdf("plots/Wimalasiri_ANOSIM_time.pdf")
dev.off()

#temperature
anosim_temperature = anosim(scounts_dist, meta$temperature)
anosim_temperature # take a look at results
summary(anosim_temperature)

plot(anosim_temperature, main = "ANOSIM for temperature", xlab = "Temperature", ylab = "ANOSIM Rank") 
pdf("plots/Wimalasiri_ANOSIM_temperature.pdf")
dev.off()


#permANOVA
#infection
adonis_infection = adonis2(scounts_dist ~ infection, meta)
adonis_infection
summary(adonis_infection)
infTab <- knitr::kable(adonis_infection, format = "pipe", digits = 3, caption = "permANOVA table for Wimalasiri infection")
kableExtra::save_kable(infTab, file = "tables/Wimalasiri_inf_peranova_table.txt")
# significance for infection

#time
adonis_time = adonis2(scounts_dist ~ time, meta)
adonis_time 
summary(adonis_time)
timeTab <- knitr::kable(adonis_time, format = "pipe", digits = 3, caption = "permANOVA table for Wimalasiri time")
kableExtra::save_kable(infTab, file = "tables/Wimalasiri_time_peranova_table.txt")

# plot(adonis_time)
# not significant for time

#temperature
adonis_temperature = adonis2(scounts_dist ~ temperature, meta)
adonis_temperature 
summary(adonis_temperature)

tempTab <- knitr::kable(adonis_temperature, format = "pipe", digits = 3, caption = "permANOVA table for Wimalasiri temperature")
kableExtra::save_kable(infTab, file = "tables/Wilmasiri_temp_peranova_table.txt")

# plot(adonis_temperature)
# approaches sig for temperature