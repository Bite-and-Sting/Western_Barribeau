library(vegan)
library(ggplot2)

x<-read.csv('totalCounts_notEuk_taxLevel1_integer_zeros.csv', sep=",")
# still some NAs
x[,6:149][is.na(x[,6:149])] <- 0
namesD <- as.data.frame(names(x[6:149]))
names(namesD)<-'sample'
rownames(x) <- x[, 1]
dim(x)

scounts <- x[,c(6:149)]
# head(x)
samples <- names(x)[6:length(names(x))]
samples

meta <- read.csv("Ferreira.metadata.csv", sep = ",", header = T)
meta <- meta[, c(1, 11, 12, 14, 16, 17)]
names(meta) <- c("sample", "temperature", "infection", "time", "treatment", "inf-temp")
table(meta$time)
head(meta)

# transpose counts
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
meta$time=factor(meta$time,labels=c("24","48"))
table(meta$temperature)
meta$temperature=factor(meta$temperature,labels=c("20","28","36"))
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, infection=meta$infection, time=meta$time, temperature=meta$temperature)
head(NMDS)
# how do samples fall in ordination space?
ggplot(NMDS, aes(x = MDS1, y = MDS2, col = infection)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    labs(title = "NMDS Plot")
# save ggplot
ggsave("plots/Ferriera_NMDS_infection.png", width = 5, height = 5, units = "in", dpi = 300)

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=time)) +
  geom_point() +
  stat_ellipse()+
  theme_bw() +
  labs(title = "NMDS Plot")
ggsave("plots/Ferriera_NMDS_time.png", width = 5, height = 5, units = "in", dpi = 300)

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=temperature)) +
  geom_point() +
  stat_ellipse()+
  theme_bw() +
  labs(title = "NMDS Plot")
ggsave("plots/Ferriera_NMDS_temperature.png", width = 5, height = 5, units = "in", dpi = 300)

#ANOSIM manually saved to a text file.
#infection
anosim_infection = anosim(scounts_dist, meta$infection)
anosim_infection # take a look at results
summary(anosim_infection)
plot(anosim_infection, main = "ANOSIM for infection", xlab = "Infection", ylab = "ANOSIM Rank")
# save plot
pdf("plots/Ferriera_ANOSIM_infection.pdf")
dev.off()

#time
anosim_time = anosim(scounts_dist, meta$time)
anosim_time # take a look at results
summary(anosim_time)

plot(anosim_time, main = "ANOSIM for time", xlab = "Time", ylab = "ANOSIM Rank")
pdf("plots/Ferriera_ANOSIM_time.pdf")
dev.off()

#temperature
anosim_temperature = anosim(scounts_dist, meta$temperature)
anosim_temperature # take a look at results
summary(anosim_temperature)
plot(anosim_temperature, main = "ANOSIM for temperature", xlab = "Temperature", ylab = "ANOSIM Rank")
pdf("plots/Ferriera_ANOSIM_temperature.pdf")
dev.off()
#all of these showed significance which I think means taht our communities were different frome achother

#permANOVA
#infection
adonis_infection = adonis2(scounts_dist ~ infection, meta)
adonis_infection
summary(adonis_infection)
infTab <- knitr::kable(adonis_infection, format = "pipe", digits = 3, caption = "permANOVA table for Ferriera infection")
kableExtra::save_kable(infTab, file = "tables/Ferriera_inf_peranova_table.txt")


# plot(adonis_infection)

# infection p 0.094
#time
adonis_time = adonis2(scounts_dist ~ time, meta)
adonis_time 
summary(adonis_time)
timeTab <- knitr::kable(adonis_time, format = "pipe", digits = 3, caption = "permANOVA table for Ferriera time")
kableExtra::save_kable(infTab, file = "tables/Ferriera_time_peranova_table.txt")


# plot(adonis_time)

#significant for time
#temperature
adonis_temperature = adonis2(scounts_dist ~ temperature, meta)
adonis_temperature 
summary(adonis_temperature)
tempTab <- knitr::kable(adonis_temperature, format = "pipe", digits = 3, caption = "permANOVA table for Ferriera temperature")
kableExtra::save_kable(infTab, file = "tables/Ferriera_temp_peranova_table.txt")
# plot(adonis_temperature)
#significant for temperature
