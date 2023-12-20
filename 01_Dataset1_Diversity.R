library(vegan)
library(ggplot2)

# changed x to df to be more stylisticly consistent
df<-read.csv('totalCounts_notEuk_taxLevel1_integer_zeros.csv', sep=",")

# make a data frame of the sample names
namesD <- as.data.frame(names(df[6:149]))
names(namesD)<-'sample'
rownames(df) <- df[,1] 
scounts <- df[,c(6:149)]

meta<-read.csv('Ferreira.metadata.csv',sep=',',header=T)
meta <- meta[,c(1,11,12,14)]
names(meta)<-c('sample', 'temperature', 'infection', 'time')
table(meta$time)


#treatment subsets
# subset the metadata sample names according to temperature and time
table(meta$temperature, meta$time)

c20_t24 <- meta$sample[meta$temperature == 20 & meta$time == 24]
c20_t48 <- meta$sample[meta$temperature == 20 & meta$time == 48]
c28_t24 <- meta$sample[meta$temperature == 28 & meta$time == 24]
c28_t48 <- meta$sample[meta$temperature == 28 & meta$time == 48]
c36_t24 <- meta$sample[meta$temperature == 36 & meta$time == 24]
c36_t48 <- meta$sample[meta$temperature == 36 & meta$time == 48]

# match the sample names to the column names in the count matrix
c20_t24 <- df[, c(match(c20_t24, colnames(df)))]
c20_t48 <- df[, c(match(c20_t48, colnames(df)))]
c28_t24 <- df[, c(match(c28_t24, colnames(df)))]
c28_t48 <- df[, c(match(c28_t48, colnames(df)))]
c36_t24 <- df[, c(match(c36_t24, colnames(df)))]
c36_t48 <- df[, c(match(c36_t48, colnames(df)))]
# head(c20_t24)

# dList <- list(c20_t24, c20_t48, c28_t24, c28_t48, c36_t24, c36_t48)
# names(dList) <- c("c20_t24", "c20_t48", "c28_t24", "c28_t48", "c36_t24", "c36_t48")
# sampleNames <- c(names(c20_t24), names(c20_t48), names(c28_t24), names(c28_t48), names(c36_t24), names(c36_t48))

# #shannon each subset
# shanLists <- lapply(dList, diversity, index="shannon")
# names(shanLists) <- c("c20_t24", "c20_t48", "c28_t24", "c28_t48", "c36_t24", "c36_t48")
# # names(c20_t24)
# # head(shanLists)

# #make a data frame of the shannon diversity values but keeping the sample names
# dfShan <- as.data.frame(shanLists)

# head(dfShan)

# make a histogram of each column in dfShan labled with that column name
# par(mfrow = c(3, 2))
# for (i in 1:6) {
#   hist(dfShan[, i],
#     main = names(dfShan)[i],
#     xlab = "Shannon's Index",
#   )
# }
# # save output as a pdf
# dev.copy2pdf(file = "Ferriera_shannonDiversity.pdf")

# make a ggplot2 boxplot of each column in single plot labled with that column name
par(mfrow = c(1, 1))

# df_long <- dfShan %>%
#   gather(key = "sample", value = "shannonDiversity")
# head(df_long) 

# ggplot(df_long, aes(x = sample, y = log(shannonDiversity))) +
#   geom_boxplot() +
#   ylab("Log Shannon's Index") +
#   xlab("Temperature and Time") +
#   theme_bw()

# # save
# ggsave("Ferriera_shannonDiversity_boxplot.pdf")

# get shannon, simpson, inverse simpson index into a single dataframe in the right format for glm
#shannon temperatures
shannondiv.T<-diversity(t(df[,6:149]), index="shannon")
simpsondiv.T <- diversity(t(df[, 6:149]), index = "simpson")
invsimpsondiv.T <- diversity(t(df[, 6:149]), index = "invsimpson")

# head(shannondiv.T)
# length(shannondiv.T)
# head(simpsondiv.T)
# length(simpsondiv.T)

# combine into a single dataframe
dfDiv <- as.data.frame(cbind(shannondiv.T, simpsondiv.T, invsimpsondiv.T))
names(dfDiv) <- c("shannon", "simpson", "invsimpson")
dfDiv$sample=rownames(dfDiv)
head(dfDiv)

# add sample metadata to the dataframe
dfDiv <- merge(dfDiv, meta, by = "sample")

# dfDiv$temperature <- factor(dfDiv$temperature, levels = c(20, 28, 36))
dfDiv$infection <- factor(dfDiv$infection, levels = c("Control", "Zika"))
# dfDiv$time <- factor(dfDiv$time, levels = c(24, 48))

shanGlm <- glm(shannon ~ temperature * infection * time, data = dfDiv)
shanAnovaTab <- anova(shanGlm, test = "F")
rownames(shanAnovaTab)[1] <- '(Intercept)'

# check assumptions of the test
par(mfrow = c(2, 2))
plot(shanGlm)

# save output as a pdf
# make plots dir
if(!dir.exists("plots/QC")){
  dir.create("plots/QC",
    recursive = TRUE
  )
  }

dev.copy2pdf(file = "plots/QC/Ferriera_shannonDiversity_glm_diagnostics.pdf")

simpGlm <- glm(simpson ~ temperature * infection * time, data = dfDiv)
simpAnovaTab <- anova(simpGlm, test = "F")
rownames(simpAnovaTab)[1] <- '(Intercept)'

# check assumptions of the test
plot(simpGlm)
dev.copy2pdf(file = "plots/QC/Ferriera_simpsonDiversity_glm_diagnostics.pdf")

invsimpGlm <- glm(invsimpson ~ temperature * infection * time, data = dfDiv)
invsimpAnovaTab <- anova(invsimpGlm, test = "F")
rownames(invsimpAnovaTab)[1] <- '(Intercept)'

# check assumptions of the test
plot(invsimpGlm)
dev.copy2pdf(file = "plots/QC/Ferriera_invsimpsonDiversity_glm_diagnostics.pdf")

# make a table of each of the results and save
if(!dir.exists("tables")){
  dir.create("tables",
    recursive = TRUE
  )
  }

# saving tables
shanTab<-knitr::kable(shanAnovaTab, format= 'pipe',digits = 3, caption = "ANOVA table for Shannon's diversity index")
kableExtra::save_kable(shanTab, file = "tables/Ferriera_shannonDiversity_glm_anova_table.txt")
# kableExtra::save_kable(shanTab, file = "tables/Ferriera_shannonDiversity_glm_anova_table.html") # would work if pandoc/latex properly installed
# kableExtra::save_kable(shanTab, file = "tables/Ferriera_shannonDiversity_glm_anova_table.pdf")

simpTab <- knitr::kable(simpAnovaTab, format = "pipe", digits = 3, caption = "ANOVA table for Simpson's diversity index")
kableExtra::save_kable(simpTab, file = "tables/Ferriera_simpsonDiversity_glm_anova_table.txt")

invSimpTab <- knitr::kable(invsimpAnovaTab, format = "pipe", digits = 3, caption = "ANOVA table for inverse Simpson's diversity index")
kableExtra::save_kable(invSimpTab, file = "tables/Ferriera_invsimpsonDiversity_glm_anova_table.txt")

# Plot effects of temperature, infection, and time on shannon diversity with ggplot in boxplots
p <- ggplot(dfDiv, aes(x = as.factor(temperature), y = shannon, fill = c(infection))) +
  scale_fill_manual("Infection", values = c("cornflowerblue", "firebrick")) +
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Shannon's Index") +
  xlab("Temperature") +
  theme_bw()
p

# ns but approaches. Save
ggsave("plots/Ferriera_shannonDiversity_boxplot_temp_infection.pdf")

p <- ggplot(dfDiv, aes(x = as.factor(temperature), y = shannon, fill = c(time))) +
  scale_fill_manual("DPI", values = c("#bdc4cc", "#87a1db")) +
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Shannon's Index") +
  xlab("Temperature") +
  theme_bw()
p
ggsave("plots/Ferriera_shannonDiversity_boxplot_temp_time.pdf")

p <- ggplot(dfDiv, aes(x = as.factor(time), y = shannon, fill = c(infection))) +
  scale_fill_manual("Infection", values = c("cornflowerblue", "firebrick")) +
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Shannon's Index") +
  xlab("DPI") +
  theme_bw()
p
ggsave("plots/Ferriera_shannonDiversity_boxplot_time_infection.pdf")

# plot and save significant interaction for simpsons
p <- ggplot(dfDiv, aes(x = as.factor(temperature), y = simpson, fill = c(infection))) +
  scale_fill_manual("Infection", values = c("cornflowerblue", "firebrick")) +
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Simpsons's Index") +
  xlab("Temperature") +
  theme_bw()
p

# ns but approaches. Save
ggsave("plots/Ferriera_simpsonEvenness_boxplot_temp_infection.pdf")

p <- ggplot(dfDiv, aes(x = as.factor(temperature), y = simpson, fill = c(time))) +
  scale_fill_manual("DPI", values = c("#bdc4cc", "#87a1db")) +
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Shannon's Index") +
  xlab("Temperature") +
  theme_bw()
p
ggsave("plots/Ferriera_simpsonEvenness_boxplot_temp_time.pdf")

p <- ggplot(dfDiv, aes(x = as.factor(time), y = simpson, fill = c(infection))) +
  scale_fill_manual("Infection", values = c("cornflowerblue", "firebrick")) +
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Shannon's Index") +
  xlab("DPI") +
  theme_bw()
p
ggsave("plots/Ferriera_simpsonEvenness_boxplot_time_infection.pdf")


# plot and save significant interaction for inverse simpsons
p <- ggplot(dfDiv, aes(x = as.factor(time), y = invsimpson, fill = c(infection))) +
  scale_fill_manual("Infection", values = c("cornflowerblue", "firebrick")) +
  geom_boxplot() +
  #   geom_dotplot(binaxis='y', dotsize=0.4)+
  ylab("Inverse Simpson's Index") +
  xlab("DPI") +
  theme_bw()
p
ggsave("plots/Ferriera_invsimpsonDiversity_boxplot_time_infection.pdf")


# ---
# indidual level plots. Each metric is sig for temperature, infection, and time
# ---

# put all three on one plot
p1 <- ggplot(dfDiv, aes(x = as.factor(temperature), y = shannon)) +
  geom_boxplot() +
  ylab("Shannon's Index") +
  xlab("Temperature") +
  theme_bw()
p1
p2 <- ggplot(dfDiv, aes(x = as.factor(temperature), y = simpson)) +
  geom_boxplot() +
  ylab("Simpson's Index") +
  xlab("Temperature") +
  theme_bw()
p2
p3 <- ggplot(dfDiv, aes(x = as.factor(temperature), y = invsimpson)) +
  geom_boxplot() +
  ylab("Inverse Simpson's Index") +
  xlab("Temperature") +
  theme_bw()
p3

p4 <- ggplot(dfDiv, aes(x = as.factor(infection), y = shannon)) +
  geom_boxplot() +
  ylab("Shannon's Index") +
  xlab("Infection") +
  theme_bw()
p4
p5 <- ggplot(dfDiv, aes(x = as.factor(infection), y = simpson)) +
  geom_boxplot() +
  ylab("Simpson's Index") +
  xlab("Infection") +
  theme_bw()
p5
p6 <- ggplot(dfDiv, aes(x = as.factor(infection), y = invsimpson)) +
  geom_boxplot() +
  ylab("Inverse Simpson's Index") +
  xlab("Infection") +
  theme_bw()
p6
p7 <- ggplot(dfDiv, aes(x = as.factor(time), y = shannon)) +
  geom_boxplot() +
  ylab("Shannon's Index") +
  xlab("DPI") +
  theme_bw()
p7
p8 <- ggplot(dfDiv, aes(x = as.factor(time), y = simpson)) +
  geom_boxplot() +
  ylab("Simpson's Index") +
  xlab("DPI") +
  theme_bw()
p8
p9 <- ggplot(dfDiv, aes(x = as.factor(time), y = invsimpson)) +
  geom_boxplot() +
  ylab("Inverse Simpson's Index") +
  xlab("DPI") +
  theme_bw()
p9

# plot together
library(gridExtra)

all <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
# save
ggsave("plots/Ferriera_diversity_boxplot_temp_infection_time.pdf", all, width = 20, height = 20, units = "cm")

# save the dataframes
write.csv(dfDiv, "Ferriera_diversityIndex.csv", row.names = FALSE)
