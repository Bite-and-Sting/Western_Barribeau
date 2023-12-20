# setwd("~/4th Year Modules/LIFE700/Wimalasiri-Yapa")
library(vegan)
library(ggplot2)

df<-read.csv('Wimalasiri_totalCounts_notEuk_taxLevel1.csv', sep=",")
# still some NAs
df[, 7:75][is.na(df[, 7:75])] <- 0

namesD <- as.data.frame(names(df[7:75]))
names(namesD)<-'sample'
rownames(df) <- df[,1] 
scounts <- df[,c(7:75)]
head(scounts)
# # Count per millions
# scountsPerMillion <- cpm(scounts)
# summary(scountsPerMillion)

meta <- read.csv('Wimalasiri-Yapa.metadata.csv', sep= ',', header=T)
meta <- meta[,c(1,11,12,13)]
names(meta)<-c('sample', 'time', 'temperature', 'infection')
meta$time <- as.factor(meta$time)
meta$temperature <- as.factor(meta$temperature)
meta$infection <- as.factor(meta$infection, levels=c("Control","Chikungunya"))
table(meta$time)
table(meta$temperature)
table(meta$infection)

table(meta$temperature, meta$infection, meta$time)

# get the sample names for each subset
c18_Control_d3 <- meta$sample[meta$temperature == "18" & meta$infection == "Control" & meta$time == "3"]
c18_Control_d7 <- meta$sample[meta$temperature == "18" & meta$infection == "Control" & meta$time == "7"]
c18_Chik_d3 <- meta$sample[meta$temperature == "18" & meta$infection == "Chikungunya" & meta$time == "3"]
c18_Chik_d7 <- meta$sample[meta$temperature == "18" & meta$infection == "Chikungunya" & meta$time == "7"]

c28_Control_d3 <- meta$sample[meta$temperature == "28" & meta$infection == "Control" & meta$time == "3"]
c28_Control_d7 <- meta$sample[meta$temperature == "28" & meta$infection == "Control" & meta$time == "7"]
c28_Chik_d3 <- meta$sample[meta$temperature == "28" & meta$infection == "Chikungunya" & meta$time == "3"]
c28_Chik_d7 <- meta$sample[meta$temperature == "28" & meta$infection == "Chikungunya" & meta$time == "7"]

c32_Control_d3 <- meta$sample[meta$temperature == "32" & meta$infection == "Control" & meta$time == "3"]
c32_Control_d7 <- meta$sample[meta$temperature == "32" & meta$infection == "Control" & meta$time == "7"]
c32_Chik_d3 <- meta$sample[meta$temperature == "32" & meta$infection == "Chikungunya" & meta$time == "3"]
c32_Chik_d7 <- meta$sample[meta$temperature == "32" & meta$infection == "Chikungunya" & meta$time == "7"]

# match the sample names to the column names in the count matrix
c18_Control_d3 <- df[,match(c18_Control_d3, colnames(df))]
c18_Control_d7 <- df[, match(c18_Control_d7, colnames(df))]
c18_Chik_d3 <- df[,match(c18_Chik_d3, colnames(df))]
c18_Chik_d7 <- df[, match(c18_Chik_d7, colnames(df))]

c28_Control_d3 <- df[, match(c28_Control_d3, colnames(df))]
c28_Control_d7 <- df[, match(c28_Control_d7, colnames(df))]
c28_Chik_d3 <- df[, match(c28_Chik_d3, colnames(df))]
c28_Chik_d7 <- df[, match(c28_Chik_d7, colnames(df))]

c32_Control_d3 <- df[, match(c32_Control_d3, colnames(df))]
c32_Control_d7 <- df[, match(c32_Control_d7, colnames(df))]
c32_Chik_d3 <- df[, match(c32_Chik_d3, colnames(df))]
c32_Chik_d7 <- df[, match(c32_Chik_d7, colnames(df))]

# dList <- list(
#   c18_Control_d3, c18_Control_d7, 
#   c18_Chik_d3, c18_Chik_d7, 
#   c28_Control_d3, c28_Control_d7, 
#   c28_Chik_d3, c28_Chik_d7, 
#   c32_Control_d3, c32_Control_d7, 
#   c32_Chik_d3, c32_Chik_d7
# )

# names(dList) <- c(
#   "c18_Control_d3", "c18_Control_d7", 
#   "c18_Chik_d3", "c18_Chik_d7", 
#   "c28_Control_d3", "c28_Control_d7", 
#   "c28_Chik_d3", "c28_Chik_d7", 
#   "c32_Control_d3", "c32_Control_d7", 
#   "c32_Chik_d3", "c32_Chik_d7"
# )

# sampleNames <- c(
#   names(c18_Control_d3), names(c18_Control_d7), 
#   names(c18_Chik_d3), names(c18_Chik_d7), 
#   names(c28_Control_d3), names(c28_Control_d7), 
#   names(c28_Chik_d3), names(c28_Chik_d7), 
#   names(c32_Control_d3), names(c32_Control_d7), 
#   names(c32_Chik_d3), names(c32_Chik_d7)
# )

# # shannon each subset
# shanLists <- lapply(dList, diversity, index = "shannon")
# names(shanLists) <- c(
#   "c18_Control_d3", "c18_Control_d7", 
#   "c18_Chik_d3", "c18_Chik_d7", 
#   "c28_Control_d3", "c28_Control_d7", 
#   "c28_Chik_d3", "c28_Chik_d7", 
#   "c32_Control_d3", "c32_Control_d7", 
#   "c32_Chik_d3", "c32_Chik_d7"
# )



# get right way around for glm
# shannon temperatures
shannondiv.T <- diversity(t(df[, 7:75]), index = "shannon")
simpsondiv.T <- diversity(t(df[, 7:75]), index = "simpson")
invsimpsondiv.T <- diversity(t(df[, 7:75]), index = "invsimpson")

# combine into a single dataframe
dfDiv <- as.data.frame(cbind(shannondiv.T, simpsondiv.T, invsimpsondiv.T))
names(dfDiv) <- c("shannon", "simpson", "invsimpson")
dfDiv$sample <- rownames(dfDiv)
head(dfDiv)
str(dfDiv)

# add sample metadata to the dataframe
dfDiv <- merge(dfDiv, meta, by = "sample")
dfDiv$infection <- factor(meta$infection, levels = c("Control", "Chikungunya"))
str(dfDiv)

shanGlm <- glm(shannon ~ temperature * infection * time, data = dfDiv)
shanAnovaTab <- anova(shanGlm, test = "F")
rownames(shanAnovaTab)[1] <- "(Intercept)"
# ns all
par(mfrow = c(2, 2))
plot(shanGlm)
dev.copy2pdf(file = "plots/QC/Wimalasiri_shannonDiversity_glm_diagnostics.pdf")

# simpson index glm
simpGlm <- glm(simpson ~ temperature * infection * time, data = dfDiv)
simpAnovaTab <- anova(simpGlm, test = "F")
plot(simpGlm)
dev.copy2pdf(file = "plots/QC/Wimalasiri_simpsonDiversity_glm_diagnostics.pdf")

# invsimpson index glm
invsimpGlm <- glm(invsimpson ~ temperature * infection * time, data = dfDiv)
invsimpAnovaTab <- anova(invsimpGlm, test = "F")
plot(invsimpGlm)
dev.copy2pdf(file = "plots/QC/Wimalasiri_invsimpsonDiversity_glm_diagnostics.pdf")


# save tables
shanTab <- knitr::kable(shanAnovaTab, format = 'pipe',  digits = 3, caption = "Shannon diversity index glm anova table")
kableExtra::save_kable(shanTab, file = "tables/Wimalasiri_shannonDiversity_glm_anovaTable.md")
simpTab <- knitr::kable(simpAnovaTab, format = 'pipe',  digits = 3, caption = "Simpson diversity index glm anova table")
kableExtra::save_kable(simpTab, file = "tables/Wimalasiri_simpsonDiversity_glm_anovaTable.md")
invsimpTab <- knitr::kable(invsimpAnovaTab, format = 'pipe',  digits = 3, caption = "Inverse Simpson diversity index glm anova table")
kableExtra::save_kable(invsimpTab, file = "tables/Wimalasiri_invsimpsonDiversity_glm_anovaTable.md")


# no significant effects
# plot individual plots on one page
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

all <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
ggsave("plots/Wimalasiri_diversityIndex_boxplots.pdf", all, width = 20, height = 20, units = "cm")

# save the dataframes
write.csv(dfDiv, "Wimalasiri_diversityIndex.csv", row.names = FALSE)

# # plot and save significant interaction for simpsons
# p1 <- ggplot(dfDiv, aes(x = as.factor(time), y = simpson)) +
#   # scale_fill_manual(values = c("cornflowerblue", "firebrick")) +
#   geom_boxplot() +
#   #   geom_dotplot(binaxis='y', dotsize=0.4)+
#   ylab("Simpsons's Index") +
#   xlab("DPI") +
#   theme_bw()
# p1

# p2 <- ggplot(dfDiv, aes(x = as.factor(temperature), y = simpson, fill = c(infection))) +
#   scale_fill_manual("Infection", values = c("cornflowerblue", "firebrick")) +
#   geom_boxplot() +
#   #   geom_dotplot(binaxis='y', dotsize=0.4)+
#   ylab("Simpsons's Index") +
#   xlab("Temperature") +
#   # legend("right", title = "Infection", c("Control", "Chikungunya")) +
#   theme_bw()
# p2


# p3 <- ggplot(dfDiv, aes(x = as.factor(time), y = simpson, fill = c(infection))) +
#   scale_fill_manual("Infection", values = c("cornflowerblue", "firebrick", "darkorange")) +
#   geom_boxplot() +
#   #   geom_dotplot(binaxis='y', dotsize=0.4)+
#   ylab("Simpsons's Index") +
#   xlab("DPI") +
#   # legend("right", title = "Temperature", c("18", "28", "32")) +
#   theme_bw()
# p3

# # inv simpson time, temperature:infection, temperature:time, infection:time

# p4 <- ggplot(dfDiv, aes(x = as.factor(time), y = invsimpson)) +
#   # scale_fill_manual(values = c("cornflowerblue", "firebrick")) +
#   geom_boxplot() +
#   #   geom_dotplot(binaxis='y', dotsize=0.4)+
#   ylab("Inverse Simpson's Index") +
#   xlab("DPI") +
#   theme_bw()
# p4

# p5 <- ggplot(dfDiv, aes(x = as.factor(temperature), y = invsimpson, fill = c(infection))) +
#   scale_fill_manual("Infection", values = c("cornflowerblue", "firebrick")) +
#   geom_boxplot() +
#   #   geom_dotplot(binaxis='y', dotsize=0.4)+
#   ylab("Inverse Simpson's Index") +
#   xlab("Temperature") +
#   # legend("right", title = "Infection", c("Control", "Chikungunya")) +
#   theme_bw()
# p5

# p6 <- ggplot(dfDiv, aes(x = as.factor(time), y = invsimpson, fill = c(temperature))) +
#   scale_fill_manual("Temperature", values = c("cornflowerblue", "darkorange", "firebrick")) +
#   geom_boxplot() +
#   #   geom_dotplot(binaxis='y', dotsize=0.4)+
#   ylab("Inverse Simpson's Index") +
#   xlab("DPI") +
#   # legend("right", title = "Temperature", c("18", "28", "32")) +
#   theme_bw()
# p6

# p7 <- ggplot(dfDiv, aes(x = as.factor(time), y = invsimpson, fill = c(infection))) +
#   scale_fill_manual("Infection", values = c("cornflowerblue", "firebrick")) +
#   geom_boxplot() +
#   #   geom_dotplot(binaxis='y', dotsize=0.4)+
#   ylab("Inverse Simpson's Index") +
#   xlab("DPI") +
#   # legend("right", title = "Temperature", c("18", "28", "32")) +
#   theme_bw()
# p7

# # save each plot
# ggsave("plots/Wimalasiri_simpsonEvenness_boxplot_time.pdf", p1, width = 10, height = 10, units = "cm")
# ggsave("plots/Wimalasiri_simpsonEvenness_boxplot_temp_infection.pdf", p2, width = 10, height = 10, units = "cm")
# ggsave("plots/Wimalasiri_simpsonEvenness_boxplot_time_infection.pdf", p3, width = 10, height = 10, units = "cm")
# ggsave("plots/Wimalasiri_invsimpsonEvenness_boxplot_time.pdf", p4, width = 10, height = 10, units = "cm")
# ggsave("plots/Wimalasiri_invsimpsonEvenness_boxplot_temp_infection.pdf", p5, width = 10, height = 10, units = "cm")
# ggsave("plots/Wimalasiri_invsimpsonEvenness_boxplot_time_temp.pdf", p6, width = 10, height = 10, units = "cm")
# ggsave("plots/Wimalasiri_invsimpsonEvenness_boxplot_time_infection.pdf", p7, width = 10, height = 10, units = "cm")

