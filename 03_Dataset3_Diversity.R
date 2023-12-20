df<-read.csv('Murrieta_nrnt_notEuk.csv', sep=",")
# still some NAs
df[, 6:186,][is.na(df[, 6:186])] <- 0
namesD <- as.data.frame(names(df[6:186]))
names(namesD)<-'sample'
rownames(df) <- df[,1] 
scounts <- df[,c(6:186)]

meta <- read.csv('Murrieta.metadata.csv', sep= ',', header=T)
meta <- meta[,c(1,12)]
names(meta)<-c('sample', 'temperature')
table(meta$temperature)
meta <- meta[meta$temperature != "25-35",]
meta$temperature <- factor(meta$temperature, levels = c("25", "28", "32", "35"))
table(meta$temperature)

keep = scounts[, names(scounts) %in% meta$sample] # here we drop all the columns that aren't in our new metadata

table(meta$temperature)

# get the sample names for each subset
c25 <- meta$sample[meta$temperature == "25"]
c28 <- meta$sample[meta$temperature == "28"]
c32 <- meta$sample[meta$temperature == "32"]
c35 <- meta$sample[meta$temperature == "35"]


#treatment subsets
c25 <- df[, match(c25, colnames(df))]
c28 <- df[, match(c28, colnames(df))]
c32 <- df[, match(c32, colnames(df))]
c35 <- df[, match(c35, colnames(df))]

dList <- list(c25, c28, c32, c35)
names(dList)  <- c("25C", "28C", "32C", "35C")
sampleNames <- c(
  names(c25), names(c28), names(c32), names(c35)
)

# shannon each subset
shanLists <- lapply(dList, diversity, index = "shannon")
names(shanLists) <- c(names(dList))

# shannon temperatures
shannondiv.T <- diversity(t(df[, 6:186]), index = "shannon")
simpsondiv.T <- diversity(t(df[, 6:186]), index = "simpson")
invsimpsondiv.T <- diversity(t(df[, 6:186]), index = "invsimpson")

# combine into a single dataframe
dfDiv <- as.data.frame(cbind(shannondiv.T, simpsondiv.T, invsimpsondiv.T))
names(dfDiv) <- c("shannon", "simpson", "invsimpson")
dfDiv$sample <- rownames(dfDiv)
head(dfDiv)

# add sample metadata to the dataframe
dfDiv <- merge(dfDiv, meta, by = "sample")
head(dfDiv)

shanGlm <- glm(shannon ~ temperature, data = dfDiv)
shanAnovaTab <- anova(shanGlm, test = "F")
rownames(shanAnovaTab)[1] <- "(Intercept)"

# check assumptions of the test
par(mfrow = c(2, 2))
plot(shanGlm)

dev.copy2pdf(file = "plots/QC/Murrieta_shannonDiversity_glm_diagnostics.pdf")

simpGlm <- glm(simpson ~ temperature, data = dfDiv)
simpAnovaTab <- anova(simpGlm, test = "F")
rownames(simpAnovaTab)[1] <- "(Intercept)"

# check assumptions of the test
plot(simpGlm)
dev.copy2pdf(file = "plots/QC/Murrieta_simpsonDiversity_glm_diagnostics.pdf")

invsimpGlm <- glm(invsimpson ~ temperature, data = dfDiv)
invsimpAnovaTab <- anova(invsimpGlm, test = "F")
rownames(invsimpAnovaTab)[1] <- "(Intercept)"

plot(invsimpGlm)
dev.copy2pdf(file = "plots/QC/Murrieta_invsimpsonDiversity_glm_diagnostics.pdf")

# save tables
shanTab <- knitr::kable(shanAnovaTab, format = "pipe", digits = 3, caption = "ANOVA table for Shannon's diversity index")
kableExtra::save_kable(shanTab, file = "tables/Murrieta_shannonDiversity_glm_anova_table.txt")

simpTab <- knitr::kable(simpAnovaTab, format = "pipe", digits = 3, caption = "ANOVA table for Simpson's diversity index")
kableExtra::save_kable(simpTab, file = "tables/Murrieta_simpsonDiversity_glm_anova_table.txt")

invSimpTab <- knitr::kable(invsimpAnovaTab, format = "pipe", digits = 3, caption = "ANOVA table for inverse Simpson's diversity index")
kableExtra::save_kable(invSimpTab, file = "tables/Murrieta_invsimpsonDiversity_glm_anova_table.txt")


# plot the diversity indices (none significant)
p1 <- ggplot(dfDiv, aes(x = temperature, y = shannon)) +
  geom_boxplot() +
  labs(x = "Temperature (C)", y = "Shannon's diversity index") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- ggplot(dfDiv, aes(x = temperature, y = simpson)) +
  geom_boxplot() +
  labs(x = "Temperature (C)", y = "Simpson's diversity index") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- ggplot(dfDiv, aes(x = temperature, y = invsimpson)) +
  geom_boxplot() +
  labs(x = "Temperature (C)", y = "Inverse Simpson's diversity index") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3

# save plots
all <- gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
all
dev.copy2pdf(file = "plots/Murrieta_diversity_indices.pdf")
