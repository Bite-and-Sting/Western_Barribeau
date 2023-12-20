# setwd("~/4th Year Modules/LIFE700/Ferreira")

library(limma)
library(edgeR)

x<-read.csv('totalCounts_notEuk_taxLevel1_integer_zeros.csv', sep=",")
# still some NAs
x[, 6:149][is.na(x[, 6:149])] <- 0

namesD <- as.data.frame(names(x[6:149]))
names(namesD)<-'sample'

rownames(x) <- x[,1] 
scounts <- x[,c(6:149)]

# Count per millions
scountsPerMillion <- cpm(scounts)
summary(scountsPerMillion)


meta <- read.csv('Ferreira.metadata.csv', sep= ',', header=T)
head(meta)
meta <- meta[,c(1,11, 12, 14)]
names(meta)<-c('sample', 'temperature', 'infection', 'time')
table(meta$time)
table(meta$temperature)
table(meta$infection)

meta$time <- factor(meta$time, levels = c('24', '48'))
meta$temperature <- factor(meta$temperature, levels = c("20", "28", "36"))

meta$multicond = paste(meta$temperature, meta$infection, meta$time, sep = '_')
multicond=meta$multicond
# group=meta$multicond
sdge <- DGEList(scounts)
#skeep <- filterByExpr(sdge)
sD <- sdge #[skeep, keep.lib.sizes=FALSE ]
dim(sD)

design = model.matrix(~meta$temperature*meta$infection*meta$time) # Create design matrix for glm

dge = estimateGLMCommonDisp(sD, design)

sdge2 <- DGEList(scounts, group = meta$multicond)
table(meta$multicond)

fit = glmFit(dge,design)

# Testing for DEGs
lrt_temp <- glmLRT(fit, coef=2)
lrt_inf <- glmLRT(fit, coef=3)
lrt_time <- glmLRT(fit, coef=4)
lrt_temp_inf <- glmLRT(fit, coef=5)
lrt_temp_time <- glmLRT(fit, coef=6)
lrt_inf_time <- glmLRT(fit, coef=7)
lrt_temp_inf_time <- glmLRT(fit, coef=8)

topTags(lrt_temp)
topTags(lrt_inf)
topTags(lrt_time)
topTags(lrt_temp_inf)
topTags(lrt_temp_time)
topTags(lrt_inf_time)
topTags(lrt_temp_inf_time)

temp.table <- topTags( lrt_temp, n=nrow(sD$counts) )$table
inf.table <- topTags( lrt_inf, n=nrow(sD$counts) )$table
time.table <- topTags( lrt_time, n=nrow(sD$counts) )$table
temp_inf.table <- topTags( lrt_temp_inf, n=nrow(sD$counts) )$table
temp_time.table <- topTags( lrt_temp_time, n=nrow(sD$counts) )$table
inf_time.table <- topTags( lrt_inf_time, n=nrow(sD$counts) )$table
temp_inf_time.table <- topTags(lrt_temp_inf_time, n=nrow(sD$counts))$table

temp.table$Gene.ID<-row.names(temp.table)
inf.table$Gene.ID<-row.names(inf.table)
time.table$Gene.ID<-row.names(time.table)
temp_inf.table$Gene.ID<-row.names(temp_inf.table)
temp_time.table$Gene.ID<-row.names(temp_time.table)
inf_time.table$Gene.ID<-row.names(inf_time.table)
temp_inf_time.table$Gene.ID<-row.names(temp_inf_time.table)

namedf<-x[,1:5]

temp.table<-merge(temp.table, namedf, by='Gene.ID')
inf.table<-merge(inf.table, namedf, by='Gene.ID')
time.table<-merge(time.table, namedf, by='Gene.ID')
temp_inf.table<-merge(temp_inf.table, namedf, by='Gene.ID')
temp_time.table<-merge(temp_time.table, namedf, by='Gene.ID')
inf_time.table<-merge(inf_time.table, namedf, by='Gene.ID')
temp_inf_time.table<-merge(temp_inf_time.table, namedf, by='Gene.ID')

# make output directory
if(!dir.exists("tables/Ferriera_DE")) {
  dir.create("tables/Ferriera_DE", recursive = TRUE)
}
write.csv(temp.table, "tables/Ferriera_DE/DEtemp.csv", quote = F, row.names = F)
write.csv(inf.table, "tables/Ferriera_DE/DEinf.csv", quote = F, row.names = F)
write.csv(time.table, "tables/Ferriera_DE/DEtime.csv", quote = F, row.names = F)
write.csv(temp_inf.table, "tables/Ferriera_DE/DEtemp_inf.csv", quote = F, row.names = F)
write.csv(temp_time.table, "tables/Ferriera_DE/DEtemp_time.csv", quote = F, row.names = F)
write.csv(inf_time.table, "tables/Ferriera_DE/DEinf_time.csv", quote = F, row.names = F)
write.csv(temp_inf_time.table, "tables/Ferriera_DE/DEtemp_inf_time.csv", quote = F, row.names = F)


# unclear what this is doing.
# #To separate logFC by temperature
design = model.matrix(~tempTreat*meta$infection*meta$time) # Create design matrix for glm
dge = estimateGLMCommonDisp(sD, design)
sdge2 <- DGEList(scounts, group = meta$multicond)
table(meta$multicond)
fit = glmFit(dge,design)
lrt_temp <- glmLRT(fit, coef=1:3)
topTags(lrt_temp)
temp.table <- topTags( lrt_temp, n=nrow(sD$counts) )$table
temp.table$Gene.ID<-row.names(temp.table)
namedf<-x[,1:5]
temp.table<-merge(temp.table, namedf, by='Gene.ID')
write.csv(temp.table, "tables/Ferriera_DE/newDEtemp.csv", quote = F, row.names = F)

head(temp.table)
