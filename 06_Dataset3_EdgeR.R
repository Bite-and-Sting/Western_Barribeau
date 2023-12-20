# setwd("~/4th Year Modules/LIFE700/Murrieta")
library(limma)
library(edgeR)

x<-read.csv('Murrieta_nrnt_notEuk.csv', sep=",")
# still some NAs
x[, 6:186][is.na(x[, 6:186])] <- 0

namesD <- as.data.frame(names(x[6:186]))
names(namesD)<-'sample'

rownames(x) <- x[,1] 
scounts <- x[,c(6:186)]

# Count per millions
scountsPerMillion <- cpm(scounts)
summary(scountsPerMillion)

meta <- read.csv('Murrieta.metadata.csv', sep= ',', header=T)

meta <- meta[,c(1,12)]
names(meta)<-c('sample', 'temperature')
table(meta$temperature)
tempTreat <- factor(meta$temperature)
#remove 25-35 treatment
meta <- meta[,c(1,2)]
names(meta)<-c('sample', 'temperature')
table(meta$temperature)
meta2=meta[meta$temperature != '25-35',]
names(meta2)<-c('sample', 'temperature')
table(meta2$temperature)
meta2$temperature=factor(meta2$temperature)
keep = scounts[, names(scounts) %in% meta2$sample] # here we drop all the columns that aren't in our new metadata


meta2$multicond = paste(meta2$temperature)
multicond=meta2$multicond

# group=meta$multicond
sdge <- DGEList(keep)
#skeep <- filterByExpr(sdge)
sD <- sdge#[skeep, keep.lib.sizes=FALSE ]
dim(sD)

design = model.matrix(~meta2$temperature) # Create design matrix for glm

dge = estimateGLMCommonDisp(sD, design)

sdge2 <- DGEList(keep, group = meta2$multicond)
table(meta2$multicond)

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

if(!dir.exists('tables/Murrieta_DE')) {
  dir.create("tables/Murrieta_DE", recursive = TRUE)
}

write.csv(temp.table, 'tables/Murrieta_DE/M.DEtemp.csv', quote=F, row.names = F)
write.csv(inf.table, 'tables/Murrieta_DE/M.DEinf.csv', quote=F, row.names = F)
write.csv(time.table, 'tables/Murrieta_DE/M.DEtime.csv', quote=F, row.names = F)
write.csv(temp_inf.table, 'tables/Murrieta_DE/M.DEtemp_inf.csv', quote=F, row.names = F)
write.csv(temp_time.table, 'tables/Murrieta_DE/M.DEtemp_time.csv', quote=F, row.names = F)
write.csv(inf_time.table, 'tables/Murrieta_DE/M.DEinf_time.csv', quote=F, row.names = F)
write.csv(temp_inf_time.table, 'tables/Murrieta_DE/M.DEtemp_inf_time.csv', quote=F, row.names = F)

#To separate logFC by temperature
design = model.matrix(~meta2$temperature) # Create design matrix for glm
dge = estimateGLMCommonDisp(sD, design)
sdge2 <- DGEList(scounts, group = meta$multicond)
table(meta$multicond)
fit = glmFit(dge,design)
lrt_temp <- glmLRT(fit, coef=1:4)
topTags(lrt_temp)
temp.table <- topTags( lrt_temp, n=nrow(sD$counts) )$table
temp.table$Gene.ID<-row.names(temp.table)
namedf<-x[,1:5]
temp.table<-merge(temp.table, namedf, by='Gene.ID')
write.csv(temp.table, "tables/Murrieta_DE/MnewDEtemp.csv", quote = F, row.names = F)
