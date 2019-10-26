##### differentially expressed genes with different SCORAD scores #####

library(ggplot2)
library(dplyr)
library(edgeR)
library(limma)
library(UpSetR)
library(flashClust)
library(GEOquery)
library(dendextend)
library(RColorBrewer)
library(knitr)
library(kableExtra)
library(FactoMineR)


## load data
alldata <- read.csv(file = "M2PHDS_19-20_OMICS_CLIN_DATA_MAARS_all_Fri_Apr_04_14h_CEST_2014.csv", header = TRUE, sep = "\t")
fulldata <- read.csv(file = "M2PHDS_19-20_OMICS_CLIN_DATA_MAARS_AD_full_20190131_12-34-49.csv", header = TRUE, sep = "\t")
omicsdata <- read.delim("~/Documents/Omics/Project/OmicsProject/M2PHDS_19-20_OMICS_TRANSC_MAARS_normTranscriptome_618samples_16042014.txt")

# remove all samples from alldata that are not used in microarray
alldata <- alldata %>% filter(sample_id %in% colnames(omicsdata))
alldata$allergy <- ifelse(rowSums(alldata[, 10:16]) > 0, TRUE, FALSE)

# separate alldata based on clinical group
alldata_ad <- alldata %>% filter(clinical_group == "AD")
alldata_ctrl <- alldata %>% filter(clinical_group == "CTRL")
alldata_pso <- alldata %>% filter(clinical_group == "PSO")


## build annotation dataset
sample_id <- colnames(omicsdata)
annotations <- data.frame(sample_id = colnames(omicsdata), MAARS_identifier = gsub('.{3}$', '', sample_id))
annotations <- dplyr::left_join(annotations, alldata, by = "sample_id") %>% select(sample_id, MAARS_identifier.x, clinical_group, lesional, CUSTOM_Age, Gender, Institution, allergy) %>% 
  rename(MAARS_identifier = MAARS_identifier.x)

annotations <- dplyr::left_join(annotations, fulldata, by = c("sample_id" = "involved.skin.biopsy.involved.skin.biopsy.MAARS.Sample.identifier..MAARS_Sample_identifier.")) %>% 
  select(sample_id, MAARS_identifier, clinical_group, lesional, CUSTOM_Age, Gender, Institution, patient.SCORAD.index.SCORAD.SCORAD.Score..SCORAD_Score., allergy) %>% 
  rename(SCORAD_Score = patient.SCORAD.index.SCORAD.SCORAD.Score..SCORAD_Score.)

annotations$SCORAD_severity[annotations$SCORAD_Score < 25] <- "mild"
annotations$SCORAD_severity[annotations$SCORAD_Score < 50 & annotations$SCORAD_Score >= 25] <- "moderate"
annotations$SCORAD_severity[annotations$SCORAD_Score >= 50] <- "severe"
annotations$SCORAD_severity[is.na(annotations$SCORAD_Score)] <- "no eczema"
annotations$SCORAD_severity <- as.factor(annotations$SCORAD_severity)


annotations$sample_group[annotations$clinical_group == "AD" & annotations$lesional == "LES"] <- "AD lesional"
annotations$sample_group[annotations$clinical_group == "AD" & annotations$lesional == "NON_LES"] <- "AD non-lesional"
annotations$sample_group[annotations$clinical_group == "PSO" & annotations$lesional == "LES"] <- "PSO lesional"
annotations$sample_group[annotations$clinical_group == "PSO" & annotations$lesional == "NON_LES"] <- "PSO non-lesional"
annotations$sample_group[annotations$clinical_group == "CTRL"] <- "CTRL"
annotations$sample_group<- as.factor(annotations$sample_group)

# relevel the dataset so that R picks the correct reference classes later on
annotations$clinical_group <- relevel(annotations$clinical_group, ref = "CTRL")
annotations$lesional <- relevel(annotations$lesional, ref = "NON_LES")

# some descriptions of SCORAD etc.
summary(annotations$SCORAD_Score)
summary(annotations$SCORAD_severity)
table(annotations$SCORAD_severity, annotations$sample_group)
summary(annotations$sample_group)
summary(annotations$clinical_group)
# there are 83 AD lesional skin samples and 81 AD non-lesional skin samples 

# remove all patientes from fulldata that are not subjected to microarray analysis
fulldata <- fulldata %>% filter(fulldata$patient.Identification.MAARS.identifier..MAARS_identifier. %in% annotations$MAARS_identifier)


### build model matrix ###
# select only AL (AD lesional) samples
annotations_AL <- filter(annotations, sample_group == "AD lesional")
annotations_AL$SCORAD_severity <- factor(annotations_AL$SCORAD_severity)
annotations_AL$SCORAD_severity <- ordered(annotations_AL$SCORAD_severity, levels=c("mild", "moderate", "severe"))

omicsdata_AL <- omicsdata[, colnames(omicsdata) %in% annotations_AL$sample_id]

DGE_AL <- DGEList(counts = omicsdata_AL, samples = annotations_AL, genes = rownames(omicsdata))

# model matrix
design_AL <- as.data.frame(model.matrix(~SCORAD_Score, data = annotations_AL))


# fit the linear model
fit <- lmFit(omicsdata_AL, design_AL)
Bayesfit <- eBayes(fit)
volcanoplot(Bayesfit, coef = 2)
AL_signif <- decideTests(Bayesfit, adjust.method = "BH", p.value = 0.05, lfc = log2(1.01))
# lfc is very small because it is only the change in 1 unit of SCORAD score, have to find optimal lfc value!
summary(AL_signif)


# get significantly upregulated genes => IT WORKS!!!!!!
up_AL <- which(AL_signif[, 2] == 1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score 
upGenes_AL <- DGE_AL$genes$genes[up_AL]
upGenes_AL

upTable <- topTable(Bayesfit, adjust.method= "BH", sort.by="p", n = Inf)
upTable$Ensembl <- rownames(upTable)
upTable <- filter(upTable, Ensembl %in% upGenes_AL)


# get significantly downregulated genes
down_AL <- which(AL_signif[, 2] == -1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score 
downGenes_AL <- DGE_AL$genes$genes[down_AL]
downGenes_AL

downTable <- topTable(Bayesfit, adjust.method= "BH", sort.by="p", n = Inf)
downTable$Ensembl <- rownames(downTable)
downTable <- filter(downTable, Ensembl %in% downGenes_AL)


# see if the most differentially expressed gene is linearly associated with SCORAD  
gene <- data.frame(t(omicsdata[upGenes_AL[1], ]))
gene$sample_id <- rownames(gene)
gene <- left_join(gene, annotations_AL, by = "sample_id") %>% na.omit(gene)
plot(gene$SCORAD_Score, gene$ENSG00000007933_at)

lm <- lm(gene$ENSG00000007933_at ~ SCORAD_Score, data = gene)
summary(lm)




### SCORAD as categorical variable ###

# model matrix
design_AL_cat <- as.data.frame(model.matrix(~SCORAD_severity, data = annotations_AL))
# this is the solution for an ordinal variable, results in a very weird design matrix, have to do more research on this if it is actually applicable
# first column should indicate linear trend, second quadratic trend

# fit the linear model
fit_cat <- lmFit(omicsdata_AL, design_AL_cat)
Bayesfit_cat <- eBayes(fit_cat)
volcanoplot(Bayesfit_cat, coef = 2)
AL_signif_cat <- decideTests(Bayesfit_cat, adjust.method = "BH", p.value = 0.05, lfc = log2(1.2))
summary(AL_signif_cat)


# get significantly upregulated genes 
up_AL_cat <- which(AL_signif_cat[, 2] == 1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score moderate, col 3 is SCORAD_Score severe
upGenes_AL_cat <- DGE_AL$genes$genes[up_AL_cat]
upGenes_AL_cat

upTable_cat <- topTable(Bayesfit_cat, adjust.method= "BH", n = Inf)
upTable_cat$Ensembl <- rownames(upTable_cat)
upTable_cat <- filter(upTable_cat, Ensembl %in% upGenes_AL_cat)


# get significantly downregulated genes
down_AL_cat <- which(AL_signif_cat[, 2] == -1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score 
downGenes_AL_cat <- DGE_AL$genes$genes[down_AL_cat]
downGenes_AL_cat

downTable_cat <- topTable(Bayesfit_cat, adjust.method= "BH", n = Inf)
downTable_cat$Ensembl <- rownames(downTable_cat)
downTable_cat <- filter(downTable_cat, Ensembl %in% downGenes_AL_cat)



### see which genes are regulated in both analyses
up_both <- intersect(upGenes_AL, upGenes_AL_cat)
up_both

down_both <- intersect(downGenes_AL, downGenes_AL_cat)
down_both
# not very satisfying results 


##### adjusted models for age, sex, and allergies #####

### continuous adjusted analysis

# model matrix
design_AL_adj <- as.data.frame(model.matrix(~SCORAD_Score + CUSTOM_Age + Gender + allergy, data = annotations_AL))


# fit the linear model
fit_adj <- lmFit(omicsdata_AL, design_AL_adj)
Bayesfit_adj <- eBayes(fit_adj)
volcanoplot(Bayesfit_adj, coef = 2)
AL_signif_adj <- decideTests(Bayesfit_adj, adjust.method = "BH", p.value = 0.05, lfc = log2(1.01))
# lfc is very small because it is only the change in 1 unit of SCORAD score, have to find optimal lfc value!
summary(AL_signif_adj)


# get significantly upregulated genes
up_AL_adj <- which(AL_signif_adj[, 2] == 1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score 
upGenes_AL_adj <- DGE_AL$genes$genes[up_AL_adj]
upGenes_AL_adj

upTable_adj <- topTable(Bayesfit_adj, adjust.method= "BH", n = Inf)
upTable_adj$Ensembl <- rownames(upTable_adj)
upTable_adj <- filter(upTable_adj, Ensembl %in% upGenes_AL_adj)
upTable_adj


# get significantly downregulated genes
down_AL_adj <- which(AL_signif_adj[, 2] == -1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score 
downGenes_AL_adj <- DGE_AL$genes$genes[down_AL_adj]
downGenes_AL_adj

downTable_adj <- topTable(Bayesfit_adj, adjust.method= "BH", n = Inf)
downTable_adj$Ensembl <- rownames(downTable_adj)
downTable_adj <- filter(downTable_adj, Ensembl %in% downGenes_AL_adj)
downTable_adj
max(downTable_adj$SCORAD_Score)

## the genes that are most up- or downregulated aren't even in the significant genes (look at volcano plot, e.g. biggest downward change = -0.06)
## should find other method to retrieve those genes as they migth be important


### categorical adjusted analysis

# model matrix
design_AL_cat_adj <- as.data.frame(model.matrix(~SCORAD_severity + CUSTOM_Age + Gender + allergy, data = annotations_AL))
# this is the solution for an ordinal variable, results in a very weird design matrix, have to do more research on this if it is actually applicable
# first column should indicate linear trend, second quadratic trend

# fit the linear model
fit_cat_adj <- lmFit(omicsdata_AL, design_AL_cat_adj)
Bayesfit_cat_adj <- eBayes(fit_cat_adj)
volcanoplot(Bayesfit_cat_adj, coef = 2)
AL_signif_cat_adj <- decideTests(Bayesfit_cat_adj, adjust.method = "BH", p.value = 0.05, lfc = log2(1.2))
summary(AL_signif_cat_adj)


# get significantly upregulated genes 
up_AL_cat_adj <- which(AL_signif_cat_adj[, 2] == 1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score moderate, col 3 is SCORAD_Score severe
upGenes_AL_cat_adj <- DGE_AL$genes$genes[up_AL_cat_adj]
upGenes_AL_cat_adj

upTable_cat_adj <- topTable(Bayesfit_cat_adj, adjust.method= "BH", n = Inf)
upTable_cat_adj$Ensembl <- rownames(upTable_cat_adj)
upTable_cat_adj <- filter(upTable_cat_adj, Ensembl %in% upGenes_AL_cat_adj)


# get significantly downregulated genes
down_AL_cat_adj <- which(AL_signif_cat_adj[, 2] == -1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score 
downGenes_AL_cat_adj <- DGE_AL$genes$genes[down_AL_cat_adj]
downGenes_AL_cat_adj

downTable_cat_adj <- topTable(Bayesfit_cat_adj, adjust.method= "BH", n = Inf)
downTable_cat_adj$Ensembl <- rownames(downTable_cat_adj)
downTable_cat_adj <- filter(downTable_cat_adj, Ensembl %in% downGenes_AL_cat_adj)


##### Comparison of adjusted and unadjusted analyses #####
up_both_adj_unadj <- intersect(upGenes_AL, upGenes_AL_adj)
up_both_adj_unadj

down_both_adj_unadj <- intersect(downGenes_AL, downGenes_AL_adj)
down_both_adj_unadj

## adjustment does not have a huge effect (most genes are similar between both analyses)

### see which genes are regulated in adjusted continuous and categorical analyses
up_both_adj <- intersect(upGenes_AL_adj, upGenes_AL_cat_adj)
up_both_adj

down_both_adj <- intersect(downGenes_AL_adj, downGenes_AL_cat_adj)
down_both_adj