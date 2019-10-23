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

# remove all patientes from fulldata that are not subjected to microarray analysis
fulldata <- fulldata %>% filter(fulldata$patient.Identification.MAARS.identifier..MAARS_identifier. %in% annotations$MAARS_identifier)

# separate alldata based on clinical group
alldata_ad <- alldata %>% filter(clinical_group == "AD")
alldata_ctrl <- alldata %>% filter(clinical_group == "CTRL")
alldata_pso <- alldata %>% filter(clinical_group == "PSO")


## build annotation dataset
sample_id <- colnames(omicsdata)
annotations <- data.frame(sample_id = colnames(omicsdata), MAARS_identifier = gsub('.{3}$', '', sample_id))
annotations <- dplyr::left_join(annotations, alldata, by = "sample_id") %>% select(sample_id, MAARS_identifier.x, clinical_group, lesional, CUSTOM_Age, Gender, Institution) %>% 
  rename(MAARS_identifier = MAARS_identifier.x)

annotations <- dplyr::left_join(annotations, fulldata, by = c("sample_id" = "involved.skin.biopsy.involved.skin.biopsy.MAARS.Sample.identifier..MAARS_Sample_identifier.")) %>% 
  select(sample_id, MAARS_identifier, clinical_group, lesional, CUSTOM_Age, Gender, Institution, patient.SCORAD.index.SCORAD.SCORAD.Score..SCORAD_Score.) %>% 
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


### build model matrix ###
# select only AL (AD lesional) samples
annotations_AL <- filter(annotations, sample_group == "AD lesional")
omicsdata_AL <- omicsdata[, colnames(omicsdata) %in% annotations_AL$sample_id]

DGE_AL <- DGEList(counts = omicsdata_AL, samples = annotations_AL, genes = rownames(omicsdata))

# model matrix
design_AL <- as.data.frame(model.matrix(~SCORAD_Score, data = annotations))


# fit the linear model
fit <- lmFit(omicsdata_AL, design_AL)
Bayesfit <- eBayes(fit)
volcanoplot(Bayesfit, coef = 2)
AL_signif <- decideTests(Bayesfit, adjust.method = "BH", p.value = 0.05, lfc = log2(1.02))
# lfc is very small because it is only the change in 1 unit of SCORAD score, have to find optimal lfc value!
summary(AL_signif)


# get significantly upregulated genes => IT WORKS!!!!!!
up_AL <- which(AL_signif[, 2] == 1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score 
upGenes_AL <- DGE_AL$genes$genes[up_AL]
upGenes_AL

# upTable <- topTable(upGenes, adjust.method= "BH", sort.by="p", n = 100) 
# questionable to do this, I don't know what the output really means and why I could put topGenes as model input


# get significantly downregulated genes
down_AL <- which(AL_signif[, 2] == -1) # 1 is upregulated, 0 not significant, -1 is downregulated
# column 1 is intercept, col 2 is SCORAD_Score 
downGenes_AL <- DGE_AL$genes$genes[down_AL]
downGenes_AL

# see if the most differentially expressed gene is linearly associated with SCORAD  
gene <- data.frame(t(omicsdata[upGenes_AL[1], ]))
gene$sample_id <- rownames(gene)
gene <- left_join(gene, annotations_AL, by = "sample_id") %>% na.omit(gene)
plot(gene$SCORAD_Score, gene$ENSG00000153993_at)

lm <- lm(gene$ENSG00000003400_at ~ SCORAD_Score, data = gene)
summary(lm)

## to do: get topTable to work, so I can look at lfc and p-value in one table
