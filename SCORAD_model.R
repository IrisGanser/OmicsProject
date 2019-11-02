##### build a model based on the DEGs and other clinical variables #####
### in order to do this, the DEG_limma.R script must be run first! ###

# install required packages


# load packages
library(mixOmics)
library(ComplexHeatmap)
library(RColorBrewer)
library(GGally)
library(caret)

# not the outcome is the SCORAD and DEGs are explanatory factors

#significantly upregulated genes in unadjusted continuous analysis
upGenes_AL
# significantly downregulated genes in unadjusted continuous analyses
downGenes_AL
# combine regulated genes in one vector
DEG_AL_unadj <- c(upGenes_AL, downGenes_AL)
DEG_AL_unadj
# remove the "_at" from the end of the IDs
DEG_AL_unadj_ID <- gsub("_at", "", DEG_AL_unadj)
# get gene names
ensembl2gene <- read.csv("~/Documents/Omics/Project/OmicsProject/ensembl2gene.txt", stringsAsFactors=FALSE)
ensembl2gene_DEG <- subset(ensembl2gene, EnsemblID %in% DEG_AL_unadj_ID)

DEG_AL_unadj_name <- ensembl2gene_DEG$GeneSymbol
list_DEG_AL <- data.frame(DEG_AL_unadj_ID, DEG_AL_unadj_name)

# write them into txt file
# write(DEG_AL_unadj_ID, "DEG_AL_unadj_ID.txt")
# write(DEG_AL_unadj_name, "DEG_AL_unadj_name.txt")
# write.table(list_DEG_AL, "list_DEG_AL.txt", quote = FALSE, row.names = FALSE)

### heatmap for differentially expressed genes

# subset expression dataframe to only contain DEG values
omicsdata_AL_DEG <- as.matrix(subset(omicsdata_AL, subset = rownames(omicsdata_AL) %in% DEG_AL_unadj))

# topTable from DEGs only
DEGtopTable <- rbind(upTable, downTable) %>% arrange(desc(logFC))

# heatmap.2(omicsdata_AL_DEG, trace = "none", density.info = "none")


row_labels <- structure(ensembl2gene_DEG$GeneSymbol, names = ensembl2gene_DEG$EnsemblID) # to do: subset df so that only DEGs are in df!
column_labels <- structure(gsub("MAARS_", "", colnames(omicsdata_AL_DEG)), names = colnames(omicsdata_AL_DEG))

# annotations_AL$SCORAD_Score[column_order(hm)]
# SCORAD score is ordered automatically, so I don't have to do this
  
row_ha <- HeatmapAnnotation("SCORAD" = anno_barplot(annotations_AL$SCORAD_Score, fill = "grey"), 
                            annotation_name_gp = gpar(fontsize = 10, fontface = "bold"))
Heatmap(omicsdata_AL_DEG, border = TRUE, column_title = "Samples", row_title = "Differentially expressed genes", 
        show_row_names = FALSE, column_labels = column_labels, row_names_gp = gpar(fontsize = 5),
        clustering_method_rows = "complete", clustering_method_columns = "complete",
        column_names_gp = gpar(fontsize = 8), bottom_annotation = row_ha, 
        heatmap_legend_param = list(title = "normalized expression level", title_position = "leftcenter-rot", legend_height = unit(5, "cm")))
  
### perform PCA descriptive analysis
# df_PCA_AL <- data.frame(t(omicsdata_AL)[, -1], subset(annotations_AL, select = c(clinical_group, lesional, sample_group, Institution, SCORAD_Score, SCORAD_severity)))
PCA_AL <- pca(t(omicsdata_AL)[, -1])
plotIndiv(PCA_AL, group = annotations_AL$SCORAD_severity, legend = TRUE, pch = annotations_AL$Institution, ind.names = FALSE, 
          title = "PCA of all genes in AL samples, comp 1 & 2", legend.title.pch = "Institution", legend.title = "SCORAD Severity")



### perform sPLS to reduce the number of dimensions
# sPLS-DA for SCORAD as categorical variable

splsda_DEG <- splsda(X = t(omicsdata_AL_DEG), Y = annotations_AL$SCORAD_severity, ncomp = 2, scale = TRUE, 
                        keepX = rep(10, 2), mode = "regression") 

plotIndiv(splsda_DEG, ind.names = annotations_AL$SCORAD_Score, legend = TRUE, ellipse = TRUE)
plotVar(splsda_DEG, overlap = TRUE, var.names = FALSE)
selectVar(splsda_DEG)
auroc(splsda_DEG) # I don't know what this really means...  


# categorical variable with SCORAD as 10-point categories
annotations_AL$SCORAD_cat10 <- cut(annotations_AL$SCORAD_Score, breaks=c(-Inf, 25, 37.5, 50, 60, 70, 80, Inf), 
                                   labels=c("<25", "25-37.5", "37.6-50", "50.1-60", "60.1-70", "70.1-80", ">80"), right = TRUE)
# col_vector <- c("#4575B4","#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027")
col_vector <- c("#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026")


# sPLS for SCORAD as continuous variable
dim(t(omicsdata_AL_DEG))
SCORAD_df <- data.frame(annotations_AL$SCORAD_Score)
dim(SCORAD_df)

spls_DEG <- spls(X = t(omicsdata_AL_DEG), Y = SCORAD_df, ncomp = 10, keepX = rep(20, 10), scale = TRUE, mode = "regression")

# spls tuning
perf.pls <- perf(spls_DEG, validation = "loo", folds = 5, progressBar = FALSE)
plot(perf.pls$Q2.total, type = "b", xlab = "Number of components", ylab = "Q2 total value", main = "Evolution the cross-validation error")
perf.pls$Q2.total
perf.pls$R2
# according to this, chose 2 components

list.keepX <- c(2:10, 15, 20)
# tuning based on MAE

tune.spls.MAE <- tune.spls(X = t(omicsdata_AL_DEG), Y = SCORAD_df, ncomp = 4, test.keepX = list.keepX, validation = "loo", 
                           folds = 5, progressBar = FALSE, measure = 'MAE')
plot(tune.spls.MAE, legend.position = 'topright')
tune.spls.MAE$choice.keepX
# keep 9 variables per component

## optimized sPLS
spls_DEG_tuned <- spls(X = t(omicsdata_AL_DEG), Y = SCORAD_df, ncomp = 2, keepX = rep(9, 2), scale = TRUE, mode = "regression")

plotIndiv(spls_DEG_tuned, rep.space = "X-variate", ind.names = FALSE, pch = 16, group = annotations_AL$SCORAD_cat10, legend = TRUE,
          col.per.group = col_vector, legend.title = "SCORAD Score", title = "sPLS of DEG with 2 comp, 9 variables per comp")
plotVar(spls_DEG_tuned, cex = c(4, 4))

var_comp1 <- selectVar(spls_DEG_tuned, comp = 1)
var_comp2 <- selectVar(spls_DEG_tuned, comp = 2)

DEG_spls <- c(var_comp1$X$name, var_comp2$X$name)


## repeat heatmap with selected 18 genes
DEG_spls_ID <- gsub("_at", "", DEG_spls)
ensembl2gene_spls_DEG <- subset(ensembl2gene, EnsemblID %in% DEG_spls_ID)
omicsdata_AL_spls_DEG <- as.matrix(subset(omicsdata_AL, subset = rownames(omicsdata_AL) %in% DEG_spls))

row_labels <- structure(ensembl2gene_spls_DEG$GeneSymbol, names = ensembl2gene_spls_DEG$EnsemblID)
column_labels <- structure(gsub("MAARS_", "", colnames(omicsdata_AL_spls_DEG)), names = colnames(omicsdata_AL_spls_DEG))


row_ha <- HeatmapAnnotation("SCORAD" = anno_barplot(annotations_AL$SCORAD_Score, fill = "grey"), 
                            annotation_name_gp = gpar(fontsize = 10, fontface = "bold"))
Heatmap(omicsdata_AL_spls_DEG, border = TRUE, column_title = "Samples", row_title = "Differentially expressed genes", 
        row_labels = row_labels, column_labels = column_labels, row_names_gp = gpar(fontsize = 8),
        clustering_method_rows = "complete", clustering_method_columns = "complete",
        column_names_gp = gpar(fontsize = 8), bottom_annotation = row_ha, 
        heatmap_legend_param = list(title = "normalized expression level", title_position = "leftcenter-rot", legend_height = unit(5, "cm")))

# spls_DEG_tuned$variates$X 
# are these the values of the components? If so, can they be used in the final model?

# extract the biggest loading values
spls_DEG_loadings <- spls_DEG_tuned$loadings$X %>% as.data.frame() %>% mutate(geneID = spls_DEG_tuned$names$colnames$X)  %>%
  filter(comp1 != 0 | comp2 != 0)

# export the 17 genes in txt file
# write.table(ensembl2gene_spls_DEG, "list_DEG_AL_spls.txt", quote = FALSE, row.names = FALSE)


##### final SCORAD prediction model #####
# build the dataframe for the model
df_DEG <- data.frame(t(omicsdata_AL_DEG))
df_DEG <- subset(df_DEG, select = colnames(df_DEG) %in% DEG_spls)
df_DEG$sample_id = rownames(df_DEG)
df_DEG <- dplyr::select(df_DEG, sample_id, everything())


sum(is.na(df_DEG)) # no missing values
names(df_DEG)


# some exploratory visualizations
ggpairs(df_DEG, columns = c(8, 5, 6, 7, 9, 10, 11))
ggpairs(df_DEG, columns = c(8, 16:30)) 
boxplot(annotations_AL$SCORAD_severity ~ annotations_AL$Institution)
# genes are all correlated individually with SCORAD (because we selected them this way!), but also correlated among each other

## do univariate analyses of SCORAD against all other variables
alldata_AL <- filter(alldata, lesional == "LES" & clinical_group == "AD")
alldata_AL <- dplyr::left_join(alldata_AL, fulldata, by = c("sample_id" = "involved.skin.biopsy.involved.skin.biopsy.MAARS.Sample.identifier..MAARS_Sample_identifier."))%>% 
  subset(select = c(1:33, 35, 51, 261)) %>% 
  rename(SCORAD_Score = patient.SCORAD.index.SCORAD.SCORAD.Score..SCORAD_Score., ethnicity = patient.Diagnostic...Phenotypic.Data.Ethnicity.Family.History.Ethnicity..Ethnicity.)
alldata_AL$SCORAD_severity <- cut(alldata_AL$SCORAD_Score, breaks=c(-Inf, 25, 50, Inf), 
                               labels=c("mild", "moderate", "severe"), right = FALSE)



df_DEG <- dplyr::left_join(alldata_AL, df_DEG,  by = "sample_id")

spls_comp <- data.frame(spls_DEG_tuned$variates$X) 
spls_comp$sample_id <- rownames(spls_comp)
df_DEG <- left_join(df_DEG, spls_comp,  by = "sample_id")


variables <- names(alldata_AL)[c(7:29, 31:35)]

fit_univariable <- list()
# build univariate models based on this dataframe
for (i in 1:length(variables)){
 fit_univariable[[i]] <- lm(as.formula(paste("SCORAD_Score", "~", variables[i], sep = " ")), data = alldata_AL)
}

# extract coefficients and p-values
lapply(fit_univariable, summary)
# variables with p<0.2: Institution, Gender, Known_Allergies_v2..House_dust_mite, Known_Allergies_v2..Food, Known_Allergies_v2..Drug_Allergy,
# Other_concurrent_chronic_diseases_v2..Asthma, CUSTOM_Fam._hist._Atopic_dermatitis, ethnicity

# check model assumptions
par(mfrow = c(2, 2))
lapply(fit_univariable, plot)



# build a model with only the genes
m1 <- lm(as.formula(paste(colnames(df_DEG)[36], "~",
                          paste(colnames(df_DEG)[38:54], collapse = "+"), sep = "")), data = df_DEG)
par(mfrow = c(2, 2))
plot(m1)
summary(m1)
# because of correlation between genes, only some of them are significant

m2 <- lm(SCORAD_Score ~ comp1 + comp2, data = df_DEG)
summary(m2)

# build a model with DEG and clinical variables
m3 <- lm(as.formula(paste(colnames(df_DEG)[36], "~",
                          paste(colnames(df_DEG)[c(9, 11, 12, 15, 28, 33, 35, 38:54)], collapse = "+"), sep = "")), data = df_DEG)
summary(m3)

# remove gender, asthma, food allergy because of p-values 
m4 <- lm(as.formula(paste(colnames(df_DEG)[36], "~",
                          paste(colnames(df_DEG)[c(11, 15, 33, 35, 38:54)], collapse = "+"), sep = "")), data = df_DEG)
summary(m4)
anova(m3, m4)

# remove ethnicity because of few participants in categories
summary(df_DEG$ethnicity)
m5 <- lm(as.formula(paste(colnames(df_DEG)[36], "~",
                          paste(colnames(df_DEG)[c(11, 15, 33, 38:54)], collapse = "+"), sep = "")), data = df_DEG)
summary(m5)
anova(m4, m5) # no significant contribution of ethnicity to the model

# remove dust mite allergy, family history, and drug allergy (all non-significant) => back to DEG only model

## now try the same without individual genes, but with the components
m6 <- lm(as.formula(paste(colnames(df_DEG)[36], "~", 
                          paste(colnames(df_DEG)[c(9, 11, 12, 15, 28, 33, 35, 55, 56)], collapse = "+"), sep = "")), data = df_DEG)
summary(m6)
# same thing as the above model...


## model calibration plots
par(mfrow = c(1, 1))
plot(m1$fitted.values, df_DEG$SCORAD_Score)
abline(0, 1) # 45° line
plot(m2$fitted.values, df_DEG$SCORAD_Score)
abline(0, 1)


##### model evaluations
### repeated cross-validation
set.seed(123)
train.control_cv <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
# Evaluate the model with all genes
m1_cv <- train(as.formula(paste(colnames(df_DEG)[36], "~", paste(colnames(df_DEG)[38:54], collapse = "+"), sep = "")), 
               data = df_DEG, method = "lm", trControl = train.control_cv)
print(m1_cv)
AIC(m1)

m1_cv_params <- round(m1_cv$results[2:4], 3)

# Evaluate the model with components only
m2_cv <- train(SCORAD_Score ~ comp1 + comp2, 
               data = df_DEG, method = "lm", trControl = train.control_cv)
print(m2_cv)
AIC(m2)

m2_cv_params <- round(m2_cv$results[2:4], 3)
# component-only model has higher R-squared and smaller RMSE and MAE than all-gene model! And also a smaller AIC

### bootstrapping
set.seed(123)
train.control_boot <- trainControl(method = "boot", number = 500)
# Evaluate the model with all genes
m1_boot <- train(as.formula(paste(colnames(df_DEG)[36], "~", paste(colnames(df_DEG)[38:54], collapse = "+"), sep = "")), 
               data = df_DEG, method = "lm", trControl = train.control_boot)
print(m1_boot)
m1_boot_params <- round(m1_boot$results[2:4], 3)

# Evaluate the model with components only
m2_boot <- train(SCORAD_Score ~ comp1 + comp2, 
               data = df_DEG, method = "lm", trControl = train.control_boot)
print(m2_boot)
m2_boot_params <- round(m2_boot$results[2:4], 3)
# component-only model still has higher R-squared and smaller RMSE and MAE than all-gene model!

# combine all performance metrics into one data frame
m1_m2_perf_params <- matrix(c(m1_cv_params, m2_cv_params, m1_boot_params, m2_boot_params), nrow = 4, byrow = F) %>% data.frame()
rownames(m1_m2_perf_params) <-  c("m1 CV", "m2 CV", "m1 boot", "m2 boot")
colnames(m1_m2_perf_params) <- c("RMSE", "Rsquared", "MAE")
m1_m2_perf_params
