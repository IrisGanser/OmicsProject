##### build a model based on the DEGs and other clinical variables #####
### in order to do this, the DEG_limma.R script must be run first! ###

# install required packages


# load packages
library(mixOmics)
library(ComplexHeatmap)
library(RColorBrewer)

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
plotIndiv(PCA_AL, group = annotations_AL$SCORAD_severity, legend = TRUE, pch = 16, ind.names = FALSE, 
          title = "PCA of all genes in AL samples, comp 1 & 2")


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
df_DEG <- dplyr::left_join(annotations_AL, df_DEG,  by = "sample_id")

sum(is.na(df_DEG)) # no missing values
names(df_DEG)

# build a model
m1 <- lm(as.formula(paste(colnames(df_DEG)[8], "~",
                          paste(colnames(df_DEG)[c(5, 6, 9, 10, 16:30)], collapse = "+"), sep = "")),
         data = df_DEG)

summary(m1)





##### cross-validation

#Randomly shuffle the data
annotations_AL_shuffled <- annotations_AL[sample(nrow(annotations_AL)),]

#Create 10 equally size folds
folds <- cut(seq(1, nrow(annotations_AL_shuffled)), breaks=10, labels=FALSE)

# initialze the list with up- and downregulated genes
up_crossval <- vector("list", 10)
down_crossval <- vector("list", 10)

#Perform 10 fold cross validation
for(i in 1:10){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i, arr.ind=TRUE)
  valData <- annotations_AL_shuffled[testIndexes, ]
  trainData <- annotations_AL_shuffled[-testIndexes, ]
  
  # design matrix
  design_AL_crossval <- as.data.frame(model.matrix(~SCORAD_Score, data = trainData))
  omicsdata_crossval <- omicsdata_AL[colnames(omicsdata_AL) %in% trainData$sample_id]
  
  # fit the linear model
  fit <- lmFit(omicsdata_crossval, design_AL_crossval)
  Bayesfit <- eBayes(fit)
  AL_signif_crossval <- decideTests(Bayesfit, adjust.method = "BH", p.value = 0.05, lfc = log2(1.01))
  # lfc is very small because it is only the change in 1 unit of SCORAD score, have to find optimal lfc value!
  
  
  # get significantly upregulated genes
  up_AL_crossval <- which(AL_signif_crossval[, 2] == 1) # 1 is upregulated, 0 not significant, -1 is downregulated
  # column 1 is intercept, col 2 is SCORAD_Score 
  up_crossval[[i]] <- DGE_AL$genes$genes[up_AL_crossval]
  
  # get significantly downregulated genes
  down_AL_crossval <- which(AL_signif_crossval[, 2] == -1) # 1 is upregulated, 0 not significant, -1 is downregulated
  # column 1 is intercept, col 2 is SCORAD_Score 
  down_crossval[[i]] <- DGE_AL$genes$genes[down_AL_crossval]
  
}

## somehow now all the genes are non-significant, and I don't know if it's due to the reduced sample size or because sth else went wrong