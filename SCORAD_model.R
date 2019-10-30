##### build a model based on the DEGs and other clinical variables #####
### in order to do this, the DEG_limma.R script must be run first! ###

# install required packages
# install.packages("ComplexHeatmap")

# load packages
library("mixOmics")
library(gplots)
library(ComplexHeatmap)

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


# write them into txt file
# write(DEG_AL_unadj_ID, "DEG_AL_unadj_ID.txt")
# write(DEG_AL_unadj_name, "DEG_AL_unadj_name.txt")


### heatmap for differentially expressed genes

# subset expression dataframe to only contain DEG values
omicsdata_AL_DEG <- as.matrix(subset(omicsdata_AL, subset = rownames(omicsdata_AL) %in% DEG_AL_unadj))

# topTable from DEGs only
DEGtopTable <- rbind(upTable, downTable) %>% arrange(desc(logFC))

heatmap.2(omicsdata_AL_DEG, trace = "none", density.info = "none")


row_labels <- structure(ensembl2gene_DEG$GeneSymbol, names = ensembl2gene_DEG$EnsemblID) # to do: subset df so that only DEGs are in df!
column_labels <- structure(gsub("MAARS_", "", colnames(omicsdata_AL_DEG)), names = colnames(omicsdata_AL_DEG))

# annotations_AL$SCORAD_Score[column_order(hm)]
# SCORAD score is ordered automatically, so I don't have to do this
  
row_ha <- HeatmapAnnotation("SCORAD" = anno_barplot(annotations_AL$SCORAD_Score, fill = "grey"), 
                            annotation_name_gp = gpar(fontsize = 10, fontface = "bold"))
Heatmap(omicsdata_AL_DEG, border = TRUE, column_title = "Samples", row_title = "Differentially expressed genes", 
        show_row_names = FALSE, column_labels = column_labels, row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8), bottom_annotation = row_ha, 
        heatmap_legend_param = list(title = "normalized expression level", title_position = "leftcenter-rot", legend_height = unit(5, "cm")))
  

# perform PLS to reduce the number of dimensions

splsda_DEG <- splsda(X = t(omicsdata_AL_DEG), Y = annotations_AL$SCORAD_Score, ncomp = 2, 
                        keepX = rep(5, 2), mode = "regression") 

plotIndiv(splsda_DEG, ind.names = annotations_AL$SCORAD_Score)
plotVar(splsda_DEG, overlap = TRUE)

selectVar(splsda_DEG)

# build the dataframe for the model
df_DEG <- data.frame(t(omicsdata_AL_DEG)) %>% mutate(sample_id = rownames(df_DEG))
df_DEG <- select(df_DEG$sample_id, everything())
df_DEG <- dplyr::left_join(df_DEG, annotations_AL, by = c(rownames(df_DEG) = "sample_id"))


# build a model
m1 <- lm()





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