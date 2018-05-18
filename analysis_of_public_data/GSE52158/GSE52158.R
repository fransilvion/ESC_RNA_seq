#Affymetrix Human Genome U133 Plus 2.0 Array
library(limma)
library(ggplot2)
library(dplyr)
library(affy)
library(hgu133plus2.db)
library(tibble)

if (file.exists("GSE52158.Rdata")) {
  # if previously downloaded
  load("GSE52158.Rdata")
} else {
  # Get geo object that contains our data and phenotype information
  geo_obj <- getGEO("GSE52158", GSEMatrix = TRUE)
  geo_obj <- geo_obj[[1]]
  save(geo_obj, file = "GSE52158.Rdata")
}

#show(geo_obj)
#APS - anteriomost primitive streak
#AFG - anterior foregut
#PFG - posterior foregut
#MHG - midgut/hindgut
geo_metadata <- pData(geo_obj)[, c("organism_ch1", "title", colnames(pData(geo_obj))[grep("characteristics", 
                                                                                          colnames(pData(geo_obj)))])]

geo_metadata <- geo_metadata[,-3]
colnames(geo_metadata) <- c("organism", "sample", "description")
#only hESC, SR1 APS and SR1 DE 
geo_metadata <- geo_metadata[1:9,]
geo_metadata$cell_type <- as.factor(c(rep("ESC", 3), rep("APS", 3), rep("DE", 3)))
geo_metadata$cell_type <- relevel(geo_metadata$cell_type, "ESC")

#data in RMA
data <- exprs(geo_obj)
data <- data[, rownames(geo_metadata)]
#filtering
keep.exprs <-rowSums(data > 20) >= 3 
data <- data[keep.exprs,]

#data is RMA-normalized signal intensity (not log transformed)
#thus we have to log normalized it 
log_data <- log(data+1)

#for plotting
hist(log_data)

#for plotting
meltedExpressionMatrix <- log_data %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  melt(id = "gene") 

#for plotting
meltedExpressionMatrix %>% 
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#all(colnames(data) == rownames(geo_metadata)) (TRUE)
designMatrix <- model.matrix(~cell_type, geo_metadata)

#linear fit 
cellTypeFit <- lmFit(log_data, designMatrix)
# run ebayes to calculate moderated t-statistics
#trend True - for accounting of low expressed probes
cellTypeFitEb <- eBayes(cellTypeFit, trend = T)

cellTypeRes <- decideTests(cellTypeFitEb, p.value = 0.05, lfc = 1)
summary(cellTypeRes)

topProbesDE <- topTable(cellTypeFitEb, coef = "cell_typeDE", p.value = 0.05, lfc = 1, number = Inf) 
downTopProbesDE <- topProbesDE %>%
  rownames_to_column("probes") %>%
  filter(logFC < 0) #616

upTopProbesDE <- topProbesDE %>% 
  rownames_to_column("probes") %>%
  filter(logFC > 0) #775

#get gene symbols
up_genes <- select(hgu133plus2.db, upTopProbesDE$probes, c("SYMBOL"))
up_genes <- up_genes$SYMBOL
up_genes <- up_genes[!is.na(up_genes)]

down_genes <- select(hgu133plus2.db, downTopProbesDE$probes, c("SYMBOL"))
down_genes <- down_genes$SYMBOL
down_genes <- down_genes[!is.na(down_genes)]

#check with golden list
up_golden_list <- scan("~/Papers/CMAP/up_genes_in_DE.grp", character(), quote = "")
up_golden_list[up_golden_list %in% up_genes]
up_golden_list[!up_golden_list %in% up_genes]

down_golden_list <- scan("~/Papers/CMAP/down_genes_in_DE.grp", character(), quote = "")
down_golden_list[down_golden_list %in% down_genes]
down_golden_list[!down_golden_list %in% down_genes]