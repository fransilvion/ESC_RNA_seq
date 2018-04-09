GSE75748 data analysis
================
German Novakovskiy
4 March 2018

Here we analyze gene expression data from [this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1033-x). Some analysis steps are made according to [this workflow article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/).

Let's read our data into R (GSE75748):

``` r
if (file.exists("GSE75748.Rdata")) {
    # if previously downloaded
    load("GSE75748.Rdata")
} else {
    # Get geo object that contains our data and phenotype information
    geo_obj <- getGEO("GSE75748", GSEMatrix = TRUE)
    #geo_obj contains 1810 entries ((1018 + 758) single cells + (19 + 15) bulk RNA-seq)
    geo_obj <- geo_obj[[1]]
    save(geo_obj, file = "GSE75748.Rdata")
}
```

Further analysis of the downloaded data:

``` r
#show data structure
show(geo_obj)
```

    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 0 features, 1810 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: GSM1964932 GSM1964933 ... GSM1966747 (1810 total)
    ##   varLabels: title geo_accession ... passage:ch1 (45 total)
    ##   varMetadata: labelDescription
    ## featureData: none
    ## experimentData: use 'experimentData(object)'
    ## Annotation: GPL16791

Let's load metadata:

``` r
#no confounding factors or batch effects (the same platform, unknown sex)
geo_metadata <- pData(geo_obj)

#due to complexity of metadata table it's more convenient to create your own (see below)
```

assayData has 0 features, it's not possible to load expression data from this dataset. That's why we will use expression data downloaded directly from the website

``` r
#cell type
bulk_cell_type_ec <- read.csv("GSE75748_bulk_cell_type_ec.csv/GSE75748_bulk_cell_type_ec.csv", header=TRUE)
#time course
bulk_time_course_ec <- read.csv("GSE75748_bulk_time_course_ec.csv/GSE75748_bulk_time_course_ec.csv", header=TRUE)

#single cell/cell type
sc_cell_type <- read.csv("GSE75748_sc_time_course_ec.csv", header=TRUE)
rownames(sc_cell_type) <- sc_cell_type$X
sc_cell_type <- sc_cell_type[,-1]
```

We are interested in 0 time point. We have only single cell data for it. Let's extract it:

``` r
zero_cols <- colnames(sc_cell_type)
zero_cols <- zero_cols[grepl("00h", zero_cols)]

zero_expr <- sc_cell_type[,zero_cols]
dim(zero_expr)
```

    ## [1] 19189    92

There are 92 cells from 0 time point, just as authors described it.

Metadata for samples of 0h:

``` r
zero_metadata <- geo_metadata %>%
  filter(grepl("for 0 hours", title))
dim(zero_metadata)
```

    ## [1] 92 45

Let's make 3 replicas out of 92 cells (31 cells in first, 31 in second, 30 in the last - clarify this moment with Sara):

``` r
zero_first_replica <- apply(zero_expr[1:31], 1, mean) #mean or median?
zero_second_replica <- apply(zero_expr[32:62], 1, mean) #mean or median?
zero_third_replica <- apply(zero_expr[63:92], 1, mean) #mean or median?

zero_expr_replicas <- data.frame(H9_0h_rep1 = zero_first_replica, H9_0h_rep2 = zero_second_replica,
                                 H9_0h_rep3 = zero_third_replica)

rownames(bulk_time_course_ec) <- bulk_time_course_ec$X
bulk_time_course_ec <- bulk_time_course_ec[,-1]

#delete extra rows from zero expr data frame
delete_rows <- setdiff(rownames(zero_expr_replicas), rownames(bulk_time_course_ec))
zero_expr_replicas <- zero_expr_replicas[!rownames(zero_expr_replicas) %in% delete_rows,]

bulk_time_course_ec <- merge(zero_expr_replicas, bulk_time_course_ec, by="row.names")
rownames(bulk_time_course_ec) <- bulk_time_course_ec$Row.names
bulk_time_course_ec <- bulk_time_course_ec[,-1]
```

Let's create metadata data frame:

``` r
sample_names <- colnames(bulk_time_course_ec) #without gene name

time_fact <- paste(c(0, 12,24,36,72,96), sep="", 'h') #factors of time
time_fact <- rep(time_fact, each=3)

time_int <- rep(c(0, 12,24,36,72,96), each=3) #time as continious variable

metadata <- data.frame(samples = sample_names, time = time_fact, age = time_int)

metadata %>% kable()
```

| samples       | time |  age|
|:--------------|:-----|----:|
| H9\_0h\_rep1  | 0h   |    0|
| H9\_0h\_rep2  | 0h   |    0|
| H9\_0h\_rep3  | 0h   |    0|
| H9\_12h\_rep1 | 12h  |   12|
| H9\_12h\_rep2 | 12h  |   12|
| H9\_12h\_rep3 | 12h  |   12|
| H9\_24h\_rep1 | 24h  |   24|
| H9\_24h\_rep2 | 24h  |   24|
| H9\_24h\_rep3 | 24h  |   24|
| H9\_36h\_rep1 | 36h  |   36|
| H9\_36h\_rep2 | 36h  |   36|
| H9\_36h\_rep3 | 36h  |   36|
| H9\_72h\_rep1 | 72h  |   72|
| H9\_72h\_rep2 | 72h  |   72|
| H9\_72h\_rep3 | 72h  |   72|
| H9\_96h\_rep1 | 96h  |   96|
| H9\_96h\_rep2 | 96h  |   96|
| H9\_96h\_rep3 | 96h  |   96|

``` r
str(metadata)
```

    ## 'data.frame':    18 obs. of  3 variables:
    ##  $ samples: Factor w/ 18 levels "H9_0h_rep1","H9_0h_rep2",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ time   : Factor w/ 6 levels "0h","12h","24h",..: 1 1 1 2 2 2 3 3 3 4 ...
    ##  $ age    : num  0 0 0 12 12 12 24 24 24 36 ...

Our primary goal here is to analyze bulk\_time\_course data (each stage with the next stage, step-by-step).

Deleting low-expressed genes (CPM &gt; 3). We don't want to filter out a gene that is highly expressed in one group, but lowly expressed in another, because that gene may be biologically interesting! So, we identify the smallest number of samples that could be part of one group (3), and make sure we don't remove genes that are expressed in more than that number:

``` r
#first let's create a edgeR DGElist object
bulk_time_course_ec <- as.matrix(bulk_time_course_ec)
rows <- rownames(bulk_time_course_ec)
#bulk_time_course_ec <- bulk_time_course_ec[,-1]
bulk_time_course_ec <- apply(bulk_time_course_ec, 2, as.double)
rownames(bulk_time_course_ec) <- rows

#for testing, delete after
#bulk_time_course_ec <- bulk_time_course_ec[,-c(1:3)]

DGE_bulk_time_course_ec <- DGEList(counts = bulk_time_course_ec) 

cpm <- cpm(DGE_bulk_time_course_ec)
keep.exprs <-rowSums(cpm > 5) >= 3

DGE_bulk_time_course_ec <- DGE_bulk_time_course_ec[keep.exprs,,]

dim(DGE_bulk_time_course_ec)
```

    ## [1] 11942    18

Very useful insight about [the filtering step](https://support.bioconductor.org/p/77178/).

Our data here is in gene expected counts that were calculated using RSEM 1.2.3 (according to GEO query). Thus, we need to perform a library size normalization using edgeR:

``` r
normalized_factors_expression <- calcNormFactors(DGE_bulk_time_course_ec, method = "TMM") #calculation of scaling factors (for library size)

normalized_factors_expression$samples$norm.factors
```

    ##  [1] 0.8428422 0.8703352 0.8118465 1.0859853 1.1004586 1.0886113 1.0490378
    ##  [8] 1.0540197 1.0514873 1.0052039 1.0037277 1.0121749 0.9378989 1.0534128
    ## [15] 1.0056914 1.0215541 1.0297469 1.0400207

For this dataset the effect of TMM-normalisation is mild, as evident in the magnitude of the scaling factors, which are all relatively close to 1 (except samples that are from single cell data, 0h time point).

Let's look at distribution of values:

``` r
#removing gene column and transforming into matrix (for hist)
data <- as.matrix(DGE_bulk_time_course_ec$counts)

hist(data, main="GSE75748", xlim = c(0,1500), xlab = "Expression",
     ylab = "Frequency", breaks = 300)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-13-1.png)

Let's now perform RNA-seq analysis with limma, using only time factor variable (time column in metadata) and let's look separately at DE gene at each stage (from 12 to 24, from 24 to 36 and etc.):

``` r
metadata_time <- metadata[,-3]

metadata_time$samples <- as.character(metadata_time$samples)
metadata_time %>% kable()
```

| samples       | time |
|:--------------|:-----|
| H9\_0h\_rep1  | 0h   |
| H9\_0h\_rep2  | 0h   |
| H9\_0h\_rep3  | 0h   |
| H9\_12h\_rep1 | 12h  |
| H9\_12h\_rep2 | 12h  |
| H9\_12h\_rep3 | 12h  |
| H9\_24h\_rep1 | 24h  |
| H9\_24h\_rep2 | 24h  |
| H9\_24h\_rep3 | 24h  |
| H9\_36h\_rep1 | 36h  |
| H9\_36h\_rep2 | 36h  |
| H9\_36h\_rep3 | 36h  |
| H9\_72h\_rep1 | 72h  |
| H9\_72h\_rep2 | 72h  |
| H9\_72h\_rep3 | 72h  |
| H9\_96h\_rep1 | 96h  |
| H9\_96h\_rep2 | 96h  |
| H9\_96h\_rep3 | 96h  |

``` r
designMatrix <- model.matrix(~0 + time, metadata_time)
head(designMatrix, 10) %>% kable()
```

|  time0h|  time12h|  time24h|  time36h|  time72h|  time96h|
|-------:|--------:|--------:|--------:|--------:|--------:|
|       1|        0|        0|        0|        0|        0|
|       1|        0|        0|        0|        0|        0|
|       1|        0|        0|        0|        0|        0|
|       0|        1|        0|        0|        0|        0|
|       0|        1|        0|        0|        0|        0|
|       0|        1|        0|        0|        0|        0|
|       0|        0|        1|        0|        0|        0|
|       0|        0|        1|        0|        0|        0|
|       0|        0|        1|        0|        0|        0|
|       0|        0|        0|        1|        0|        0|

We have estimated counts data here. So we can apply voom (which usually takes count data as an input and transforms them to logCPM) that estimates the mean-variance relationship and uses this to compute appropriate observation-level weights. The data are then ready for linear modelling.

``` r
after_voom_cpm <- voom(normalized_factors_expression, designMatrix, plot=TRUE)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
hist(after_voom_cpm$E, main="cleaned GSE75748 - log2 transformed CPM", xlab = "Expression",
     ylab = "Frequency")
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-17-1.png) Let's build boxplots and explore the data:

``` r
cleaned_log_cpm_df <- as.data.frame(after_voom_cpm$E)

cleaned_log_cpm_df <- cleaned_log_cpm_df %>% rownames_to_column("gene")

meltedLogedBultTimeCourseEc <- melt(cleaned_log_cpm_df, id='gene')

meltedLogedBultTimeCourseEc %>%
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-18-1.png)

Looks much better now!

``` r
plotMDS(cleaned_log_cpm_df[,-1], cex=1.5)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-19-1.png)

It's clear that time ponts 0h, 12h, 24h and 36h are well (mostly) separated from each other and 72h with 96h, however the last group (72 and 96) is clustered closer together, thus we don't expect to see much DE genes in these groups.

Let's create a contrast matrix because we are interested in DE genes across different time points (12h compared to 0h, 24h compared to 12h and so on):

``` r
rownames(cleaned_log_cpm_df) <- cleaned_log_cpm_df$gene
cleaned_log_cpm_df <- cleaned_log_cpm_df[,-1]

# construct the contrast matrix
contrastMatrix <- makeContrasts(
  v12v0 = time12h - time0h,
  v24v12 = time24h - time12h,
  v36v24 = time36h - time24h,
  v72v36 = time72h - time36h,
  v96v72 = time96h - time72h,
  levels = designMatrix
)

contrastMatrix %>% kable()
```

|         |  v12v0|  v24v12|  v36v24|  v72v36|  v96v72|
|---------|------:|-------:|-------:|-------:|-------:|
| time0h  |     -1|       0|       0|       0|       0|
| time12h |      1|      -1|       0|       0|       0|
| time24h |      0|       1|      -1|       0|       0|
| time36h |      0|       0|       1|      -1|       0|
| time72h |      0|       0|       0|       1|      -1|
| time96h |      0|       0|       0|       0|       1|

``` r
# keep the fit around as we will need to it for looking at other contrasts later 
time_course_Fit <- lmFit(after_voom_cpm, designMatrix)

# fit the contrast using the original fitted model
contrastFit <- contrasts.fit(time_course_Fit, contrastMatrix)

# apply eBayes() for moderated statistics
contrastFitEb <- eBayes(contrastFit)

contrastGenes <- topTable(contrastFitEb, number = Inf, p.value = 0.05)

plotSA(contrastFitEb)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-21-1.png)

``` r
#contrastGenes %>% kable()
```

``` r
cutoff <- 5e-02 #0.05 p value
#adjust method by default is BH (equivalent to fdr)
time_course_res <- decideTests(contrastFitEb, p.value = cutoff, lfc = 1)
summary(time_course_res)
```

    ##        v12v0 v24v12 v36v24 v72v36 v96v72
    ## Down    2625    148    271   1732    129
    ## NotSig  6515  11563  11293   8542  11657
    ## Up      2802    231    378   1668    156

We see there are different number of genes up and down regulated at each stage.

Here are the genes that upregulated from 24 hours to 36 hours.

``` r
hits3 <- time_course_res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(v36v24 > 0)

head(hits3) %>% kable()
```

| gene   |  v12v0|  v24v12|  v36v24|  v72v36|  v96v72|
|:-------|------:|-------:|-------:|-------:|-------:|
| AADAT  |      0|       0|       1|       0|       0|
| ABHD6  |      0|       0|       1|      -1|       0|
| ABLIM1 |      0|       0|       1|      -1|       0|
| ABTB2  |      0|       0|       1|       0|       0|
| ACSS3  |      0|       1|       1|       0|       0|
| ADAM12 |      0|       0|       1|       0|       1|

``` r
#function for plotting genes
plotGenes <- function(genes, expressionMatrix, metadata) {
  
  expressionDataForGenes <- expressionMatrix %>%
    rownames_to_column("gene") %>%
    filter(gene %in% genes) %>%
    melt()
  
  colnames(expressionDataForGenes) <- c("gene", "samples", "expression")
  expressionDataForGenes <- expressionDataForGenes %>%
    left_join(metadata, id="time")
  
  expressionDataForGenes %>% 
    ggplot(aes(x = time, y = expression, color=gene)) +
    geom_point() +
    geom_jitter() +
    stat_summary(aes(y = expression, group=1), fun.y = "mean", geom="line", size=2) +
    facet_wrap(~gene)
}
```

Let's plot 4 random genes that are upregulated at each stage.

``` r
hits1 <- time_course_res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(v12v0 > 0)

#stage 24v12
hits2 <- time_course_res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(v24v12 > 0)

#stage 72v36
hits4 <- time_course_res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(v72v36 > 0)

#stage 96v72
hits5 <- time_course_res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(v96v72 > 0)
```

``` r
sample_genes_12v0 <- sample(hits1$gene,4)
plotGenes(sample_genes_12v0, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-26-1.png)

``` r
sample_genes_24v12 <- sample(hits2$gene,4)
plotGenes(sample_genes_24v12, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
sample_genes_36v24 <- sample(hits3$gene,4)
plotGenes(sample_genes_36v24, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-28-1.png)

``` r
sample_genes_72v36 <- sample(hits4$gene,4)
plotGenes(sample_genes_72v36, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-29-1.png)

``` r
sample_genes_96v72 <- sample(hits5$gene,4)
plotGenes(sample_genes_96v72, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-30-1.png)

let's build those genes that according to the paper are DE in 24v12 transition:

``` r
sample_genes <- c("T", "CDX1", "MSX2")
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-31-1.png)

COMMENTS
--------

-   paper reports that CDX1, MSX2 and T are over expressed from 12 to 24h transition. I see only T and CDX1 (p-value cutoff 5e-02, FDR)
-   paper reports that CER1 and GATA4 are over expressed from 24 to 36h transition.I see both (p-value cutoff 5e-02, FDR)
-   paper reports that DKK4 and MYCT1 are over expressed from 36 to 72. I see both (p-value cutoff 5e-02, FDR)

Let's look at the expression of these genes EOMES, CER1, GATA4, PRDM1, and POU2AF1 at 96h because they are expected to be highly expressed during 96h stage, according to the paper (and KLF8 - "KLF8 may play a specific role during the transition from mesendoderm toward DE cells"):

``` r
sample_genes <- c("EOMES", "CER1", "GATA4", "PRDM1", "POU2AF1", "KLF8")
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-32-1.png) Let's look at the expression of pluripotency genes POU5F1, NANOG, and SOX2:

``` r
sample_genes <- c("POU5F1", "NANOG", "SOX2")
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-33-1.png) Their expression is going down, just as expected.

Now let's look at key DE markers CXCR4, SOX17, HNF1B, KIT, and KRT19:

``` r
#sample_genes <- c("FOXA2")
#HNF1B is filtered as low expressed gene (cpm < 5)
sample_genes <- c("CXCR4", "SOX17", "HNF1B", "KIT", "KRT19")
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-34-1.png) It's clear that CXCR4, SOX17 and KIT are upregulated at 96h, expression of KRT19 and HNF1B did not change so much, but KRT19 remained highly expressed, while HNF1B remained to be down-regulated.

Interesting genes: - FOXA2, which is regulated by long non-coding RNA DEANR1 (is low filtered here); - GSC is controlled by DIGIT lncRNA; - EOMES, MIXL1, SOX17 are DE markers (MIXL1 is mesodendoderm marker)

``` r
sample_genes <- c("FOXA2", "GSC", "EOMES", "MIXL1", "SOX17")
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-35-1.png)

Comparisons with the paper
==========================

"A total of 3247 differentially expressed genes were identified in the scRNA-seq time course experiment listed in Additional file 4: Table S3."

In my work I found this number of DE genes:

``` r
#list all DE genes
allDEGenes <- topTable(contrastFitEb, number = Inf, p.value = 0.05, lfc = 1)
nrow(allDEGenes)
```

    ## [1] 7886

According to Decide tests we have this number of DE genes (at different time stages):

``` r
de_genes_rows <- apply(time_course_res, 1, function(x) any(x != 0))
de_genes_at_stages <- time_course_res[de_genes_rows,]
nrow(de_genes_at_stages)
```

    ## [1] 7480

We will check both numbers.

``` r
all_de_genes <- rownames(allDEGenes)
stages_de_genes <- rownames(de_genes_at_stages)
```

Get list of DE genes from the study:

``` r
paper_de_genes <- read_excel("13059_2016_1033_MOESM4_ESM.xlsx")
head(paper_de_genes) %>% kable()
```

| GeneID | most likely pattern |  PP pattern|
|:-------|:--------------------|-----------:|
| A2ML1  | Down-NC-NC-NC-NC    |   0.2976393|
| AAK1   | Down-NC-Up-NC-NC    |   0.4062348|
| AARS   | Down-NC-NC-NC-NC    |   0.4685744|
| AARS2  | Down-Up-NC-NC-NC    |   0.2714841|
| AASS   | Down-NC-Up-NC-NC    |   0.5190453|
| AATF   | NC-NC-Up-NC-NC      |   0.3152373|

``` r
paper_de_genes <- paper_de_genes$GeneID
length(paper_de_genes)
```

    ## [1] 3247

So there is in fact 3247 genes in the result of this study.

Comparison of DE genes from decideTest with DE genes from the paper:

Comparison of DE genes from topTable with DE genes from the paper:

``` r
temp <- venn.diagram(list(My_TopTable = all_de_genes, Paper_genes = paper_de_genes),fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4, lty =2, fontfamily =3, filename = NULL, main = "Comparison of paper DE genes with my all DE genes (from TopTable)", category.names = c("My topTable", "Paper genes"))

grid::grid.newpage()
grid::grid.draw(temp)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-41-1.png)

``` r
temp2 <- venn.diagram(list(My_TopTable = stages_de_genes, Paper_genes = paper_de_genes),fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4, lty =2, fontfamily =3, filename = NULL, main = "Comparison of paper DE genes with my stage DE genes (from decideTests)", category.names = c("Stage DE genes", "Paper genes"))

grid::grid.newpage()
grid::grid.draw(temp2)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-42-1.png)

Performing analysis without contrast matrix, using 0 as a reference
===================================================================

``` r
#0 hours as a reference
designMatrixReference <- model.matrix(~time, metadata_time)
head(designMatrixReference, 10) %>% kable()
```

|  (Intercept)|  time12h|  time24h|  time36h|  time72h|  time96h|
|------------:|--------:|--------:|--------:|--------:|--------:|
|            1|        0|        0|        0|        0|        0|
|            1|        0|        0|        0|        0|        0|
|            1|        0|        0|        0|        0|        0|
|            1|        1|        0|        0|        0|        0|
|            1|        1|        0|        0|        0|        0|
|            1|        1|        0|        0|        0|        0|
|            1|        0|        1|        0|        0|        0|
|            1|        0|        1|        0|        0|        0|
|            1|        0|        1|        0|        0|        0|
|            1|        0|        0|        1|        0|        0|

``` r
v <- voom (normalized_factors_expression, designMatrixReference, plot = FALSE)

# keep the fit around as we will need to it for looking at other contrasts later 
time_course_Fit_Reference <- lmFit(v, designMatrixReference)

# apply eBayes() for moderated statistics
time_course_Fit_Reference_Ebayes <- eBayes(time_course_Fit_Reference)

genesReference <- topTable(time_course_Fit_Reference_Ebayes, number = Inf, p.value = 0.05, lfc = 1)
```

    ## Removing intercept from test coefficients

``` r
dim(genesReference)
```

    ## [1] 8510    9

``` r
genesReference75748 <- genesReference
save(genesReference75748, file="GSE75748_topGenes.Rdata")
```

Comparisons of different papers
===============================

``` r
load("../GSE109658/GSE109658_topGenes.Rdata")

genes109658 <- rownames(genesReference109658)
genes75748 <- rownames(genesReference75748)

commonGenes <- intersect(genes109658, genes75748)
length(commonGenes)
```

    ## [1] 3321

``` r
commonGenes <- genesReference75748[commonGenes,]
commonGenes <- commonGenes[with(commonGenes, order(adj.P.Val)),]


sample_genes <- rownames(commonGenes)[1:6]
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-47-1.png)

``` r
all(rownames(commonGenes) %in% rownames(time_course_res))
```

    ## [1] TRUE

``` r
commonGenesTrends <- time_course_res[rownames(commonGenes),]
head(commonGenesTrends) %>% kable()
```

|         |  v12v0|  v24v12|  v36v24|  v72v36|  v96v72|
|---------|------:|-------:|-------:|-------:|-------:|
| ANP32E  |     -1|       0|       0|      -1|       0|
| TERF1   |     -1|      -1|       0|       0|       0|
| SRSF6   |     -1|       0|       0|      -1|       0|
| AK4     |     -1|      -1|      -1|       0|       0|
| PPP3CA  |     -1|       0|       0|      -1|       0|
| RHOBTB3 |      0|       1|       1|       1|       0|

``` r
upRegulatedToDE <- commonGenesTrends %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(`v72v36` == 1)

#79 are upregulated to DE
sample_genes <- upRegulatedToDE$gene[1:6]
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-50-1.png)

``` r
dim(upRegulatedToDE)
```

    ## [1] 816   6

Using age as continious
=======================

``` r
metadata_cont_time <- metadata[,-2]
metadata_cont_time %>% kable()
```

| samples       |  age|
|:--------------|----:|
| H9\_0h\_rep1  |    0|
| H9\_0h\_rep2  |    0|
| H9\_0h\_rep3  |    0|
| H9\_12h\_rep1 |   12|
| H9\_12h\_rep2 |   12|
| H9\_12h\_rep3 |   12|
| H9\_24h\_rep1 |   24|
| H9\_24h\_rep2 |   24|
| H9\_24h\_rep3 |   24|
| H9\_36h\_rep1 |   36|
| H9\_36h\_rep2 |   36|
| H9\_36h\_rep3 |   36|
| H9\_72h\_rep1 |   72|
| H9\_72h\_rep2 |   72|
| H9\_72h\_rep3 |   72|
| H9\_96h\_rep1 |   96|
| H9\_96h\_rep2 |   96|
| H9\_96h\_rep3 |   96|

``` r
designMatrix <- model.matrix(~age, metadata_cont_time)
head(designMatrix, 10) %>% kable()
```

|  (Intercept)|  age|
|------------:|----:|
|            1|    0|
|            1|    0|
|            1|    0|
|            1|   12|
|            1|   12|
|            1|   12|
|            1|   24|
|            1|   24|
|            1|   24|
|            1|   36|

``` r
expressionFit_age <- lmFit(cleaned_log_cpm_df, designMatrix)
expressionFitBayes_age <- eBayes(expressionFit_age)

topGenesAge <- topTable(expressionFitBayes_age, number = Inf, p.value = 0.05)#, lfc = 1)
```

    ## Removing intercept from test coefficients

``` r
nrow(topGenesAge)
```

    ## [1] 5326

``` r
head(topGenesAge) %>% kable()
```

|         |       logFC|   AveExpr|          t|  P.Value|  adj.P.Val|         B|
|---------|-----------:|---------:|----------:|--------:|----------:|---------:|
| CYP26A1 |   0.0995670|  6.154112|   20.25512|        0|          0|  22.66588|
| AP1M2   |  -0.0917140|  2.899708|  -19.61839|        0|          0|  22.08642|
| CRABP1  |  -0.0739464|  2.894697|  -19.44200|        0|          0|  21.92265|
| CDH3    |  -0.0599511|  4.131228|  -18.15966|        0|          0|  20.68725|
| COL3A1  |   0.1273994|  2.570885|   17.34050|        0|          0|  19.85395|
| ESRP1   |  -0.0940553|  4.839029|  -16.65435|        0|          0|  19.12717|

``` r
sample_genes <- rownames(topGenesAge)[1:6]
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-56-1.png)

``` r
temp <- venn.diagram(list(My_TopTable = rownames(topGenesAge), Paper_genes = paper_de_genes),fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4, lty =2, fontfamily =3, filename = NULL, main = "Comparison of paper DE genes with my all DE genes (from TopTable)", category.names = c("My topTable", "Paper genes"))

grid::grid.newpage()
grid::grid.draw(temp)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-57-1.png) \# Finding discrepancies

``` r
paper_de_genes <- read_excel("13059_2016_1033_MOESM4_ESM.xlsx")
foo <- data.frame(do.call('rbind', strsplit(as.character(paper_de_genes$`most likely pattern`),'-',fixed=TRUE)))
df <- data.frame(Gene = paper_de_genes$GeneID)
dfx <- cbind(df, foo)
paper_de_genes <- dfx
colnames(paper_de_genes) <- c("Gene", "12v0", "24v12", "36v24", "72v36", "96v72")
head(paper_de_genes) %>% kable()
```

| Gene  | 12v0 | 24v12 | 36v24 | 72v36 | 96v72 |
|:------|:-----|:------|:------|:------|:------|
| A2ML1 | Down | NC    | NC    | NC    | NC    |
| AAK1  | Down | NC    | Up    | NC    | NC    |
| AARS  | Down | NC    | NC    | NC    | NC    |
| AARS2 | Down | Up    | NC    | NC    | NC    |
| AASS  | Down | NC    | Up    | NC    | NC    |
| AATF  | NC   | NC    | Up    | NC    | NC    |

``` r
#deleting column 12v0, because I didn't  have that comparison
#paper_de_genes <- paper_de_genes[,-2]

paper_de_genes <- paper_de_genes %>%
  filter('12v0' != "NC" | `24v12` != "NC" | `36v24` != "NC" | `72v36` != "NC" | `96v72` != "NC")
dim(paper_de_genes)
```

    ## [1] 3247    6

``` r
temp2 <- venn.diagram(list(My_TopTable = stages_de_genes, Paper_genes = paper_de_genes$Gene),fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4, lty =2, fontfamily =3, filename = NULL, main = "Comparison of paper DE genes with my stage DE genes (from decideTests)", category.names = c("Stage DE genes", "Paper genes"))

grid::grid.newpage()
grid::grid.draw(temp2)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-60-1.png)

``` r
#reminder
summary(time_course_res)
```

    ##        v12v0 v24v12 v36v24 v72v36 v96v72
    ## Down    2625    148    271   1732    129
    ## NotSig  6515  11563  11293   8542  11657
    ## Up      2802    231    378   1668    156

``` r
paper_res <- list(`12v0` = table(paper_de_genes$`12v0`), `24v12` = table(paper_de_genes$`24v12`), `36v24` = table(paper_de_genes$`36v24`),
                  `72v36` = table(paper_de_genes$`72v36`), `96v72` = table(paper_de_genes$`96v72`))

paper_res
```

    ## $`12v0`
    ## 
    ## Down   NC   Up 
    ## 2051 1023  173 
    ## 
    ## $`24v12`
    ## 
    ## Down   NC   Up 
    ##   53 2417  777 
    ## 
    ## $`36v24`
    ## 
    ## Down   NC   Up 
    ##  328 1930  989 
    ## 
    ## $`72v36`
    ## 
    ## Down   NC   Up 
    ##  235 2781  231 
    ## 
    ## $`96v72`
    ## 
    ## Down   NC   Up 
    ##   10 3236    1

``` r
#hits at 36v24
head(hits3) %>% kable()
```

| gene   |  v12v0|  v24v12|  v36v24|  v72v36|  v96v72|
|:-------|------:|-------:|-------:|-------:|-------:|
| AADAT  |      0|       0|       1|       0|       0|
| ABHD6  |      0|       0|       1|      -1|       0|
| ABLIM1 |      0|       0|       1|      -1|       0|
| ABTB2  |      0|       0|       1|       0|       0|
| ACSS3  |      0|       1|       1|       0|       0|
| ADAM12 |      0|       0|       1|       0|       1|

``` r
paper_up_36v24 <- paper_de_genes %>%
  filter(`36v24` == "Up")
nrow(paper_up_36v24)
```

    ## [1] 989

``` r
#in paper but not in my research
s <- setdiff(paper_up_36v24$Gene, hits2$gene)
length(s)
```

    ## [1] 946

``` r
#plot genes that are DE in paper at 36 to 24 but not in my results
set.seed(123)
sample_genes <- sample(s, 4)
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-66-1.png)
