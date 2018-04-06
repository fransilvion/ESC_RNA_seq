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

Deleting low-expressed genes (CPM &gt; 1). We don't want to filter out a gene that is highly expressed in one group, but lowly expressed in another, because that gene may be biologically interesting! So, we identify the smallest number of samples that could be part of one group (3), and make sure we don't remove genes that are expressed in more than that number:

``` r
#first let's create a edgeR DGElist object
bulk_time_course_ec <- as.matrix(bulk_time_course_ec)
rows <- rownames(bulk_time_course_ec)
#bulk_time_course_ec <- bulk_time_course_ec[,-1]
bulk_time_course_ec <- apply(bulk_time_course_ec, 2, as.double)
rownames(bulk_time_course_ec) <- rows

DGE_bulk_time_course_ec <- DGEList(counts = bulk_time_course_ec) 

cpm <- cpm(DGE_bulk_time_course_ec)
keep.exprs <-rowSums(cpm > 1) >= 3

DGE_bulk_time_course_ec <- DGE_bulk_time_course_ec[keep.exprs,,]

dim(DGE_bulk_time_course_ec)
```

    ## [1] 14333    18

Our data here is in gene expected counts that were calculated using RSEM 1.2.3 (according to GEO query). Thus, we need to perform a library size normalization using edgeR:

``` r
normalized_factors_expression <- calcNormFactors(DGE_bulk_time_course_ec, method = "TMM") #calculation of scaling factors (for library size)

normalized_factors_expression$samples$norm.factors
```

    ##  [1] 0.8442666 0.8716866 0.8044489 1.0858852 1.1035056 1.0820767 1.0468473
    ##  [8] 1.0447433 1.0476810 1.0021610 1.0003328 1.0034298 0.9439126 1.0600265
    ## [15] 1.0124795 1.0272051 1.0362534 1.0480844

For this dataset the effect of TMM-normalisation is mild, as evident in the magnitude of the scaling factors, which are all relatively close to 1.

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
time_course_Fit <- lmFit(cleaned_log_cpm_df, designMatrix)

# fit the contrast using the original fitted model
contrastFit <- contrasts.fit(time_course_Fit, contrastMatrix)

# apply eBayes() for moderated statistics
contrastFitEb <- eBayes(contrastFit)

contrastGenes <- topTable(contrastFitEb)

contrastGenes %>% kable()
```

|          |       v12v0|      v24v12|      v36v24|      v72v36|      v96v72|   AveExpr|         F|  P.Value|  adj.P.Val|
|----------|-----------:|-----------:|-----------:|-----------:|-----------:|---------:|---------:|--------:|----------:|
| GRM4     |   2.2089230|   1.1734871|  -2.4766781|  -6.0614660|   0.0292857|  2.485032|  977.3394|        0|          0|
| EOMES    |   8.8302459|   0.6724103|   0.1751554|   0.7446285|  -0.2371836|  7.096943|  842.7476|        0|          0|
| RSPO3    |  -0.8314077|   5.7547933|   2.2045844|   0.4417324|   0.0925639|  3.268075|  765.7165|        0|          0|
| WLS      |   6.4954736|   2.2393209|   0.1029077|   0.6515993|  -0.0697174|  6.297953|  707.2975|        0|          0|
| HTRA1    |   2.3781906|  -0.7974383|  -0.2921808|   6.8048875|   0.0139365|  4.578772|  687.9387|        0|          0|
| SERPINE2 |   3.5305349|   2.4714634|   0.5893767|   0.0998826|  -0.1530383|  8.917734|  663.6921|        0|          0|
| TBX3     |   5.6716824|   2.6807795|   1.3694251|  -2.8127483|  -0.2857683|  5.072513|  623.4500|        0|          0|
| DKK1     |   3.7353694|   3.9606430|   0.9531446|  -0.1500933|  -0.3256291|  4.984929|  596.5328|        0|          0|
| ZEB2     |   4.2841122|   1.5376152|   1.9141175|  -0.1566615|  -0.3736645|  4.436634|  554.0512|        0|          0|
| MMP14    |   5.4107705|   1.2799127|   1.6240900|   0.9572361|  -0.8326844|  5.685395|  528.9915|        0|          0|

``` r
cutoff <- 5e-02 #0.05 p value
#adjust method by default is BH (equivalent to fdr)
time_course_res <- decideTests(contrastFitEb, p.value = cutoff, lfc = 1)
summary(time_course_res)
```

    ##        v12v0 v24v12 v36v24 v72v36 v96v72
    ## Down    2894     94    253   1840    116
    ## NotSig  8417  14019  13729  10733  14090
    ## Up      3022    220    351   1760    127

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
| ABHD6  |      0|       1|       1|      -1|       0|
| ABLIM1 |      0|       0|       1|      -1|       0|
| ABTB2  |      0|       0|       1|       0|       0|
| ACSS3  |      0|       1|       1|       0|       0|
| ACVRL1 |     -1|       0|       1|       1|       0|
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
sample_genes <- c("CXCR4", "SOX17", "HNF1B", "KIT", "KRT19")
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-34-1.png) It's clear that CXCR4, SOX17 and KIT are upregulated at 96h, expression of KRT19 and HNF1B did not change so much, but KRT19 remained highly expressed, while HNF1B remained to be down-regulated.

Comparisons with the paper
==========================

"A total of 3247 differentially expressed genes were identified in the scRNA-seq time course experiment listed in Additional file 4: Table S3."

In my work I found this number of DE genes:

``` r
#list all DE genes
allDEGenes <- topTable(contrastFitEb, number = Inf, p.value = 0.05, lfc = 1)
nrow(allDEGenes)
```

    ## [1] 8817

According to Decide tests we have this number of DE genes (at different time stages):

``` r
de_genes_rows <- apply(time_course_res, 1, function(x) any(x != 0))
de_genes_at_stages <- time_course_res[de_genes_rows,]
nrow(de_genes_at_stages)
```

    ## [1] 8090

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

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-40-1.png)

``` r
temp2 <- venn.diagram(list(My_TopTable = stages_de_genes, Paper_genes = paper_de_genes),fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4, lty =2, fontfamily =3, filename = NULL, main = "Comparison of paper DE genes with my stage DE genes (from decideTests)", category.names = c("Stage DE genes", "Paper genes"))

grid::grid.newpage()
grid::grid.draw(temp2)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-41-1.png)

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

genesReference <- topTable(contrastFitEb, number = Inf, p.value = 0.05, lfc = 1)

dim(genesReference)
```

    ## [1] 8817    9

``` r
all_de_genes_ref <- rownames(genesReference)
all(all_de_genes == all_de_genes_ref)
```

    ## [1] TRUE

We are getting the same result as a covariat matrix top table for all genes.

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

    ## [1] 5968

``` r
head(topGenesAge) %>% kable()
```

|         |       logFC|   AveExpr|          t|  P.Value|  adj.P.Val|         B|
|---------|-----------:|---------:|----------:|--------:|----------:|---------:|
| CYP26A1 |   0.0994021|  6.154112|   19.99919|        0|          0|  22.44215|
| AP1M2   |  -0.0918790|  2.899708|  -19.37828|        0|          0|  21.87479|
| CRABP1  |  -0.0741114|  2.894697|  -19.18767|        0|          0|  21.69703|
| CDH3    |  -0.0601161|  4.131228|  -17.73228|        0|          0|  20.28048|
| COL3A1  |   0.1272345|  2.570885|   17.27234|        0|          0|  19.80961|
| ESRP1   |  -0.0942202|  4.839029|  -16.55618|        0|          0|  19.05243|

``` r
sample_genes <- rownames(topGenesAge)[1:6]
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-49-1.png)

``` r
temp <- venn.diagram(list(My_TopTable = rownames(topGenesAge), Paper_genes = paper_de_genes),fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4, lty =2, fontfamily =3, filename = NULL, main = "Comparison of paper DE genes with my all DE genes (from TopTable)", category.names = c("My topTable", "Paper genes"))

grid::grid.newpage()
grid::grid.draw(temp)
```

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-50-1.png) \# Finding discrepancies

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

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-53-1.png)

``` r
#reminder
summary(time_course_res)
```

    ##        v12v0 v24v12 v36v24 v72v36 v96v72
    ## Down    2894     94    253   1840    116
    ## NotSig  8417  14019  13729  10733  14090
    ## Up      3022    220    351   1760    127

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
| ABHD6  |      0|       1|       1|      -1|       0|
| ABLIM1 |      0|       0|       1|      -1|       0|
| ABTB2  |      0|       0|       1|       0|       0|
| ACSS3  |      0|       1|       1|       0|       0|
| ACVRL1 |     -1|       0|       1|       1|       0|
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

    ## [1] 941

``` r
#plot genes that are DE in paper at 36 to 24 but not in my results
set.seed(123)
sample_genes <- sample(s, 4)
plotGenes(sample_genes, cleaned_log_cpm_df, metadata)
```

    ## Using gene as id variables

    ## Joining, by = "samples"

![](GSE75748_data_analysis_files/figure-markdown_github/unnamed-chunk-59-1.png)
