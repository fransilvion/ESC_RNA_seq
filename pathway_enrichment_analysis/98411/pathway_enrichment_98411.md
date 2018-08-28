Pathway enrichment analysis, GSE98411
================
German Novakovskiy
August 23, 2018

Analysis of first day of differentiation 0h VS 24 h (formation of anterior Primitive streak, mesodendoderm)
-----------------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE98411/DEgenes_0h_24h_98411.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
#absolute Log fold changes for ermineR
ermineInputGeneScores <- DEgenes_0h_24h_98411 %>% 
  rownames_to_column("gene") %>%
  mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

#exact log fold changes for fgsea
# scores forfgsea
scoresFGSEADF <- DEgenes_0h_24h_98411 %>% 
  rownames_to_column("gene") %>%
  #mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(logFC)) %>% 
  column_to_rownames("gene")


scoresFGSEA <- scoresFGSEADF$logFC
names(scoresFGSEA) <- rownames(scoresFGSEADF)

#for WikiPathways
eg <-  bitr(rownames(ermineInputGeneScores), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
ermineInputGeneScoresWiki <- ermineInputGeneScores %>% rownames_to_column("SYMBOL") %>% 
  filter(SYMBOL %in% eg$SYMBOL)
ermineInputGeneScoresWiki <- merge(ermineInputGeneScoresWiki, eg, sort=FALSE)
rownames(ermineInputGeneScoresWiki) <- ermineInputGeneScoresWiki$ENTREZID
ermineInputGeneScoresWiki <- ermineInputGeneScoresWiki %>% dplyr::select("absolute_logFC")

scoresFGSEADFWiki <- scoresFGSEADF %>% rownames_to_column("SYMBOL") %>% 
  filter(SYMBOL %in% eg$SYMBOL)
scoresFGSEADFWiki <- merge(scoresFGSEADFWiki, eg, sort=FALSE)
rownames(scoresFGSEADFWiki) <- scoresFGSEADFWiki$ENTREZID
scoresFGSEADFWiki <- scoresFGSEADFWiki %>% dplyr::select("logFC")

scoresFGSEAWiki <- scoresFGSEADFWiki$logFC
names(scoresFGSEAWiki) <- rownames(scoresFGSEADFWiki)
```

### Reactome

#### ErmineR Reactome pathways

``` r
enrichmentResultReactome <- precRecall(scores = ermineInputGeneScores,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/Human_Reactome_August_01_2018_symbol.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

#enrichmentResultReactome$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

``` r
enrichmentResultReactome$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  head(10) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                    | Corrected p-value | Corrected MF p-value |
|:-------------------------------------------|:------------------|:---------------------|
| TCF dependent signaling in response to WNT | 0.0000000         | 0.0000000            |
| Neuronal System                            | 0.0000000         | 0.0000000            |
| Extracellular matrix organization          | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)       | 0.0000000         | 0.0000000            |
| G alpha (s) signalling events              | 0.0247714         | 0.0674333            |
| Class B 2 (Secretin family receptors)      | 0.0289000         | 0.0693600            |
| G alpha (i) signalling events              | 0.0000000         | 0.0908286            |
| PI3K Cascade                               | 0.2968818         | 0.3757000            |
| IRS-mediated signalling                    | 0.2578769         | 0.3861040            |
| IRS-related events triggered by IGF1R      | 0.2526074         | 0.3912615            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 4 enriched Reactome pathways with ermineR

#### FGSEA Reactome pathways

``` r
pathwaysReactome <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/Human_Reactome_August_01_2018_symbol.gmt")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysReactome, scoresFGSEA, minSize=15, maxSize=300, nperm=10000)

#up-regulated pathways
upPathwaysReactome <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))

#down-regulated pathways
downPathwaysReactome <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 0 up-regulated and 0 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome98411_024 <- upPathwaysReactome
save(upPathwaysReactome98411_024, file="upPathwaysReactome98411_024.Rdata")
```

### KEGG pathways

#### ErmineR KEGG pathways

``` r
enrichmentResultKEGG <- precRecall(scores = ermineInputGeneScores,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/KeggPathways.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

#enrichmentResultKEGG$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

``` r
enrichmentResultKEGG$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  head(10) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                      | Corrected p-value | Corrected MF p-value |
|:---------------------------------------------|:------------------|:---------------------|
| Kegg wnt signaling pathway                   | 0.0000000         | 0.0000000            |
| Kegg pathways in cancer                      | 0.0051000         | 0.0000000            |
| Kegg focal adhesion                          | 0.0714000         | 0.0000000            |
| Kegg mapk signaling pathway                  | 0.0918000         | 0.0038250            |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.0061200            |
| Kegg melanogenesis                           | 0.0707625         | 0.0127500            |
| Kegg calcium signaling pathway               | 0.0981750         | 0.0374000            |
| Kegg ecm receptor interaction                | 0.0790500         | 0.0401625            |
| Kegg basal cell carcinoma                    | 0.0890182         | 0.0415286            |
| Kegg cytokine cytokine receptor interaction  | 0.0918000         | 0.0443700            |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 14 enriched KEGG pathways with ermineR.

#### FGSEA KEGG pathways

``` r
pathwaysKEGG <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/KeggPathways.gmt")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysKEGG, scoresFGSEA, minSize=15, maxSize=300, nperm=10000)

#up-regulated pathways
upPathwaysKEGG <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))

#down-regulated pathways
downPathwaysKEGG <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 5 up-regulated and 0 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                             |      padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------|---------:|----------:|---------:|-------------:|
| KEGG\_WNT\_SIGNALING\_PATHWAY       |  0.007429|  0.5319154|  1.828706|             0|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY  |  0.007429|  0.7310908|  2.031668|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY |  0.007429|  0.6162020|  1.924140|             0|
| KEGG\_PATHWAYS\_IN\_CANCER          |  0.007429|  0.4818222|  1.820277|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA        |  0.007429|  0.7290053|  2.057050|             0|

``` r
upPathwaysKEGG98411_024 <- upPathwaysKEGG
save(upPathwaysKEGG98411_024, file="upPathwaysKEGG98411_024.Rdata")
```

### WikiPathways

#### ErmineR Wiki pathways

``` r
enrichmentResultWiki <- precRecall(scores = ermineInputGeneScoresWiki,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/Wikipathways_changed_column.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

#enrichmentResultReactome$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

``` r
enrichmentResultWiki$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  filter(CorrectedMFPvalue <= 0.05) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                                                                                        | Corrected p-value | Corrected MF p-value |
|:---------------------------------------------------------------------------------------------------------------|:------------------|:---------------------|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens | 0.000000          | 0.000000             |
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       | 0.000000          | 0.000000             |
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            | 0.000000          | 0.000000             |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  | 0.000000          | 0.000000             |
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                | 0.000000          | 0.003463             |
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                            | 0.049860          | 0.003957             |
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                             | 0.003957          | 0.004617             |
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens        | 0.000000          | 0.005540             |

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 8 enriched Wiki pathways with ermineR

#### FGSEA Wiki pathways

``` r
pathwaysWiki <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/Wikipathways.gmt")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysWiki, scoresFGSEAWiki, minSize=15, maxSize=300, nperm=10000)

#up-regulated pathways
upPathwaysWiki <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))

#down-regulated pathways
downPathwaysWiki <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 12 up-regulated and 0 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                        |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                             |  0.0090490|  0.8388928|  2.195693|             0|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                           |  0.0090490|  0.5915900|  1.941045|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       |  0.0090490|  0.5915214|  2.046894|             0|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            |  0.0090490|  0.5857901|  2.011356|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                   |  0.0090490|  0.7598814|  2.020384|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  |  0.0090490|  0.7089267|  1.947049|             0|
| Breast cancer pathway%WikiPathways\_20180810%WP4262%Homo sapiens                                               |  0.0090490|  0.5268981|  1.796892|             0|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens |  0.0152098|  0.6085278|  1.927874|             1|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                     |  0.0199942|  0.7574625|  1.873863|             2|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                      |  0.0433329|  0.4263492|  1.545713|             6|
| Dopaminergic Neurogenesis%WikiPathways\_20180810%WP2855%Homo sapiens                                           |  0.0433329|  0.7650124|  1.855146|             7|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                |  0.0456593|  0.5300696|  1.739193|             8|

``` r
upPathwaysWiki98411_024 <- upPathwaysWiki
save(upPathwaysWiki98411_024, file="upPathwaysWiki98411_024.Rdata")
```

Analysis of second-fourth days of differentiation 24h VS 96h (formation of definitive endoderm from APS)
--------------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE98411/DEgenes_24h_72h_98411.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
#absolute Log fold changes for ermineR
ermineInputGeneScores <- DEgenes_24h_72h_98411 %>% 
  rownames_to_column("gene") %>%
  mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

#exact log fold changes for fgsea
# scores forfgsea
scoresFGSEADF <- DEgenes_24h_72h_98411 %>% 
  rownames_to_column("gene") %>%
  #mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(logFC)) %>% 
  column_to_rownames("gene")


scoresFGSEA <- scoresFGSEADF$logFC
names(scoresFGSEA) <- rownames(scoresFGSEADF)

#for WikiPathways
eg <-  bitr(rownames(ermineInputGeneScores), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
ermineInputGeneScoresWiki <- ermineInputGeneScores %>% rownames_to_column("SYMBOL") %>% 
  filter(SYMBOL %in% eg$SYMBOL)
ermineInputGeneScoresWiki <- merge(ermineInputGeneScoresWiki, eg, sort=FALSE)
rownames(ermineInputGeneScoresWiki) <- ermineInputGeneScoresWiki$ENTREZID
ermineInputGeneScoresWiki <- ermineInputGeneScoresWiki %>% dplyr::select("absolute_logFC")

scoresFGSEADFWiki <- scoresFGSEADF %>% rownames_to_column("SYMBOL") %>% 
  filter(SYMBOL %in% eg$SYMBOL)
scoresFGSEADFWiki <- merge(scoresFGSEADFWiki, eg, sort=FALSE)
rownames(scoresFGSEADFWiki) <- scoresFGSEADFWiki$ENTREZID
scoresFGSEADFWiki <- scoresFGSEADFWiki %>% dplyr::select("logFC")

scoresFGSEAWiki <- scoresFGSEADFWiki$logFC
names(scoresFGSEAWiki) <- rownames(scoresFGSEADFWiki)
```

### Reactome

#### ErmineR Reactome pathways

``` r
enrichmentResultReactome <- precRecall(scores = ermineInputGeneScores,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/Human_Reactome_August_01_2018_symbol.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

#enrichmentResultReactome$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

``` r
enrichmentResultReactome$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  head(10) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                                           | Corrected p-value | Corrected MF p-value |
|:------------------------------------------------------------------|:------------------|:---------------------|
| G alpha (i) signalling events                                     | 0.0000000         | 0.0000000            |
| Extracellular matrix organization                                 | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)                              | 0.0000000         | 0.0000000            |
| Post-translational protein phosphorylation                        | 0.0433500         | 0.0433500            |
| Collagen formation                                                | 0.0693600         | 0.0924800            |
| Signaling by FGFR4                                                | 0.0954956         | 0.1045905            |
| Signaling by FGFR2 in disease                                     | 0.0998364         | 0.1098200            |
| IRS-related events triggered by IGF1R                             | 0.1045905         | 0.1156000            |
| Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R) | 0.1045905         | 0.1156000            |
| IRS-mediated signalling                                           | 0.1098200         | 0.1220222            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 4 enriched Reactome pathways with ermineR

#### FGSEA Reactome pathways

``` r
pathwaysReactome <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/Human_Reactome_August_01_2018_symbol.gmt")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysReactome, scoresFGSEA, minSize=15, maxSize=300, nperm=10000)

#up-regulated pathways
upPathwaysReactome <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))

#down-regulated pathways
downPathwaysReactome <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 41 up-regulated and 230 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                              |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                 |  0.0072082|  0.6831081|  1.895106|             3|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                |  0.0072082|  0.4963027|  1.671817|             0|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                        |  0.0072082|  0.7137058|  2.057333|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                            |  0.0072082|  0.5782104|  1.971336|             0|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                  |  0.0072082|  0.6030993|  1.900734|             0|
| G ALPHA (I) SIGNALLING EVENTS%REACTOME%R-HSA-418594.5                                                |  0.0072082|  0.4870838|  1.643245|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                          |  0.0072082|  0.6005742|  1.852461|             0|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                   |  0.0072082|  0.6404478|  1.990870|             0|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                        |  0.0072082|  0.6221719|  1.948357|             0|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                           |  0.0072082|  0.6815479|  1.926785|             2|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090 |  0.0072082|  0.6869231|  1.971729|             0|
| PLASMA LIPOPROTEIN REMODELING%REACTOME DATABASE ID RELEASE 65%8963899                                |  0.0072082|  0.8505300|  1.997870|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                          |  0.0072082|  0.6356000|  1.990408|             0|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                           |  0.0072301|  0.6394826|  1.835557|             5|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                  |  0.0077697|  0.5934265|  1.775288|             7|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                      |  0.0081142|  0.5905612|  1.752989|             8|
| MUSCLE CONTRACTION%REACTOME%R-HSA-397014.2                                                           |  0.0092935|  0.4738146|  1.555706|            13|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                             |  0.0093857|  0.5952607|  1.750482|            12|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                      |  0.0097460|  0.6478909|  1.777586|            12|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                                 |  0.0109394|  0.6707494|  1.752868|            14|
| STRIATED MUSCLE CONTRACTION%REACTOME DATABASE ID RELEASE 65%390522                                   |  0.0111198|  0.6993129|  1.769505|            14|
| CELL JUNCTION ORGANIZATION%REACTOME DATABASE ID RELEASE 65%446728                                    |  0.0111626|  0.5570686|  1.693579|            17|
| VISUAL PHOTOTRANSDUCTION%REACTOME DATABASE ID RELEASE 65%2187338                                     |  0.0139487|  0.5728466|  1.689446|            22|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                      |  0.0154924|  0.5332283|  1.641671|            28|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                         |  0.0163553|  0.5155537|  1.619802|            31|
| DEATH RECEPTOR SIGNALLING%REACTOME%R-HSA-73887.3                                                     |  0.0202579|  0.4690446|  1.530915|            45|
| CELL-CELL COMMUNICATION%REACTOME%R-HSA-1500931.3                                                     |  0.0212224|  0.4939098|  1.565365|            46|
| RETINOID METABOLISM AND TRANSPORT%REACTOME%R-HSA-975634.2                                            |  0.0212598|  0.6348848|  1.682931|            39|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                               |  0.0251573|  0.5127598|  1.589975|            55|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                         |  0.0257980|  0.5570028|  1.627462|            54|
| NON-INTEGRIN MEMBRANE-ECM INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000171                       |  0.0272188|  0.5919711|  1.651709|            55|
| P75 NTR RECEPTOR-MEDIATED SIGNALLING%REACTOME%R-HSA-193704.1                                         |  0.0359507|  0.4838781|  1.521893|            86|
| METABOLISM OF FAT-SOLUBLE VITAMINS%REACTOME DATABASE ID RELEASE 65%6806667                           |  0.0362071|  0.5997627|  1.622078|            75|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                            |  0.0372134|  0.5765545|  1.617105|            80|
| ACTIVATION OF MATRIX METALLOPROTEINASES%REACTOME%R-HSA-1592389.1                                     |  0.0374910|  0.7073713|  1.661595|            72|
| GAP JUNCTION TRAFFICKING AND REGULATION%REACTOME%R-HSA-157858.1                                      |  0.0377319|  0.7055499|  1.657316|            73|
| CELL-CELL JUNCTION ORGANIZATION%REACTOME%R-HSA-421270.4                                              |  0.0381253|  0.5617993|  1.601146|            85|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                    |  0.0386599|  0.6378608|  1.640827|            80|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                           |  0.0417048|  0.5920108|  1.601112|            90|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                           |  0.0457639|  0.4349069|  1.434776|           124|
| NRAGE SIGNALS DEATH THROUGH JNK%REACTOME DATABASE ID RELEASE 65%193648                               |  0.0490006|  0.5251015|  1.548636|           118|

``` r
upPathwaysReactome98411_2496 <- upPathwaysReactome
save(upPathwaysReactome98411_2496, file="upPathwaysReactome98411_2496.Rdata")
```

### KEGG pathways

#### ErmineR KEGG pathways

``` r
enrichmentResultKEGG <- precRecall(scores = ermineInputGeneScores,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/KeggPathways.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

#enrichmentResultKEGG$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

``` r
enrichmentResultKEGG$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  head(10) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                           | Corrected p-value | Corrected MF p-value |
|:--------------------------------------------------|:------------------|:---------------------|
| Kegg axon guidance                                | 0.0214200         | 0.0153000            |
| Kegg ecm receptor interaction                     | 0.0153000         | 0.0229500            |
| Kegg neuroactive ligand receptor interaction      | 0.0280500         | 0.0255000            |
| Kegg cytokine cytokine receptor interaction       | 0.0204000         | 0.0382500            |
| Kegg intestinal immune network for iga production | 0.0655714         | 0.1499400            |
| Kegg endocytosis                                  | 0.4125333         | 0.1708500            |
| Kegg complement and coagulation cascades          | 0.1239300         | 0.1945286            |
| Kegg cell adhesion molecules cams                 | 0.1088000         | 0.2161125            |
| Kegg basal cell carcinoma                         | 0.1772250         | 0.2340900            |
| Kegg focal adhesion                               | 0.0153000         | 0.2448000            |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 4 enriched KEGG pathways with ermineR.

#### FGSEA KEGG pathways

``` r
pathwaysKEGG <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/KeggPathways.gmt")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysKEGG, scoresFGSEA, minSize=15, maxSize=300, nperm=10000)

#up-regulated pathways
upPathwaysKEGG <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))

#down-regulated pathways
downPathwaysKEGG <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 10 up-regulated and 14 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                         |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_AXON\_GUIDANCE                            |  0.0083836|  0.5525302|  1.801816|             0|
| KEGG\_FOCAL\_ADHESION                           |  0.0083836|  0.5099686|  1.713044|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                |  0.0155789|  0.5948298|  1.798833|             5|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES    |  0.0157084|  0.7010039|  1.839894|             5|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION |  0.0280379|  0.5188667|  1.645433|            20|
| KEGG\_PPAR\_SIGNALING\_PATHWAY                  |  0.0344620|  0.5762875|  1.665854|            33|
| KEGG\_LEUKOCYTE\_TRANSENDOTHELIAL\_MIGRATION    |  0.0344620|  0.5133445|  1.611694|            36|
| KEGG\_DILATED\_CARDIOMYOPATHY                   |  0.0344620|  0.5429973|  1.639791|            35|
| KEGG\_HYPERTROPHIC\_CARDIOMYOPATHY\_HCM         |  0.0350122|  0.5495000|  1.637502|            38|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                  |  0.0463155|  0.4155230|  1.414208|            69|

``` r
upPathwaysKEGG98411_2496 <- upPathwaysKEGG
save(upPathwaysKEGG98411_2496, file="upPathwaysKEGG98411_2496.Rdata")
```

### WikiPathways

#### ErmineR Wiki pathways

``` r
enrichmentResultWiki <- precRecall(scores = ermineInputGeneScoresWiki,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/Wikipathways_changed_column.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

#enrichmentResultReactome$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

``` r
enrichmentResultWiki$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  filter(CorrectedMFPvalue <= 0.05) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                                                            | Corrected p-value | Corrected MF p-value |
|:-----------------------------------------------------------------------------------|:------------------|:---------------------|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                | 0.0672714         | 0.0000000            |
| Peptide GPCRs%WikiPathways\_20180810%WP24%Homo sapiens                             | 0.0461667         | 0.0346250            |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens      | 0.0277000         | 0.0369333            |
| Hematopoietic Stem Cell Differentiation%WikiPathways\_20180810%WP2849%Homo sapiens | 0.0277000         | 0.0415500            |

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 4 enriched Wiki pathways with ermineR

#### FGSEA Wiki pathways

``` r
pathwaysWiki <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/Wikipathways.gmt")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysWiki, scoresFGSEAWiki, minSize=15, maxSize=300, nperm=10000)

#up-regulated pathways
upPathwaysWiki <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))

#down-regulated pathways
downPathwaysWiki <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 21 up-regulated and 16 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                                                                     |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| GPCRs, Class A Rhodopsin-like%WikiPathways\_20180810%WP455%Homo sapiens                                                                                     |  0.0115793|  0.6165773|  1.810599|             2|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                                                                    |  0.0115793|  0.4972543|  1.667416|             2|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                                                                  |  0.0115793|  0.5971515|  1.795074|             1|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                                                                   |  0.0115793|  0.4637730|  1.567497|             6|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                                                                         |  0.0115793|  0.6663857|  2.169645|             0|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                                                                |  0.0115793|  0.5522147|  1.726074|             1|
| Peptide GPCRs%WikiPathways\_20180810%WP24%Homo sapiens                                                                                                      |  0.0115793|  0.8203847|  1.972412|             0|
| Hypothesized Pathways in Pathogenesis of Cardiovascular Disease%WikiPathways\_20180810%WP3668%Homo sapiens                                                  |  0.0115793|  0.7269469|  1.860057|             3|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                                                                       |  0.0115793|  0.5285761|  1.678666|             5|
| VEGFA-VEGFR2 Signaling Pathway%WikiPathways\_20180810%WP3888%Homo sapiens                                                                                   |  0.0115793|  0.4667619|  1.588535|             3|
| Myometrial Relaxation and Contraction Pathways%WikiPathways\_20180810%WP289%Homo sapiens                                                                    |  0.0119371|  0.5005985|  1.634172|             7|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                                                                      |  0.0186073|  0.5134752|  1.632640|            14|
| Complement and Coagulation Cascades%WikiPathways\_20180810%WP558%Homo sapiens                                                                               |  0.0249955|  0.6840384|  1.764332|            17|
| Oligodendrocyte Specification and differentiation(including remyelination), leading to Myelin Components for CNS%WikiPathways\_20180810%WP4304%Homo sapiens |  0.0254715|  0.7313979|  1.758466|            17|
| Striated Muscle Contraction%WikiPathways\_20180810%WP383%Homo sapiens                                                                                       |  0.0259989|  0.6719445|  1.762297|            20|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                                                                   |  0.0259989|  0.4293101|  1.472537|            26|
| Statin Pathway%WikiPathways\_20180810%WP430%Homo sapiens                                                                                                    |  0.0268216|  0.7328210|  1.717216|            21|
| Canonical and Non-Canonical TGF-B signaling%WikiPathways\_20180810%WP3874%Homo sapiens                                                                      |  0.0268216|  0.7362772|  1.725315|            21|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                                                                                         |  0.0437253|  0.5608112|  1.638129|            43|
| Vitamin B12 Metabolism%WikiPathways\_20180810%WP1533%Homo sapiens                                                                                           |  0.0446761|  0.6305813|  1.688125|            43|
| Photodynamic therapy-induced unfolded protein response%WikiPathways\_20180810%WP3613%Homo sapiens                                                           |  0.0473382|  0.6431187|  1.673816|            46|

``` r
upPathwaysWiki98411_2496 <- upPathwaysWiki
save(upPathwaysWiki98411_2496, file="upPathwaysWiki98411_2496.Rdata")
```

Analysis of whole differentiation process 0h VS 96 h (formation of Definitive endoderm from hESC)
-------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE98411/DEgenes_0h_72h_98411.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
#absolute Log fold changes for ermineR
ermineInputGeneScores <- DEgenes_0h_72h_98411 %>% 
  rownames_to_column("gene") %>%
  mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

#exact log fold changes for fgsea
# scores forfgsea
scoresFGSEADF <- DEgenes_0h_72h_98411 %>% 
  rownames_to_column("gene") %>%
  #mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(logFC)) %>% 
  column_to_rownames("gene")


scoresFGSEA <- scoresFGSEADF$logFC
names(scoresFGSEA) <- rownames(scoresFGSEADF)

#for WikiPathways
eg <-  bitr(rownames(ermineInputGeneScores), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
ermineInputGeneScoresWiki <- ermineInputGeneScores %>% rownames_to_column("SYMBOL") %>% 
  filter(SYMBOL %in% eg$SYMBOL)
ermineInputGeneScoresWiki <- merge(ermineInputGeneScoresWiki, eg, sort=FALSE)
rownames(ermineInputGeneScoresWiki) <- ermineInputGeneScoresWiki$ENTREZID
ermineInputGeneScoresWiki <- ermineInputGeneScoresWiki %>% dplyr::select("absolute_logFC")

scoresFGSEADFWiki <- scoresFGSEADF %>% rownames_to_column("SYMBOL") %>% 
  filter(SYMBOL %in% eg$SYMBOL)
scoresFGSEADFWiki <- merge(scoresFGSEADFWiki, eg, sort=FALSE)
rownames(scoresFGSEADFWiki) <- scoresFGSEADFWiki$ENTREZID
scoresFGSEADFWiki <- scoresFGSEADFWiki %>% dplyr::select("logFC")

scoresFGSEAWiki <- scoresFGSEADFWiki$logFC
names(scoresFGSEAWiki) <- rownames(scoresFGSEADFWiki)
```

### Reactome

#### ErmineR Reactome pathways

``` r
enrichmentResultReactome <- precRecall(scores = ermineInputGeneScores,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/Human_Reactome_August_01_2018_symbol.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

#enrichmentResultReactome$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

``` r
enrichmentResultReactome$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  head(10) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                        | Corrected p-value | Corrected MF p-value |
|:-----------------------------------------------|:------------------|:---------------------|
| SLC-mediated transmembrane transport           | 0.0289000         | 0.0000000            |
| G alpha (i) signalling events                  | 0.0000000         | 0.0000000            |
| Extracellular matrix organization              | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)           | 0.0000000         | 0.0000000            |
| Neuronal System                                | 0.1300500         | 0.0192667            |
| Collagen formation                             | 0.0346800         | 0.0231200            |
| Platelet activation, signaling and aggregation | 0.1156000         | 0.0330286            |
| Post-translational protein phosphorylation     | 0.0990857         | 0.0939250            |
| Signaling by Interleukins                      | 0.9130441         | 0.1412889            |
| Iron uptake and transport                      | 0.3236800         | 0.2489846            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 7 enriched Reactome pathways with ermineR

#### FGSEA Reactome pathways

``` r
pathwaysReactome <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/Human_Reactome_August_01_2018_symbol.gmt")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysReactome, scoresFGSEA, minSize=15, maxSize=300, nperm=10000)

#up-regulated pathways
upPathwaysReactome <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))

#down-regulated pathways
downPathwaysReactome <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 33 up-regulated and 163 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                          |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                             |  0.0106160|  0.6906764|  1.908196|             2|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                                     |  0.0106160|  0.6582874|  1.923287|             1|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                              |  0.0106160|  0.6095565|  1.827407|             1|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                        |  0.0106160|  0.5530032|  1.918271|             0|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                                       |  0.0106160|  0.7017027|  1.892675|             3|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                                      |  0.0106160|  0.5937514|  1.844573|             1|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                                    |  0.0106160|  0.5818321|  1.833946|             0|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                  |  0.0106160|  0.7193802|  1.962679|             2|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                                         |  0.0106160|  0.6007229|  1.767216|             4|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                       |  0.0106160|  0.5028360|  1.677825|             2|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090             |  0.0106160|  0.6493038|  1.857627|             2|
| INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING%REACTOME DATABASE ID RELEASE 65%6785807                               |  0.0106761|  0.5562357|  1.713633|             5|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                              |  0.0109899|  0.5414102|  1.719442|             6|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                                      |  0.0117332|  0.5419548|  1.708252|             7|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                       |  0.0118007|  0.6202001|  1.774363|             7|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                                |  0.0118007|  0.7010078|  1.794661|             6|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                                     |  0.0118007|  0.5324235|  1.698424|             8|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                               |  0.0135460|  0.5408211|  1.690260|            11|
| TRANSPORT OF BILE SALTS AND ORGANIC ACIDS, METAL IONS AND AMINE COMPOUNDS%REACTOME DATABASE ID RELEASE 65%425366 |  0.0135460|  0.6137404|  1.762604|            10|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                                         |  0.0145247|  0.6956923|  1.781053|            10|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                                    |  0.0151323|  0.6072272|  1.743899|            12|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                        |  0.0181340|  0.6339416|  1.766483|            16|
| IRE1ALPHA ACTIVATES CHAPERONES%REACTOME DATABASE ID RELEASE 65%381070                                            |  0.0240544|  0.5578257|  1.663922|            26|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                                     |  0.0246427|  0.5141659|  1.625320|            29|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                            |  0.0322354|  0.4337733|  1.482052|            49|
| SIGNALING BY WNT IN CANCER%REACTOME DATABASE ID RELEASE 65%4791275                                               |  0.0329729|  0.6290613|  1.696742|            40|
| SIGNALING BY INTERLEUKINS%REACTOME%R-HSA-449147.10                                                               |  0.0377206|  0.3884707|  1.379406|            66|
| XBP1(S) ACTIVATES CHAPERONE GENES%REACTOME%R-HSA-381038.2                                                        |  0.0394428|  0.5462921|  1.620926|            57|
| SIGNALING BY FGFR1%REACTOME DATABASE ID RELEASE 65%5654736                                                       |  0.0408245|  0.6004543|  1.658931|            56|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                                       |  0.0420987|  0.5821210|  1.638884|            62|
| DISEASES OF SIGNAL TRANSDUCTION%REACTOME DATABASE ID RELEASE 65%5663202                                          |  0.0437593|  0.3828059|  1.359705|            87|
| RUNX1 INTERACTS WITH CO-FACTORS WHOSE PRECISE EFFECT ON RUNX1 TARGETS IS NOT KNOWN%REACTOME%R-HSA-8939243.1      |  0.0444753|  0.6012725|  1.651107|            66|
| RHO GTPASE CYCLE%REACTOME DATABASE ID RELEASE 65%194840                                                          |  0.0487773|  0.4496612|  1.487066|            92|

``` r
upPathwaysReactome98411_096 <- upPathwaysReactome
save(upPathwaysReactome98411_096, file="upPathwaysReactome98411_096.Rdata")
```

### KEGG pathways

#### ErmineR KEGG pathways

``` r
enrichmentResultKEGG <- precRecall(scores = ermineInputGeneScores,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/KeggPathways.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

#enrichmentResultKEGG$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

``` r
enrichmentResultKEGG$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  head(10) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                      | Corrected p-value | Corrected MF p-value |
|:---------------------------------------------|:------------------|:---------------------|
| Kegg focal adhesion                          | 0.0229500         | 0.0000000            |
| Kegg calcium signaling pathway               | 0.0382500         | 0.0076500            |
| Kegg pathways in cancer                      | 0.1355143         | 0.0076500            |
| Kegg ecm receptor interaction                | 0.0306000         | 0.0102000            |
| Kegg cytokine cytokine receptor interaction  | 0.0255000         | 0.0122400            |
| Kegg neuroactive ligand receptor interaction | 0.0306000         | 0.0612000            |
| Kegg wnt signaling pathway                   | 0.1281375         | 0.0612000            |
| Kegg tgf beta signaling pathway              | 0.1096500         | 0.0631125            |
| Kegg mapk signaling pathway                  | 0.1338750         | 0.0748000            |
| Kegg cell adhesion molecules cams            | 0.1224000         | 0.1168364            |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 5 enriched KEGG pathways with ermineR.

#### FGSEA KEGG pathways

``` r
pathwaysKEGG <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/KeggPathways.gmt")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysKEGG, scoresFGSEA, minSize=15, maxSize=300, nperm=10000)

#up-regulated pathways
upPathwaysKEGG <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))

#down-regulated pathways
downPathwaysKEGG <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 14 up-regulated and 12 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                         |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION |  0.0047795|  0.6292265|  2.013032|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY             |  0.0047795|  0.6570327|  2.018304|             0|
| KEGG\_AXON\_GUIDANCE                            |  0.0047795|  0.5635419|  1.863018|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                      |  0.0047795|  0.4809057|  1.693448|             0|
| KEGG\_FOCAL\_ADHESION                           |  0.0096081|  0.4930680|  1.680894|             2|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                |  0.0113890|  0.5720013|  1.736038|             5|
| KEGG\_SMALL\_CELL\_LUNG\_CANCER                 |  0.0226872|  0.5395996|  1.687651|            16|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                  |  0.0271277|  0.4328452|  1.499464|            24|
| KEGG\_WNT\_SIGNALING\_PATHWAY                   |  0.0325661|  0.4734545|  1.569656|            29|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES    |  0.0412680|  0.6656616|  1.740144|            33|
| KEGG\_BASAL\_CELL\_CARCINOMA                    |  0.0412680|  0.6020553|  1.689423|            34|
| KEGG\_REGULATION\_OF\_ACTIN\_CYTOSKELETON       |  0.0459234|  0.4397563|  1.495214|            54|
| KEGG\_PANCREATIC\_CANCER                        |  0.0489942|  0.5231049|  1.596138|            55|
| KEGG\_CHRONIC\_MYELOID\_LEUKEMIA                |  0.0496962|  0.5090933|  1.575676|            67|

``` r
upPathwaysKEGG98411_096 <- upPathwaysKEGG
save(upPathwaysKEGG98411_096, file="upPathwaysKEGG98411_096.Rdata")
```

### WikiPathways

#### ErmineR Wiki pathways

``` r
enrichmentResultWiki <- precRecall(scores = ermineInputGeneScoresWiki,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/Wikipathways_changed_column.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

#enrichmentResultReactome$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

``` r
enrichmentResultWiki$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  filter(CorrectedMFPvalue <= 0.05) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                                                       | Corrected p-value | Corrected MF p-value |
|:------------------------------------------------------------------------------|:------------------|:---------------------|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens           | 0.006925          | 0.000000             |
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens      | 0.000000          | 0.000000             |
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens           | 0.000000          | 0.000000             |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens | 0.000000          | 0.020775             |

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 4 enriched Wiki pathways with ermineR

#### FGSEA Wiki pathways

``` r
pathwaysWiki <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/Wikipathways.gmt")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysWiki, scoresFGSEAWiki, minSize=15, maxSize=300, nperm=10000)

#up-regulated pathways
upPathwaysWiki <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))

#down-regulated pathways
downPathwaysWiki <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 38 up-regulated and 13 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                                                                     |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                                                                    |  0.0042290|  0.4896696|  1.664664|             0|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                                                                          |  0.0042290|  0.8160812|  2.128388|             0|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                                                                   |  0.0042290|  0.5131921|  1.766090|             0|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                                                                      |  0.0042290|  0.5985506|  1.920877|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                                                                    |  0.0042290|  0.5754341|  1.912614|             0|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                                                                         |  0.0042290|  0.6917237|  2.278921|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                                                                |  0.0042290|  0.7991357|  2.116948|             0|
| Peptide GPCRs%WikiPathways\_20180810%WP24%Homo sapiens                                                                                                      |  0.0042290|  0.8320267|  1.998107|             0|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                                                                   |  0.0042290|  0.4962045|  1.736088|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens                                                        |  0.0042290|  0.5683104|  1.872826|             0|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                                                                                       |  0.0052825|  0.4586217|  1.606298|             1|
| GPCRs, Class A Rhodopsin-like%WikiPathways\_20180810%WP455%Homo sapiens                                                                                     |  0.0054559|  0.6254128|  1.839924|             1|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                                                                  |  0.0054559|  0.6396532|  1.934970|             1|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                                                               |  0.0054559|  0.7183660|  1.961242|             1|
| Canonical and Non-Canonical TGF-B signaling%WikiPathways\_20180810%WP3874%Homo sapiens                                                                      |  0.0057005|  0.8359753|  1.958198|             1|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                                                                        |  0.0057774|  0.5776968|  1.840484|             2|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                                                                       |  0.0057774|  0.5612577|  1.798091|             2|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                                                                         |  0.0066178|  0.5238820|  1.731541|             3|
| Spinal Cord Injury%WikiPathways\_20180810%WP2431%Homo sapiens                                                                                               |  0.0075508|  0.5592672|  1.742050|             4|
| Oncostatin M Signaling Pathway%WikiPathways\_20180810%WP2374%Homo sapiens                                                                                   |  0.0085744|  0.5865768|  1.761556|             5|
| TGF-B Signaling in Thyroid Cells for Epithelial-Mesenchymal Transition%WikiPathways\_20180810%WP3859%Homo sapiens                                           |  0.0106884|  0.7828937|  1.833859|             6|
| VEGFA-VEGFR2 Signaling Pathway%WikiPathways\_20180810%WP3888%Homo sapiens                                                                                   |  0.0111349|  0.4432686|  1.538087|            10|
| TGF-beta Receptor Signaling%WikiPathways\_20180810%WP560%Homo sapiens                                                                                       |  0.0119309|  0.5943945|  1.743237|             9|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                                                             |  0.0168653|  0.5165918|  1.645809|            16|
| Dopaminergic Neurogenesis%WikiPathways\_20180810%WP2855%Homo sapiens                                                                                        |  0.0171905|  0.7356367|  1.789514|            13|
| Hypothesized Pathways in Pathogenesis of Cardiovascular Disease%WikiPathways\_20180810%WP3668%Homo sapiens                                                  |  0.0200114|  0.6886972|  1.754184|            17|
| MAPK Signaling Pathway%WikiPathways\_20180810%WP382%Homo sapiens                                                                                            |  0.0210678|  0.4375281|  1.508667|            26|
| NRF2 pathway%WikiPathways\_20180810%WP2884%Homo sapiens                                                                                                     |  0.0215883|  0.5310219|  1.655963|            24|
| Simplified Interaction Map Between LOXL4 and Oxidative Stress Pathway%WikiPathways\_20180810%WP3670%Homo sapiens                                            |  0.0224456|  0.7641441|  1.789940|            20|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens                                              |  0.0228128|  0.5289301|  1.639693|            28|
| Chromosomal and microsatellite instability in colorectal cancer %WikiPathways\_20180810%WP4216%Homo sapiens                                                 |  0.0301070|  0.5240649|  1.620190|            38|
| Oligodendrocyte Specification and differentiation(including remyelination), leading to Myelin Components for CNS%WikiPathways\_20180810%WP4304%Homo sapiens |  0.0310127|  0.7303100|  1.753835|            32|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                                                                  |  0.0321763|  0.6932492|  1.712649|            35|
| ErbB Signaling Pathway%WikiPathways\_20180810%WP673%Homo sapiens                                                                                            |  0.0321763|  0.5692805|  1.660380|            41|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                                                                |  0.0345932|  0.4990197|  1.572651|            49|
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens                                                     |  0.0365232|  0.5597345|  1.641587|            49|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                                                                                         |  0.0423778|  0.5593290|  1.635614|            58|
| Aryl Hydrocarbon Receptor Pathway%WikiPathways\_20180810%WP2873%Homo sapiens                                                                                |  0.0477879|  0.6542665|  1.681899|            60|

``` r
upPathwaysWiki98411_096 <- upPathwaysWiki
save(upPathwaysWiki98411_096, file="upPathwaysWiki98411_096.Rdata")
```
