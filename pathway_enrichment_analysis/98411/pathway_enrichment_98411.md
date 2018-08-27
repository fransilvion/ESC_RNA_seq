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

| Pathway                                                           | Corrected p-value | Corrected MF p-value |
|:------------------------------------------------------------------|:------------------|:---------------------|
| TCF dependent signaling in response to WNT                        | 0.0000000         | 0.0000000            |
| Neuronal System                                                   | 0.0000000         | 0.0000000            |
| Extracellular matrix organization                                 | 0.0000000         | 0.0000000            |
| Class B 2 (Secretin family receptors)                             | 0.0115600         | 0.0289000            |
| Class A 1 (Rhodopsin-like receptors)                              | 0.0192667         | 0.0346800            |
| G alpha (i) signalling events                                     | 0.0000000         | 0.0481667            |
| G alpha (s) signalling events                                     | 0.0247714         | 0.0660571            |
| IRS-related events triggered by IGF1R                             | 0.2450720         | 0.3579154            |
| Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R) | 0.2450720         | 0.3579154            |
| Collagen formation                                                | 0.3468000         | 0.3620105            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 6 enriched Reactome pathways with ermineR

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
| Kegg pathways in cancer                      | 0.0000000         | 0.0000000            |
| Kegg mapk signaling pathway                  | 0.0726750         | 0.0051000            |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.0061200            |
| Kegg focal adhesion                          | 0.0627300         | 0.0076500            |
| Kegg melanogenesis                           | 0.0663000         | 0.0102000            |
| Kegg calcium signaling pathway               | 0.1160250         | 0.0367200            |
| Kegg cell adhesion molecules cams            | 0.0703800         | 0.0374000            |
| Kegg basal cell carcinoma                    | 0.1210091         | 0.0382500            |
| Kegg ecm receptor interaction                | 0.0590143         | 0.0415286            |

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

| pathway                             |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY  |  0.0122068|  0.7310908|  2.030777|             0|
| KEGG\_PATHWAYS\_IN\_CANCER          |  0.0122068|  0.4818222|  1.820212|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA        |  0.0122068|  0.7290053|  2.051424|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY       |  0.0138305|  0.5319154|  1.825567|             1|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY |  0.0138305|  0.6162020|  1.922657|             1|

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
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                | 0.000000          | 0.000000             |
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens | 0.000000          | 0.000000             |
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       | 0.000000          | 0.000000             |
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            | 0.000000          | 0.000000             |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  | 0.000000          | 0.000000             |
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                             | 0.000000          | 0.004617             |
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                            | 0.055400          | 0.006925             |
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens        | 0.003957          | 0.007914             |

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

There are 14 up-regulated and 0 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                        |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                             |  0.0104237|  0.8388928|  2.189173|             0|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                           |  0.0104237|  0.5915900|  1.936939|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       |  0.0104237|  0.5915214|  2.051230|             0|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            |  0.0104237|  0.5857901|  2.014432|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                   |  0.0104237|  0.7598814|  2.018717|             0|
| Breast cancer pathway%WikiPathways\_20180810%WP4262%Homo sapiens                                               |  0.0104237|  0.5268981|  1.795414|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  |  0.0153582|  0.7089267|  1.945526|             1|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens |  0.0153582|  0.6085278|  1.929377|             1|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                     |  0.0197801|  0.7574625|  1.883922|             2|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                |  0.0246004|  0.5300696|  1.735513|             3|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                      |  0.0266490|  0.4263492|  1.548492|             4|
| Dopaminergic Neurogenesis%WikiPathways\_20180810%WP2855%Homo sapiens                                           |  0.0266490|  0.7650124|  1.854734|             4|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens           |  0.0336766|  0.4681823|  1.604651|             6|
| DNA Damage Response (only ATM dependent)%WikiPathways\_20180810%WP710%Homo sapiens                             |  0.0350345|  0.5061212|  1.670647|             7|

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

| Pathway                                    | Corrected p-value | Corrected MF p-value |
|:-------------------------------------------|:------------------|:---------------------|
| G alpha (i) signalling events              | 0.0000000         | 0.0000000            |
| Extracellular matrix organization          | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)       | 0.0000000         | 0.0000000            |
| Post-translational protein phosphorylation | 0.0144500         | 0.0289000            |
| Collagen formation                         | 0.0578000         | 0.0924800            |
| Signaling by FGFR1                         | 0.1156000         | 0.1206261            |
| Insulin receptor signalling cascade        | 0.1206261         | 0.1261091            |
| Signaling by FGFR2 in disease              | 0.1156000         | 0.1271600            |
| Signaling by FGFR4                         | 0.1156000         | 0.1277684            |
| FGFR1 mutant receptor activation           | 0.1109760         | 0.1321143            |

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

There are 39 up-regulated and 234 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                              |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                 |  0.0065376|  0.6831081|  1.877221|             1|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                |  0.0065376|  0.4963027|  1.667731|             0|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                        |  0.0065376|  0.7137058|  2.042020|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                  |  0.0065376|  0.5934265|  1.763427|             3|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                            |  0.0065376|  0.5782104|  1.969083|             0|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                  |  0.0065376|  0.6030993|  1.896274|             0|
| G ALPHA (I) SIGNALLING EVENTS%REACTOME%R-HSA-418594.5                                                |  0.0065376|  0.4870838|  1.640676|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                          |  0.0065376|  0.6005742|  1.850388|             0|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                   |  0.0065376|  0.6404478|  1.980675|             0|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                        |  0.0065376|  0.6221719|  1.937707|             0|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                           |  0.0065376|  0.6815479|  1.912033|             0|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                           |  0.0065376|  0.6394826|  1.822284|             1|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090 |  0.0065376|  0.6869231|  1.957472|             0|
| PLASMA LIPOPROTEIN REMODELING%REACTOME DATABASE ID RELEASE 65%8963899                                |  0.0065376|  0.8505300|  1.995750|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                          |  0.0065376|  0.6356000|  1.979528|             0|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                      |  0.0070447|  0.5905612|  1.745242|             7|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                             |  0.0073830|  0.5952607|  1.743306|             8|
| CELL JUNCTION ORGANIZATION%REACTOME DATABASE ID RELEASE 65%446728                                    |  0.0089511|  0.5570686|  1.687657|            12|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                      |  0.0096059|  0.6478909|  1.761470|            13|
| STRIATED MUSCLE CONTRACTION%REACTOME DATABASE ID RELEASE 65%390522                                   |  0.0104036|  0.6993129|  1.759538|            14|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                         |  0.0121773|  0.5155537|  1.610787|            22|
| VISUAL PHOTOTRANSDUCTION%REACTOME DATABASE ID RELEASE 65%2187338                                     |  0.0123446|  0.5728466|  1.682139|            22|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                                 |  0.0128942|  0.6707494|  1.749022|            21|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                      |  0.0131353|  0.5332283|  1.640395|            25|
| MUSCLE CONTRACTION%REACTOME%R-HSA-397014.2                                                           |  0.0135897|  0.4738146|  1.550432|            28|
| RETINOID METABOLISM AND TRANSPORT%REACTOME%R-HSA-975634.2                                            |  0.0208907|  0.6348848|  1.679286|            39|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                               |  0.0209039|  0.5127598|  1.582592|            45|
| DEATH RECEPTOR SIGNALLING%REACTOME%R-HSA-73887.3                                                     |  0.0209768|  0.4690446|  1.523459|            48|
| CELL-CELL COMMUNICATION%REACTOME%R-HSA-1500931.3                                                     |  0.0224254|  0.4939098|  1.559381|            50|
| NON-INTEGRIN MEMBRANE-ECM INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000171                       |  0.0262117|  0.5919711|  1.638923|            54|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                         |  0.0269740|  0.5570028|  1.621303|            59|
| ACTIVATION OF MATRIX METALLOPROTEINASES%REACTOME%R-HSA-1592389.1                                     |  0.0304630|  0.7073713|  1.659832|            57|
| GAP JUNCTION TRAFFICKING AND REGULATION%REACTOME%R-HSA-157858.1                                      |  0.0318867|  0.7055499|  1.655558|            60|
| METABOLISM OF FAT-SOLUBLE VITAMINS%REACTOME DATABASE ID RELEASE 65%6806667                           |  0.0353182|  0.5997627|  1.617145|            75|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                            |  0.0362479|  0.5765545|  1.603183|            79|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                    |  0.0373633|  0.6378608|  1.637377|            77|
| P75 NTR RECEPTOR-MEDIATED SIGNALLING%REACTOME%R-HSA-193704.1                                         |  0.0382492|  0.4838781|  1.514589|            95|
| CELL-CELL JUNCTION ORGANIZATION%REACTOME%R-HSA-421270.4                                              |  0.0382492|  0.5617993|  1.586390|            86|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                           |  0.0424564|  0.5920108|  1.596243|            95|

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
| Kegg neuroactive ligand receptor interaction      | 0.0178500         | 0.0051000            |
| Kegg cytokine cytokine receptor interaction       | 0.0076500         | 0.0076500            |
| Kegg ecm receptor interaction                     | 0.0000000         | 0.0153000            |
| Kegg axon guidance                                | 0.0214200         | 0.0153000            |
| Kegg endocytosis                                  | 0.4109143         | 0.1377000            |
| Kegg intestinal immune network for iga production | 0.0808714         | 0.1591200            |
| Kegg complement and coagulation cascades          | 0.1292000         | 0.1792286            |
| Kegg focal adhesion                               | 0.0102000         | 0.2218500            |
| Kegg basal cell carcinoma                         | 0.1479000         | 0.2261000            |
| Kegg calcium signaling pathway                    | 0.0975375         | 0.2336727            |

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

There are 11 up-regulated and 14 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                         |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_FOCAL\_ADHESION                           |  0.0160562|  0.5099686|  1.712414|             0|
| KEGG\_AXON\_GUIDANCE                            |  0.0167177|  0.5525302|  1.798639|             1|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                |  0.0188656|  0.5948298|  1.799818|             5|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES    |  0.0188656|  0.7010039|  1.822021|             6|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION |  0.0238512|  0.5188667|  1.643097|            17|
| KEGG\_PPAR\_SIGNALING\_PATHWAY                  |  0.0345567|  0.5762875|  1.658438|            30|
| KEGG\_LEUKOCYTE\_TRANSENDOTHELIAL\_MIGRATION    |  0.0345567|  0.5133445|  1.606738|            40|
| KEGG\_HYPERTROPHIC\_CARDIOMYOPATHY\_HCM         |  0.0345567|  0.5495000|  1.637990|            42|
| KEGG\_DILATED\_CARDIOMYOPATHY                   |  0.0345567|  0.5429973|  1.636977|            38|
| KEGG\_LYSOSOME                                  |  0.0465328|  0.4705680|  1.513350|            65|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                  |  0.0494981|  0.4155230|  1.411700|            77|

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
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                | 0.08310           | 0.00000              |
| Peptide GPCRs%WikiPathways\_20180810%WP24%Homo sapiens                             | 0.00000           | 0.00000              |
| Hematopoietic Stem Cell Differentiation%WikiPathways\_20180810%WP2849%Homo sapiens | 0.01385           | 0.02770              |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens      | 0.02770           | 0.04155              |

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

There are 23 up-regulated and 17 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                                                                     |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| GPCRs, Class A Rhodopsin-like%WikiPathways\_20180810%WP455%Homo sapiens                                                                                     |  0.0102708|  0.6165773|  1.822945|             1|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                                                                    |  0.0102708|  0.4972543|  1.668972|             2|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                                                                  |  0.0102708|  0.5971515|  1.803762|             0|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                                                                         |  0.0102708|  0.6663857|  2.176782|             0|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                                                                |  0.0102708|  0.5522147|  1.733318|             2|
| Peptide GPCRs%WikiPathways\_20180810%WP24%Homo sapiens                                                                                                      |  0.0102708|  0.8203847|  1.979421|             0|
| Hypothesized Pathways in Pathogenesis of Cardiovascular Disease%WikiPathways\_20180810%WP3668%Homo sapiens                                                  |  0.0102708|  0.7269469|  1.868881|             2|
| VEGFA-VEGFR2 Signaling Pathway%WikiPathways\_20180810%WP3888%Homo sapiens                                                                                   |  0.0102708|  0.4667619|  1.591127|             4|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                                                                   |  0.0105693|  0.4637730|  1.571802|             5|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                                                                       |  0.0121747|  0.5285761|  1.687050|             8|
| Myometrial Relaxation and Contraction Pathways%WikiPathways\_20180810%WP289%Homo sapiens                                                                    |  0.0132666|  0.5005985|  1.641101|            10|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                                                                      |  0.0143438|  0.5134752|  1.640737|            11|
| Oligodendrocyte Specification and differentiation(including remyelination), leading to Myelin Components for CNS%WikiPathways\_20180810%WP4304%Homo sapiens |  0.0205891|  0.7313979|  1.764714|            13|
| Complement and Coagulation Cascades%WikiPathways\_20180810%WP558%Homo sapiens                                                                               |  0.0234207|  0.6840384|  1.774264|            17|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                                                                   |  0.0237410|  0.4293101|  1.474463|            24|
| Striated Muscle Contraction%WikiPathways\_20180810%WP383%Homo sapiens                                                                                       |  0.0241109|  0.6719445|  1.769603|            19|
| Canonical and Non-Canonical TGF-B signaling%WikiPathways\_20180810%WP3874%Homo sapiens                                                                      |  0.0366388|  0.7362772|  1.731970|            29|
| Vitamin B12 Metabolism%WikiPathways\_20180810%WP1533%Homo sapiens                                                                                           |  0.0369136|  0.6305813|  1.697110|            33|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                                                                                         |  0.0373038|  0.5608112|  1.647834|            37|
| Statin Pathway%WikiPathways\_20180810%WP430%Homo sapiens                                                                                                    |  0.0380636|  0.7328210|  1.723840|            33|
| Selenium Micronutrient Network%WikiPathways\_20180810%WP15%Homo sapiens                                                                                     |  0.0427162|  0.5514807|  1.620418|            45|
| Physiological and Pathological Hypertrophy of the Heart%WikiPathways\_20180810%WP1528%Homo sapiens                                                          |  0.0443787|  0.6669471|  1.677430|            44|
| Photodynamic therapy-induced unfolded protein response%WikiPathways\_20180810%WP3613%Homo sapiens                                                           |  0.0443787|  0.6431187|  1.683134|            44|

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
| G alpha (i) signalling events                  | 0.0000000         | 0.0000000            |
| Extracellular matrix organization              | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)           | 0.0000000         | 0.0000000            |
| SLC-mediated transmembrane transport           | 0.0144500         | 0.0144500            |
| Neuronal System                                | 0.1156000         | 0.0231200            |
| Platelet activation, signaling and aggregation | 0.0825714         | 0.0385333            |
| Collagen formation                             | 0.0963333         | 0.0660571            |
| Post-translational protein phosphorylation     | 0.1517250         | 0.1541333            |
| Signaling by Interleukins                      | 0.9375356         | 0.1661750            |
| Collagen biosynthesis and modifying enzymes    | 0.3146889         | 0.2543200            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 6 enriched Reactome pathways with ermineR

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

There are 30 up-regulated and 149 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                          |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                             |  0.0107666|  0.6906764|  1.899279|             3|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                                     |  0.0107666|  0.6582874|  1.919521|             2|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                              |  0.0107666|  0.6095565|  1.821493|             4|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                        |  0.0107666|  0.5530032|  1.914216|             0|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                                       |  0.0107666|  0.7017027|  1.885387|             5|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                                      |  0.0107666|  0.5937514|  1.835589|             1|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                                    |  0.0107666|  0.5818321|  1.825623|             1|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                  |  0.0107666|  0.7193802|  1.956638|             3|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                       |  0.0107666|  0.5028360|  1.669751|             2|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090             |  0.0107666|  0.6493038|  1.858531|             3|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                              |  0.0136619|  0.5414102|  1.713964|            10|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                                      |  0.0136877|  0.5419548|  1.700500|            10|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                                         |  0.0142855|  0.6007229|  1.762912|            10|
| TRANSPORT OF BILE SALTS AND ORGANIC ACIDS, METAL IONS AND AMINE COMPOUNDS%REACTOME DATABASE ID RELEASE 65%425366 |  0.0152151|  0.6137404|  1.765577|            11|
| INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING%REACTOME DATABASE ID RELEASE 65%6785807                               |  0.0153533|  0.5562357|  1.704690|            12|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                       |  0.0160598|  0.6202001|  1.775226|            12|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                                     |  0.0162402|  0.5324235|  1.695007|            14|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                        |  0.0170541|  0.6339416|  1.760789|            13|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                                    |  0.0192177|  0.6072272|  1.746841|            17|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                               |  0.0192177|  0.5408211|  1.681164|            18|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                                |  0.0204299|  0.7010078|  1.784715|            17|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                                         |  0.0212479|  0.6956923|  1.771182|            18|
| IRE1ALPHA ACTIVATES CHAPERONES%REACTOME DATABASE ID RELEASE 65%381070                                            |  0.0298316|  0.5578257|  1.659253|            33|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                                     |  0.0341014|  0.5141659|  1.620357|            42|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                            |  0.0385249|  0.4337733|  1.480071|            58|
| SIGNALING BY WNT IN CANCER%REACTOME DATABASE ID RELEASE 65%4791275                                               |  0.0445346|  0.6290613|  1.690209|            55|
| XBP1(S) ACTIVATES CHAPERONE GENES%REACTOME%R-HSA-381038.2                                                        |  0.0481070|  0.5462921|  1.614861|            68|
| SIGNALING BY INTERLEUKINS%REACTOME%R-HSA-449147.10                                                               |  0.0482551|  0.3884707|  1.373126|            88|
| SIGNALING BY FGFR1%REACTOME DATABASE ID RELEASE 65%5654736                                                       |  0.0488311|  0.6004543|  1.651179|            68|
| RUNX1 INTERACTS WITH CO-FACTORS WHOSE PRECISE EFFECT ON RUNX1 TARGETS IS NOT KNOWN%REACTOME%R-HSA-8939243.1      |  0.0488311|  0.6012725|  1.644533|            68|

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
| Kegg focal adhesion                          | 0.0153000         | 0.0000000            |
| Kegg pathways in cancer                      | 0.1748571         | 0.0051000            |
| Kegg cytokine cytokine receptor interaction  | 0.0153000         | 0.0076500            |
| Kegg ecm receptor interaction                | 0.0306000         | 0.0153000            |
| Kegg calcium signaling pathway               | 0.0336600         | 0.0191250            |
| Kegg neuroactive ligand receptor interaction | 0.0229500         | 0.0561000            |
| Kegg tgf beta signaling pathway              | 0.1096500         | 0.0786857            |
| Kegg wnt signaling pathway                   | 0.1361700         | 0.0841500            |
| Kegg mapk signaling pathway                  | 0.1440750         | 0.0867000            |
| Kegg complement and coagulation cascades     | 0.1513000         | 0.1210091            |

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

There are 15 up-regulated and 11 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                         |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION |  0.0047170|  0.6292265|  2.010945|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY             |  0.0047170|  0.6570327|  2.021503|             0|
| KEGG\_AXON\_GUIDANCE                            |  0.0047170|  0.5635419|  1.856192|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                      |  0.0047170|  0.4809057|  1.694392|             0|
| KEGG\_FOCAL\_ADHESION                           |  0.0066995|  0.4930680|  1.677180|             1|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                |  0.0106986|  0.5720013|  1.733687|             4|
| KEGG\_SMALL\_CELL\_LUNG\_CANCER                 |  0.0150967|  0.5395996|  1.690950|             8|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                  |  0.0250829|  0.4328452|  1.499178|            22|
| KEGG\_WNT\_SIGNALING\_PATHWAY                   |  0.0281984|  0.4734545|  1.562072|            25|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES    |  0.0295630|  0.6656616|  1.732482|            22|
| KEGG\_REGULATION\_OF\_ACTIN\_CYTOSKELETON       |  0.0357221|  0.4397563|  1.493791|            38|
| KEGG\_BASAL\_CELL\_CARCINOMA                    |  0.0357221|  0.6020553|  1.679916|            32|
| KEGG\_PANCREATIC\_CANCER                        |  0.0423282|  0.5231049|  1.596391|            50|
| KEGG\_CHRONIC\_MYELOID\_LEUKEMIA                |  0.0456493|  0.5090933|  1.575532|            60|
| KEGG\_JAK\_STAT\_SIGNALING\_PATHWAY             |  0.0482578|  0.4937263|  1.547197|            67|

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
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens           | 0.0138500         | 0.000000             |
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens           | 0.0184667         | 0.000000             |
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens      | 0.0000000         | 0.009233             |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens | 0.0277000         | 0.027700             |

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

There are 39 up-regulated and 13 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                                                                     |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                                                                    |  0.0048549|  0.4896696|  1.665172|             1|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                                                                  |  0.0048549|  0.6396532|  1.939548|             1|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                                                                          |  0.0048549|  0.8160812|  2.155541|             0|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                                                                   |  0.0048549|  0.5131921|  1.766858|             0|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                                                                      |  0.0048549|  0.5985506|  1.925896|             1|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                                                                        |  0.0048549|  0.5776968|  1.847032|             1|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                                                                    |  0.0048549|  0.5754341|  1.917070|             0|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                                                                         |  0.0048549|  0.6917237|  2.285408|             0|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                                                                         |  0.0048549|  0.5238820|  1.735133|             1|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                                                                |  0.0048549|  0.7991357|  2.132811|             1|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                                                               |  0.0048549|  0.7183660|  1.974883|             1|
| Canonical and Non-Canonical TGF-B signaling%WikiPathways\_20180810%WP3874%Homo sapiens                                                                      |  0.0048549|  0.8359753|  1.949959|             0|
| Peptide GPCRs%WikiPathways\_20180810%WP24%Homo sapiens                                                                                                      |  0.0048549|  0.8320267|  1.999267|             0|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                                                                       |  0.0048549|  0.5612577|  1.802527|             1|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                                                                   |  0.0048549|  0.4962045|  1.735540|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens                                                        |  0.0048549|  0.5683104|  1.878165|             0|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                                                                                       |  0.0051401|  0.4586217|  1.606694|             2|
| Spinal Cord Injury%WikiPathways\_20180810%WP2431%Homo sapiens                                                                                               |  0.0055939|  0.5592672|  1.749923|             2|
| GPCRs, Class A Rhodopsin-like%WikiPathways\_20180810%WP455%Homo sapiens                                                                                     |  0.0078716|  0.6254128|  1.853147|             5|
| Oncostatin M Signaling Pathway%WikiPathways\_20180810%WP2374%Homo sapiens                                                                                   |  0.0078716|  0.5865768|  1.768154|             5|
| TGF-B Signaling in Thyroid Cells for Epithelial-Mesenchymal Transition%WikiPathways\_20180810%WP3859%Homo sapiens                                           |  0.0102549|  0.7828937|  1.826144|             6|
| TGF-beta Receptor Signaling%WikiPathways\_20180810%WP560%Homo sapiens                                                                                       |  0.0114616|  0.5943945|  1.756724|             9|
| VEGFA-VEGFR2 Signaling Pathway%WikiPathways\_20180810%WP3888%Homo sapiens                                                                                   |  0.0114616|  0.4432686|  1.537589|            11|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                                                             |  0.0158957|  0.5165918|  1.651665|            15|
| NRF2 pathway%WikiPathways\_20180810%WP2884%Homo sapiens                                                                                                     |  0.0160753|  0.5310219|  1.663525|            16|
| Simplified Interaction Map Between LOXL4 and Oxidative Stress Pathway%WikiPathways\_20180810%WP3670%Homo sapiens                                            |  0.0160753|  0.7641441|  1.782409|            13|
| MAPK Signaling Pathway%WikiPathways\_20180810%WP382%Homo sapiens                                                                                            |  0.0179910|  0.4375281|  1.508321|            22|
| Dopaminergic Neurogenesis%WikiPathways\_20180810%WP2855%Homo sapiens                                                                                        |  0.0200221|  0.7356367|  1.792640|            18|
| Hypothesized Pathways in Pathogenesis of Cardiovascular Disease%WikiPathways\_20180810%WP3668%Homo sapiens                                                  |  0.0241775|  0.6886972|  1.769343|            25|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens                                              |  0.0241775|  0.5289301|  1.642538|            28|
| Chromosomal and microsatellite instability in colorectal cancer %WikiPathways\_20180810%WP4216%Homo sapiens                                                 |  0.0241775|  0.5240649|  1.622813|            29|
| Oligodendrocyte Specification and differentiation(including remyelination), leading to Myelin Components for CNS%WikiPathways\_20180810%WP4304%Homo sapiens |  0.0241775|  0.7303100|  1.754853|            24|
| ErbB Signaling Pathway%WikiPathways\_20180810%WP673%Homo sapiens                                                                                            |  0.0270120|  0.5692805|  1.670941|            33|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                                                                  |  0.0309557|  0.6932492|  1.728234|            35|
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens                                                     |  0.0333104|  0.5597345|  1.654287|            44|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                                                                |  0.0333133|  0.4990197|  1.577426|            48|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                                                                                         |  0.0348786|  0.5593290|  1.648068|            48|
| Aryl Hydrocarbon Receptor Pathway%WikiPathways\_20180810%WP2873%Homo sapiens                                                                                |  0.0364106|  0.6542665|  1.698411|            46|
| Thymic Stromal LymphoPoietin (TSLP) Signaling Pathway%WikiPathways\_20180810%WP2203%Homo sapiens                                                            |  0.0431330|  0.5884621|  1.643068|            59|
