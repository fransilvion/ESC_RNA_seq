Pathway enrichment analysis, GEO109658
================
German Novakovskiy
August 17, 2018

Analysis of first day of differentiation 0h VS 24 h (formation of anterior Primitive streak, mesodendoderm)
-----------------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE109658/DEgenes_0h_24h_109658.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
#absolute Log fold changes for ermineR
ermineInputGeneScores <- DEgenes_0h_24h_109658 %>% 
  rownames_to_column("gene") %>%
  mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

#exact log fold changes for fgsea
# scores forfgsea
scoresFGSEADF <- DEgenes_0h_24h_109658 %>% 
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

    ## Warning in bitr(rownames(ermineInputGeneScores), fromType = "SYMBOL",
    ## toType = "ENTREZID", : 8.91% of input gene IDs are fail to map...

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

| Pathway                                                                              | Corrected p-value | Corrected MF p-value |
|:-------------------------------------------------------------------------------------|:------------------|:---------------------|
| Neuronal System                                                                      | 0.0000000         | 0.0000000            |
| G alpha (i) signalling events                                                        | 0.0000000         | 0.0000000            |
| Extracellular matrix organization                                                    | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)                                                 | 0.0000000         | 0.0000000            |
| TCF dependent signaling in response to WNT                                           | 0.0000000         | 0.0122800            |
| Class B 2 (Secretin family receptors)                                                | 0.0307000         | 0.0409333            |
| Activation of anterior HOX genes in hindbrain development during early embryogenesis | 0.1052571         | 0.1315714            |
| Activation of HOX genes during differentiation                                       | 0.1052571         | 0.1315714            |
| WNT ligand biogenesis and trafficking                                                | 0.1432667         | 0.1611750            |
| Incretin synthesis, secretion, and inactivation                                      | 0.2701600         | 0.3206444            |

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

There are 12 up-regulated and 26 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                       |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                     |  0.0248403|  0.4847848|  1.780218|             0|
| ACTIVATION OF HOX GENES DURING DIFFERENTIATION%REACTOME DATABASE ID RELEASE 65%5619507                        |  0.0248403|  0.5440262|  1.742851|             2|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                                |  0.0248403|  0.9108422|  2.195390|             0|
| ACTIVATION OF ANTERIOR HOX GENES IN HINDBRAIN DEVELOPMENT DURING EARLY EMBRYOGENESIS%REACTOME%R-HSA-5617472.2 |  0.0248403|  0.5440262|  1.742851|             2|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                             |  0.0248403|  0.7189226|  1.908512|             1|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                    |  0.0248403|  0.4782794|  1.706794|             1|
| PEPTIDE HORMONE METABOLISM%REACTOME DATABASE ID RELEASE 65%2980736                                            |  0.0248403|  0.5758822|  1.773739|             2|
| SIGNALING BY NODAL%REACTOME DATABASE ID RELEASE 65%1181150                                                    |  0.0281331|  0.7644339|  1.842504|             3|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                               |  0.0281331|  0.6777577|  1.923094|             3|
| INCRETIN SYNTHESIS, SECRETION, AND INACTIVATION%REACTOME%R-HSA-400508.1                                       |  0.0353002|  0.7653867|  1.792000|             6|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                          |  0.0445294|  0.6162814|  1.761959|            10|
| CLASS B 2 (SECRETIN FAMILY RECEPTORS)%REACTOME DATABASE ID RELEASE 65%373080                                  |  0.0445294|  0.5761229|  1.748651|            10|

``` r
upPathwaysReactome109658_024 <- upPathwaysReactome
save(upPathwaysReactome109658_024, file="upPathwaysReactome109658_024.Rdata")
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

| Pathway                                                   | Corrected p-value | Corrected MF p-value |
|:----------------------------------------------------------|:------------------|:---------------------|
| Kegg wnt signaling pathway                                | 0.0000000         | 0.000000             |
| Kegg pathways in cancer                                   | 0.0000000         | 0.000000             |
| Kegg neuroactive ligand receptor interaction              | 0.0000000         | 0.000000             |
| Kegg melanogenesis                                        | 0.0027330         | 0.000000             |
| Kegg focal adhesion                                       | 0.1164400         | 0.000000             |
| Kegg cytokine cytokine receptor interaction               | 0.0000000         | 0.000000             |
| Kegg basal cell carcinoma                                 | 0.0000000         | 0.000000             |
| Kegg hedgehog signaling pathway                           | 0.0046860         | 0.002050             |
| Kegg mapk signaling pathway                               | 0.1020444         | 0.007289             |
| Kegg arrhythmogenic right ventricular cardiomyopathy arvc | 0.2121067         | 0.008200             |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 22 enriched KEGG pathways with ermineR.

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

There are 7 up-regulated and 1 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                             |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY  |  0.0062671|  0.6907024|  2.019305|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY |  0.0062671|  0.6045600|  1.940383|             0|
| KEGG\_PATHWAYS\_IN\_CANCER          |  0.0062671|  0.4536458|  1.699764|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA        |  0.0062671|  0.6672720|  1.966077|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY       |  0.0086050|  0.4915260|  1.708527|             1|
| KEGG\_ECM\_RECEPTOR\_INTERACTION    |  0.0271482|  0.5443352|  1.727487|             7|
| KEGG\_MELANOGENESIS                 |  0.0290084|  0.5240401|  1.697418|             9|

``` r
upPathwaysKEGG109658_024 <- upPathwaysKEGG
save(upPathwaysKEGG109658_024, file="upPathwaysKEGG109658_024.Rdata")
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
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                | 0.0000000         | 0.0000000            |
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens | 0.0000000         | 0.0000000            |
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       | 0.0000000         | 0.0000000            |
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            | 0.0000000         | 0.0000000            |
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                             | 0.0000000         | 0.0000000            |
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens        | 0.0000000         | 0.0000000            |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  | 0.0000000         | 0.0000000            |
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                        | 0.0000000         | 0.0000000            |
| GPCRs, Class A Rhodopsin-like%WikiPathways\_20180810%WP455%Homo sapiens                                        | 0.0107636         | 0.0164444            |
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                            | 0.0318769         | 0.0177600            |
| Wnt Signaling Pathway and Pluripotency%WikiPathways\_20180810%WP399%Homo sapiens                               | 0.0059200         | 0.0242182            |

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 11 enriched Wiki pathways with ermineR

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

There are 17 up-regulated and 4 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                           |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                           |  0.0058039|  0.6050129|  1.967408|             0|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                              |  0.0058039|  0.5598027|  1.862355|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                          |  0.0058039|  0.5507020|  1.915016|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                      |  0.0058039|  0.7626768|  2.152499|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                     |  0.0058039|  0.6719177|  1.968151|             0|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens    |  0.0058039|  0.5696587|  1.852442|             0|
| Breast cancer pathway%WikiPathways\_20180810%WP4262%Homo sapiens                                                  |  0.0058039|  0.5349072|  1.840416|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens              |  0.0058039|  0.5702915|  1.978250|             0|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                            |  0.0081341|  0.5251821|  1.767138|             1|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                                |  0.0116554|  0.6859226|  1.899441|             2|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                   |  0.0187762|  0.5027597|  1.678039|             5|
| TGF-B Signaling in Thyroid Cells for Epithelial-Mesenchymal Transition%WikiPathways\_20180810%WP3859%Homo sapiens |  0.0232744|  0.7757247|  1.871421|             6|
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens           |  0.0305984|  0.5482734|  1.698407|            11|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                               |  0.0305984|  0.4659379|  1.606639|            13|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                        |  0.0305984|  0.6835662|  1.845111|            10|
| T-Cell antigen Receptor (TCR) Signaling Pathway%WikiPathways\_20180810%WP69%Homo sapiens                          |  0.0469841|  0.5110336|  1.627143|            21|
| DNA Damage Response (only ATM dependent)%WikiPathways\_20180810%WP710%Homo sapiens                                |  0.0487455|  0.4841329|  1.607717|            24|

``` r
upPathwaysWiki109658_024 <- upPathwaysWiki
save(upPathwaysWiki109658_024, file="upPathwaysWiki109658_024.Rdata")
```

Analysis of second-fourth days of differentiation 24h VS 96h (formation of definitive endoderm from APS)
--------------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE109658/DEgenes_24h_96h_109658.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
#absolute Log fold changes for ermineR
ermineInputGeneScores <- DEgenes_24h_96h_109658 %>% 
  rownames_to_column("gene") %>%
  mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

#exact log fold changes for fgsea
# scores forfgsea
scoresFGSEADF <- DEgenes_24h_96h_109658 %>% 
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

    ## Warning in bitr(rownames(ermineInputGeneScores), fromType = "SYMBOL",
    ## toType = "ENTREZID", : 8.91% of input gene IDs are fail to map...

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

| Pathway                                                                              | Corrected p-value | Corrected MF p-value |
|:-------------------------------------------------------------------------------------|:------------------|:---------------------|
| Extracellular matrix organization                                                    | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)                                                 | 0.0000000         | 0.0000000            |
| Activation of anterior HOX genes in hindbrain development during early embryogenesis | 0.0000000         | 0.0000000            |
| Activation of HOX genes during differentiation                                       | 0.0000000         | 0.0000000            |
| Neuronal System                                                                      | 0.0122800         | 0.0122800            |
| G alpha (i) signalling events                                                        | 0.0153500         | 0.0153500            |
| Class B 2 (Secretin family receptors)                                                | 0.0307000         | 0.0409333            |
| Incretin synthesis, secretion, and inactivation                                      | 0.1978444         | 0.2543714            |
| SLC-mediated transmembrane transport                                                 | 0.2072250         | 0.2609500            |
| Diseases associated with O-glycosylation of proteins                                 | 0.2795660         | 0.3003261            |

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

There are 30 up-regulated and 191 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                             |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                |  0.0077130|  0.6303752|  1.851646|             1|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244           |  0.0077130|  0.5005085|  1.900650|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                           |  0.0077130|  0.6312265|  1.903997|             0|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275  |  0.0080451|  0.5369792|  1.788835|             2|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814 |  0.0086329|  0.5624466|  1.775616|             4|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                         |  0.0096020|  0.5127644|  1.728042|             6|
| G ALPHA (I) SIGNALLING EVENTS%REACTOME%R-HSA-418594.5                               |  0.0096917|  0.4089427|  1.534720|             7|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                       |  0.0100657|  0.5745773|  1.775497|             6|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                         |  0.0113304|  0.5184447|  1.702125|             8|
| DOPAMINE NEUROTRANSMITTER RELEASE CYCLE%REACTOME%R-HSA-212676.2                     |  0.0147255|  0.7125226|  1.809653|            10|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                       |  0.0151530|  0.4834313|  1.639343|            13|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                     |  0.0159333|  0.5962055|  1.733138|            13|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                          |  0.0196947|  0.5511119|  1.683328|            18|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                          |  0.0201601|  0.5342714|  1.668401|            19|
| OTHER SEMAPHORIN INTERACTIONS%REACTOME DATABASE ID RELEASE 65%416700                |  0.0233367|  0.7196216|  1.753019|            21|
| SURFACTANT METABOLISM%REACTOME DATABASE ID RELEASE 65%5683826                       |  0.0241164|  0.7092872|  1.753009|            23|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2               |  0.0256653|  0.3932485|  1.472504|            34|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                 |  0.0275053|  0.4605667|  1.573127|            34|
| REGULATION OF COMPLEMENT CASCADE%REACTOME%R-HSA-977606.5                            |  0.0288625|  0.7146985|  1.718085|            30|
| COMPLEMENT CASCADE%REACTOME DATABASE ID RELEASE 65%166658                           |  0.0315265|  0.6993303|  1.703589|            34|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                |  0.0317112|  0.6136479|  1.694497|            36|
| ION CHANNEL TRANSPORT%REACTOME%R-HSA-983712.2                                       |  0.0325573|  0.4166411|  1.489823|            45|
| EFFECTS OF PIP2 HYDROLYSIS%REACTOME%R-HSA-114508.2                                  |  0.0325573|  0.6320716|  1.702772|            37|
| KERATAN SULFATE KERATIN METABOLISM%REACTOME%R-HSA-1638074.1                         |  0.0388389|  0.6146049|  1.667571|            48|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                            |  0.0388389|  0.6159588|  1.671244|            47|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                       |  0.0388389|  0.6142315|  1.666557|            48|
| PLASMA LIPOPROTEIN REMODELING%REACTOME DATABASE ID RELEASE 65%8963899               |  0.0412998|  0.6633448|  1.660677|            52|
| CELL-CELL JUNCTION ORGANIZATION%REACTOME%R-HSA-421270.4                             |  0.0416325|  0.5041505|  1.580054|            60|
| MUSCLE CONTRACTION%REACTOME%R-HSA-397014.2                                          |  0.0458723|  0.3955380|  1.434875|            77|
| DISEASES OF GLYCOSYLATION%REACTOME DATABASE ID RELEASE 65%3781865                   |  0.0487045|  0.4851009|  1.552078|            76|

``` r
upPathwaysReactome109658_2496 <- upPathwaysReactome
save(upPathwaysReactome109658_2496, file="upPathwaysReactome109658_2496.Rdata")
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
| Kegg pathways in cancer                      | 0.0164000         | 0.0000000            |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.0000000            |
| Kegg cytokine cytokine receptor interaction  | 0.0000000         | 0.0000000            |
| Kegg wnt signaling pathway                   | 0.0273333         | 0.0041000            |
| Kegg hematopoietic cell lineage              | 0.0246000         | 0.0131200            |
| Kegg cell adhesion molecules cams            | 0.0262400         | 0.0164000            |
| Kegg ecm receptor interaction                | 0.0492000         | 0.0351429            |
| Kegg melanogenesis                           | 0.1822222         | 0.0451000            |
| Kegg focal adhesion                          | 0.4829800         | 0.0730546            |
| Kegg maturity onset diabetes of the young    | 0.0861000         | 0.0765333            |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 8 enriched KEGG pathways with ermineR.

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

There are 6 up-regulated and 4 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                                             |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_AXON\_GUIDANCE                                                |  0.0112782|  0.4619294|  1.646794|             2|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                                    |  0.0112782|  0.5578373|  1.840238|             1|
| KEGG\_LEUKOCYTE\_TRANSENDOTHELIAL\_MIGRATION                        |  0.0112782|  0.5034342|  1.734094|             1|
| KEGG\_CELL\_ADHESION\_MOLECULES\_CAMS                               |  0.0201064|  0.4921250|  1.695139|             5|
| KEGG\_GLYCOSPHINGOLIPID\_BIOSYNTHESIS\_LACTO\_AND\_NEOLACTO\_SERIES |  0.0247335|  0.6982301|  1.815233|             8|
| KEGG\_FOCAL\_ADHESION                                               |  0.0247335|  0.4212441|  1.563146|            10|

``` r
upPathwaysKEGG109658_2496 <- upPathwaysKEGG
save(upPathwaysKEGG109658_2496, file="upPathwaysKEGG109658_2496.Rdata")
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

| Pathway                                                                 | Corrected p-value | Corrected MF p-value |
|:------------------------------------------------------------------------|:------------------|:---------------------|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens     | 0.0740000         | 0.0444000            |
| GPCRs, Class A Rhodopsin-like%WikiPathways\_20180810%WP455%Homo sapiens | 0.0493333         | 0.0493333            |

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 2 enriched Wiki pathways with ermineR

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

There are 4 up-regulated and 8 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                       |  0.0311153|  0.6217300|  2.202875|             0|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens |  0.0393722|  0.3960488|  1.516293|             6|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                     |  0.0393722|  0.3972533|  1.529666|             4|
| Phosphodiesterases in neuronal function%WikiPathways\_20180810%WP4222%Homo sapiens        |  0.0396987|  0.5983603|  1.786129|             7|

``` r
upPathwaysWiki109658_2496 <- upPathwaysWiki
save(upPathwaysWiki109658_2496, file="upPathwaysWiki109658_2496.Rdata")
```

Analysis of whole differentiation process 0h VS 96 h (formation of Definitive endoderm from hESC)
-------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE109658/DEgenes_0h_96h_109658.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
#absolute Log fold changes for ermineR
ermineInputGeneScores <- DEgenes_0h_96h_109658 %>% 
  rownames_to_column("gene") %>%
  mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

#exact log fold changes for fgsea
# scores forfgsea
scoresFGSEADF <- DEgenes_0h_96h_109658 %>% 
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

    ## Warning in bitr(rownames(ermineInputGeneScores), fromType = "SYMBOL",
    ## toType = "ENTREZID", : 8.91% of input gene IDs are fail to map...

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

| Pathway                                         | Corrected p-value | Corrected MF p-value |
|:------------------------------------------------|:------------------|:---------------------|
| Incretin synthesis, secretion, and inactivation | 0.0000000         | 0.0000000            |
| Extracellular matrix organization               | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)            | 0.0000000         | 0.0000000            |
| G alpha (i) signalling events                   | 0.0000000         | 0.0153500            |
| Neuronal System                                 | 0.0122800         | 0.1228000            |
| G alpha (q) signalling events                   | 0.1842000         | 0.2105143            |
| Amine ligand-binding receptors                  | 0.2017429         | 0.2353667            |
| Class B 2 (Secretin family receptors)           | 0.3795636         | 0.3963091            |
| Other interleukin signaling                     | 0.4009059         | 0.4045176            |
| Signaling by BMP                                | 0.4042167         | 0.4144500            |

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

There are 34 up-regulated and 221 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                         |       padj|         ES|       NES|  nMoreExtreme|
|:----------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                            |  0.0049687|  0.6778319|  2.006600|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                             |  0.0049687|  0.6026327|  1.896839|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                       |  0.0049687|  0.5498580|  2.047014|             0|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                                  |  0.0049687|  0.8103541|  2.027593|             0|
| CLASS B 2 (SECRETIN FAMILY RECEPTORS)%REACTOME DATABASE ID RELEASE 65%373080                                    |  0.0049687|  0.6183446|  1.935982|             0|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                                        |  0.0049687|  0.6841620|  1.876828|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                                     |  0.0049687|  0.5818259|  1.902946|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                       |  0.0049687|  0.6573691|  1.995181|             0|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                                   |  0.0049687|  0.6976477|  1.913823|             0|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                 |  0.0049687|  0.7379812|  2.166014|             0|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                      |  0.0049687|  0.6334096|  1.979759|             0|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                               |  0.0049687|  0.7410476|  2.032880|             0|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                      |  0.0049687|  0.4563494|  1.661867|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090            |  0.0049687|  0.6042254|  1.859026|             3|
| KERATAN SULFATE KERATIN METABOLISM%REACTOME%R-HSA-1638074.1                                                     |  0.0056314|  0.6564317|  1.800757|             4|
| KERATAN SULFATE BIOSYNTHESIS%REACTOME%R-HSA-2022854.1                                                           |  0.0056314|  0.6994094|  1.838213|             4|
| PEPTIDE HORMONE METABOLISM%REACTOME DATABASE ID RELEASE 65%2980736                                              |  0.0131458|  0.5259866|  1.669968|            18|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                                          |  0.0140804|  0.4788939|  1.615572|            21|
| SIGNALING BY NODAL%REACTOME DATABASE ID RELEASE 65%1181150                                                      |  0.0161681|  0.6946138|  1.737998|            20|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                                   |  0.0176236|  0.6232640|  1.719890|            24|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                                   |  0.0200507|  0.4696021|  1.587591|            34|
| INCRETIN SYNTHESIS, SECRETION, AND INACTIVATION%REACTOME%R-HSA-400508.1                                         |  0.0211855|  0.7041486|  1.706743|            29|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                              |  0.0220606|  0.4785204|  1.588952|            38|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                                     |  0.0260329|  0.4638452|  1.559026|            47|
| SURFACTANT METABOLISM%REACTOME DATABASE ID RELEASE 65%5683826                                                   |  0.0272277|  0.6749192|  1.688720|            41|
| INTERLEUKIN-6 FAMILY SIGNALING%REACTOME DATABASE ID RELEASE 65%6783589                                          |  0.0295148|  0.6508816|  1.693446|            46|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                               |  0.0301882|  0.6499894|  1.691125|            48|
| SIGNALING BY PDGF%REACTOME DATABASE ID RELEASE 65%186797                                                        |  0.0315933|  0.5212166|  1.603633|            57|
| IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON-LYMPHOID CELL%REACTOME DATABASE ID RELEASE 65%198933 |  0.0323126|  0.5408741|  1.616068|            57|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                                 |  0.0383916|  0.4603672|  1.521909|            77|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                                    |  0.0410993|  0.4381799|  1.501316|            87|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                                         |  0.0416903|  0.6753657|  1.636978|            69|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                                      |  0.0468029|  0.5103004|  1.564948|            90|
| ACTIVATION OF MATRIX METALLOPROTEINASES%REACTOME%R-HSA-1592389.1                                                |  0.0469925|  0.6523676|  1.632294|            80|

``` r
upPathwaysReactome109658_096 <- upPathwaysReactome
save(upPathwaysReactome109658_096, file="upPathwaysReactome109658_096.Rdata")
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
| Kegg pathways in cancer                      | 0.0304571         | 0.0000000            |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.0000000            |
| Kegg focal adhesion                          | 0.0355333         | 0.0000000            |
| Kegg cytokine cytokine receptor interaction  | 0.0000000         | 0.0000000            |
| Kegg ecm receptor interaction                | 0.0273333         | 0.0098400            |
| Kegg hematopoietic cell lineage              | 0.0287000         | 0.0164000            |
| Kegg cell adhesion molecules cams            | 0.0393600         | 0.0164000            |
| Kegg tgf beta signaling pathway              | 0.0307500         | 0.0205000            |
| Kegg calcium signaling pathway               | 0.2427200         | 0.0874667            |
| Kegg systemic lupus erythematosus            | 0.2186667         | 0.1872333            |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 8 enriched KEGG pathways with ermineR.

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

There are 10 up-regulated and 7 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                         |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY              |  0.0032638|  0.6493254|  1.948780|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY             |  0.0032638|  0.5752288|  1.902979|             0|
| KEGG\_AXON\_GUIDANCE                            |  0.0032638|  0.5025458|  1.769772|             0|
| KEGG\_FOCAL\_ADHESION                           |  0.0032638|  0.4824189|  1.764018|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                |  0.0032638|  0.6018477|  1.967144|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                      |  0.0032638|  0.4193480|  1.589378|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA                    |  0.0032638|  0.6480634|  1.965679|             0|
| KEGG\_HEMATOPOIETIC\_CELL\_LINEAGE              |  0.0056717|  0.6173100|  1.863430|             1|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION |  0.0168991|  0.4471618|  1.574731|            12|
| KEGG\_WNT\_SIGNALING\_PATHWAY                   |  0.0241539|  0.4376114|  1.553223|            19|

``` r
upPathwaysKEGG109658_096 <- upPathwaysKEGG
save(upPathwaysKEGG109658_096, file="upPathwaysKEGG109658_096.Rdata")
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

| Pathway                                                             | Corrected p-value | Corrected MF p-value |
|:--------------------------------------------------------------------|:------------------|:---------------------|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens | 0.1036            | 0.0493333            |

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 1 enriched Wiki pathways with ermineR

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

There are 25 up-regulated and 13 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                        |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                     |  0.0050250|  0.5738558|  1.850280|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       |  0.0050250|  0.5103239|  1.827561|             0|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                            |  0.0050250|  0.5535169|  1.950735|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                   |  0.0050250|  0.7855874|  2.313349|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  |  0.0050250|  0.6736308|  2.055530|             0|
| Peptide GPCRs%WikiPathways\_20180810%WP24%Homo sapiens                                                         |  0.0050250|  0.7201486|  1.907934|             0|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                      |  0.0050250|  0.4274579|  1.615726|             0|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                                          |  0.0050250|  0.4303336|  1.632475|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens           |  0.0050250|  0.5806698|  2.075207|             0|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                       |  0.0067135|  0.4637370|  1.701601|             1|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                         |  0.0070049|  0.4970351|  1.723048|             2|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                           |  0.0070049|  0.5284407|  1.812161|             2|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                     |  0.0070049|  0.6828241|  1.922083|             1|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens |  0.0070049|  0.5188105|  1.741162|             2|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                   |  0.0092535|  0.5067811|  1.712171|             4|
| Hematopoietic Stem Cell Differentiation%WikiPathways\_20180810%WP2849%Homo sapiens                             |  0.0098872|  0.6112717|  1.831103|             4|
| TGF-beta Receptor Signaling%WikiPathways\_20180810%WP560%Homo sapiens                                          |  0.0109401|  0.5792501|  1.802880|             5|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                      |  0.0129867|  0.4031603|  1.509313|             9|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                |  0.0168665|  0.4742329|  1.630994|            12|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                             |  0.0176234|  0.6106621|  1.763162|            11|
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                        |  0.0253081|  0.4812308|  1.615042|            20|
| miRNA targets in ECM and membrane receptors%WikiPathways\_20180810%WP2911%Homo sapiens                         |  0.0315207|  0.6587497|  1.745266|            22|
| BMP Signaling Pathway in Eyelid Development%WikiPathways\_20180810%WP3927%Homo sapiens                         |  0.0326369|  0.7084519|  1.776322|            23|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                          |  0.0344858|  0.4408707|  1.543211|            33|
| Aryl Hydrocarbon Receptor Pathway%WikiPathways\_20180810%WP2873%Homo sapiens                                   |  0.0431487|  0.6007809|  1.691140|            37|

``` r
upPathwaysWiki109658_096 <- upPathwaysWiki
save(upPathwaysWiki109658_096, file="upPathwaysWiki109658_096.Rdata")
```
