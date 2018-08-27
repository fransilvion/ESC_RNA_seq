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
| Class B 2 (Secretin family receptors)                                                | 0.0263143         | 0.0409333            |
| Activation of anterior HOX genes in hindbrain development during early embryogenesis | 0.0307000         | 0.0701714            |
| Activation of HOX genes during differentiation                                       | 0.0307000         | 0.0701714            |
| WNT ligand biogenesis and trafficking                                                | 0.1432667         | 0.1611750            |
| Incretin synthesis, secretion, and inactivation                                      | 0.3315600         | 0.3888667            |

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

There are 13 up-regulated and 26 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                         |       padj|         ES|       NES|  nMoreExtreme|
|:----------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                       |  0.0194656|  0.4847848|  1.776369|             0|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                                  |  0.0194656|  0.9108422|  2.186788|             0|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                 |  0.0194656|  0.6777577|  1.920416|             0|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                               |  0.0194656|  0.7189226|  1.888199|             1|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                      |  0.0194656|  0.4782794|  1.704356|             1|
| ACTIVATION OF HOX GENES DURING DIFFERENTIATION%REACTOME DATABASE ID RELEASE 65%5619507                          |  0.0248743|  0.5440262|  1.744191|             4|
| ACTIVATION OF ANTERIOR HOX GENES IN HINDBRAIN DEVELOPMENT DURING EARLY EMBRYOGENESIS%REACTOME%R-HSA-5617472.2   |  0.0248743|  0.5440262|  1.744191|             4|
| PEPTIDE HORMONE METABOLISM%REACTOME DATABASE ID RELEASE 65%2980736                                              |  0.0248743|  0.5758822|  1.774636|             4|
| CLASS B 2 (SECRETIN FAMILY RECEPTORS)%REACTOME DATABASE ID RELEASE 65%373080                                    |  0.0289047|  0.5761229|  1.747605|             5|
| SIGNALING BY NODAL%REACTOME DATABASE ID RELEASE 65%1181150                                                      |  0.0320405|  0.7644339|  1.835285|             5|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                            |  0.0386402|  0.6162814|  1.762655|             8|
| IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON-LYMPHOID CELL%REACTOME DATABASE ID RELEASE 65%198933 |  0.0494146|  0.6037601|  1.748189|            12|
| INCRETIN SYNTHESIS, SECRETION, AND INACTIVATION%REACTOME%R-HSA-400508.1                                         |  0.0496362|  0.7653867|  1.787272|            11|

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
| Kegg wnt signaling pathway                   | 0.0000000         | 0.000000             |
| Kegg pathways in cancer                      | 0.0000000         | 0.000000             |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.000000             |
| Kegg melanogenesis                           | 0.0027330         | 0.000000             |
| Kegg focal adhesion                          | 0.1180800         | 0.000000             |
| Kegg cytokine cytokine receptor interaction  | 0.0000000         | 0.000000             |
| Kegg basal cell carcinoma                    | 0.0000000         | 0.000000             |
| Kegg mapk signaling pathway                  | 0.1056889         | 0.001822             |
| Kegg hedgehog signaling pathway              | 0.0070290         | 0.002050             |
| Kegg cell adhesion molecules cams            | 0.0184500         | 0.013120             |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 23 enriched KEGG pathways with ermineR.

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
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY  |  0.0062396|  0.6907024|  2.014565|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY |  0.0062396|  0.6045600|  1.938675|             0|
| KEGG\_PATHWAYS\_IN\_CANCER          |  0.0062396|  0.4536458|  1.696470|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA        |  0.0062396|  0.6672720|  1.961995|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY       |  0.0108467|  0.4915260|  1.701823|             2|
| KEGG\_ECM\_RECEPTOR\_INTERACTION    |  0.0237103|  0.5443352|  1.721987|             6|
| KEGG\_MELANOGENESIS                 |  0.0261296|  0.5240401|  1.693749|             8|

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
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  | 0.0000000         | 0.0000000            |
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                        | 0.0000000         | 0.0000000            |
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens        | 0.0000000         | 0.0037000            |
| GPCRs, Class A Rhodopsin-like%WikiPathways\_20180810%WP455%Homo sapiens                                        | 0.0059200         | 0.0131556            |
| Wnt Signaling Pathway and Pluripotency%WikiPathways\_20180810%WP399%Homo sapiens                               | 0.0065780         | 0.0207200            |
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                            | 0.0500923         | 0.0269091            |

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

There are 15 up-regulated and 3 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                           |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                           |  0.0057409|  0.6050129|  1.965483|             0|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                            |  0.0057409|  0.5251821|  1.769434|             0|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                              |  0.0057409|  0.5598027|  1.862854|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                          |  0.0057409|  0.5507020|  1.916223|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                      |  0.0057409|  0.7626768|  2.170951|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                     |  0.0057409|  0.6719177|  1.976613|             0|
| Breast cancer pathway%WikiPathways\_20180810%WP4262%Homo sapiens                                                  |  0.0057409|  0.5349072|  1.842839|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens              |  0.0057409|  0.5702915|  1.980160|             0|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens    |  0.0091971|  0.5696587|  1.850629|             1|
| TGF-B Signaling in Thyroid Cells for Epithelial-Mesenchymal Transition%WikiPathways\_20180810%WP3859%Homo sapiens |  0.0095018|  0.7757247|  1.885462|             1|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                                |  0.0101372|  0.6859226|  1.914298|             2|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                        |  0.0101372|  0.6835662|  1.858062|             2|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                   |  0.0152431|  0.5027597|  1.675641|             5|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                               |  0.0205840|  0.4659379|  1.609632|             8|
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens           |  0.0407813|  0.5482734|  1.699790|            16|

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
| Neuronal System                                                                      | 0.0000000         | 0.0122800            |
| G alpha (i) signalling events                                                        | 0.0000000         | 0.0153500            |
| Class B 2 (Secretin family receptors)                                                | 0.0921000         | 0.0921000            |
| Incretin synthesis, secretion, and inactivation                                      | 0.1315714         | 0.1315714            |
| SLC-mediated transmembrane transport                                                 | 0.1773778         | 0.2763000            |
| Keratan sulfate keratin metabolism                                                   | 0.2532750         | 0.2824400            |

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

There are 30 up-regulated and 188 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                             |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244           |  0.0050692|  0.5005085|  1.897086|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                           |  0.0050692|  0.6312265|  1.905791|             1|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                |  0.0068358|  0.6303752|  1.860367|             3|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275  |  0.0075077|  0.5369792|  1.791875|             4|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                       |  0.0095508|  0.5745773|  1.779483|             6|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                         |  0.0099326|  0.5127644|  1.734311|             7|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814 |  0.0109425|  0.5624466|  1.774667|             8|
| G ALPHA (I) SIGNALLING EVENTS%REACTOME%R-HSA-418594.5                               |  0.0109515|  0.4089427|  1.533419|            10|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                         |  0.0135778|  0.5184447|  1.702630|            11|
| DOPAMINE NEUROTRANSMITTER RELEASE CYCLE%REACTOME%R-HSA-212676.2                     |  0.0141073|  0.7125226|  1.809762|            10|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                       |  0.0178566|  0.4834313|  1.646091|            17|
| OTHER SEMAPHORIN INTERACTIONS%REACTOME DATABASE ID RELEASE 65%416700                |  0.0206760|  0.7196216|  1.750987|            18|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                          |  0.0214026|  0.5511119|  1.688575|            21|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                     |  0.0215693|  0.5962055|  1.744505|            21|
| SURFACTANT METABOLISM%REACTOME DATABASE ID RELEASE 65%5683826                       |  0.0218078|  0.7092872|  1.752413|            20|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                          |  0.0246910|  0.5342714|  1.672355|            27|
| EFFECTS OF PIP2 HYDROLYSIS%REACTOME%R-HSA-114508.2                                  |  0.0260590|  0.6320716|  1.710793|            27|
| REGULATION OF COMPLEMENT CASCADE%REACTOME%R-HSA-977606.5                            |  0.0279000|  0.7146985|  1.710884|            29|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                 |  0.0279785|  0.4605667|  1.576998|            36|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2               |  0.0286555|  0.3932485|  1.472080|            42|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                |  0.0300728|  0.6136479|  1.700634|            34|
| COMPLEMENT CASCADE%REACTOME DATABASE ID RELEASE 65%166658                           |  0.0315754|  0.6993303|  1.701614|            36|
| ION CHANNEL TRANSPORT%REACTOME%R-HSA-983712.2                                       |  0.0354979|  0.4166411|  1.489751|            53|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                            |  0.0354979|  0.6159588|  1.681914|            43|
| KERATAN SULFATE KERATIN METABOLISM%REACTOME%R-HSA-1638074.1                         |  0.0369730|  0.6146049|  1.678217|            46|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                       |  0.0369730|  0.6142315|  1.677198|            46|
| MUSCLE CONTRACTION%REACTOME%R-HSA-397014.2                                          |  0.0456231|  0.3955380|  1.438880|            76|
| CELL-CELL JUNCTION ORGANIZATION%REACTOME%R-HSA-421270.4                             |  0.0470851|  0.5041505|  1.583028|            69|
| PLASMA LIPOPROTEIN REMODELING%REACTOME DATABASE ID RELEASE 65%8963899               |  0.0470851|  0.6633448|  1.664950|            62|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1       |  0.0496513|  0.5935540|  1.634261|            68|

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
| Kegg pathways in cancer                      | 0.0246000         | 0.0000000            |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.0000000            |
| Kegg cytokine cytokine receptor interaction  | 0.0000000         | 0.0000000            |
| Kegg hematopoietic cell lineage              | 0.0109333         | 0.0041000            |
| Kegg wnt signaling pathway                   | 0.0229600         | 0.0098400            |
| Kegg cell adhesion molecules cams            | 0.0300667         | 0.0191333            |
| Kegg melanogenesis                           | 0.1590800         | 0.0374857            |
| Kegg ecm receptor interaction                | 0.0632571         | 0.0451000            |
| Kegg basal cell carcinoma                    | 0.1658222         | 0.0606800            |
| Kegg focal adhesion                          | 0.4842316         | 0.0656000            |

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
| KEGG\_ECM\_RECEPTOR\_INTERACTION                                    |  0.0125332|  0.5578373|  1.825807|             1|
| KEGG\_LEUKOCYTE\_TRANSENDOTHELIAL\_MIGRATION                        |  0.0125332|  0.5034342|  1.730487|             1|
| KEGG\_CELL\_ADHESION\_MOLECULES\_CAMS                               |  0.0153674|  0.4921250|  1.691613|             3|
| KEGG\_AXON\_GUIDANCE                                                |  0.0167388|  0.4619294|  1.644947|             5|
| KEGG\_FOCAL\_ADHESION                                               |  0.0167388|  0.4212441|  1.564264|             5|
| KEGG\_GLYCOSPHINGOLIPID\_BIOSYNTHESIS\_LACTO\_AND\_NEOLACTO\_SERIES |  0.0350834|  0.6982301|  1.817136|            12|

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
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens | 0.0000000         | 0.0000000            |
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens       | 0.0444000         | 0.0296000            |
| GPCRs, Class A Rhodopsin-like%WikiPathways\_20180810%WP455%Homo sapiens       | 0.0394667         | 0.0394667            |
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens      | 0.0473600         | 0.0473600            |

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

There are 5 up-regulated and 9 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                       |  0.0244547|  0.6217300|  2.205633|             0|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens |  0.0250635|  0.3960488|  1.515990|             5|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                     |  0.0250635|  0.3972533|  1.530507|             5|
| Phosphodiesterases in neuronal function%WikiPathways\_20180810%WP4222%Homo sapiens        |  0.0250635|  0.5983603|  1.798489|             4|
| Hematopoietic Stem Cell Differentiation%WikiPathways\_20180810%WP2849%Homo sapiens        |  0.0462067|  0.5863630|  1.734946|            13|

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

| Pathway                                                                              | Corrected p-value | Corrected MF p-value |
|:-------------------------------------------------------------------------------------|:------------------|:---------------------|
| G alpha (i) signalling events                                                        | 0.0000000         | 0.0000000            |
| Extracellular matrix organization                                                    | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)                                                 | 0.0000000         | 0.0000000            |
| Incretin synthesis, secretion, and inactivation                                      | 0.0153500         | 0.0153500            |
| Neuronal System                                                                      | 0.0245600         | 0.1228000            |
| G alpha (q) signalling events                                                        | 0.1330333         | 0.1842000            |
| Amine ligand-binding receptors                                                       | 0.1578857         | 0.1842000            |
| Other interleukin signaling                                                          | 0.3903286         | 0.4166429            |
| Activation of anterior HOX genes in hindbrain development during early embryogenesis | 0.4298000         | 0.4239524            |
| Activation of HOX genes during differentiation                                       | 0.4298000         | 0.4239524            |

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

There are 32 up-regulated and 220 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                         |       padj|         ES|       NES|  nMoreExtreme|
|:----------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                            |  0.0054916|  0.6778319|  1.998935|             1|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                             |  0.0054916|  0.6026327|  1.888871|             1|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                       |  0.0054916|  0.5498580|  2.047200|             0|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                                  |  0.0054916|  0.8103541|  2.014676|             0|
| CLASS B 2 (SECRETIN FAMILY RECEPTORS)%REACTOME DATABASE ID RELEASE 65%373080                                    |  0.0054916|  0.6183446|  1.926980|             1|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                                     |  0.0054916|  0.5818259|  1.896239|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                       |  0.0054916|  0.6573691|  1.981918|             1|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                                   |  0.0054916|  0.6976477|  1.898277|             2|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                 |  0.0054916|  0.7379812|  2.152629|             0|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                      |  0.0054916|  0.6334096|  1.970161|             0|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                               |  0.0054916|  0.7410476|  2.016367|             0|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                      |  0.0054916|  0.4563494|  1.658089|             2|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                                        |  0.0057687|  0.6841620|  1.861583|             4|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090            |  0.0057687|  0.6042254|  1.850693|             4|
| KERATAN SULFATE BIOSYNTHESIS%REACTOME%R-HSA-2022854.1                                                           |  0.0088884|  0.6994094|  1.824733|             8|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                                          |  0.0114068|  0.4788939|  1.611305|            15|
| PEPTIDE HORMONE METABOLISM%REACTOME DATABASE ID RELEASE 65%2980736                                              |  0.0119377|  0.5259866|  1.663535|            15|
| KERATAN SULFATE KERATIN METABOLISM%REACTOME%R-HSA-1638074.1                                                     |  0.0128762|  0.6564317|  1.786130|            15|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                                   |  0.0178986|  0.4696021|  1.581449|            28|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                                   |  0.0208761|  0.6232640|  1.709158|            29|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                              |  0.0208761|  0.4785204|  1.587013|            34|
| SIGNALING BY NODAL%REACTOME DATABASE ID RELEASE 65%1181150                                                      |  0.0216531|  0.6946138|  1.726926|            29|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                                     |  0.0248095|  0.4638452|  1.554841|            44|
| INCRETIN SYNTHESIS, SECRETION, AND INACTIVATION%REACTOME%R-HSA-400508.1                                         |  0.0261303|  0.7041486|  1.700235|            37|
| SIGNALING BY PDGF%REACTOME DATABASE ID RELEASE 65%186797                                                        |  0.0315309|  0.5212166|  1.596444|            54|
| IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON-LYMPHOID CELL%REACTOME DATABASE ID RELEASE 65%198933 |  0.0315309|  0.5408741|  1.610866|            53|
| SURFACTANT METABOLISM%REACTOME DATABASE ID RELEASE 65%5683826                                                   |  0.0324320|  0.6749192|  1.677962|            49|
| INTERLEUKIN-6 FAMILY SIGNALING%REACTOME DATABASE ID RELEASE 65%6783589                                          |  0.0347914|  0.6508816|  1.680380|            56|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                               |  0.0352524|  0.6499894|  1.678077|            57|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                                    |  0.0422890|  0.4381799|  1.494026|            88|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                                 |  0.0425836|  0.4603672|  1.519702|            86|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                                         |  0.0495101|  0.6753657|  1.630736|            83|

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
| Kegg pathways in cancer                      | 0.0348500         | 0.0000000            |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.0000000            |
| Kegg focal adhesion                          | 0.0351429         | 0.0041000            |
| Kegg cytokine cytokine receptor interaction  | 0.0000000         | 0.0054670            |
| Kegg hematopoietic cell lineage              | 0.0205000         | 0.0082000            |
| Kegg ecm receptor interaction                | 0.0218667         | 0.0098400            |
| Kegg tgf beta signaling pathway              | 0.0355333         | 0.0187429            |
| Kegg cell adhesion molecules cams            | 0.0328000         | 0.0205000            |
| Kegg calcium signaling pathway               | 0.2214000         | 0.0856444            |
| Kegg systemic lupus erythematosus            | 0.2087273         | 0.1804000            |

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

There are 11 up-regulated and 7 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                         |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY              |  0.0028488|  0.6493254|  1.962970|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY             |  0.0028488|  0.5752288|  1.908204|             0|
| KEGG\_AXON\_GUIDANCE                            |  0.0028488|  0.5025458|  1.776486|             0|
| KEGG\_FOCAL\_ADHESION                           |  0.0028488|  0.4824189|  1.767004|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                |  0.0028488|  0.6018477|  1.973276|             0|
| KEGG\_HEMATOPOIETIC\_CELL\_LINEAGE              |  0.0028488|  0.6173100|  1.873857|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                      |  0.0028488|  0.4193480|  1.594976|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA                    |  0.0028488|  0.6480634|  1.974066|             0|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION |  0.0234486|  0.4471618|  1.580705|            17|
| KEGG\_WNT\_SIGNALING\_PATHWAY                   |  0.0242847|  0.4376114|  1.557108|            19|
| KEGG\_LYSOSOME                                  |  0.0498851|  0.4384520|  1.530303|            44|

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
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens      | 0.0296000         | 0.0296000            |
| GPCRs, Class A Rhodopsin-like%WikiPathways\_20180810%WP455%Homo sapiens       | 0.0148000         | 0.0296000            |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens | 0.0394667         | 0.0493333            |

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 3 enriched Wiki pathways with ermineR

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

There are 24 up-regulated and 13 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                        |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                     |  0.0059937|  0.5738558|  1.841895|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       |  0.0059937|  0.5103239|  1.825835|             0|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                            |  0.0059937|  0.5535169|  1.950233|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                   |  0.0059937|  0.7855874|  2.295872|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  |  0.0059937|  0.6736308|  2.040369|             0|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                                          |  0.0059937|  0.4303336|  1.629774|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens           |  0.0059937|  0.5806698|  2.073771|             0|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                       |  0.0066352|  0.4637370|  1.697689|             1|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                         |  0.0066352|  0.4970351|  1.724609|             1|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                           |  0.0066352|  0.5284407|  1.810255|             1|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                      |  0.0066352|  0.4274579|  1.613788|             1|
| Peptide GPCRs%WikiPathways\_20180810%WP24%Homo sapiens                                                         |  0.0074725|  0.7201486|  1.898644|             1|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                     |  0.0100806|  0.6828241|  1.909568|             3|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens |  0.0124563|  0.5188105|  1.738063|             5|
| TGF-beta Receptor Signaling%WikiPathways\_20180810%WP560%Homo sapiens                                          |  0.0125975|  0.5792501|  1.794713|             5|
| Hematopoietic Stem Cell Differentiation%WikiPathways\_20180810%WP2849%Homo sapiens                             |  0.0137649|  0.6112717|  1.813765|             6|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                   |  0.0161091|  0.5067811|  1.709391|             9|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                      |  0.0181306|  0.4031603|  1.507891|            14|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                |  0.0200750|  0.4742329|  1.629488|            15|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                             |  0.0245811|  0.6106621|  1.752753|            17|
| BMP Signaling Pathway in Eyelid Development%WikiPathways\_20180810%WP3927%Homo sapiens                         |  0.0273626|  0.7084519|  1.768914|            18|
| miRNA targets in ECM and membrane receptors%WikiPathways\_20180810%WP2911%Homo sapiens                         |  0.0298899|  0.6587497|  1.736768|            21|
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                        |  0.0307745|  0.4812308|  1.612167|            27|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                          |  0.0364513|  0.4408707|  1.546007|            36|
