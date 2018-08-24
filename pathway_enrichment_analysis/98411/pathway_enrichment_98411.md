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
| Neuronal System                            | 0.0000000         | 0.0000000            |
| Extracellular matrix organization          | 0.0000000         | 0.0000000            |
| TCF dependent signaling in response to WNT | 0.0000000         | 0.0144500            |
| Class B 2 (Secretin family receptors)      | 0.0000000         | 0.0192667            |
| G alpha (i) signalling events              | 0.0000000         | 0.0578000            |
| Class A 1 (Rhodopsin-like receptors)       | 0.0330286         | 0.0578000            |
| G alpha (s) signalling events              | 0.0578000         | 0.0990857            |
| PI3K Cascade                               | 0.3146889         | 0.3468000            |
| FGFR2 mutant receptor activation           | 0.2526074         | 0.3543391            |
| Negative regulation of FGFR3 signaling     | 0.2312000         | 0.3651909            |

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
| Kegg focal adhesion                          | 0.0633857         | 0.0000000            |
| Kegg mapk signaling pathway                  | 0.0765000         | 0.0076500            |
| Kegg melanogenesis                           | 0.0595000         | 0.0091800            |
| Kegg neuroactive ligand receptor interaction | 0.0051000         | 0.0102000            |
| Kegg basal cell carcinoma                    | 0.1029273         | 0.0284143            |
| Kegg ecm receptor interaction                | 0.0573750         | 0.0401625            |
| Kegg regulation of actin cytoskeleton        | 0.1972895         | 0.0417273            |
| Kegg calcium signaling pathway               | 0.1045500         | 0.0425000            |

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
| KEGG\_WNT\_SIGNALING\_PATHWAY       |  0.0071562|  0.5319154|  1.826301|             0|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY  |  0.0071562|  0.7310908|  2.020332|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY |  0.0071562|  0.6162020|  1.924486|             0|
| KEGG\_PATHWAYS\_IN\_CANCER          |  0.0071562|  0.4818222|  1.817863|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA        |  0.0071562|  0.7290053|  2.048621|             0|

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
| Post-translational protein phosphorylation | 0.0000000         | 0.0000000            |
| G alpha (i) signalling events              | 0.0000000         | 0.0000000            |
| Extracellular matrix organization          | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)       | 0.0000000         | 0.0000000            |
| Collagen formation                         | 0.0578000         | 0.1040400            |
| Signaling by FGFR3 in disease              | 0.1194533         | 0.1194533            |
| Negative regulation of FGFR3 signaling     | 0.1197286         | 0.1197286            |
| Negative regulation of FGFR2 signaling     | 0.1523818         | 0.1289385            |
| Collagen chain trimerization               | 0.1222692         | 0.1387200            |
| FGFR1 mutant receptor activation           | 0.1247263         | 0.1394000            |

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

There are 43 up-regulated and 243 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                              |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                 |  0.0063102|  0.6831081|  1.885341|             0|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                |  0.0063102|  0.4963027|  1.670513|             0|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                        |  0.0063102|  0.7137058|  2.053218|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                  |  0.0063102|  0.5934265|  1.772728|             3|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                            |  0.0063102|  0.5782104|  1.971507|             0|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                  |  0.0063102|  0.6030993|  1.899199|             0|
| G ALPHA (I) SIGNALLING EVENTS%REACTOME%R-HSA-418594.5                                                |  0.0063102|  0.4870838|  1.642524|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                          |  0.0063102|  0.6005742|  1.853498|             0|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                   |  0.0063102|  0.6404478|  1.987522|             0|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                        |  0.0063102|  0.6221719|  1.943392|             0|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                           |  0.0063102|  0.6815479|  1.916315|             1|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                           |  0.0063102|  0.6394826|  1.832695|             2|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090 |  0.0063102|  0.6869231|  1.968655|             0|
| PLASMA LIPOPROTEIN REMODELING%REACTOME DATABASE ID RELEASE 65%8963899                                |  0.0063102|  0.8505300|  2.007580|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                          |  0.0063102|  0.6356000|  1.985335|             0|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                      |  0.0071208|  0.5905612|  1.749456|             8|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                      |  0.0071208|  0.6478909|  1.768344|             7|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                             |  0.0071208|  0.5952607|  1.748914|             8|
| CELL JUNCTION ORGANIZATION%REACTOME DATABASE ID RELEASE 65%446728                                    |  0.0071507|  0.5570686|  1.693026|             9|
| VISUAL PHOTOTRANSDUCTION%REACTOME DATABASE ID RELEASE 65%2187338                                     |  0.0087816|  0.5728466|  1.686473|            14|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                      |  0.0092105|  0.5332283|  1.643172|            16|
| STRIATED MUSCLE CONTRACTION%REACTOME DATABASE ID RELEASE 65%390522                                   |  0.0103593|  0.6993129|  1.766069|            15|
| MUSCLE CONTRACTION%REACTOME%R-HSA-397014.2                                                           |  0.0108719|  0.4738146|  1.551778|            21|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                         |  0.0126313|  0.5155537|  1.615185|            25|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                                 |  0.0127643|  0.6707494|  1.749141|            21|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                               |  0.0156617|  0.5127598|  1.588646|            33|
| CELL-CELL COMMUNICATION%REACTOME%R-HSA-1500931.3                                                     |  0.0192704|  0.4939098|  1.559680|            44|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                         |  0.0192704|  0.5570028|  1.624477|            40|
| DEATH RECEPTOR SIGNALLING%REACTOME%R-HSA-73887.3                                                     |  0.0205632|  0.4690446|  1.524868|            49|
| RETINOID METABOLISM AND TRANSPORT%REACTOME%R-HSA-975634.2                                            |  0.0206512|  0.6348848|  1.683572|            40|
| NON-INTEGRIN MEMBRANE-ECM INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000171                       |  0.0221017|  0.5919711|  1.640565|            45|
| ACTIVATION OF MATRIX METALLOPROTEINASES%REACTOME%R-HSA-1592389.1                                     |  0.0283806|  0.7073713|  1.669670|            54|
| GAP JUNCTION TRAFFICKING AND REGULATION%REACTOME%R-HSA-157858.1                                      |  0.0306006|  0.7055499|  1.665371|            59|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                            |  0.0311303|  0.5765545|  1.606399|            68|
| CELL-CELL JUNCTION ORGANIZATION%REACTOME%R-HSA-421270.4                                              |  0.0316029|  0.5617993|  1.595171|            71|
| METABOLISM OF FAT-SOLUBLE VITAMINS%REACTOME DATABASE ID RELEASE 65%6806667                           |  0.0328990|  0.5997627|  1.621046|            72|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                    |  0.0347688|  0.6378608|  1.640440|            73|
| P75 NTR RECEPTOR-MEDIATED SIGNALLING%REACTOME%R-HSA-193704.1                                         |  0.0364908|  0.4838781|  1.517076|            94|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                           |  0.0392329|  0.5920108|  1.600094|            89|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                           |  0.0464977|  0.4349069|  1.433381|           131|
| NRAGE SIGNALS DEATH THROUGH JNK%REACTOME DATABASE ID RELEASE 65%193648                               |  0.0466464|  0.5251015|  1.545910|           116|
| CELL DEATH SIGNALLING VIA NRAGE, NRIF AND NADE%REACTOME DATABASE ID RELEASE 65%204998                |  0.0486054|  0.4969968|  1.513192|           127|
| SHC1 EVENTS IN ERBB2 SIGNALING%REACTOME DATABASE ID RELEASE 65%1250196                               |  0.0490205|  0.6730409|  1.612745|           106|

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
| Kegg neuroactive ligand receptor interaction      | 0.0229500         | 0.0051000            |
| Kegg cytokine cytokine receptor interaction       | 0.0000000         | 0.0076500            |
| Kegg axon guidance                                | 0.0275400         | 0.0153000            |
| Kegg ecm receptor interaction                     | 0.0229500         | 0.0382500            |
| Kegg intestinal immune network for iga production | 0.0677571         | 0.1499400            |
| Kegg endocytosis                                  | 0.4371429         | 0.1938000            |
| Kegg complement and coagulation cascades          | 0.1411000         | 0.2032714            |
| Kegg basal cell carcinoma                         | 0.1517250         | 0.2142000            |
| Kegg focal adhesion                               | 0.0255000         | 0.2312000            |
| Kegg calcium signaling pathway                    | 0.0822375         | 0.2475818            |

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

There are 12 up-regulated and 13 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                         |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_AXON\_GUIDANCE                            |  0.0114179|  0.5525302|  1.799373|             0|
| KEGG\_FOCAL\_ADHESION                           |  0.0114179|  0.5099686|  1.709489|             1|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                |  0.0114179|  0.5948298|  1.804096|             4|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES    |  0.0114179|  0.7010039|  1.824223|             1|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION |  0.0207598|  0.5188667|  1.644826|            11|
| KEGG\_LEUKOCYTE\_TRANSENDOTHELIAL\_MIGRATION    |  0.0329174|  0.5133445|  1.609967|            29|
| KEGG\_PPAR\_SIGNALING\_PATHWAY                  |  0.0340909|  0.5762875|  1.664645|            38|
| KEGG\_HYPERTROPHIC\_CARDIOMYOPATHY\_HCM         |  0.0340909|  0.5495000|  1.638802|            39|
| KEGG\_DILATED\_CARDIOMYOPATHY                   |  0.0340909|  0.5429973|  1.642353|            34|
| KEGG\_LYSOSOME                                  |  0.0413867|  0.4705680|  1.519597|            55|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                  |  0.0431102|  0.4155230|  1.411902|            64|
| KEGG\_CELL\_ADHESION\_MOLECULES\_CAMS           |  0.0493087|  0.5031673|  1.555181|            68|

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
| Neuronal System                                | 0.0939250         | 0.0000000            |
| G alpha (i) signalling events                  | 0.0000000         | 0.0000000            |
| Extracellular matrix organization              | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)           | 0.0000000         | 0.0000000            |
| SLC-mediated transmembrane transport           | 0.0144500         | 0.0115600            |
| Platelet activation, signaling and aggregation | 0.0990857         | 0.0192667            |
| Collagen formation                             | 0.0346800         | 0.0247714            |
| Post-translational protein phosphorylation     | 0.1059667         | 0.1011500            |
| Signaling by Interleukins                      | 0.8959000         | 0.1541333            |
| Collagen chain trimerization                   | 0.2825778         | 0.2543200            |

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

There are 32 up-regulated and 144 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                          |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                             |  0.0089004|  0.6906764|  1.910710|             2|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                                     |  0.0089004|  0.6582874|  1.923559|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                              |  0.0089004|  0.6095565|  1.827085|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                        |  0.0089004|  0.5530032|  1.919082|             0|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                              |  0.0089004|  0.5414102|  1.721309|             3|
| INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING%REACTOME DATABASE ID RELEASE 65%6785807                               |  0.0089004|  0.5562357|  1.715465|             3|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                                       |  0.0089004|  0.7017027|  1.894088|             2|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                                      |  0.0089004|  0.5937514|  1.847148|             0|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                                    |  0.0089004|  0.5818321|  1.833213|             1|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                  |  0.0089004|  0.7193802|  1.966847|             0|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                                         |  0.0089004|  0.6007229|  1.770771|             3|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                       |  0.0089004|  0.5028360|  1.684362|             2|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                                     |  0.0089004|  0.5324235|  1.705389|             3|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090             |  0.0089004|  0.6493038|  1.865956|             1|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                                      |  0.0089004|  0.5419548|  1.707570|             3|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                               |  0.0091033|  0.5408211|  1.689659|             4|
| TRANSPORT OF BILE SALTS AND ORGANIC ACIDS, METAL IONS AND AMINE COMPOUNDS%REACTOME DATABASE ID RELEASE 65%425366 |  0.0092293|  0.6137404|  1.774005|             4|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                       |  0.0100163|  0.6202001|  1.782319|             6|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                                    |  0.0107006|  0.6072272|  1.755179|             8|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                        |  0.0119109|  0.6339416|  1.774015|             9|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                                |  0.0124841|  0.7010078|  1.793804|             9|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                                         |  0.0162813|  0.6956923|  1.780202|            13|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                                     |  0.0179532|  0.5141659|  1.627829|            19|
| IRE1ALPHA ACTIVATES CHAPERONES%REACTOME DATABASE ID RELEASE 65%381070                                            |  0.0263473|  0.5578257|  1.660629|            29|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                            |  0.0283267|  0.4337733|  1.486368|            38|
| SIGNALING BY WNT IN CANCER%REACTOME DATABASE ID RELEASE 65%4791275                                               |  0.0326994|  0.6290613|  1.698009|            36|
| SIGNALING BY FGFR1%REACTOME DATABASE ID RELEASE 65%5654736                                                       |  0.0394884|  0.6004543|  1.661117|            50|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                                       |  0.0404024|  0.5821210|  1.642850|            53|
| XBP1(S) ACTIVATES CHAPERONE GENES%REACTOME%R-HSA-381038.2                                                        |  0.0426866|  0.5462921|  1.618779|            60|
| RUNX1 INTERACTS WITH CO-FACTORS WHOSE PRECISE EFFECT ON RUNX1 TARGETS IS NOT KNOWN%REACTOME%R-HSA-8939243.1      |  0.0431552|  0.6012725|  1.652724|            57|
| SIGNALING BY INTERLEUKINS%REACTOME%R-HSA-449147.10                                                               |  0.0436052|  0.3884707|  1.379264|            78|
| RHO GTPASE CYCLE%REACTOME DATABASE ID RELEASE 65%194840                                                          |  0.0448704|  0.4496612|  1.493905|            75|

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
| Kegg ecm receptor interaction                | 0.0000000         | 0.0000000            |
| Kegg calcium signaling pathway               | 0.0382500         | 0.0051000            |
| Kegg pathways in cancer                      | 0.1514700         | 0.0076500            |
| Kegg cytokine cytokine receptor interaction  | 0.0306000         | 0.0091800            |
| Kegg wnt signaling pathway                   | 0.1702125         | 0.0655714            |
| Kegg neuroactive ligand receptor interaction | 0.0397800         | 0.0707625            |
| Kegg tgf beta signaling pathway              | 0.1122000         | 0.0765000            |
| Kegg mapk signaling pathway                  | 0.1447615         | 0.0884000            |
| Kegg cell adhesion molecules cams            | 0.1639286         | 0.0994500            |

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

There are 12 up-regulated and 8 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                         |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION |  0.0037890|  0.6292265|  2.010726|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY             |  0.0037890|  0.6570327|  2.016665|             0|
| KEGG\_AXON\_GUIDANCE                            |  0.0037890|  0.5635419|  1.858536|             0|
| KEGG\_FOCAL\_ADHESION                           |  0.0037890|  0.4930680|  1.684289|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                      |  0.0037890|  0.4809057|  1.689648|             0|
| KEGG\_SMALL\_CELL\_LUNG\_CANCER                 |  0.0156638|  0.5395996|  1.690922|            10|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                |  0.0164175|  0.5720013|  1.738258|            11|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                  |  0.0239156|  0.4328452|  1.498346|            21|
| KEGG\_WNT\_SIGNALING\_PATHWAY                   |  0.0348337|  0.4734545|  1.567926|            32|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES    |  0.0348337|  0.6656616|  1.735588|            26|
| KEGG\_REGULATION\_OF\_ACTIN\_CYTOSKELETON       |  0.0359806|  0.4397563|  1.497467|            40|
| KEGG\_BASAL\_CELL\_CARCINOMA                    |  0.0359806|  0.6020553|  1.695145|            31|
