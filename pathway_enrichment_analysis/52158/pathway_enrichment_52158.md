Pathway enrichment analysis, GEO52158
================
German Novakovskiy
August 21, 2018

This is an analysis of microarray data, which mostly will not be considered.

Analysis of first day of differentiation 0h VS 24 h (formation of anterior Primitive streak, mesodendoderm)
-----------------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE52158/DEgenes_0h_24h_52158.Rdata")
```

Converting probe ids to gene names:

``` r
x <- hgu133plus2SYMBOL
# Get the probe identifiers - gene symbol mappings
mapped_probes <- mappedkeys(x)
# Convert to a dataframe
genesym.probeid <- as.data.frame(x[mapped_probes])
head(genesym.probeid)
```

    ##    probe_id symbol
    ## 1   1053_at   RFC2
    ## 2    117_at  HSPA6
    ## 3    121_at   PAX8
    ## 4 1255_g_at GUCA1A
    ## 5   1316_at   THRA
    ## 6   1320_at PTPN21

Sorted log Fold Changes give us a ranked list:

``` r
probes_to_genes <- genesym.probeid %>%
  filter(probe_id %in% rownames(DEgenes_0h_24h_52158))

topProbes <- DEgenes_0h_24h_52158 %>%
  rownames_to_column("probes") %>%
  filter(probes %in% probes_to_genes$probe_id)

probes_to_genes <- probes_to_genes %>% column_to_rownames('probe_id')
symbs <- probes_to_genes[topProbes$probes,]
topProbes$Symbol <- symbs

#for ermineR
ermineInputProbeScores <- topProbes %>% 
  #as.data.frame() %>%
  mutate(absolute_logFC = logFC) %>% 
  dplyr::select(Symbol, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) 

#randomly delete repeats (we will delete the lowerst ones)
repeats <- which(duplicated(ermineInputProbeScores$Symbol))
ermineInputProbeScores <- ermineInputProbeScores[-repeats,]

rownames(ermineInputProbeScores) <- NULL
ermineInputProbeScores <- ermineInputProbeScores %>% 
  column_to_rownames("Symbol")

#for fgsea
ermineInputProbeFGSEA <- topProbes %>% 
  #as.data.frame() %>%
  #mutate(absolute_logFC = logFC) %>% 
  dplyr::select(Symbol, logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(logFC)) 

#randomly delete repeats (we will delete the lowerst ones)
repeats <- which(duplicated(ermineInputProbeFGSEA$Symbol))
ermineInputProbeFGSEA <- ermineInputProbeFGSEA[-repeats,]

rownames(ermineInputProbeFGSEA) <- NULL
ermineInputProbeFGSEA <- ermineInputProbeFGSEA %>% 
  column_to_rownames("Symbol")

scoresFGSEA <- ermineInputProbeFGSEA$logFC
names(scoresFGSEA) <- rownames(ermineInputProbeFGSEA)
```

### Reactome

#### ErmineR Reactome pathways

``` r
enrichmentResultReactome <- precRecall(scores = ermineInputProbeScores,
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
| TCF dependent signaling in response to WNT                                           | 0.0000000         | 0.5409000            |
| SIRT1 negatively regulates rRNA expression                                           | 0.4337652         | 0.6942586            |
| PI5P, PP2A and IER3 Regulate PI3K AKT Signaling                                      | 0.3828593         | 0.7083214            |
| Regulation of beta-cell development                                                  | 0.3090857         | 0.7235115            |
| Signaling by FGFR1                                                                   | 0.3411081         | 0.7278778            |
| Negative regulation of the PI3K AKT network                                          | 0.3218689         | 0.7281884            |
| Activation of anterior HOX genes in hindbrain development during early embryogenesis | 0.3730345         | 0.7455262            |
| Activation of HOX genes during differentiation                                       | 0.3730345         | 0.7455262            |
| Signaling by FGFR2 in disease                                                        | 0.3409121         | 0.7457864            |
| Pre-NOTCH Transcription and Translation                                              | 0.3975846         | 0.7519178            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 0 enriched Reactome pathways with ermineR

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

There are 15 up-regulated and 4 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                          |       padj|         ES|       NES|  nMoreExtreme|
|:-------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| PROCESSING OF CAPPED INTRON-CONTAINING PRE-MRNA%REACTOME DATABASE ID RELEASE 65%72203            |  0.0289904|  0.4349910|  1.647490|             1|
| PRE-MRNA SPLICING%REACTOME DATABASE ID RELEASE 65%72163                                          |  0.0289904|  0.4592533|  1.690582|             1|
| RESOLUTION OF SISTER CHROMATID COHESION%REACTOME DATABASE ID RELEASE 65%2500257                  |  0.0289904|  0.5089663|  1.734457|             1|
| FORMATION OF THE BETA-CATENIN:TCF TRANSACTIVATING COMPLEX%REACTOME DATABASE ID RELEASE 65%201722 |  0.0289904|  0.5962274|  1.825655|             1|
| HATS ACETYLATE HISTONES%REACTOME%R-HSA-3214847.1                                                 |  0.0289904|  0.5295682|  1.793410|             1|
| SIGNALING BY WNT IN CANCER%REACTOME DATABASE ID RELEASE 65%4791275                               |  0.0289904|  0.6771217|  1.888521|             1|
| CHROMATIN ORGANIZATION%REACTOME%R-HSA-4839726.2                                                  |  0.0289904|  0.5246601|  1.976443|             0|
| CHROMATIN MODIFYING ENZYMES%REACTOME%R-HSA-3247509.4                                             |  0.0289904|  0.5246601|  1.976443|             0|
| MRNA SPLICING%REACTOME%R-HSA-72172.3                                                             |  0.0345004|  0.4567411|  1.690196|             2|
| ACTIVATION OF GENE EXPRESSION BY SREBF (SREBP)%REACTOME%R-HSA-2426168.2                          |  0.0362397|  0.6490888|  1.887083|             2|
| CONDENSATION OF PROPHASE CHROMOSOMES%REACTOME DATABASE ID RELEASE 65%2299718                     |  0.0403859|  0.6681044|  1.885595|             3|
| RNA POLYMERASE I PROMOTER OPENING%REACTOME%R-HSA-73728.2                                         |  0.0403859|  0.7140352|  1.902126|             3|
| B-WICH COMPLEX POSITIVELY REGULATES RRNA EXPRESSION%REACTOME%R-HSA-5250924.2                     |  0.0403859|  0.5853037|  1.792207|             3|
| MITOTIC PROMETAPHASE%REACTOME%R-HSA-68877.3                                                      |  0.0409109|  0.4378091|  1.613843|             4|
| SIRT1 NEGATIVELY REGULATES RRNA EXPRESSION%REACTOME%R-HSA-427359.2                               |  0.0441396|  0.6623209|  1.830636|             4|

### KEGG pathways

#### ErmineR KEGG pathways

``` r
enrichmentResultKEGG <- precRecall(scores = ermineInputProbeScores,
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

| Pathway                                | Corrected p-value | Corrected MF p-value |
|:---------------------------------------|:------------------|:---------------------|
| Kegg wnt signaling pathway             | 0.0159000         | 0.000000             |
| Kegg prostate cancer                   | 0.5024400         | 0.000000             |
| Kegg pathways in cancer                | 0.6479250         | 0.000000             |
| Kegg neurotrophin signaling pathway    | 0.6042000         | 0.002650             |
| Kegg erbb signaling pathway            | 0.5807684         | 0.003180             |
| Kegg colorectal cancer                 | 0.7504800         | 0.003975             |
| Kegg t cell receptor signaling pathway | 0.7415923         | 0.004543             |
| Kegg endometrial cancer                | 0.5724000         | 0.005300             |
| Kegg chronic myeloid leukemia          | 0.7068273         | 0.005963             |
| Kegg renal cell carcinoma              | 0.6024333         | 0.008480             |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 25 enriched KEGG pathways with ermineR.

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

There are 0 up-regulated and 0 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

pathway padj ES NES nMoreExtreme -------- ----- --- ---- -------------

Analysis of second-fourth days of differentiation 24h VS 96h (formation of definitive endoderm from APS)
--------------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE52158/DEgenes_24h_96h_52158.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
probes_to_genes <- genesym.probeid %>%
  filter(probe_id %in% rownames(DEgenes_24h_96h_52158))

topProbes <- DEgenes_24h_96h_52158 %>%
  rownames_to_column("probes") %>%
  filter(probes %in% probes_to_genes$probe_id)

probes_to_genes <- probes_to_genes %>% column_to_rownames('probe_id')
symbs <- probes_to_genes[topProbes$probes,]
topProbes$Symbol <- symbs

#for ermineR
ermineInputProbeScores <- topProbes %>% 
  #as.data.frame() %>%
  mutate(absolute_logFC = logFC) %>% 
  dplyr::select(Symbol, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) 

#randomly delete repeats (we will delete the lowerst ones)
repeats <- which(duplicated(ermineInputProbeScores$Symbol))
ermineInputProbeScores <- ermineInputProbeScores[-repeats,]

rownames(ermineInputProbeScores) <- NULL
ermineInputProbeScores <- ermineInputProbeScores %>% 
  column_to_rownames("Symbol")

#for fgsea
ermineInputProbeFGSEA <- topProbes %>% 
  #as.data.frame() %>%
  #mutate(absolute_logFC = logFC) %>% 
  dplyr::select(Symbol, logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(logFC)) 

#randomly delete repeats (we will delete the lowerst ones)
repeats <- which(duplicated(ermineInputProbeFGSEA$Symbol))
ermineInputProbeFGSEA <- ermineInputProbeFGSEA[-repeats,]

rownames(ermineInputProbeFGSEA) <- NULL
ermineInputProbeFGSEA <- ermineInputProbeFGSEA %>% 
  column_to_rownames("Symbol")

scoresFGSEA <- ermineInputProbeFGSEA$logFC
names(scoresFGSEA) <- rownames(ermineInputProbeFGSEA)
```

### Reactome

#### ErmineR Reactome pathways

``` r
enrichmentResultReactome <- precRecall(scores = ermineInputProbeScores,
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

| Pathway                                      | Corrected p-value | Corrected MF p-value |
|:---------------------------------------------|:------------------|:---------------------|
| Extracellular matrix organization            | 0.0000000         | 0.0000000            |
| Post-translational protein phosphorylation   | 0.0120200         | 0.0150250            |
| Lipoprotein metabolism                       | 0.0200333         | 0.0200333            |
| G alpha (i) signalling events                | 0.0000000         | 0.0300500            |
| Retinoid metabolism and transport            | 0.0400667         | 0.0343429            |
| Metabolism of fat-soluble vitamins           | 0.0150250         | 0.0400667            |
| Class A 1 (Rhodopsin-like receptors)         | 0.0429286         | 0.0480800            |
| Visual phototransduction                     | 0.0375625         | 0.0601000            |
| Response to elevated platelet cytosolic Ca2+ | 0.0701167         | 0.1001667            |
| Collagen formation                           | 0.0734556         | 0.1021700            |

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

There are 63 up-regulated and 203 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                              |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| ION CHANNEL TRANSPORT%REACTOME%R-HSA-983712.2                                                        |  0.0055461|  0.5201462|  1.803811|             0|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                 |  0.0055461|  0.7074541|  1.970932|             0|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                |  0.0055461|  0.4828371|  1.762239|             1|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                        |  0.0055461|  0.7175620|  2.169155|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                  |  0.0055461|  0.6277923|  1.929329|             2|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                            |  0.0055461|  0.5479179|  2.029385|             0|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                  |  0.0055461|  0.6206102|  2.100146|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                          |  0.0055461|  0.6227896|  2.004444|             0|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                   |  0.0055461|  0.6000689|  1.984038|             0|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                        |  0.0055461|  0.6293615|  2.113525|             0|
| OTHER SEMAPHORIN INTERACTIONS%REACTOME DATABASE ID RELEASE 65%416700                                 |  0.0055461|  0.7976115|  1.947158|             0|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                      |  0.0055461|  0.5714703|  1.889481|             0|
| ASPARAGINE N-LINKED GLYCOSYLATION%REACTOME%R-HSA-446203.4                                            |  0.0055461|  0.4407060|  1.655061|             0|
| UNFOLDED PROTEIN RESPONSE (UPR)%REACTOME%R-HSA-381119.2                                              |  0.0055461|  0.5726784|  1.903260|             0|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                             |  0.0055461|  0.5567924|  1.834257|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090 |  0.0055461|  0.6662090|  1.978666|             1|
| CHEMOKINE RECEPTORS BIND CHEMOKINES%REACTOME DATABASE ID RELEASE 65%380108                           |  0.0055461|  0.7870083|  1.947258|             0|
| PLASMA LIPOPROTEIN REMODELING%REACTOME DATABASE ID RELEASE 65%8963899                                |  0.0055461|  0.7841155|  1.940101|             1|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                          |  0.0055461|  0.5833557|  1.962734|             0|
| MUSCLE CONTRACTION%REACTOME%R-HSA-397014.2                                                           |  0.0055461|  0.5551219|  1.953554|             0|
| IRE1ALPHA ACTIVATES CHAPERONES%REACTOME DATABASE ID RELEASE 65%381070                                |  0.0055461|  0.6288664|  1.926380|             1|
| G ALPHA (I) SIGNALLING EVENTS%REACTOME%R-HSA-418594.5                                                |  0.0055726|  0.4431098|  1.628608|             3|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                                 |  0.0055726|  0.7055468|  1.910916|             2|
| CARDIAC CONDUCTION%REACTOME%R-HSA-5576891.2                                                          |  0.0057070|  0.5206248|  1.725309|             4|
| XBP1(S) ACTIVATES CHAPERONE GENES%REACTOME%R-HSA-381038.2                                            |  0.0064883|  0.6077636|  1.845439|             6|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                         |  0.0070965|  0.4829511|  1.676919|             8|
| TRANSPORT TO THE GOLGI AND SUBSEQUENT MODIFICATION%REACTOME DATABASE ID RELEASE 65%948021            |  0.0078771|  0.4459903|  1.595314|            10|
| COPI-MEDIATED ANTEROGRADE TRANSPORT%REACTOME%R-HSA-6807878.1                                         |  0.0084267|  0.5150550|  1.691331|            10|
| SLC-MEDIATED TRANSMEMBRANE TRANSPORT%REACTOME DATABASE ID RELEASE 65%425407                          |  0.0087229|  0.4449974|  1.594880|            12|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                           |  0.0090437|  0.4385312|  1.593253|            13|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                           |  0.0093042|  0.6077477|  1.799306|            11|
| VISUAL PHOTOTRANSDUCTION%REACTOME DATABASE ID RELEASE 65%2187338                                     |  0.0095335|  0.5667778|  1.777668|            12|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                           |  0.0097173|  0.6027405|  1.798111|            12|
| CELL JUNCTION ORGANIZATION%REACTOME DATABASE ID RELEASE 65%446728                                    |  0.0105986|  0.5416545|  1.731363|            14|
| RETINOID METABOLISM AND TRANSPORT%REACTOME%R-HSA-975634.2                                            |  0.0147154|  0.6228380|  1.769864|            19|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                      |  0.0154053|  0.5385592|  1.696954|            22|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                    |  0.0176494|  0.6866765|  1.757879|            22|
| ACTIVATION OF MATRIX METALLOPROTEINASES%REACTOME%R-HSA-1592389.1                                     |  0.0179813|  0.6786020|  1.769177|            23|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                        |  0.0188158|  0.6464456|  1.750845|            25|
| FORMATION OF FIBRIN CLOT (CLOTTING CASCADE)%REACTOME DATABASE ID RELEASE 65%140877                   |  0.0191330|  0.6767461|  1.752909|            25|
| ROS, RNS PRODUCTION IN PHAGOCYTES%REACTOME%R-HSA-1222556.7                                           |  0.0204284|  0.6732692|  1.755274|            27|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                           |  0.0204284|  0.6421321|  1.754313|            28|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                        |  0.0222150|  0.6352033|  1.735383|            31|
| SMOOTH MUSCLE CONTRACTION%REACTOME%R-HSA-445355.2                                                    |  0.0227784|  0.6134645|  1.729350|            33|
| METABOLISM OF FAT-SOLUBLE VITAMINS%REACTOME DATABASE ID RELEASE 65%6806667                           |  0.0238566|  0.5944045|  1.717292|            36|
| FATTY ACYL-COA BIOSYNTHESIS%REACTOME%R-HSA-75105.6                                                   |  0.0247304|  0.6307471|  1.723209|            36|
| ER TO GOLGI ANTEROGRADE TRANSPORT%REACTOME%R-HSA-199977.3                                            |  0.0248166|  0.4435518|  1.550071|            45|
| G ALPHA (Q) SIGNALLING EVENTS%REACTOME DATABASE ID RELEASE 65%416476                                 |  0.0290125|  0.4448053|  1.542537|            54|
| STRIATED MUSCLE CONTRACTION%REACTOME DATABASE ID RELEASE 65%390522                                   |  0.0290346|  0.6479715|  1.705468|            43|
| GASTRIN-CREB SIGNALLING PATHWAY VIA PKC AND MAPK%REACTOME DATABASE ID RELEASE 65%881907              |  0.0291867|  0.7069934|  1.725938|            42|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                              |  0.0313701|  0.7181541|  1.726732|            46|
| KERATINIZATION%REACTOME DATABASE ID RELEASE 65%6805567                                               |  0.0342130|  0.5667515|  1.670432|            58|
| FATTY ACID METABOLISM%REACTOME%R-HSA-8978868.4                                                       |  0.0360223|  0.4366320|  1.517136|            72|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                               |  0.0388493|  0.4622392|  1.542423|            76|
| CELL-CELL JUNCTION ORGANIZATION%REACTOME%R-HSA-421270.4                                              |  0.0416069|  0.5478410|  1.634333|            76|
| CELL-CELL COMMUNICATION%REACTOME%R-HSA-1500931.3                                                     |  0.0431677|  0.4575350|  1.532817|            89|
| AMINO ACID SYNTHESIS AND INTERCONVERSION (TRANSAMINATION)%REACTOME%R-HSA-70614.5                     |  0.0433857|  0.6241597|  1.660032|            75|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                      |  0.0433857|  0.5653624|  1.643097|            80|
| CRMPS IN SEMA3A SIGNALING%REACTOME DATABASE ID RELEASE 65%399956                                     |  0.0450334|  0.6996935|  1.682345|            74|
| SEMA3A PAK DEPENDENT AXON REPULSION%REACTOME%R-HSA-399954.1                                          |  0.0451245|  0.6863326|  1.675500|            75|
| STIMULI-SENSING CHANNELS%REACTOME DATABASE ID RELEASE 65%2672351                                     |  0.0454857|  0.5017315|  1.585090|            91|
| METABOLISM OF VITAMINS AND COFACTORS%REACTOME DATABASE ID RELEASE 65%196854                          |  0.0464411|  0.4102320|  1.461705|           106|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                            |  0.0471511|  0.5511011|  1.624305|            90|

### KEGG pathways

#### ErmineR KEGG pathways

``` r
enrichmentResultKEGG <- precRecall(scores = ermineInputProbeScores,
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
| Kegg ppar signaling pathway                  | 0.0318000         | 0.0265000            |
| Kegg focal adhesion                          | 0.0477000         | 0.0318000            |
| Kegg cytokine cytokine receptor interaction  | 0.0318000         | 0.0318000            |
| Kegg calcium signaling pathway               | 0.0954000         | 0.0914250            |
| Kegg axon guidance                           | 0.0993750         | 0.0927500            |
| Kegg ecm receptor interaction                | 0.1144800         | 0.0954000            |
| Kegg abc transporters                        | 0.1431000         | 0.1635429            |
| Kegg cell adhesion molecules cams            | 0.2420333         | 0.2208333            |
| Kegg leukocyte transendothelial migration    | 0.3564250         | 0.2226000            |
| Kegg neuroactive ligand receptor interaction | 0.1249286         | 0.2245875            |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 3 enriched KEGG pathways with ermineR.

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

There are 13 up-regulated and 7 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                           |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_BIOSYNTHESIS\_OF\_UNSATURATED\_FATTY\_ACIDS |  0.0052666|  0.8003185|  1.927542|             0|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION   |  0.0052666|  0.5215385|  1.828232|             0|
| KEGG\_ENDOCYTOSIS                                 |  0.0052666|  0.4596183|  1.645157|             0|
| KEGG\_FOCAL\_ADHESION                             |  0.0052666|  0.4848644|  1.754456|             0|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES      |  0.0052666|  0.7012962|  1.967556|             0|
| KEGG\_PPAR\_SIGNALING\_PATHWAY                    |  0.0065317|  0.6403443|  1.951384|             1|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                  |  0.0078784|  0.5758007|  1.838946|             3|
| KEGG\_LYSOSOME                                    |  0.0220970|  0.4912767|  1.681704|            14|
| KEGG\_AXON\_GUIDANCE                              |  0.0220970|  0.4830445|  1.661048|            14|
| KEGG\_LEUKOCYTE\_TRANSENDOTHELIAL\_MIGRATION      |  0.0220970|  0.5001548|  1.677134|            14|
| KEGG\_CALCIUM\_SIGNALING\_PATHWAY                 |  0.0314149|  0.4702313|  1.611592|            24|
| KEGG\_CELL\_ADHESION\_MOLECULES\_CAMS             |  0.0344575|  0.4799147|  1.605792|            29|
| KEGG\_VEGF\_SIGNALING\_PATHWAY                    |  0.0385645|  0.5198372|  1.637176|            33|

Analysis of whole differentiation process 0h VS 96 h (formation of Definitive endoderm from hESC)
-------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE52158/DEgenes_0h_96h_52158.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
probes_to_genes <- genesym.probeid %>%
  filter(probe_id %in% rownames(DEgenes_0h_96h_52158))

topProbes <- DEgenes_0h_96h_52158 %>%
  rownames_to_column("probes") %>%
  filter(probes %in% probes_to_genes$probe_id)

probes_to_genes <- probes_to_genes %>% column_to_rownames('probe_id')
symbs <- probes_to_genes[topProbes$probes,]
topProbes$Symbol <- symbs

#for ermineR
ermineInputProbeScores <- topProbes %>% 
  #as.data.frame() %>%
  mutate(absolute_logFC = logFC) %>% 
  dplyr::select(Symbol, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) 

#randomly delete repeats (we will delete the lowerst ones)
repeats <- which(duplicated(ermineInputProbeScores$Symbol))
ermineInputProbeScores <- ermineInputProbeScores[-repeats,]

rownames(ermineInputProbeScores) <- NULL
ermineInputProbeScores <- ermineInputProbeScores %>% 
  column_to_rownames("Symbol")

#for fgsea
ermineInputProbeFGSEA <- topProbes %>% 
  #as.data.frame() %>%
  #mutate(absolute_logFC = logFC) %>% 
  dplyr::select(Symbol, logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(logFC)) 

#randomly delete repeats (we will delete the lowerst ones)
repeats <- which(duplicated(ermineInputProbeFGSEA$Symbol))
ermineInputProbeFGSEA <- ermineInputProbeFGSEA[-repeats,]

rownames(ermineInputProbeFGSEA) <- NULL
ermineInputProbeFGSEA <- ermineInputProbeFGSEA %>% 
  column_to_rownames("Symbol")

scoresFGSEA <- ermineInputProbeFGSEA$logFC
names(scoresFGSEA) <- rownames(ermineInputProbeFGSEA)
```

### Reactome

#### ErmineR Reactome pathways

``` r
enrichmentResultReactome <- precRecall(scores = ermineInputProbeScores,
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

| Pathway                                                                                        | Corrected p-value | Corrected MF p-value |
|:-----------------------------------------------------------------------------------------------|:------------------|:---------------------|
| Extracellular matrix organization                                                              | 0.0000000         | 0.0000000            |
| Lipoprotein metabolism                                                                         | 0.0601000         | 0.0601000            |
| Metabolism of fat-soluble vitamins                                                             | 0.0801333         | 0.0801333            |
| Retinoid metabolism and transport                                                              | 0.1352250         | 0.1442400            |
| G alpha (i) signalling events                                                                  | 0.1442400         | 0.1459571            |
| Visual phototransduction                                                                       | 0.1116143         | 0.1602667            |
| Post-translational protein phosphorylation                                                     | 0.1202000         | 0.1803000            |
| Regulation of lipid metabolism by Peroxisome proliferator-activated receptor alpha (PPARalpha) | 0.1669444         | 0.1936556            |
| PPARA activates gene expression                                                                | 0.1878125         | 0.2178625            |
| ER to Golgi Anterograde Transport                                                              | 0.1966909         | 0.2404000            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 1 enriched Reactome pathways with ermineR

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

There are 18 up-regulated and 68 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| ER TO GOLGI ANTEROGRADE TRANSPORT%REACTOME%R-HSA-199977.3                                 |  0.0130739|  0.5217251|  1.753500|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                 |  0.0130739|  0.4519773|  1.613263|             2|
| XBP1(S) ACTIVATES CHAPERONE GENES%REACTOME%R-HSA-381038.2                                 |  0.0130739|  0.6528118|  1.908696|             0|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                             |  0.0130739|  0.5446356|  1.753682|             1|
| SIGNALING BY WNT IN CANCER%REACTOME DATABASE ID RELEASE 65%4791275                        |  0.0130739|  0.7041005|  1.899124|             2|
| OTHER SEMAPHORIN INTERACTIONS%REACTOME DATABASE ID RELEASE 65%416700                      |  0.0130739|  0.8076681|  1.885740|             2|
| ASPARAGINE N-LINKED GLYCOSYLATION%REACTOME%R-HSA-446203.4                                 |  0.0130739|  0.5025789|  1.823495|             0|
| UNFOLDED PROTEIN RESPONSE (UPR)%REACTOME%R-HSA-381119.2                                   |  0.0130739|  0.5978927|  1.904389|             1|
| ACTIVATION OF GENE EXPRESSION BY SREBF (SREBP)%REACTOME%R-HSA-2426168.2                   |  0.0130739|  0.6543364|  1.842286|             2|
| IRE1ALPHA ACTIVATES CHAPERONES%REACTOME DATABASE ID RELEASE 65%381070                     |  0.0130739|  0.6659925|  1.962536|             0|
| COPI-MEDIATED ANTEROGRADE TRANSPORT%REACTOME%R-HSA-6807878.1                              |  0.0143913|  0.5632953|  1.776232|             4|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                       |  0.0144425|  0.5221375|  1.693266|             5|
| TRANSPORT TO THE GOLGI AND SUBSEQUENT MODIFICATION%REACTOME DATABASE ID RELEASE 65%948021 |  0.0144425|  0.4763411|  1.642212|             5|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                             |  0.0148337|  0.6207620|  1.810275|             5|
| REGULATION OF CHOLESTEROL BIOSYNTHESIS BY SREBP (SREBF)%REACTOME%R-HSA-1655829.2          |  0.0208620|  0.5881657|  1.737277|             9|
| PEPTIDE HORMONE METABOLISM%REACTOME DATABASE ID RELEASE 65%2980736                        |  0.0283329|  0.5757659|  1.700651|            14|
| SIALIC ACID METABOLISM%REACTOME%R-HSA-4085001.1                                           |  0.0374209|  0.6632057|  1.751731|            19|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                           |  0.0392830|  0.5538747|  1.671497|            22|

### KEGG pathways

#### ErmineR KEGG pathways

``` r
enrichmentResultKEGG <- precRecall(scores = ermineInputProbeScores,
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

| Pathway                                | Corrected p-value | Corrected MF p-value |
|:---------------------------------------|:------------------|:---------------------|
| Kegg pathways in cancer                | 0.8702095         | 0.0000000            |
| Kegg focal adhesion                    | 0.5883000         | 0.0159000            |
| Kegg ppar signaling pathway            | 0.1033500         | 0.0212000            |
| Kegg wnt signaling pathway             | 0.4229400         | 0.0437250            |
| Kegg erbb signaling pathway            | 0.5992312         | 0.0454286            |
| Kegg axon guidance                     | 0.0795000         | 0.0530000            |
| Kegg neurotrophin signaling pathway    | 0.8101429         | 0.0604200            |
| Kegg chemokine signaling pathway       | 0.8720538         | 0.1033500            |
| Kegg t cell receptor signaling pathway | 0.8893400         | 0.1077667            |
| Kegg vegf signaling pathway            | 0.7506474         | 0.1176600            |

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

There are 2 up-regulated and 1 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway              |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------|----------:|----------:|---------:|-------------:|
| KEGG\_AXON\_GUIDANCE |  0.0210194|  0.5169020|  1.708050|             0|
| KEGG\_LYSOSOME       |  0.0212427|  0.5333283|  1.750172|             1|
