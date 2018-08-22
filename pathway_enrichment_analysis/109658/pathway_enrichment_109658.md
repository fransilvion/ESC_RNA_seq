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
| TCF dependent signaling in response to WNT                                           | 0.0000000         | 0.0000000            |
| Neuronal System                                                                      | 0.0000000         | 0.0000000            |
| G alpha (i) signalling events                                                        | 0.0000000         | 0.0000000            |
| Extracellular matrix organization                                                    | 0.0000000         | 0.0000000            |
| Class A 1 (Rhodopsin-like receptors)                                                 | 0.0000000         | 0.0000000            |
| Class B 2 (Secretin family receptors)                                                | 0.0102333         | 0.0102333            |
| WNT ligand biogenesis and trafficking                                                | 0.1159778         | 0.1304750            |
| Activation of anterior HOX genes in hindbrain development during early embryogenesis | 0.0789429         | 0.1403429            |
| Activation of HOX genes during differentiation                                       | 0.0789429         | 0.1403429            |
| FCERI mediated Ca+2 mobilization                                                     | 0.3670356         | 0.4284046            |

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

There are 8 up-regulated and 18 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                       |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                     |  0.0299798|  0.4847848|  1.777842|             0|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                                |  0.0299798|  0.9108422|  2.186894|             0|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                    |  0.0299798|  0.4782794|  1.708569|             2|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                               |  0.0332424|  0.6777577|  1.918210|             3|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                             |  0.0332424|  0.7189226|  1.891883|             2|
| ACTIVATION OF HOX GENES DURING DIFFERENTIATION%REACTOME DATABASE ID RELEASE 65%5619507                        |  0.0430515|  0.5440262|  1.740487|             7|
| ACTIVATION OF ANTERIOR HOX GENES IN HINDBRAIN DEVELOPMENT DURING EARLY EMBRYOGENESIS%REACTOME%R-HSA-5617472.2 |  0.0430515|  0.5440262|  1.740487|             7|
| SIGNALING BY NODAL%REACTOME DATABASE ID RELEASE 65%1181150                                                    |  0.0430515|  0.7644339|  1.835374|             6|

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
| Kegg melanogenesis                           | 0.0046860         | 0.000000             |
| Kegg focal adhesion                          | 0.1082400         | 0.000000             |
| Kegg cytokine cytokine receptor interaction  | 0.0000000         | 0.000000             |
| Kegg basal cell carcinoma                    | 0.0000000         | 0.000000             |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.002050             |
| Kegg hedgehog signaling pathway              | 0.0027330         | 0.002343             |
| Kegg melanoma                                | 0.2202286         | 0.008200             |
| Kegg cell adhesion molecules cams            | 0.0123000         | 0.008945             |

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
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY  |  0.0062604|  0.6907024|  2.017783|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY |  0.0062604|  0.6045600|  1.942306|             0|
| KEGG\_PATHWAYS\_IN\_CANCER          |  0.0062604|  0.4536458|  1.697974|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA        |  0.0062604|  0.6672720|  1.966419|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY       |  0.0108839|  0.4915260|  1.696861|             2|
| KEGG\_ECM\_RECEPTOR\_INTERACTION    |  0.0259253|  0.5443352|  1.727387|             7|
| KEGG\_MELANOGENESIS                 |  0.0259253|  0.5240401|  1.696434|             8|

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
| G alpha (i) signalling events                                                        | 0.0000000         | 0.0153500            |
| Neuronal System                                                                      | 0.0122800         | 0.0245600            |
| Class B 2 (Secretin family receptors)                                                | 0.0409333         | 0.0511667            |
| Incretin synthesis, secretion, and inactivation                                      | 0.1578857         | 0.1578857            |
| Insulin receptor signalling cascade                                                  | 0.2631429         | 0.2798698            |
| RA biosynthesis pathway                                                              | 0.2674311         | 0.2832773            |

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

There are 29 up-regulated and 194 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                              |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                            |  0.0047428|  0.5005085|  1.897386|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                            |  0.0047428|  0.6312265|  1.907763|             1|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                   |  0.0047428|  0.5369792|  1.796349|             1|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                          |  0.0048598|  0.5127644|  1.739873|             2|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                 |  0.0059739|  0.6303752|  1.860607|             3|
| DOPAMINE NEUROTRANSMITTER RELEASE CYCLE%REACTOME%R-HSA-212676.2                                      |  0.0088423|  0.7125226|  1.800558|             5|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                  |  0.0090735|  0.5624466|  1.774681|             6|
| G ALPHA (I) SIGNALLING EVENTS%REACTOME%R-HSA-418594.5                                                |  0.0094940|  0.4089427|  1.531389|             8|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                        |  0.0098017|  0.5745773|  1.779860|             7|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                          |  0.0115527|  0.5184447|  1.699251|             9|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                        |  0.0145730|  0.4834313|  1.650374|            13|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                      |  0.0192116|  0.5962055|  1.745156|            17|
| SURFACTANT METABOLISM%REACTOME DATABASE ID RELEASE 65%5683826                                        |  0.0215691|  0.7092872|  1.749333|            19|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                           |  0.0231696|  0.5511119|  1.683307|            23|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                  |  0.0246759|  0.4605667|  1.587514|            28|
| OTHER SEMAPHORIN INTERACTIONS%REACTOME DATABASE ID RELEASE 65%416700                                 |  0.0246759|  0.7196216|  1.754938|            24|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                           |  0.0246759|  0.5342714|  1.670249|            26|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                |  0.0248400|  0.3932485|  1.470440|            34|
| REGULATION OF COMPLEMENT CASCADE%REACTOME%R-HSA-977606.5                                             |  0.0298537|  0.7146985|  1.716180|            32|
| EFFECTS OF PIP2 HYDROLYSIS%REACTOME%R-HSA-114508.2                                                   |  0.0315358|  0.6320716|  1.696014|            37|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                                 |  0.0331608|  0.6136479|  1.688105|            40|
| COMPLEMENT CASCADE%REACTOME DATABASE ID RELEASE 65%166658                                            |  0.0372545|  0.6993303|  1.705453|            44|
| ION CHANNEL TRANSPORT%REACTOME%R-HSA-983712.2                                                        |  0.0380146|  0.4166411|  1.497065|            58|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                             |  0.0387692|  0.6159588|  1.666500|            49|
| KERATAN SULFATE KERATIN METABOLISM%REACTOME%R-HSA-1638074.1                                          |  0.0412040|  0.6146049|  1.662837|            54|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                        |  0.0412040|  0.6142315|  1.661826|            54|
| PLASMA LIPOPROTEIN REMODELING%REACTOME DATABASE ID RELEASE 65%8963899                                |  0.0414479|  0.6633448|  1.656492|            53|
| CELL-CELL JUNCTION ORGANIZATION%REACTOME%R-HSA-421270.4                                              |  0.0433506|  0.5041505|  1.581430|            63|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090 |  0.0473947|  0.5122667|  1.570330|            70|

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
| Kegg pathways in cancer                      | 0.0205000         | 0.0000000            |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.0000000            |
| Kegg cytokine cytokine receptor interaction  | 0.0000000         | 0.0000000            |
| Kegg hematopoietic cell lineage              | 0.0273333         | 0.0123000            |
| Kegg wnt signaling pathway                   | 0.0300667         | 0.0229600            |
| Kegg cell adhesion molecules cams            | 0.0328000         | 0.0273333            |
| Kegg ecm receptor interaction                | 0.0585714         | 0.0374857            |
| Kegg melanogenesis                           | 0.1731111         | 0.0615000            |
| Kegg focal adhesion                          | 0.5014947         | 0.0639600            |
| Kegg maturity onset diabetes of the young    | 0.0820000         | 0.0710667            |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 7 enriched KEGG pathways with ermineR.

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
| KEGG\_ECM\_RECEPTOR\_INTERACTION                                    |  0.0208967|  0.5578373|  1.824031|             1|
| KEGG\_LEUKOCYTE\_TRANSENDOTHELIAL\_MIGRATION                        |  0.0208967|  0.5034342|  1.727928|             0|
| KEGG\_AXON\_GUIDANCE                                                |  0.0230061|  0.4619294|  1.638762|             4|
| KEGG\_CELL\_ADHESION\_MOLECULES\_CAMS                               |  0.0230061|  0.4921250|  1.689111|             4|
| KEGG\_GLYCOSPHINGOLIPID\_BIOSYNTHESIS\_LACTO\_AND\_NEOLACTO\_SERIES |  0.0233728|  0.6982301|  1.809485|             5|
| KEGG\_FOCAL\_ADHESION                                               |  0.0294118|  0.4212441|  1.559089|            13|

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
| Neuronal System                                 | 0.0368400         | 0.0736800            |
| G alpha (q) signalling events                   | 0.1330333         | 0.2046667            |
| Amine ligand-binding receptors                  | 0.2543714         | 0.2543714            |
| Post-translational protein phosphorylation      | 0.3792353         | 0.4029375            |
| Amyloid fiber formation                         | 0.3929600         | 0.4045176            |
| Signaling by BMP                                | 0.3815571         | 0.4156308            |

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

There are 36 up-regulated and 221 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                         |       padj|         ES|       NES|  nMoreExtreme|
|:----------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                            |  0.0046050|  0.6778319|  2.015432|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                             |  0.0046050|  0.6026327|  1.903707|             1|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                       |  0.0046050|  0.5498580|  2.053307|             0|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                                  |  0.0046050|  0.8103541|  2.036106|             0|
| CLASS B 2 (SECRETIN FAMILY RECEPTORS)%REACTOME DATABASE ID RELEASE 65%373080                                    |  0.0046050|  0.6183446|  1.939815|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                                     |  0.0046050|  0.5818259|  1.903435|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                       |  0.0046050|  0.6573691|  1.991448|             0|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                                   |  0.0046050|  0.6976477|  1.916749|             2|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                 |  0.0046050|  0.7379812|  2.171155|             0|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                      |  0.0046050|  0.6334096|  1.979161|             0|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                               |  0.0046050|  0.7410476|  2.035988|             0|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                      |  0.0046050|  0.4563494|  1.667382|             2|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090            |  0.0046050|  0.6042254|  1.856775|             1|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                                        |  0.0052759|  0.6841620|  1.879698|             4|
| KERATAN SULFATE BIOSYNTHESIS%REACTOME%R-HSA-2022854.1                                                           |  0.0058308|  0.6994094|  1.835025|             5|
| KERATAN SULFATE KERATIN METABOLISM%REACTOME%R-HSA-1638074.1                                                     |  0.0077595|  0.6564317|  1.803510|             8|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                                          |  0.0098561|  0.4788939|  1.619831|            14|
| PEPTIDE HORMONE METABOLISM%REACTOME DATABASE ID RELEASE 65%2980736                                              |  0.0102178|  0.5259866|  1.672962|            14|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                                   |  0.0129543|  0.4696021|  1.589870|            20|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                                   |  0.0142999|  0.6232640|  1.724699|            19|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                              |  0.0159948|  0.4785204|  1.590700|            26|
| SIGNALING BY NODAL%REACTOME DATABASE ID RELEASE 65%1181150                                                      |  0.0183391|  0.6946138|  1.745295|            25|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                                     |  0.0197497|  0.4638452|  1.561378|            34|
| INCRETIN SYNTHESIS, SECRETION, AND INACTIVATION%REACTOME%R-HSA-400508.1                                         |  0.0218723|  0.7041486|  1.712514|            31|
| INTERLEUKIN-6 FAMILY SIGNALING%REACTOME DATABASE ID RELEASE 65%6783589                                          |  0.0264570|  0.6508816|  1.695255|            41|
| SURFACTANT METABOLISM%REACTOME DATABASE ID RELEASE 65%5683826                                                   |  0.0264570|  0.6749192|  1.695810|            40|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                               |  0.0281024|  0.6499894|  1.692931|            44|
| SIGNALING BY PDGF%REACTOME DATABASE ID RELEASE 65%186797                                                        |  0.0281024|  0.5212166|  1.601690|            49|
| IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON-LYMPHOID CELL%REACTOME DATABASE ID RELEASE 65%198933 |  0.0289037|  0.5408741|  1.622170|            50|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                                    |  0.0292596|  0.4381799|  1.498477|            58|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                                 |  0.0339794|  0.4603672|  1.527078|            67|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                                         |  0.0388746|  0.6753657|  1.642513|            64|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                             |  0.0449926|  0.4313929|  1.470635|            96|
| ACTIVATION OF MATRIX METALLOPROTEINASES%REACTOME%R-HSA-1592389.1                                                |  0.0452303|  0.6523676|  1.639147|            78|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                                      |  0.0474127|  0.5103004|  1.563589|            92|
| REGULATION OF BETA-CELL DEVELOPMENT%REACTOME DATABASE ID RELEASE 65%186712                                      |  0.0479852|  0.5738164|  1.587867|            87|

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
| Kegg pathways in cancer                      | 0.0295200         | 0.0000000            |
| Kegg neuroactive ligand receptor interaction | 0.0000000         | 0.0000000            |
| Kegg hematopoietic cell lineage              | 0.0109333         | 0.0000000            |
| Kegg focal adhesion                          | 0.0445143         | 0.0000000            |
| Kegg cytokine cytokine receptor interaction  | 0.0000000         | 0.0000000            |
| Kegg ecm receptor interaction                | 0.0123000         | 0.0054670            |
| Kegg cell adhesion molecules cams            | 0.0273333         | 0.0234286            |
| Kegg tgf beta signaling pathway              | 0.0430500         | 0.0389500            |
| Kegg calcium signaling pathway               | 0.2678667         | 0.0728889            |
| Kegg systemic lupus erythematosus            | 0.2200333         | 0.1968000            |

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
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY              |  0.0037811|  0.6493254|  1.954695|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY             |  0.0037811|  0.5752288|  1.909323|             0|
| KEGG\_AXON\_GUIDANCE                            |  0.0037811|  0.5025458|  1.768550|             0|
| KEGG\_FOCAL\_ADHESION                           |  0.0037811|  0.4824189|  1.765112|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                |  0.0037811|  0.6018477|  1.978659|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA                    |  0.0037811|  0.6480634|  1.969615|             0|
| KEGG\_HEMATOPOIETIC\_CELL\_LINEAGE              |  0.0067771|  0.6173100|  1.867366|             2|
| KEGG\_PATHWAYS\_IN\_CANCER                      |  0.0067771|  0.4193480|  1.591134|             2|
| KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION |  0.0194278|  0.4471618|  1.573644|            14|
| KEGG\_WNT\_SIGNALING\_PATHWAY                   |  0.0216472|  0.4376114|  1.552012|            17|
| KEGG\_LYSOSOME                                  |  0.0417066|  0.4384520|  1.521735|            37|
