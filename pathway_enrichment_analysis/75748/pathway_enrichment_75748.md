Pathway enrichment analysis, GEO75748
================
German Novakovskiy
August 21, 2018

Analysis of first day of differentiation 0h VS 24 h (formation of anterior Primitive streak, mesodendoderm)
-----------------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE75748/DEgenes_0h_24h_75748.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
#absolute Log fold changes for ermineR
ermineInputGeneScores <- DEgenes_0h_24h_75748 %>% 
  rownames_to_column("gene") %>%
  mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

#exact log fold changes for fgsea
# scores forfgsea
scoresFGSEADF <- DEgenes_0h_24h_75748 %>% 
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
    ## toType = "ENTREZID", : 6.27% of input gene IDs are fail to map...

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
| TCF dependent signaling in response to WNT                                           | 0.0000000         | 0.1238000            |
| Processing of DNA double-strand break ends                                           | 0.1680143         | 0.2259350            |
| Meiosis                                                                              | 0.1711353         | 0.2285538            |
| Neurotransmitter release cycle                                                       | 0.1669860         | 0.2291973            |
| Transcriptional regulation by small RNAs                                             | 0.1658920         | 0.2294829            |
| Activation of anterior HOX genes in hindbrain development during early embryogenesis | 0.1709619         | 0.2328619            |
| Activation of HOX genes during differentiation                                       | 0.1709619         | 0.2328619            |
| NCAM1 interactions                                                                   | 0.1776261         | 0.2329395            |
| RUNX1 regulates transcription of genes involved in differentiation of HSCs           | 0.1987316         | 0.2338444            |

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

There are 87 up-regulated and 51 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                 |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                         |  0.0066989|  0.5294765|  2.084714|             0|
| FACTORS INVOLVED IN MEGAKARYOCYTE DEVELOPMENT AND PLATELET PRODUCTION%REACTOME%R-HSA-983231.2           |  0.0066989|  0.3952818|  1.758705|             0|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                    |  0.0066989|  0.7081415|  2.535219|             0|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                            |  0.0066989|  0.5435469|  2.030855|             0|
| NRAGE SIGNALS DEATH THROUGH JNK%REACTOME DATABASE ID RELEASE 65%193648                                  |  0.0066989|  0.5352238|  2.049029|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                     |  0.0066989|  0.6786716|  2.617607|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                               |  0.0066989|  0.4685203|  2.271531|             0|
| SIGNALING BY NOTCH3%REACTOME DATABASE ID RELEASE 65%9012852                                             |  0.0066989|  0.5592139|  2.082210|             0|
| EPH-EPHRIN MEDIATED REPULSION OF CELLS%REACTOME%R-HSA-3928665.3                                         |  0.0066989|  0.5557962|  2.087876|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                             |  0.0066989|  0.6177198|  2.491573|             0|
| EPH-EPHRIN SIGNALING%REACTOME%R-HSA-2682334.1                                                           |  0.0066989|  0.4808375|  2.029657|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                               |  0.0066989|  0.5624801|  2.069257|             0|
| SIGNALING BY PDGF%REACTOME DATABASE ID RELEASE 65%186797                                                |  0.0066989|  0.6334639|  2.391827|             0|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                      |  0.0066989|  0.4564114|  1.905061|             0|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                              |  0.0066989|  0.5225696|  2.000584|             0|
| RHO GTPASE CYCLE%REACTOME DATABASE ID RELEASE 65%194840                                                 |  0.0066989|  0.4909874|  2.200047|             0|
| SIGNALING BY NTRK1 (TRKA)%REACTOME%R-HSA-187037.2                                                       |  0.0066989|  0.4430739|  1.817281|             0|
| SIGNALING BY RAS MUTANTS%REACTOME%R-HSA-6802949.1                                                       |  0.0066989|  0.5381687|  2.003849|             0|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                            |  0.0066989|  0.4468272|  1.929967|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090    |  0.0066989|  0.6219693|  2.315878|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                             |  0.0066989|  0.4546404|  1.924573|             0|
| CELLULAR SENESCENCE%REACTOME%R-HSA-2559583.2                                                            |  0.0098541|  0.3550414|  1.668533|             1|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                              |  0.0098541|  0.5877192|  1.985894|             1|
| SEMA4D IN SEMAPHORIN SIGNALING%REACTOME%R-HSA-400685.2                                                  |  0.0098541|  0.6456456|  2.039723|             1|
| ONCOGENIC MAPK SIGNALING%REACTOME%R-HSA-6802957.2                                                       |  0.0098541|  0.4635441|  1.869705|             1|
| SIGNALING BY WNT%REACTOME DATABASE ID RELEASE 65%195721                                                 |  0.0098541|  0.3292186|  1.639566|             1|
| NON-INTEGRIN MEMBRANE-ECM INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000171                          |  0.0098541|  0.5647645|  2.021915|             1|
| SIGNALING BY FGFR1%REACTOME DATABASE ID RELEASE 65%5654736                                              |  0.0121085|  0.5538595|  1.954075|             2|
| SIGNALING BY NTRKS%REACTOME DATABASE ID RELEASE 65%166520                                               |  0.0121085|  0.4226017|  1.807945|             2|
| NOTCH3 INTRACELLULAR DOMAIN REGULATES TRANSCRIPTION%REACTOME%R-HSA-9013508.1                            |  0.0121085|  0.6446422|  2.012513|             2|
| DOWNSTREAM SIGNALING OF ACTIVATED FGFR1%REACTOME%R-HSA-5654687.2                                        |  0.0121085|  0.6662548|  2.005703|             2|
| HS-GAG BIOSYNTHESIS%REACTOME%R-HSA-2022928.1                                                            |  0.0121085|  0.5863204|  1.964617|             2|
| O-LINKED GLYCOSYLATION%REACTOME%R-HSA-5173105.3                                                         |  0.0134666|  0.4268245|  1.777131|             3|
| GLYCOSAMINOGLYCAN METABOLISM%REACTOME DATABASE ID RELEASE 65%1630316                                    |  0.0134666|  0.3918191|  1.717710|             3|
| SEMA4D INDUCED CELL MIGRATION AND GROWTH-CONE COLLAPSE%REACTOME%R-HSA-416572.1                          |  0.0134666|  0.6602435|  1.987606|             3|
| DOWNSTREAM SIGNAL TRANSDUCTION%REACTOME DATABASE ID RELEASE 65%186763                                   |  0.0134666|  0.5913993|  1.942946|             3|
| MET PROMOTES CELL MOTILITY%REACTOME DATABASE ID RELEASE 65%8875878                                      |  0.0134666|  0.6095724|  2.002650|             3|
| PRE-NOTCH PROCESSING IN GOLGI%REACTOME%R-HSA-1912420.2                                                  |  0.0134666|  0.6905442|  1.984814|             3|
| PRE-NOTCH EXPRESSION AND PROCESSING%REACTOME DATABASE ID RELEASE 65%1912422                             |  0.0180566|  0.4081292|  1.742498|             5|
| G ALPHA (12 13) SIGNALLING EVENTS%REACTOME%R-HSA-416482.3                                               |  0.0181121|  0.4293278|  1.766800|             5|
| CHONDROITIN SULFATE DERMATAN SULFATE METABOLISM%REACTOME DATABASE ID RELEASE 65%1793185                 |  0.0189174|  0.5015270|  1.856281|             6|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                                    |  0.0189174|  0.5829550|  1.915203|             6|
| SIGNALING BY VEGF%REACTOME DATABASE ID RELEASE 65%194138                                                |  0.0189174|  0.4062985|  1.734682|             6|
| N-GLYCAN TRIMMING IN THE ER AND CALNEXIN CALRETICULIN CYCLE%REACTOME%R-HSA-532668.2                     |  0.0189174|  0.5402788|  1.892492|             6|
| LDL CLEARANCE%REACTOME%R-HSA-8964038.1                                                                  |  0.0189174|  0.6550377|  1.946366|             6|
| INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING%REACTOME DATABASE ID RELEASE 65%6785807                      |  0.0205386|  0.4048856|  1.703058|             7|
| SIGNALING BY TGF-BETA RECEPTOR COMPLEX%REACTOME DATABASE ID RELEASE 65%170834                           |  0.0205386|  0.4227119|  1.726559|             7|
| TCF DEPENDENT SIGNALING IN RESPONSE TO WNT%REACTOME%R-HSA-201681.1                                      |  0.0206049|  0.3211413|  1.539132|             8|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                          |  0.0206049|  0.6638013|  1.907947|             7|
| CELL-CELL COMMUNICATION%REACTOME%R-HSA-1500931.3                                                        |  0.0209176|  0.3838493|  1.671360|             8|
| PARADOXICAL ACTIVATION OF RAF SIGNALING BY KINASE INACTIVE BRAF%REACTOME DATABASE ID RELEASE 65%6802955 |  0.0209966|  0.5317474|  1.849929|             8|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                       |  0.0209966|  0.6333236|  1.906566|             8|
| SIGNALING BY MODERATE KINASE ACTIVITY BRAF MUTANTS%REACTOME DATABASE ID RELEASE 65%6802946              |  0.0209966|  0.5317474|  1.849929|             8|
| TRANSCRIPTIONAL REGULATION BY THE AP-2 (TFAP2) FAMILY OF TRANSCRIPTION FACTORS%REACTOME%R-HSA-8864260.3 |  0.0209966|  0.5613401|  1.880914|             8|
| ACTIVATION OF BH3-ONLY PROTEINS%REACTOME DATABASE ID RELEASE 65%114452                                  |  0.0240370|  0.5633874|  1.873514|            10|
| CHROMATIN ORGANIZATION%REACTOME%R-HSA-4839726.2                                                         |  0.0240370|  0.3022272|  1.475508|            11|
| CHROMATIN MODIFYING ENZYMES%REACTOME%R-HSA-3247509.4                                                    |  0.0240370|  0.3022272|  1.475508|            11|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                         |  0.0259521|  0.5166638|  1.849709|            11|
| VEGFA-VEGFR2 PATHWAY%REACTOME%R-HSA-4420097.3                                                           |  0.0271232|  0.3945066|  1.656764|            12|
| INTRACELLULAR SIGNALING BY SECOND MESSENGERS%REACTOME DATABASE ID RELEASE 65%9006925                    |  0.0271232|  0.2992943|  1.460553|            13|
| CLATHRIN-MEDIATED ENDOCYTOSIS%REACTOME%R-HSA-8856828.3                                                  |  0.0278848|  0.3553193|  1.592137|            13|
| L1CAM INTERACTIONS%REACTOME%R-HSA-373760.2                                                              |  0.0296804|  0.3913974|  1.649032|            14|
| DISASSEMBLY OF THE DESTRUCTION COMPLEX AND RECRUITMENT OF AXIN TO THE MEMBRANE%REACTOME%R-HSA-4641262.3 |  0.0299310|  0.5488056|  1.825023|            14|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                           |  0.0311943|  0.5607128|  1.825402|            15|
| LAMININ INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000157                                            |  0.0323928|  0.5939494|  1.832468|            16|
| HEPARAN SULFATE HEPARIN (HS-GAG) METABOLISM%REACTOME%R-HSA-1638091.1                                    |  0.0323928|  0.4529951|  1.717350|            17|
| REGULATION OF PTEN GENE TRANSCRIPTION%REACTOME DATABASE ID RELEASE 65%8943724                           |  0.0323928|  0.4358255|  1.686612|            17|
| INTRINSIC PATHWAY FOR APOPTOSIS%REACTOME%R-HSA-109606.2                                                 |  0.0325456|  0.4877898|  1.777215|            17|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                         |  0.0340023|  0.3887661|  1.613920|            19|
| TP53 REGULATES TRANSCRIPTION OF CELL CYCLE GENES%REACTOME DATABASE ID RELEASE 65%6791312                |  0.0354101|  0.4637610|  1.732751|            20|
| TGF-BETA RECEPTOR SIGNALING ACTIVATES SMADS%REACTOME DATABASE ID RELEASE 65%2173789                     |  0.0369865|  0.5199507|  1.769636|            21|
| PRE-NOTCH TRANSCRIPTION AND TRANSLATION%REACTOME%R-HSA-1912408.3                                        |  0.0376205|  0.3876789|  1.604907|            22|
| DAG AND IP3 SIGNALING%REACTOME DATABASE ID RELEASE 65%1489509                                           |  0.0379126|  0.5378068|  1.788447|            22|
| CALMODULIN INDUCED EVENTS%REACTOME DATABASE ID RELEASE 65%111933                                        |  0.0386051|  0.5727439|  1.788053|            23|
| CAM PATHWAY%REACTOME DATABASE ID RELEASE 65%111997                                                      |  0.0386051|  0.5727439|  1.788053|            23|
| ASPARAGINE N-LINKED GLYCOSYLATION%REACTOME%R-HSA-446203.4                                               |  0.0386051|  0.2836056|  1.400766|            25|
| SIGNALING BY BRAF AND RAF FUSIONS%REACTOME%R-HSA-6802952.1                                              |  0.0419388|  0.4273780|  1.648378|            26|
| REGULATION OF FZD BY UBIQUITINATION%REACTOME%R-HSA-4641263.2                                            |  0.0421770|  0.6037169|  1.817438|            27|
| SIGNALING BY FGFR3 IN DISEASE%REACTOME DATABASE ID RELEASE 65%5655332                                   |  0.0421770|  0.6282799|  1.805849|            26|
| SIGNALING BY FGFR3 POINT MUTANTS IN CANCER%REACTOME%R-HSA-8853338.2                                     |  0.0421770|  0.6282799|  1.805849|            26|
| O-GLYCOSYLATION OF TSR DOMAIN-CONTAINING PROTEINS%REACTOME%R-HSA-5173214.1                              |  0.0450314|  0.5090161|  1.732420|            31|
| HDACS DEACETYLATE HISTONES%REACTOME%R-HSA-3214815.2                                                     |  0.0450314|  0.3765636|  1.571776|            31|
| CA-DEPENDENT EVENTS%REACTOME DATABASE ID RELEASE 65%111996                                              |  0.0450314|  0.5519660|  1.770250|            31|
| CELL DEATH SIGNALLING VIA NRAGE, NRIF AND NADE%REACTOME DATABASE ID RELEASE 65%204998                   |  0.0450314|  0.3994354|  1.615986|            31|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                              |  0.0450314|  0.4581408|  1.695698|            31|
| ESR-MEDIATED SIGNALING%REACTOME DATABASE ID RELEASE 65%8939211                                          |  0.0461225|  0.3294227|  1.491255|            34|
| INTRA-GOLGI AND RETROGRADE GOLGI-TO-ER TRAFFIC%REACTOME DATABASE ID RELEASE 65%6811442                  |  0.0472566|  0.3121420|  1.465565|            36|

``` r
upPathwaysReactome75748_024 <- upPathwaysReactome
save(upPathwaysReactome75748_024, file="upPathwaysReactome75748_024.Rdata")
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
| Kegg pathways in cancer                                   | 0.0364500         | 0.000000             |
| Kegg mapk signaling pathway                               | 0.0081000         | 0.000000             |
| Kegg focal adhesion                                       | 0.0324000         | 0.000000             |
| Kegg basal cell carcinoma                                 | 0.0356400         | 0.006480             |
| Kegg melanogenesis                                        | 0.1332000         | 0.043200             |
| Kegg alzheimers disease                                   | 0.1634727         | 0.048600             |
| Kegg arrhythmogenic right ventricular cardiomyopathy arvc | 0.2106000         | 0.075600             |
| Kegg ecm receptor interaction                             | 0.1593000         | 0.076140             |
| Kegg parkinsons disease                                   | 0.1411714         | 0.078975             |

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

There are 37 up-regulated and 11 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                                        |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                                 |  0.0017565|  0.4056368|  1.975072|             0|
| KEGG\_ENDOCYTOSIS                                              |  0.0017565|  0.4037747|  1.888157|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY                                  |  0.0017565|  0.4299674|  1.938473|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY                            |  0.0017565|  0.4547765|  1.874634|             0|
| KEGG\_AXON\_GUIDANCE                                           |  0.0017565|  0.4493756|  2.005834|             0|
| KEGG\_FOCAL\_ADHESION                                          |  0.0017565|  0.5544194|  2.607923|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                               |  0.0017565|  0.6108126|  2.453503|             0|
| KEGG\_NEUROTROPHIN\_SIGNALING\_PATHWAY                         |  0.0017565|  0.4347978|  1.926153|             0|
| KEGG\_REGULATION\_OF\_ACTIN\_CYTOSKELETON                      |  0.0017565|  0.3542808|  1.666754|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                                     |  0.0017565|  0.3939500|  1.973200|             0|
| KEGG\_COLORECTAL\_CANCER                                       |  0.0017565|  0.5037259|  1.986889|             0|
| KEGG\_RENAL\_CELL\_CARCINOMA                                   |  0.0017565|  0.4860469|  1.958704|             0|
| KEGG\_PANCREATIC\_CANCER                                       |  0.0017565|  0.4877522|  1.972585|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA                                   |  0.0017565|  0.5864775|  2.147247|             0|
| KEGG\_CHRONIC\_MYELOID\_LEUKEMIA                               |  0.0017565|  0.4665465|  1.900659|             0|
| KEGG\_SMALL\_CELL\_LUNG\_CANCER                                |  0.0017565|  0.4679651|  1.955322|             0|
| KEGG\_ARRHYTHMOGENIC\_RIGHT\_VENTRICULAR\_CARDIOMYOPATHY\_ARVC |  0.0017565|  0.5006561|  1.970281|             0|
| KEGG\_ERBB\_SIGNALING\_PATHWAY                                 |  0.0024384|  0.4317931|  1.799550|             1|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY                             |  0.0024384|  0.5263141|  1.915953|             1|
| KEGG\_INSULIN\_SIGNALING\_PATHWAY                              |  0.0024384|  0.4063390|  1.813736|             1|
| KEGG\_ENDOMETRIAL\_CANCER                                      |  0.0024384|  0.4891885|  1.859069|             1|
| KEGG\_NON\_SMALL\_CELL\_LUNG\_CANCER                           |  0.0024384|  0.4732919|  1.823101|             1|
| KEGG\_DILATED\_CARDIOMYOPATHY                                  |  0.0024384|  0.4775598|  1.924502|             1|
| KEGG\_CHEMOKINE\_SIGNALING\_PATHWAY                            |  0.0043852|  0.3725686|  1.671777|             3|
| KEGG\_MELANOGENESIS                                            |  0.0043852|  0.4276664|  1.779571|             3|
| KEGG\_BLADDER\_CANCER                                          |  0.0043852|  0.5453271|  1.944315|             3|
| KEGG\_NOTCH\_SIGNALING\_PATHWAY                                |  0.0061188|  0.4997174|  1.829595|             5|
| KEGG\_GNRH\_SIGNALING\_PATHWAY                                 |  0.0081563|  0.4131290|  1.719079|             8|
| KEGG\_HYPERTROPHIC\_CARDIOMYOPATHY\_HCM                        |  0.0089098|  0.4366952|  1.735859|             9|
| KEGG\_GAP\_JUNCTION                                            |  0.0142923|  0.3969463|  1.640064|            16|
| KEGG\_ACUTE\_MYELOID\_LEUKEMIA                                 |  0.0153985|  0.4398111|  1.684903|            18|
| KEGG\_TYPE\_II\_DIABETES\_MELLITUS                             |  0.0190680|  0.4709110|  1.703613|            23|
| KEGG\_GLYCOSAMINOGLYCAN\_BIOSYNTHESIS\_CHONDROITIN\_SULFATE    |  0.0224392|  0.5813376|  1.771199|            28|
| KEGG\_VASOPRESSIN\_REGULATED\_WATER\_REABSORPTION              |  0.0257368|  0.4731124|  1.686840|            33|
| KEGG\_B\_CELL\_RECEPTOR\_SIGNALING\_PATHWAY                    |  0.0269205|  0.3984707|  1.600571|            36|
| KEGG\_ADHERENS\_JUNCTION                                       |  0.0290692|  0.3871865|  1.577355|            40|
| KEGG\_PROSTATE\_CANCER                                         |  0.0316893|  0.3672851|  1.540999|            45|

``` r
upPathwaysKEGG75748_024 <- upPathwaysKEGG
save(upPathwaysKEGG75748_024, file="upPathwaysKEGG75748_024.Rdata")
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
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       | 0.000000          | 0.000000             |
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            | 0.000000          | 0.000000             |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  | 0.009900          | 0.009900             |
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens | 0.009900          | 0.017820             |
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                | 0.007425          | 0.019800             |
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens        | 0.011880          | 0.022275             |

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 6 enriched Wiki pathways with ermineR

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

There are 107 up-regulated and 6 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                                         |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| B Cell Receptor Signaling Pathway%WikiPathways\_20180810%WP23%Homo sapiens                                                      |  0.0011922|  0.4320285|  1.812734|             0|
| EGF/EGFR Signaling Pathway%WikiPathways\_20180810%WP437%Homo sapiens                                                            |  0.0011922|  0.3903432|  1.794846|             0|
| Angiopoietin Like Protein 8 Regulatory Pathway%WikiPathways\_20180810%WP3915%Homo sapiens                                       |  0.0011922|  0.4055949|  1.802509|             0|
| TGF-beta Receptor Signaling%WikiPathways\_20180810%WP560%Homo sapiens                                                           |  0.0011922|  0.5486288|  2.064472|             0|
| Integrated Breast Cancer Pathway%WikiPathways\_20180810%WP1984%Homo sapiens                                                     |  0.0011922|  0.4025444|  1.821527|             0|
| TGF-B Signaling in Thyroid Cells for Epithelial-Mesenchymal Transition%WikiPathways\_20180810%WP3859%Homo sapiens               |  0.0011922|  0.7395393|  2.139843|             0|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                                        |  0.0011922|  0.5363231|  2.502556|             0|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                                      |  0.0011922|  0.5407072|  2.138815|             0|
| Wnt/beta-catenin Signaling Pathway in Leukemia%WikiPathways\_20180810%WP3658%Homo sapiens                                       |  0.0011922|  0.6818005|  2.116480|             0|
| Hepatitis C and Hepatocellular Carcinoma%WikiPathways\_20180810%WP3646%Homo sapiens                                             |  0.0011922|  0.5156964|  1.904841|             0|
| miR-509-3p alteration of YAP1/ECM axis%WikiPathways\_20180810%WP3967%Homo sapiens                                               |  0.0011922|  0.7271337|  2.103947|             0|
| TGF-beta Signaling Pathway%WikiPathways\_20180810%WP366%Homo sapiens                                                            |  0.0011922|  0.4403416|  1.963861|             0|
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                                         |  0.0011922|  0.5582691|  2.337051|             0|
| Arrhythmogenic Right Ventricular Cardiomyopathy%WikiPathways\_20180810%WP2118%Homo sapiens                                      |  0.0011922|  0.4991551|  1.954944|             0|
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens                         |  0.0011922|  0.5114115|  2.022933|             0|
| Leptin signaling pathway%WikiPathways\_20180810%WP2034%Homo sapiens                                                             |  0.0011922|  0.4841693|  1.962009|             0|
| DNA Damage Response (only ATM dependent)%WikiPathways\_20180810%WP710%Homo sapiens                                              |  0.0011922|  0.4218610|  1.807293|             0|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                                          |  0.0011922|  0.5048453|  2.204919|             0|
| Bladder Cancer%WikiPathways\_20180810%WP2828%Homo sapiens                                                                       |  0.0011922|  0.6401758|  2.115825|             0|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                                            |  0.0011922|  0.4872567|  2.093074|             0|
| Notch Signaling Pathway%WikiPathways\_20180810%WP61%Homo sapiens                                                                |  0.0011922|  0.4922462|  1.919337|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                                        |  0.0011922|  0.4580451|  2.076677|             0|
| miRNA targets in ECM and membrane receptors%WikiPathways\_20180810%WP2911%Homo sapiens                                          |  0.0011922|  0.7770215|  2.382859|             0|
| Integrin-mediated Cell Adhesion%WikiPathways\_20180810%WP185%Homo sapiens                                                       |  0.0011922|  0.4920657|  2.064642|             0|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                                             |  0.0011922|  0.4208167|  1.894092|             0|
| MAPK Signaling Pathway%WikiPathways\_20180810%WP382%Homo sapiens                                                                |  0.0011922|  0.4004338|  1.919980|             0|
| Prolactin Signaling Pathway%WikiPathways\_20180810%WP2037%Homo sapiens                                                          |  0.0011922|  0.4912431|  1.998855|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                                    |  0.0011922|  0.6785853|  2.456140|             0|
| Endometrial cancer%WikiPathways\_20180810%WP4155%Homo sapiens                                                                   |  0.0011922|  0.5009958|  1.976177|             0|
| Association Between Physico-Chemical Features and Toxicity Associated Pathways%WikiPathways\_20180810%WP3680%Homo sapiens       |  0.0011922|  0.5023847|  1.976034|             0|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                                                             |  0.0011922|  0.4987816|  1.933995|             0|
| IL-6 signaling pathway%WikiPathways\_20180810%WP364%Homo sapiens                                                                |  0.0011922|  0.5407779|  1.920349|             0|
| RANKL/RANK (Receptor activator of NFKB (ligand)) Signaling Pathway%WikiPathways\_20180810%WP2018%Homo sapiens                   |  0.0011922|  0.5097552|  1.953866|             0|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                                 |  0.0011922|  0.4641251|  1.990121|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                                   |  0.0011922|  0.5874748|  2.198565|             0|
| Ebola Virus Pathway on Host%WikiPathways\_20180810%WP4217%Homo sapiens                                                          |  0.0011922|  0.4229853|  1.847395|             0|
| Androgen receptor signaling pathway%WikiPathways\_20180810%WP138%Homo sapiens                                                   |  0.0011922|  0.4416392|  1.839847|             0|
| Splicing factor NOVA regulated synaptic proteins%WikiPathways\_20180810%WP4148%Homo sapiens                                     |  0.0011922|  0.6007186|  2.133203|             0|
| Interleukin-11 Signaling Pathway%WikiPathways\_20180810%WP2332%Homo sapiens                                                     |  0.0011922|  0.5155563|  1.866056|             0|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                                       |  0.0011922|  0.3913055|  1.917798|             0|
| IL-3 Signaling Pathway%WikiPathways\_20180810%WP286%Homo sapiens                                                                |  0.0011922|  0.5810400|  2.063323|             0|
| VEGFA-VEGFR2 Signaling Pathway%WikiPathways\_20180810%WP3888%Homo sapiens                                                       |  0.0011922|  0.4380901|  2.109806|             0|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens                  |  0.0011922|  0.5432686|  2.247370|             0|
| Pancreatic adenocarcinoma pathway%WikiPathways\_20180810%WP4263%Homo sapiens                                                    |  0.0011922|  0.4608888|  1.933828|             0|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                                                           |  0.0011922|  0.3672727|  1.807004|             0|
| Breast cancer pathway%WikiPathways\_20180810%WP4262%Homo sapiens                                                                |  0.0011922|  0.4285451|  1.919135|             0|
| MET in type 1 papillary renal cell carcinoma%WikiPathways\_20180810%WP4205%Homo sapiens                                         |  0.0011922|  0.5322077|  2.054222|             0|
| Ras Signaling%WikiPathways\_20180810%WP4223%Homo sapiens                                                                        |  0.0011922|  0.3782239|  1.735253|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens                            |  0.0011922|  0.4366771|  1.974864|             0|
| Wnt Signaling Pathway and Pluripotency%WikiPathways\_20180810%WP399%Homo sapiens                                                |  0.0018603|  0.4098185|  1.726532|             1|
| Corticotropin-releasing hormone signaling pathway%WikiPathways\_20180810%WP2355%Homo sapiens                                    |  0.0018603|  0.4538540|  1.871675|             1|
| Insulin Signaling%WikiPathways\_20180810%WP481%Homo sapiens                                                                     |  0.0018603|  0.3701844|  1.698345|             1|
| PDGF Pathway%WikiPathways\_20180810%WP2526%Homo sapiens                                                                         |  0.0018603|  0.5480423|  1.934937|             1|
| MicroRNAs in cardiomyocyte hypertrophy%WikiPathways\_20180810%WP1544%Homo sapiens                                               |  0.0018603|  0.4283846|  1.779621|             1|
| Signaling Pathways in Glioblastoma%WikiPathways\_20180810%WP2261%Homo sapiens                                                   |  0.0018603|  0.4426197|  1.838758|             1|
| Apoptosis-related network due to altered Notch3 in ovarian cancer%WikiPathways\_20180810%WP2864%Homo sapiens                    |  0.0018603|  0.5095398|  1.906901|             1|
| IL-4 Signaling Pathway%WikiPathways\_20180810%WP395%Homo sapiens                                                                |  0.0018603|  0.4950647|  1.881905|             1|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP363%Homo sapiens                                                                 |  0.0018603|  0.5084747|  1.902915|             1|
| Chromosomal and microsatellite instability in colorectal cancer %WikiPathways\_20180810%WP4216%Homo sapiens                     |  0.0018603|  0.4620539|  1.869541|             1|
| Regulation of Actin Cytoskeleton%WikiPathways\_20180810%WP51%Homo sapiens                                                       |  0.0025533|  0.3692273|  1.650696|             2|
| Human Thyroid Stimulating Hormone (TSH) signaling pathway%WikiPathways\_20180810%WP2032%Homo sapiens                            |  0.0035290|  0.4648067|  1.828228|             3|
| Alpha 6 Beta 4 signaling pathway%WikiPathways\_20180810%WP244%Homo sapiens                                                      |  0.0035441|  0.5713602|  1.936809|             3|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                                              |  0.0043723|  0.5585430|  1.935634|             4|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                                    |  0.0049000|  0.3942640|  1.672344|             5|
| IL-5 Signaling Pathway%WikiPathways\_20180810%WP127%Homo sapiens                                                                |  0.0049936|  0.5352106|  1.889633|             5|
| IL17 signaling pathway%WikiPathways\_20180810%WP2112%Homo sapiens                                                               |  0.0049936|  0.5926364|  1.901373|             5|
| Estrogen signaling pathway%WikiPathways\_20180810%WP712%Homo sapiens                                                            |  0.0049936|  0.6564416|  1.986938|             5|
| T-Cell antigen Receptor (TCR) Signaling Pathway%WikiPathways\_20180810%WP69%Homo sapiens                                        |  0.0062571|  0.4262051|  1.727120|             7|
| Non-small cell lung cancer%WikiPathways\_20180810%WP4255%Homo sapiens                                                           |  0.0062571|  0.4373975|  1.744379|             7|
| Chemokine signaling pathway%WikiPathways\_20180810%WP3929%Homo sapiens                                                          |  0.0073159|  0.3648770|  1.608872|             9|
| Myometrial Relaxation and Contraction Pathways%WikiPathways\_20180810%WP289%Homo sapiens                                        |  0.0073159|  0.3563647|  1.595893|             9|
| Hedgehog Signaling Pathway%WikiPathways\_20180810%WP4249%Homo sapiens                                                           |  0.0076893|  0.5100525|  1.800809|             9|
| Oncostatin M Signaling Pathway%WikiPathways\_20180810%WP2374%Homo sapiens                                                       |  0.0089153|  0.4384284|  1.724474|            11|
| Kit receptor signaling pathway%WikiPathways\_20180810%WP304%Homo sapiens                                                        |  0.0102779|  0.4389567|  1.711553|            13|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                                      |  0.0103792|  0.5589618|  1.847407|            13|
| TNF related weak inducer of apoptosis (TWEAK) Signaling Pathway%WikiPathways\_20180810%WP2036%Homo sapiens                      |  0.0122812|  0.4906609|  1.764910|            16|
| Nanoparticle-mediated activation of receptor signaling%WikiPathways\_20180810%WP2643%Homo sapiens                               |  0.0124066|  0.5658621|  1.852351|            16|
| Spinal Cord Injury%WikiPathways\_20180810%WP2431%Homo sapiens                                                                   |  0.0124790|  0.3801514|  1.604727|            18|
| Notch Signaling Pathway%WikiPathways\_20180810%WP268%Homo sapiens                                                               |  0.0124790|  0.4822076|  1.745351|            17|
| PDGFR-beta pathway%WikiPathways\_20180810%WP3972%Homo sapiens                                                                   |  0.0124790|  0.5466798|  1.821663|            17|
| Lipid Metabolism Pathway%WikiPathways\_20180810%WP3965%Homo sapiens                                                             |  0.0124790|  0.5683300|  1.823390|            17|
| Photodynamic therapy-induced AP-1 survival signaling.%WikiPathways\_20180810%WP3611%Homo sapiens                                |  0.0133797|  0.4655896|  1.719760|            19|
| Leptin Insulin Overlap%WikiPathways\_20180810%WP3935%Homo sapiens                                                               |  0.0134718|  0.6531517|  1.831104|            19|
| Viral Acute Myocarditis%WikiPathways\_20180810%WP4298%Homo sapiens                                                              |  0.0171999|  0.4020623|  1.616573|            26|
| Signaling of Hepatocyte Growth Factor Receptor%WikiPathways\_20180810%WP313%Homo sapiens                                        |  0.0177930|  0.5027303|  1.742214|            27|
| Aryl Hydrocarbon Receptor Pathway%WikiPathways\_20180810%WP2873%Homo sapiens                                                    |  0.0177930|  0.5236745|  1.745004|            27|
| Brain-Derived Neurotrophic Factor (BDNF) signaling pathway%WikiPathways\_20180810%WP2380%Homo sapiens                           |  0.0177930|  0.3383520|  1.522919|            29|
| IL-1 signaling pathway%WikiPathways\_20180810%WP195%Homo sapiens                                                                |  0.0184121|  0.4333018|  1.660824|            29|
| Aryl Hydrocarbon Receptor%WikiPathways\_20180810%WP2586%Homo sapiens                                                            |  0.0202512|  0.4709270|  1.693927|            32|
| Extracellular vesicle-mediated signaling in recipient cells%WikiPathways\_20180810%WP2870%Homo sapiens                          |  0.0203010|  0.5183525|  1.727270|            32|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                                             |  0.0232503|  0.3392901|  1.513186|            40|
| Thymic Stromal LymphoPoietin (TSLP) Signaling Pathway%WikiPathways\_20180810%WP2203%Homo sapiens                                |  0.0238256|  0.4544407|  1.651366|            39|
| BMP Signaling Pathway in Eyelid Development%WikiPathways\_20180810%WP3927%Homo sapiens                                          |  0.0288288|  0.5936611|  1.744556|            47|
| miRs in Muscle Cell Differentiation%WikiPathways\_20180810%WP2012%Homo sapiens                                                  |  0.0329154|  0.5339054|  1.696691|            55|
| ErbB Signaling Pathway%WikiPathways\_20180810%WP673%Homo sapiens                                                                |  0.0337820|  0.4194829|  1.594594|            59|
| Simplified Interaction Map Between LOXL4 and Oxidative Stress Pathway%WikiPathways\_20180810%WP3670%Homo sapiens                |  0.0363392|  0.5737966|  1.686181|            62|
| AGE/RAGE pathway%WikiPathways\_20180810%WP2324%Homo sapiens                                                                     |  0.0367943|  0.3985061|  1.560752|            65|
| The effect of progerin on the involved genes in Hutchinson-Gilford Progeria Syndrome%WikiPathways\_20180810%WP4320%Homo sapiens |  0.0397839|  0.4691903|  1.625981|            70|
| White fat cell differentiation%WikiPathways\_20180810%WP4149%Homo sapiens                                                       |  0.0398584|  0.5005813|  1.638654|            70|
| TNF alpha Signaling Pathway%WikiPathways\_20180810%WP231%Homo sapiens                                                           |  0.0408420|  0.3570467|  1.498120|            77|
| Serotonin Receptor 4/6/7 and NR3C Signaling%WikiPathways\_20180810%WP734%Homo sapiens                                           |  0.0417766|  0.5812404|  1.681808|            76|
| MAPK Cascade%WikiPathways\_20180810%WP422%Homo sapiens                                                                          |  0.0417766|  0.4896546|  1.631642|            76|
| Calcium Regulation in the Cardiac Cell%WikiPathways\_20180810%WP536%Homo sapiens                                                |  0.0425155|  0.3269138|  1.444239|            84|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                                       |  0.0437586|  0.2794731|  1.352347|            93|
| Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds%WikiPathways\_20180810%WP3664%Homo sapiens                    |  0.0437586|  0.5969873|  1.673648|            82|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                                           |  0.0452597|  0.3253178|  1.437188|            92|
| Preimplantation Embryo%WikiPathways\_20180810%WP3527%Homo sapiens                                                               |  0.0489906|  0.4265397|  1.558349|            96|

``` r
upPathwaysWiki75748_024 <- upPathwaysWiki
save(upPathwaysWiki75748_024, file="upPathwaysWiki75748_024.Rdata")
```

Analysis of second-fourth days of differentiation 24h VS 96h (formation of definitive endoderm from APS)
--------------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE75748/DEgenes_24h_96h_75748.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
#absolute Log fold changes for ermineR
ermineInputGeneScores <- DEgenes_24h_96h_75748 %>% 
  rownames_to_column("gene") %>%
  mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

#exact log fold changes for fgsea
# scores forfgsea
scoresFGSEADF <- DEgenes_24h_96h_75748 %>% 
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
    ## toType = "ENTREZID", : 6.27% of input gene IDs are fail to map...

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

| Pathway                                              | Corrected p-value | Corrected MF p-value |
|:-----------------------------------------------------|:------------------|:---------------------|
| Extracellular matrix organization                    | 0.0000000         | 0.0000000            |
| Class B 2 (Secretin family receptors)                | 0.0000000         | 0.0000000            |
| Neuronal System                                      | 0.0247600         | 0.0309500            |
| Class A 1 (Rhodopsin-like receptors)                 | 0.0309500         | 0.0412667            |
| G alpha (i) signalling events                        | 0.0000000         | 0.0495200            |
| Collagen formation                                   | 0.1326429         | 0.1591714            |
| Platelet activation, signaling and aggregation       | 0.1134833         | 0.1650667            |
| Collagen chain trimerization                         | 0.2957444         | 0.2957444            |
| ECM proteoglycans                                    | 0.3017625         | 0.3017625            |
| Binding and Uptake of Ligands by Scavenger Receptors | 0.3095000         | 0.3218800            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 5 enriched Reactome pathways with ermineR

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

There are 22 up-regulated and 60 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                              |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                 |  0.0056376|  0.6411296|  2.013795|             0|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                |  0.0056376|  0.4229049|  1.736420|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                  |  0.0056376|  0.6223854|  2.081821|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                            |  0.0056376|  0.5076242|  2.111313|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                          |  0.0056376|  0.5986648|  2.099407|             0|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                        |  0.0056376|  0.4983576|  1.842697|             0|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                           |  0.0056376|  0.4408896|  1.795852|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090 |  0.0056376|  0.6139364|  1.993003|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                          |  0.0056376|  0.5305453|  1.949592|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                            |  0.0073997|  0.6009169|  1.936648|             1|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                      |  0.0073997|  0.6215197|  1.952201|             1|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                  |  0.0088720|  0.4773152|  1.780324|             2|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                   |  0.0088884|  0.5061851|  1.832215|             2|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                           |  0.0144843|  0.5639546|  1.823764|             4|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                              |  0.0151555|  0.7521806|  1.895901|             4|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                         |  0.0153527|  0.4328159|  1.675405|             5|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                             |  0.0198040|  0.4911359|  1.736448|             7|
| COMPLEMENT CASCADE%REACTOME DATABASE ID RELEASE 65%166658                                            |  0.0211466|  0.6803816|  1.870540|             7|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                    |  0.0280678|  0.6871840|  1.844578|            10|
| DISEASES OF GLYCOSYLATION%REACTOME DATABASE ID RELEASE 65%3781865                                    |  0.0310532|  0.4987830|  1.725498|            13|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                        |  0.0337389|  0.6020824|  1.779428|            14|
| G ALPHA (Q) SIGNALLING EVENTS%REACTOME DATABASE ID RELEASE 65%416476                                 |  0.0376486|  0.4091554|  1.587525|            18|

``` r
upPathwaysReactome75748_2496 <- upPathwaysReactome
save(upPathwaysReactome75748_2496, file="upPathwaysReactome75748_2496.Rdata")
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
| Kegg neuroactive ligand receptor interaction      | 0.0000000         | 0.0000000            |
| Kegg ecm receptor interaction                     | 0.0000000         | 0.0000000            |
| Kegg cell adhesion molecules cams                 | 0.0054000         | 0.0054000            |
| Kegg cytokine cytokine receptor interaction       | 0.0135000         | 0.0283500            |
| Kegg calcium signaling pathway                    | 0.0097200         | 0.0324000            |
| Kegg focal adhesion                               | 0.0081000         | 0.0486000            |
| Kegg taste transduction                           | 0.1494000         | 0.1944000            |
| Kegg hematopoietic cell lineage                   | 0.1620000         | 0.2369250            |
| Kegg systemic lupus erythematosus                 | 0.1907182         | 0.2503636            |
| Kegg intestinal immune network for iga production | 0.1984500         | 0.2624400            |

``` r
#based on corrected MFPvalues 
sizeErmineKEGG <- enrichmentResultKEGG$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 6 enriched KEGG pathways with ermineR.

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

There are 4 up-regulated and 3 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                      |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_FOCAL\_ADHESION                        |  0.0134156|  0.4368927|  1.772686|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION             |  0.0134156|  0.5627785|  1.971726|             0|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES |  0.0303934|  0.6288321|  1.934071|             3|
| KEGG\_HEMATOPOIETIC\_CELL\_LINEAGE           |  0.0390409|  0.5695954|  1.821477|             8|

``` r
upPathwaysKEGG75748_2496 <- upPathwaysKEGG
save(upPathwaysKEGG75748_2496, file="upPathwaysKEGG75748_2496.Rdata")
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

Pathway Corrected p-value Corrected MF p-value -------- ------------------ ---------------------

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 0 enriched Wiki pathways with ermineR

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

There are 10 up-regulated and 2 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Hematopoietic Stem Cell Differentiation%WikiPathways\_20180810%WP2849%Homo sapiens        |  0.0182893|  0.6673470|  2.087039|             0|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens              |  0.0182893|  0.4983099|  1.826037|             0|
| Endothelin Pathways%WikiPathways\_20180810%WP2197%Homo sapiens                            |  0.0182893|  0.6915424|  1.962654|             0|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                       |  0.0197789|  0.4502829|  1.734859|             1|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                              |  0.0197789|  0.5958494|  1.894809|             2|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                           |  0.0197789|  0.4600371|  1.703385|             2|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                     |  0.0197789|  0.3826395|  1.617099|             2|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens |  0.0203296|  0.3673466|  1.543849|             3|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                       |  0.0208971|  0.5285724|  1.782414|             3|
| Complement and Coagulation Cascades%WikiPathways\_20180810%WP558%Homo sapiens             |  0.0294011|  0.6264701|  1.892451|             5|

``` r
upPathwaysWiki75748_2496 <- upPathwaysWiki
save(upPathwaysWiki75748_2496, file="upPathwaysWiki75748_2496.Rdata")
```

Analysis of whole differentiation process 0h VS 96 h (formation of Definitive endoderm from hESC)
-------------------------------------------------------------------------------------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE75748/DEgenes_0h_96h_75748.Rdata")
```

Sorted log Fold Changes give us a ranked list:

``` r
#absolute Log fold changes for ermineR
ermineInputGeneScores <- DEgenes_0h_96h_75748 %>% 
  rownames_to_column("gene") %>%
  mutate(absolute_logFC = abs(logFC)) %>% 
  dplyr::select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC)) %>% 
  column_to_rownames("gene")

#exact log fold changes for fgsea
# scores forfgsea
scoresFGSEADF <- DEgenes_0h_96h_75748 %>% 
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
    ## toType = "ENTREZID", : 6.27% of input gene IDs are fail to map...

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

| Pathway                                              | Corrected p-value | Corrected MF p-value |
|:-----------------------------------------------------|:------------------|:---------------------|
| Extracellular matrix organization                    | 0.0000000         | 0.0000000            |
| ECM proteoglycans                                    | 0.0000000         | 0.0000000            |
| Collagen formation                                   | 0.0000000         | 0.0000000            |
| Collagen chain trimerization                         | 0.0000000         | 0.0000000            |
| Collagen biosynthesis and modifying enzymes          | 0.0000000         | 0.0000000            |
| G alpha (i) signalling events                        | 0.0309500         | 0.1444333            |
| Signaling by NODAL                                   | 0.2785500         | 0.2862875            |
| Signaling by BMP                                     | 0.2888667         | 0.3026222            |
| Binding and Uptake of Ligands by Scavenger Receptors | 0.3183429         | 0.3183429            |
| SLC-mediated transmembrane transport                 | 0.3095000         | 0.3837800            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 5 enriched Reactome pathways with ermineR

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

There are 126 up-regulated and 116 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| LAMININ INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000157                                                                              |  0.0049800|  0.7271213|  2.200165|             0|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                                                      |  0.0049800|  0.7634609|  2.658488|             0|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                                                     |  0.0049800|  0.3726225|  1.730304|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                                                       |  0.0049800|  0.7369268|  2.739950|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                                                 |  0.0049800|  0.5727275|  2.700714|             0|
| GLYCOSAMINOGLYCAN METABOLISM%REACTOME DATABASE ID RELEASE 65%1630316                                                                      |  0.0049800|  0.4614803|  1.959137|             0|
| EPH-EPHRIN MEDIATED REPULSION OF CELLS%REACTOME%R-HSA-3928665.3                                                                           |  0.0049800|  0.5303114|  1.934427|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                                                               |  0.0049800|  0.6796685|  2.655697|             0|
| EPH-EPHRIN SIGNALING%REACTOME%R-HSA-2682334.1                                                                                             |  0.0049800|  0.4717783|  1.936137|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                                                 |  0.0049800|  0.6656855|  2.375160|             0|
| SIGNALING BY PDGF%REACTOME DATABASE ID RELEASE 65%186797                                                                                  |  0.0049800|  0.5737726|  2.097503|             0|
| HS-GAG BIOSYNTHESIS%REACTOME%R-HSA-2022928.1                                                                                              |  0.0049800|  0.6647125|  2.165515|             0|
| SIGNALING BY VEGF%REACTOME DATABASE ID RELEASE 65%194138                                                                                  |  0.0049800|  0.4638137|  1.922169|             0|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                                                        |  0.0049800|  0.5364088|  2.176827|             0|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                                                                |  0.0049800|  0.5966669|  2.138085|             0|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                                           |  0.0049800|  0.7033912|  2.449316|             0|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                                                |  0.0049800|  0.6060797|  2.239085|             0|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                                                         |  0.0049800|  0.6600929|  2.118793|             0|
| RHO GTPASE CYCLE%REACTOME DATABASE ID RELEASE 65%194840                                                                                   |  0.0049800|  0.5021043|  2.191589|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090                                      |  0.0049800|  0.7323880|  2.642019|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                                                               |  0.0049800|  0.5596483|  2.304079|             0|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                                                              |  0.0072257|  0.5364065|  1.945252|             1|
| REGULATION OF FZD BY UBIQUITINATION%REACTOME%R-HSA-4641263.2                                                                              |  0.0072257|  0.6978179|  2.055875|             1|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                                                         |  0.0072257|  0.6910985|  2.036078|             1|
| SIGNALING BY NOTCH3%REACTOME DATABASE ID RELEASE 65%9012852                                                                               |  0.0072257|  0.5394003|  1.945834|             1|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                                                                |  0.0072257|  0.6000930|  1.975750|             1|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                                                             |  0.0072257|  0.6276843|  1.996373|             1|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                                                                   |  0.0072257|  0.7379233|  2.052832|             1|
| HEPARAN SULFATE HEPARIN (HS-GAG) METABOLISM%REACTOME%R-HSA-1638091.1                                                                      |  0.0093426|  0.5241214|  1.923898|             2|
| SIGNALING BY NTRK1 (TRKA)%REACTOME%R-HSA-187037.2                                                                                         |  0.0093426|  0.4692721|  1.873978|             2|
| NRAGE SIGNALS DEATH THROUGH JNK%REACTOME DATABASE ID RELEASE 65%193648                                                                    |  0.0102636|  0.5086439|  1.879121|             3|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                                                             |  0.0102636|  0.5940740|  1.926734|             3|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                                                             |  0.0102636|  0.4303519|  1.782064|             3|
| MET PROMOTES CELL MOTILITY%REACTOME DATABASE ID RELEASE 65%8875878                                                                        |  0.0102636|  0.6098729|  1.957595|             3|
| SIGNALING BY NTRKS%REACTOME DATABASE ID RELEASE 65%166520                                                                                 |  0.0106615|  0.4218281|  1.757071|             4|
| DEFECTIVE B3GALTL CAUSES PETERS-PLUS SYNDROME (PPS)%REACTOME%R-HSA-5083635.1                                                              |  0.0106615|  0.5731837|  1.887154|             4|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                                                              |  0.0106615|  0.4235370|  1.776752|             4|
| SIGNALING BY WNT%REACTOME DATABASE ID RELEASE 65%195721                                                                                   |  0.0106615|  0.3056332|  1.478164|             4|
| NON-INTEGRIN MEMBRANE-ECM INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000171                                                            |  0.0106615|  0.5349328|  1.862718|             4|
| CELL-EXTRACELLULAR MATRIX INTERACTIONS%REACTOME%R-HSA-446353.1                                                                            |  0.0106615|  0.7135915|  1.985143|             4|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                                                           |  0.0119022|  0.4428251|  1.789220|             5|
| O-GLYCOSYLATION OF TSR DOMAIN-CONTAINING PROTEINS%REACTOME%R-HSA-5173214.1                                                                |  0.0119387|  0.5821632|  1.928416|             5|
| INTRA-GOLGI AND RETROGRADE GOLGI-TO-ER TRAFFIC%REACTOME DATABASE ID RELEASE 65%6811442                                                    |  0.0126160|  0.3474726|  1.584038|             6|
| VEGFA-VEGFR2 PATHWAY%REACTOME%R-HSA-4420097.3                                                                                             |  0.0127307|  0.4313524|  1.761343|             6|
| SIGNALING BY EGFR%REACTOME DATABASE ID RELEASE 65%177929                                                                                  |  0.0127307|  0.5160705|  1.829773|             6|
| RUNX2 REGULATES BONE DEVELOPMENT%REACTOME%R-HSA-8941326.1                                                                                 |  0.0127307|  0.5949921|  1.909830|             6|
| DISEASES OF GLYCOSYLATION%REACTOME DATABASE ID RELEASE 65%3781865                                                                         |  0.0127307|  0.4907506|  1.887855|             6|
| CHONDROITIN SULFATE DERMATAN SULFATE METABOLISM%REACTOME DATABASE ID RELEASE 65%1793185                                                   |  0.0155651|  0.5070917|  1.817103|             8|
| L1CAM INTERACTIONS%REACTOME%R-HSA-373760.2                                                                                                |  0.0155651|  0.4252722|  1.739798|             8|
| TP53 REGULATES TRANSCRIPTION OF CELL CYCLE GENES%REACTOME DATABASE ID RELEASE 65%6791312                                                  |  0.0168387|  0.4974348|  1.803923|             9|
| N-GLYCAN TRIMMING IN THE ER AND CALNEXIN CALRETICULIN CYCLE%REACTOME%R-HSA-532668.2                                                       |  0.0168856|  0.5501566|  1.876357|             9|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                                                       |  0.0182808|  0.3986534|  1.666131|            11|
| G ALPHA (12 13) SIGNALLING EVENTS%REACTOME%R-HSA-416482.3                                                                                 |  0.0182808|  0.4359478|  1.743793|            11|
| DISEASES ASSOCIATED WITH O-GLYCOSYLATION OF PROTEINS%REACTOME%R-HSA-3906995.2                                                             |  0.0198007|  0.5061652|  1.778390|            12|
| CLATHRIN-MEDIATED ENDOCYTOSIS%REACTOME%R-HSA-8856828.3                                                                                    |  0.0200563|  0.3596602|  1.569848|            13|
| TOLL LIKE RECEPTOR TLR6:TLR2 CASCADE%REACTOME%R-HSA-168188.1                                                                              |  0.0200840|  0.4101515|  1.683226|            13|
| MYD88:MAL CASCADE INITIATED ON PLASMA MEMBRANE%REACTOME%R-HSA-166058.2                                                                    |  0.0200840|  0.4101515|  1.683226|            13|
| TOLL LIKE RECEPTOR 2 (TLR2) CASCADE%REACTOME DATABASE ID RELEASE 65%181438                                                                |  0.0200840|  0.4101515|  1.683226|            13|
| TOLL LIKE RECEPTOR TLR1:TLR2 CASCADE%REACTOME DATABASE ID RELEASE 65%168179                                                               |  0.0200840|  0.4101515|  1.683226|            13|
| CLASS B 2 (SECRETIN FAMILY RECEPTORS)%REACTOME DATABASE ID RELEASE 65%373080                                                              |  0.0201530|  0.4816090|  1.771817|            13|
| TRAF6 MEDIATED INDUCTION OF NFKB AND MAP KINASES UPON TLR7 8 OR 9 ACTIVATION%REACTOME%R-HSA-975138.1                                      |  0.0209302|  0.4133133|  1.677287|            14|
| ASPARAGINE N-LINKED GLYCOSYLATION%REACTOME%R-HSA-446203.4                                                                                 |  0.0209626|  0.3045354|  1.463395|            15|
| SIGNALLING TO RAS%REACTOME DATABASE ID RELEASE 65%167044                                                                                  |  0.0228453|  0.6459158|  1.841258|            15|
| MET ACTIVATES PTK2 SIGNALING%REACTOME DATABASE ID RELEASE 65%8874081                                                                      |  0.0228453|  0.6485362|  1.848728|            15|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                                                                  |  0.0243778|  0.5650694|  1.813783|            17|
| TOLL LIKE RECEPTOR 7 8 (TLR7 8) CASCADE%REACTOME DATABASE ID RELEASE 65%168181                                                            |  0.0270582|  0.4047850|  1.646778|            20|
| MYD88 DEPENDENT CASCADE INITIATED ON ENDOSOME%REACTOME%R-HSA-975155.1                                                                     |  0.0270582|  0.4047850|  1.646778|            20|
| ACTIVATION OF MATRIX METALLOPROTEINASES%REACTOME%R-HSA-1592389.1                                                                          |  0.0270587|  0.6350531|  1.810292|            19|
| TOLL LIKE RECEPTOR 9 (TLR9) CASCADE%REACTOME DATABASE ID RELEASE 65%168138                                                                |  0.0277424|  0.4011339|  1.639233|            21|
| INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING%REACTOME DATABASE ID RELEASE 65%6785807                                                        |  0.0277424|  0.3982062|  1.627269|            21|
| SIGNALING BY WNT IN CANCER%REACTOME DATABASE ID RELEASE 65%4791275                                                                        |  0.0289597|  0.5181554|  1.740790|            22|
| O-LINKED GLYCOSYLATION%REACTOME%R-HSA-5173105.3                                                                                           |  0.0300987|  0.4037059|  1.633970|            24|
| CALNEXIN CALRETICULIN CYCLE%REACTOME%R-HSA-901042.2                                                                                       |  0.0300987|  0.5631319|  1.791062|            23|
| NOTCH3 ACTIVATION AND TRANSMISSION OF SIGNAL TO THE NUCLEUS%REACTOME DATABASE ID RELEASE 65%9013507                                       |  0.0300987|  0.5749506|  1.795427|            23|
| MYD88 CASCADE INITIATED ON PLASMA MEMBRANE%REACTOME%R-HSA-975871.1                                                                        |  0.0310441|  0.4014290|  1.624754|            26|
| SIGNALLING TO ERKS%REACTOME%R-HSA-187687.1                                                                                                |  0.0310441|  0.5403635|  1.752537|            25|
| TOLL LIKE RECEPTOR 10 (TLR10) CASCADE%REACTOME DATABASE ID RELEASE 65%168142                                                              |  0.0310441|  0.4014290|  1.624754|            26|
| TOLL LIKE RECEPTOR 5 (TLR5) CASCADE%REACTOME%R-HSA-168176.1                                                                               |  0.0310441|  0.4014290|  1.624754|            26|
| TRIF(TICAM1)-MEDIATED TLR4 SIGNALING%REACTOME%R-HSA-937061.2                                                                              |  0.0310441|  0.3917665|  1.612907|            26|
| MYD88-INDEPENDENT TLR4 CASCADE%REACTOME%R-HSA-166166.2                                                                                    |  0.0310441|  0.3917665|  1.612907|            26|
| CASPASE ACTIVATION VIA EXTRINSIC APOPTOTIC SIGNALLING PATHWAY%REACTOME%R-HSA-5357769.2                                                    |  0.0313155|  0.6251075|  1.781941|            25|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                                                             |  0.0313704|  0.4638759|  1.695761|            27|
| TOLL LIKE RECEPTOR 3 (TLR3) CASCADE%REACTOME DATABASE ID RELEASE 65%168164                                                                |  0.0313704|  0.3905721|  1.602874|            27|
| MAP KINASE ACTIVATION%REACTOME%R-HSA-450294.3                                                                                             |  0.0313704|  0.4285091|  1.648419|            27|
| DEATH RECEPTOR SIGNALLING%REACTOME%R-HSA-73887.3                                                                                          |  0.0317335|  0.3487413|  1.520458|            29|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                                                                    |  0.0317335|  0.3752864|  1.570398|            28|
| SMOOTH MUSCLE CONTRACTION%REACTOME%R-HSA-445355.2                                                                                         |  0.0323274|  0.5225519|  1.730953|            29|
| MAPK TARGETS NUCLEAR EVENTS MEDIATED BY MAP KINASES%REACTOME DATABASE ID RELEASE 65%450282                                                |  0.0323274|  0.5219748|  1.729042|            29|
| PEPTIDE HORMONE METABOLISM%REACTOME DATABASE ID RELEASE 65%2980736                                                                        |  0.0360558|  0.4425518|  1.661998|            33|
| G-PROTEIN BETA:GAMMA SIGNALLING%REACTOME DATABASE ID RELEASE 65%397795                                                                    |  0.0362171|  0.5633308|  1.740614|            33|
| TOLL LIKE RECEPTOR 4 (TLR4) CASCADE%REACTOME%R-HSA-166016.2                                                                               |  0.0372376|  0.3554246|  1.522014|            36|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                                                |  0.0382943|  0.3126855|  1.439759|            38|
| TGF-BETA RECEPTOR SIGNALING ACTIVATES SMADS%REACTOME DATABASE ID RELEASE 65%2173789                                                       |  0.0391774|  0.5140958|  1.702943|            37|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                                                            |  0.0395288|  0.6271626|  1.771262|            37|
| INTERLEUKIN-17 SIGNALING%REACTOME DATABASE ID RELEASE 65%448424                                                                           |  0.0395288|  0.4152766|  1.612291|            38|
| NOTCH3 INTRACELLULAR DOMAIN REGULATES TRANSCRIPTION%REACTOME%R-HSA-9013508.1                                                              |  0.0397511|  0.5615138|  1.716509|            38|
| NOTCH4 INTRACELLULAR DOMAIN REGULATES TRANSCRIPTION%REACTOME DATABASE ID RELEASE 65%9013695                                               |  0.0397511|  0.5904751|  1.739627|            38|
| FCERI MEDIATED MAPK ACTIVATION%REACTOME%R-HSA-2871796.2                                                                                   |  0.0410311|  0.5523948|  1.724991|            40|
| NEURODEGENERATIVE DISEASES%REACTOME DATABASE ID RELEASE 65%8863678                                                                        |  0.0418206|  0.5811884|  1.734988|            42|
| DEREGULATED CDK5 TRIGGERS MULTIPLE NEURODEGENERATIVE PATHWAYS IN ALZHEIMER'S DISEASE MODELS%REACTOME%R-HSA-8862803.2                      |  0.0418206|  0.5811884|  1.734988|            42|
| DOWNSTREAM SIGNAL TRANSDUCTION%REACTOME DATABASE ID RELEASE 65%186763                                                                     |  0.0420112|  0.5307263|  1.703547|            43|
| TP53 REGULATES TRANSCRIPTION OF ADDITIONAL CELL CYCLE GENES WHOSE EXACT ROLE IN THE P53 PATHWAY REMAIN UNCERTAIN%REACTOME%R-HSA-6804115.1 |  0.0420963|  0.5888055|  1.734708|            43|
| SIGNALING BY RAS MUTANTS%REACTOME%R-HSA-6802949.1                                                                                         |  0.0420963|  0.4587096|  1.654751|            44|
| COPI-MEDIATED ANTEROGRADE TRANSPORT%REACTOME%R-HSA-6807878.1                                                                              |  0.0427077|  0.3899728|  1.565100|            46|
| INTRINSIC PATHWAY FOR APOPTOSIS%REACTOME%R-HSA-109606.2                                                                                   |  0.0427077|  0.4682451|  1.653956|            45|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                                                           |  0.0434145|  0.4176251|  1.589320|            47|
| PRE-NOTCH PROCESSING IN GOLGI%REACTOME%R-HSA-1912420.2                                                                                    |  0.0441780|  0.6137732|  1.733447|            47|
| EPHRIN SIGNALING%REACTOME DATABASE ID RELEASE 65%3928664                                                                                  |  0.0455419|  0.5976574|  1.703692|            49|
| REGULATION OF PTEN GENE TRANSCRIPTION%REACTOME DATABASE ID RELEASE 65%8943724                                                             |  0.0455419|  0.4355678|  1.628492|            52|
| RUNX2 REGULATES OSTEOBLAST DIFFERENTIATION%REACTOME DATABASE ID RELEASE 65%8940973                                                        |  0.0455419|  0.5738155|  1.712978|            51|
| GAMMA CARBOXYLATION, HYPUSINE FORMATION AND ARYLSULFATASE ACTIVATION%REACTOME DATABASE ID RELEASE 65%163841                               |  0.0471895|  0.5140355|  1.667148|            55|
| CELL DEATH SIGNALLING VIA NRAGE, NRIF AND NADE%REACTOME DATABASE ID RELEASE 65%204998                                                     |  0.0474351|  0.4017010|  1.575063|            56|
| SYNTHESIS OF PIPS AT THE PLASMA MEMBRANE%REACTOME DATABASE ID RELEASE 65%1660499                                                          |  0.0474869|  0.4358078|  1.599725|            56|
| P75 NTR RECEPTOR-MEDIATED SIGNALLING%REACTOME%R-HSA-193704.1                                                                              |  0.0479851|  0.3745126|  1.532139|            57|
| ACTIVATION OF BH3-ONLY PROTEINS%REACTOME DATABASE ID RELEASE 65%114452                                                                    |  0.0480211|  0.5133241|  1.664841|            57|
| SIGNALING BY NOTCH1 PEST DOMAIN MUTANTS IN CANCER%REACTOME DATABASE ID RELEASE 65%2644602                                                 |  0.0493758|  0.4275683|  1.600238|            61|
| CONSTITUTIVE SIGNALING BY NOTCH1 HD+PEST DOMAIN MUTANTS%REACTOME DATABASE ID RELEASE 65%2894862                                           |  0.0493758|  0.4275683|  1.600238|            61|
| INTEGRIN SIGNALING%REACTOME DATABASE ID RELEASE 65%9006921                                                                                |  0.0493758|  0.5512139|  1.685023|            59|
| INTEGRIN ALPHAIIB BETA3 SIGNALING%REACTOME DATABASE ID RELEASE 65%354192                                                                  |  0.0493758|  0.5512139|  1.685023|            59|
| SIGNALING BY NOTCH1 IN CANCER%REACTOME DATABASE ID RELEASE 65%2644603                                                                     |  0.0493758|  0.4275683|  1.600238|            61|
| SIGNALING BY NOTCH1 HD+PEST DOMAIN MUTANTS IN CANCER%REACTOME DATABASE ID RELEASE 65%2894858                                              |  0.0493758|  0.4275683|  1.600238|            61|
| SIGNALING BY LIGAND-RESPONSIVE EGFR VARIANTS IN CANCER%REACTOME DATABASE ID RELEASE 65%5637815                                            |  0.0493758|  0.5804440|  1.684686|            58|
| G ALPHA (Q) SIGNALLING EVENTS%REACTOME DATABASE ID RELEASE 65%416476                                                                      |  0.0493758|  0.3331424|  1.456782|            63|
| CONSTITUTIVE SIGNALING BY LIGAND-RESPONSIVE EGFR CANCER VARIANTS%REACTOME DATABASE ID RELEASE 65%1236382                                  |  0.0493758|  0.5804440|  1.684686|            58|
| CONSTITUTIVE SIGNALING BY NOTCH1 PEST DOMAIN MUTANTS%REACTOME%R-HSA-2644606.1                                                             |  0.0493758|  0.4275683|  1.600238|            61|
| SIGNALING BY EGFR IN CANCER%REACTOME%R-HSA-1643713.1                                                                                      |  0.0493758|  0.5804440|  1.684686|            58|

``` r
upPathwaysReactome75748_096 <- upPathwaysReactome
save(upPathwaysReactome75748_096, file="upPathwaysReactome75748_096.Rdata")
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

| Pathway                                     | Corrected p-value | Corrected MF p-value |
|:--------------------------------------------|:------------------|:---------------------|
| Kegg pathways in cancer                     | 0.0226800         | 0.0000000            |
| Kegg focal adhesion                         | 0.0000000         | 0.0000000            |
| Kegg ecm receptor interaction               | 0.0000000         | 0.0000000            |
| Kegg wnt signaling pathway                  | 0.0108000         | 0.0040500            |
| Kegg mapk signaling pathway                 | 0.0370286         | 0.0162000            |
| Kegg cell adhesion molecules cams           | 0.0297000         | 0.0189000            |
| Kegg axon guidance                          | 0.0121500         | 0.0194400            |
| Kegg dilated cardiomyopathy                 | 0.2138400         | 0.0627750            |
| Kegg erbb signaling pathway                 | 0.2246400         | 0.0736364            |
| Kegg cytokine cytokine receptor interaction | 0.0789750         | 0.0738000            |

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

There are 48 up-regulated and 15 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                                        |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                                 |  0.0018830|  0.3649497|  1.730410|             0|
| KEGG\_ERBB\_SIGNALING\_PATHWAY                                 |  0.0018830|  0.4471755|  1.813840|             0|
| KEGG\_ENDOCYTOSIS                                              |  0.0018830|  0.3778583|  1.718132|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY                                  |  0.0018830|  0.4658496|  2.044840|             0|
| KEGG\_AXON\_GUIDANCE                                           |  0.0018830|  0.4539690|  1.970230|             0|
| KEGG\_FOCAL\_ADHESION                                          |  0.0018830|  0.5517971|  2.527500|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                               |  0.0018830|  0.6864207|  2.704448|             0|
| KEGG\_NEUROTROPHIN\_SIGNALING\_PATHWAY                         |  0.0018830|  0.5058686|  2.176880|             0|
| KEGG\_REGULATION\_OF\_ACTIN\_CYTOSKELETON                      |  0.0018830|  0.3536441|  1.620820|             0|
| KEGG\_MELANOGENESIS                                            |  0.0018830|  0.4873528|  1.972072|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                                     |  0.0018830|  0.4080179|  1.982873|             0|
| KEGG\_COLORECTAL\_CANCER                                       |  0.0018830|  0.5189777|  2.010172|             0|
| KEGG\_PANCREATIC\_CANCER                                       |  0.0018830|  0.5206546|  2.056040|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA                                   |  0.0018830|  0.6165568|  2.218047|             0|
| KEGG\_CHRONIC\_MYELOID\_LEUKEMIA                               |  0.0018830|  0.5024482|  1.997538|             0|
| KEGG\_SMALL\_CELL\_LUNG\_CANCER                                |  0.0018830|  0.5145797|  2.096118|             0|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY                             |  0.0029389|  0.5382415|  1.927247|             1|
| KEGG\_VASOPRESSIN\_REGULATED\_WATER\_REABSORPTION              |  0.0029389|  0.5789300|  2.025111|             1|
| KEGG\_RENAL\_CELL\_CARCINOMA                                   |  0.0029389|  0.4768893|  1.882926|             1|
| KEGG\_HYPERTROPHIC\_CARDIOMYOPATHY\_HCM                        |  0.0029389|  0.4956741|  1.932364|             1|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY                            |  0.0038065|  0.4493235|  1.802434|             2|
| KEGG\_INSULIN\_SIGNALING\_PATHWAY                              |  0.0038065|  0.3818684|  1.657312|             2|
| KEGG\_DILATED\_CARDIOMYOPATHY                                  |  0.0038065|  0.4723732|  1.865095|             2|
| KEGG\_CHEMOKINE\_SIGNALING\_PATHWAY                            |  0.0077313|  0.3639969|  1.589976|             6|
| KEGG\_GNRH\_SIGNALING\_PATHWAY                                 |  0.0077313|  0.4203573|  1.700975|             6|
| KEGG\_ACUTE\_MYELOID\_LEUKEMIA                                 |  0.0084352|  0.4792117|  1.804243|             7|
| KEGG\_MTOR\_SIGNALING\_PATHWAY                                 |  0.0113699|  0.4843856|  1.788512|            10|
| KEGG\_GLYCOSAMINOGLYCAN\_BIOSYNTHESIS\_HEPARAN\_SULFATE        |  0.0128604|  0.5890263|  1.847293|            12|
| KEGG\_LYSOSOME                                                 |  0.0128604|  0.3674954|  1.571546|            12|
| KEGG\_ARRHYTHMOGENIC\_RIGHT\_VENTRICULAR\_CARDIOMYOPATHY\_ARVC |  0.0136147|  0.4456180|  1.718914|            14|
| KEGG\_TOLL\_LIKE\_RECEPTOR\_SIGNALING\_PATHWAY                 |  0.0162782|  0.4177962|  1.649606|            18|
| KEGG\_PROGESTERONE\_MEDIATED\_OOCYTE\_MATURATION               |  0.0185509|  0.3950158|  1.591300|            21|
| KEGG\_VEGF\_SIGNALING\_PATHWAY                                 |  0.0203203|  0.4112727|  1.623849|            24|
| KEGG\_FC\_EPSILON\_RI\_SIGNALING\_PATHWAY                      |  0.0204460|  0.4195961|  1.625234|            25|
| KEGG\_PRION\_DISEASES                                          |  0.0204460|  0.5582828|  1.784036|            25|
| KEGG\_NON\_SMALL\_CELL\_LUNG\_CANCER                           |  0.0204460|  0.4341424|  1.639508|            26|
| KEGG\_APOPTOSIS                                                |  0.0231049|  0.3916074|  1.570909|            30|
| KEGG\_PHOSPHATIDYLINOSITOL\_SIGNALING\_SYSTEM                  |  0.0249590|  0.3982205|  1.574439|            35|
| KEGG\_VASCULAR\_SMOOTH\_MUSCLE\_CONTRACTION                    |  0.0258031|  0.3676556|  1.522348|            40|
| KEGG\_T\_CELL\_RECEPTOR\_SIGNALING\_PATHWAY                    |  0.0258031|  0.3747279|  1.543308|            40|
| KEGG\_ENDOMETRIAL\_CANCER                                      |  0.0258031|  0.4409380|  1.641719|            38|
| KEGG\_INOSITOL\_PHOSPHATE\_METABOLISM                          |  0.0262640|  0.4388364|  1.629462|            41|
| KEGG\_PROSTATE\_CANCER                                         |  0.0275423|  0.3724379|  1.523363|            46|
| KEGG\_ADIPOCYTOKINE\_SIGNALING\_PATHWAY                        |  0.0367229|  0.3998956|  1.539599|            64|
| KEGG\_HEMATOPOIETIC\_CELL\_LINEAGE                             |  0.0390802|  0.4481770|  1.604758|            68|
| KEGG\_GAP\_JUNCTION                                            |  0.0391628|  0.3715653|  1.494107|            70|
| KEGG\_BLADDER\_CANCER                                          |  0.0395526|  0.4578896|  1.601709|            72|
| KEGG\_B\_CELL\_RECEPTOR\_SIGNALING\_PATHWAY                    |  0.0422827|  0.3788605|  1.492684|            83|

``` r
upPathwaysKEGG75748_096 <- upPathwaysKEGG
save(upPathwaysKEGG75748_096, file="upPathwaysReactome75748_096.Rdata")
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
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            | 0.0000000         | 0.00000              |
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       | 0.0099000         | 0.01485              |
| miRNA targets in ECM and membrane receptors%WikiPathways\_20180810%WP2911%Homo sapiens                         | 0.0118800         | 0.01980              |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  | 0.0148500         | 0.02970              |
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens | 0.0148500         | 0.02970              |
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                | 0.0169714         | 0.03465              |

``` r
#based on corrected MFPvalues 
sizeErmineWiki <- enrichmentResultWiki$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 6 enriched Wiki pathways with ermineR

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

There are 107 up-regulated and 7 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| EGF/EGFR Signaling Pathway%WikiPathways\_20180810%WP437%Homo sapiens                                                      |  0.0016180|  0.3989701|  1.790294|             0|
| TGF-beta Receptor Signaling%WikiPathways\_20180810%WP560%Homo sapiens                                                     |  0.0016180|  0.6176381|  2.262197|             0|
| Hematopoietic Stem Cell Differentiation%WikiPathways\_20180810%WP2849%Homo sapiens                                        |  0.0016180|  0.5777307|  2.007356|             0|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                                  |  0.0016180|  0.5097105|  2.315118|             0|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                                |  0.0016180|  0.6442048|  2.484615|             0|
| Wnt/beta-catenin Signaling Pathway in Leukemia%WikiPathways\_20180810%WP3658%Homo sapiens                                 |  0.0016180|  0.6734948|  2.061697|             0|
| TGF-beta Signaling Pathway%WikiPathways\_20180810%WP366%Homo sapiens                                                      |  0.0016180|  0.4500663|  1.961037|             0|
| Corticotropin-releasing hormone signaling pathway%WikiPathways\_20180810%WP2355%Homo sapiens                              |  0.0016180|  0.5046611|  2.023612|             0|
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                                   |  0.0016180|  0.5158416|  2.110404|             0|
| Myometrial Relaxation and Contraction Pathways%WikiPathways\_20180810%WP289%Homo sapiens                                  |  0.0016180|  0.4348730|  1.899397|             0|
| DNA Damage Response (only ATM dependent)%WikiPathways\_20180810%WP710%Homo sapiens                                        |  0.0016180|  0.4507485|  1.877629|             0|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                                    |  0.0016180|  0.5203722|  2.217944|             0|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                                      |  0.0016180|  0.4622360|  1.932371|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                                  |  0.0016180|  0.4089892|  1.808893|             0|
| miRNA targets in ECM and membrane receptors%WikiPathways\_20180810%WP2911%Homo sapiens                                    |  0.0016180|  0.8041657|  2.427921|             0|
| Nanoparticle-mediated activation of receptor signaling%WikiPathways\_20180810%WP2643%Homo sapiens                         |  0.0016180|  0.6690084|  2.146250|             0|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                                       |  0.0016180|  0.4296359|  1.872017|             0|
| Integrin-mediated Cell Adhesion%WikiPathways\_20180810%WP185%Homo sapiens                                                 |  0.0016180|  0.4787230|  1.958358|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                              |  0.0016180|  0.7275293|  2.567705|             0|
| Aryl Hydrocarbon Receptor Pathway%WikiPathways\_20180810%WP2873%Homo sapiens                                              |  0.0016180|  0.6066411|  1.979869|             0|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                              |  0.0016180|  0.5388780|  2.229630|             0|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                                                       |  0.0016180|  0.6271042|  2.355461|             0|
| RANKL/RANK (Receptor activator of NFKB (ligand)) Signaling Pathway%WikiPathways\_20180810%WP2018%Homo sapiens             |  0.0016180|  0.5374950|  2.001902|             0|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                           |  0.0016180|  0.5332985|  2.226291|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                             |  0.0016180|  0.6085853|  2.219721|             0|
| Apoptosis-related network due to altered Notch3 in ovarian cancer%WikiPathways\_20180810%WP2864%Homo sapiens              |  0.0016180|  0.5388780|  1.965474|             0|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                                |  0.0016180|  0.6437870|  2.083805|             0|
| Hypothesized Pathways in Pathogenesis of Cardiovascular Disease%WikiPathways\_20180810%WP3668%Homo sapiens                |  0.0016180|  0.6625199|  2.069017|             0|
| ErbB Signaling Pathway%WikiPathways\_20180810%WP673%Homo sapiens                                                          |  0.0016180|  0.5158083|  1.905061|             0|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                                 |  0.0016180|  0.4246874|  2.026805|             0|
| VEGFA-VEGFR2 Signaling Pathway%WikiPathways\_20180810%WP3888%Homo sapiens                                                 |  0.0016180|  0.4759047|  2.226352|             0|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens            |  0.0016180|  0.5700028|  2.290323|             0|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                                                     |  0.0016180|  0.4132502|  1.978287|             0|
| Breast cancer pathway%WikiPathways\_20180810%WP4262%Homo sapiens                                                          |  0.0016180|  0.4106147|  1.793444|             0|
| Ras Signaling%WikiPathways\_20180810%WP4223%Homo sapiens                                                                  |  0.0016180|  0.3830946|  1.714445|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens                      |  0.0016180|  0.4829492|  2.130409|             0|
| miR-509-3p alteration of YAP1/ECM axis%WikiPathways\_20180810%WP3967%Homo sapiens                                         |  0.0023027|  0.6964823|  1.990981|             1|
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens                   |  0.0023027|  0.5036266|  1.942423|             1|
| Prolactin Signaling Pathway%WikiPathways\_20180810%WP2037%Homo sapiens                                                    |  0.0023027|  0.4870824|  1.927017|             1|
| MicroRNAs in cardiomyocyte hypertrophy%WikiPathways\_20180810%WP1544%Homo sapiens                                         |  0.0023027|  0.4495802|  1.818231|             1|
| Association Between Physico-Chemical Features and Toxicity Associated Pathways%WikiPathways\_20180810%WP3680%Homo sapiens |  0.0023027|  0.5103733|  1.955828|             1|
| IL-6 signaling pathway%WikiPathways\_20180810%WP364%Homo sapiens                                                          |  0.0023027|  0.5356893|  1.861281|             1|
| Human Thyroid Stimulating Hormone (TSH) signaling pathway%WikiPathways\_20180810%WP2032%Homo sapiens                      |  0.0023027|  0.5014604|  1.921673|             1|
| Canonical and Non-Canonical TGF-B signaling%WikiPathways\_20180810%WP3874%Homo sapiens                                    |  0.0023027|  0.7019356|  2.006570|             1|
| Splicing factor NOVA regulated synaptic proteins%WikiPathways\_20180810%WP4148%Homo sapiens                               |  0.0023027|  0.5451452|  1.894136|             1|
| Interleukin-11 Signaling Pathway%WikiPathways\_20180810%WP2332%Homo sapiens                                               |  0.0023027|  0.5417170|  1.911909|             1|
| IL-3 Signaling Pathway%WikiPathways\_20180810%WP286%Homo sapiens                                                          |  0.0023027|  0.5416330|  1.881933|             1|
| Chromosomal and microsatellite instability in colorectal cancer %WikiPathways\_20180810%WP4216%Homo sapiens               |  0.0023027|  0.4741848|  1.866070|             1|
| MAPK Signaling Pathway%WikiPathways\_20180810%WP382%Homo sapiens                                                          |  0.0030968|  0.3575754|  1.666351|             2|
| Photodynamic therapy-induced AP-1 survival signaling.%WikiPathways\_20180810%WP3611%Homo sapiens                          |  0.0031652|  0.5379339|  1.931921|             2|
| T-Cell antigen Receptor (TCR) Signaling Pathway%WikiPathways\_20180810%WP69%Homo sapiens                                  |  0.0031652|  0.4604169|  1.815810|             2|
| Simplified Interaction Map Between LOXL4 and Oxidative Stress Pathway%WikiPathways\_20180810%WP3670%Homo sapiens          |  0.0032028|  0.6747657|  1.955323|             2|
| Angiopoietin Like Protein 8 Regulatory Pathway%WikiPathways\_20180810%WP3915%Homo sapiens                                 |  0.0038004|  0.3909708|  1.696635|             3|
| Leptin signaling pathway%WikiPathways\_20180810%WP2034%Homo sapiens                                                       |  0.0038004|  0.4436301|  1.749605|             3|
| Regulation of toll-like receptor signaling pathway%WikiPathways\_20180810%WP1449%Homo sapiens                             |  0.0038004|  0.4109383|  1.723356|             3|
| MET in type 1 papillary renal cell carcinoma%WikiPathways\_20180810%WP4205%Homo sapiens                                   |  0.0038004|  0.4875766|  1.824142|             3|
| TGF-B Signaling in Thyroid Cells for Epithelial-Mesenchymal Transition%WikiPathways\_20180810%WP3859%Homo sapiens         |  0.0038504|  0.6823482|  1.950577|             3|
| Calcium Regulation in the Cardiac Cell%WikiPathways\_20180810%WP536%Homo sapiens                                          |  0.0043211|  0.3874143|  1.673046|             4|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                                     |  0.0043211|  0.3887124|  1.678652|             4|
| MAPK Cascade%WikiPathways\_20180810%WP422%Homo sapiens                                                                    |  0.0053306|  0.5826328|  1.901514|             5|
| Aryl Hydrocarbon Receptor%WikiPathways\_20180810%WP2586%Homo sapiens                                                      |  0.0053306|  0.5267824|  1.845973|             5|
| Angiogenesis%WikiPathways\_20180810%WP1539%Homo sapiens                                                                   |  0.0061487|  0.5959847|  1.861231|             6|
| Toll-like Receptor Signaling Pathway%WikiPathways\_20180810%WP75%Homo sapiens                                             |  0.0064183|  0.4319551|  1.685682|             7|
| Oncostatin M Signaling Pathway%WikiPathways\_20180810%WP2374%Homo sapiens                                                 |  0.0064183|  0.4495244|  1.722646|             7|
| Signaling Pathways in Glioblastoma%WikiPathways\_20180810%WP2261%Homo sapiens                                             |  0.0064183|  0.4090615|  1.654362|             7|
| Hedgehog Signaling Pathway%WikiPathways\_20180810%WP4249%Homo sapiens                                                     |  0.0064183|  0.5251672|  1.813799|             7|
| Pancreatic adenocarcinoma pathway%WikiPathways\_20180810%WP4263%Homo sapiens                                              |  0.0064183|  0.4176090|  1.708353|             7|
| Chemokine signaling pathway%WikiPathways\_20180810%WP3929%Homo sapiens                                                    |  0.0067322|  0.3765263|  1.623218|             8|
| Insulin Signaling%WikiPathways\_20180810%WP481%Homo sapiens                                                               |  0.0067322|  0.3456344|  1.548477|             8|
| Hepatitis C and Hepatocellular Carcinoma%WikiPathways\_20180810%WP3646%Homo sapiens                                       |  0.0076714|  0.4909576|  1.763212|             9|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                                        |  0.0076714|  0.5333665|  1.805131|             9|
| B Cell Receptor Signaling Pathway%WikiPathways\_20180810%WP23%Homo sapiens                                                |  0.0080560|  0.4121499|  1.686020|            10|
| BMP Signaling Pathway in Eyelid Development%WikiPathways\_20180810%WP3927%Homo sapiens                                    |  0.0083884|  0.6462598|  1.872719|            10|
| Spinal Cord Injury%WikiPathways\_20180810%WP2431%Homo sapiens                                                             |  0.0085566|  0.4004925|  1.649952|            11|
| Alpha 6 Beta 4 signaling pathway%WikiPathways\_20180810%WP244%Homo sapiens                                                |  0.0087464|  0.5398652|  1.791270|            11|
| Structural Pathway of Interleukin 1 (IL-1)%WikiPathways\_20180810%WP2637%Homo sapiens                                     |  0.0091246|  0.4708771|  1.724662|            12|
| IL-1 signaling pathway%WikiPathways\_20180810%WP195%Homo sapiens                                                          |  0.0096733|  0.4556654|  1.697127|            13|
| Arrhythmogenic Right Ventricular Cardiomyopathy%WikiPathways\_20180810%WP2118%Homo sapiens                                |  0.0101465|  0.4427742|  1.688542|            14|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP363%Homo sapiens                                                           |  0.0101465|  0.4774436|  1.741402|            14|
| Signal Transduction of S1P Receptor%WikiPathways\_20180810%WP26%Homo sapiens                                              |  0.0109619|  0.5631428|  1.771879|            15|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                                       |  0.0119760|  0.3531775|  1.552604|            18|
| Kit receptor signaling pathway%WikiPathways\_20180810%WP304%Homo sapiens                                                  |  0.0130273|  0.4431447|  1.679734|            19|
| Estrogen signaling pathway%WikiPathways\_20180810%WP712%Homo sapiens                                                      |  0.0131913|  0.6054016|  1.804734|            19|
| Canonical and Non-canonical Notch signaling%WikiPathways\_20180810%WP3845%Homo sapiens                                    |  0.0157039|  0.5497633|  1.729782|            23|
| TNF related weak inducer of apoptosis (TWEAK) Signaling Pathway%WikiPathways\_20180810%WP2036%Homo sapiens                |  0.0157364|  0.4870180|  1.706629|            24|
| Regulation of Actin Cytoskeleton%WikiPathways\_20180810%WP51%Homo sapiens                                                 |  0.0157364|  0.3496077|  1.525229|            25|
| Ebola Virus Pathway on Host%WikiPathways\_20180810%WP4217%Homo sapiens                                                    |  0.0179414|  0.3590918|  1.530531|            29|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                                 |  0.0199882|  0.3040116|  1.429692|            34|
| Hypertrophy Model%WikiPathways\_20180810%WP516%Homo sapiens                                                               |  0.0199894|  0.6215551|  1.776792|            31|
| Lipid Metabolism Pathway%WikiPathways\_20180810%WP3965%Homo sapiens                                                       |  0.0208421|  0.5429885|  1.708465|            33|
| TNF alpha Signaling Pathway%WikiPathways\_20180810%WP231%Homo sapiens                                                     |  0.0208724|  0.3741306|  1.530492|            35|
| IL-5 Signaling Pathway%WikiPathways\_20180810%WP127%Homo sapiens                                                          |  0.0217278|  0.4821209|  1.665128|            36|
| Inflammatory Response Pathway%WikiPathways\_20180810%WP453%Homo sapiens                                                   |  0.0231073|  0.5757371|  1.738254|            38|
| Notch Signaling Pathway%WikiPathways\_20180810%WP61%Homo sapiens                                                          |  0.0251860|  0.4203454|  1.593313|            43|
| Brain-Derived Neurotrophic Factor (BDNF) signaling pathway%WikiPathways\_20180810%WP2380%Homo sapiens                     |  0.0262984|  0.3362748|  1.478297|            48|
| Endometrial cancer%WikiPathways\_20180810%WP4155%Homo sapiens                                                             |  0.0286385|  0.4077715|  1.568353|            51|
| Extracellular vesicle-mediated signaling in recipient cells%WikiPathways\_20180810%WP2870%Homo sapiens                    |  0.0311300|  0.5064609|  1.652915|            55|
| IL-4 Signaling Pathway%WikiPathways\_20180810%WP395%Homo sapiens                                                          |  0.0311300|  0.4267597|  1.576173|            56|
| Thymic Stromal LymphoPoietin (TSLP) Signaling Pathway%WikiPathways\_20180810%WP2203%Homo sapiens                          |  0.0336476|  0.4539192|  1.610102|            61|
| NOTCH1 regulation of human endothelial cell calcification%WikiPathways\_20180810%WP3413%Homo sapiens                      |  0.0360853|  0.6045801|  1.699446|            64|
| PDGF Pathway%WikiPathways\_20180810%WP2526%Homo sapiens                                                                   |  0.0367325|  0.4616787|  1.594525|            68|
| RAC1/PAK1/p38/MMP2 Pathway%WikiPathways\_20180810%WP3303%Homo sapiens                                                     |  0.0379465|  0.3955655|  1.529476|            72|
| Serotonin Receptor 2 and ELK-SRF/GATA4 signaling%WikiPathways\_20180810%WP732%Homo sapiens                                |  0.0379465|  0.5713155|  1.678535|            70|
| EPO Receptor Signaling%WikiPathways\_20180810%WP581%Homo sapiens                                                          |  0.0379465|  0.5277873|  1.648254|            71|
| Signaling of Hepatocyte Growth Factor Receptor%WikiPathways\_20180810%WP313%Homo sapiens                                  |  0.0426382|  0.4700848|  1.590960|            82|
| Imatinib and Chronic Myeloid Leukemia%WikiPathways\_20180810%WP3640%Homo sapiens                                          |  0.0434775|  0.5643316|  1.658017|            83|
| Physiological and Pathological Hypertrophy of the Heart%WikiPathways\_20180810%WP1528%Homo sapiens                        |  0.0447517|  0.5178686|  1.617278|            87|

``` r
upPathwaysWiki75748_096 <- upPathwaysWiki
save(upPathwaysWiki75748_096, file="upPathwaysWiki75748_096.Rdata")
```
