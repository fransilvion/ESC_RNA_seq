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

| Pathway                                    | Corrected p-value | Corrected MF p-value |
|:-------------------------------------------|:------------------|:---------------------|
| Extracellular matrix organization          | 0.0000000         | 0.0000000            |
| TCF dependent signaling in response to WNT | 0.0000000         | 0.0309500            |
| Dopamine Neurotransmitter Release Cycle    | 0.1583212         | 0.2174436            |
| Transcriptional regulation by small RNAs   | 0.1555049         | 0.2174865            |
| ECM proteoglycans                          | 0.1606766         | 0.2215368            |
| Meiosis                                    | 0.1857000         | 0.2218083            |
| Neurotransmitter release cycle             | 0.1535120         | 0.2228400            |
| Reproduction                               | 0.1623156         | 0.2238936            |
| Class A 1 (Rhodopsin-like receptors)       | 0.1599083         | 0.2257529            |
| NCAM1 interactions                         | 0.1645977         | 0.2281457            |

``` r
#based on corrected MFPvalues 
sizeErmineReactome <- enrichmentResultReactome$results %>% filter(CorrectedMFPvalue <= 0.05)
```

There are 2 enriched Reactome pathways with ermineR

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

There are 86 up-regulated and 46 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                 |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                         |  0.0083810|  0.5294765|  2.074165|             0|
| FACTORS INVOLVED IN MEGAKARYOCYTE DEVELOPMENT AND PLATELET PRODUCTION%REACTOME%R-HSA-983231.2           |  0.0083810|  0.3952818|  1.748819|             0|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                    |  0.0083810|  0.7081415|  2.525339|             0|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                            |  0.0083810|  0.5435469|  2.017068|             0|
| NRAGE SIGNALS DEATH THROUGH JNK%REACTOME DATABASE ID RELEASE 65%193648                                  |  0.0083810|  0.5352238|  2.029486|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                     |  0.0083810|  0.6786716|  2.595703|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                               |  0.0083810|  0.4685203|  2.268934|             0|
| SIGNALING BY NOTCH3%REACTOME DATABASE ID RELEASE 65%9012852                                             |  0.0083810|  0.5592139|  2.064059|             0|
| EPH-EPHRIN MEDIATED REPULSION OF CELLS%REACTOME%R-HSA-3928665.3                                         |  0.0083810|  0.5557962|  2.071194|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                             |  0.0083810|  0.6177198|  2.480593|             0|
| EPH-EPHRIN SIGNALING%REACTOME%R-HSA-2682334.1                                                           |  0.0083810|  0.4808375|  2.025460|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                               |  0.0083810|  0.5624801|  2.059117|             0|
| SIGNALING BY PDGF%REACTOME DATABASE ID RELEASE 65%186797                                                |  0.0083810|  0.6334639|  2.373283|             0|
| RHO GTPASE CYCLE%REACTOME DATABASE ID RELEASE 65%194840                                                 |  0.0083810|  0.4909874|  2.196273|             0|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                            |  0.0083810|  0.4468272|  1.921162|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090    |  0.0083810|  0.6219693|  2.295689|             0|
| SIGNALING BY WNT%REACTOME DATABASE ID RELEASE 65%195721                                                 |  0.0083810|  0.3292186|  1.639749|             0|
| DOWNSTREAM SIGNALING OF ACTIVATED FGFR1%REACTOME%R-HSA-5654687.2                                        |  0.0121110|  0.6662548|  1.999059|             1|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                              |  0.0121110|  0.5877192|  1.982181|             1|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                      |  0.0121110|  0.4564114|  1.899185|             1|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                              |  0.0121110|  0.5225696|  1.981503|             1|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                             |  0.0121110|  0.4546404|  1.919740|             1|
| SIGNALING BY NTRKS%REACTOME DATABASE ID RELEASE 65%166520                                               |  0.0142229|  0.4226017|  1.806720|             2|
| CELLULAR SENESCENCE%REACTOME%R-HSA-2559583.2                                                            |  0.0142229|  0.3550414|  1.665656|             2|
| SEMA4D IN SEMAPHORIN SIGNALING%REACTOME%R-HSA-400685.2                                                  |  0.0142229|  0.6456456|  2.038662|             2|
| SIGNALING BY NTRK1 (TRKA)%REACTOME%R-HSA-187037.2                                                       |  0.0142229|  0.4430739|  1.814763|             2|
| SIGNALING BY RAS MUTANTS%REACTOME%R-HSA-6802949.1                                                       |  0.0142229|  0.5381687|  1.986381|             2|
| NON-INTEGRIN MEMBRANE-ECM INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000171                          |  0.0142229|  0.5647645|  2.014035|             2|
| SIGNALING BY FGFR1%REACTOME DATABASE ID RELEASE 65%5654736                                              |  0.0168467|  0.5538595|  1.951702|             3|
| NOTCH3 INTRACELLULAR DOMAIN REGULATES TRANSCRIPTION%REACTOME%R-HSA-9013508.1                            |  0.0168467|  0.6446422|  2.011185|             3|
| SEMA4D INDUCED CELL MIGRATION AND GROWTH-CONE COLLAPSE%REACTOME%R-HSA-416572.1                          |  0.0168467|  0.6602435|  1.981022|             3|
| ONCOGENIC MAPK SIGNALING%REACTOME%R-HSA-6802957.2                                                       |  0.0168467|  0.4635441|  1.861466|             3|
| O-LINKED GLYCOSYLATION%REACTOME%R-HSA-5173105.3                                                         |  0.0185265|  0.4268245|  1.771336|             4|
| G ALPHA (12 13) SIGNALLING EVENTS%REACTOME%R-HSA-416482.3                                               |  0.0185265|  0.4293278|  1.761312|             4|
| MET PROMOTES CELL MOTILITY%REACTOME DATABASE ID RELEASE 65%8875878                                      |  0.0187234|  0.6095724|  2.002316|             4|
| TCF DEPENDENT SIGNALING IN RESPONSE TO WNT%REACTOME%R-HSA-201681.1                                      |  0.0192207|  0.3211413|  1.539724|             5|
| PRE-NOTCH PROCESSING IN GOLGI%REACTOME%R-HSA-1912420.2                                                  |  0.0199896|  0.6905442|  1.974675|             5|
| LDL CLEARANCE%REACTOME%R-HSA-8964038.1                                                                  |  0.0199896|  0.6550377|  1.936789|             5|
| PRE-NOTCH EXPRESSION AND PROCESSING%REACTOME DATABASE ID RELEASE 65%1912422                             |  0.0209637|  0.4081292|  1.737474|             6|
| GLYCOSAMINOGLYCAN METABOLISM%REACTOME DATABASE ID RELEASE 65%1630316                                    |  0.0209637|  0.3918191|  1.709541|             6|
| SIGNALING BY VEGF%REACTOME DATABASE ID RELEASE 65%194138                                                |  0.0209637|  0.4062985|  1.729680|             6|
| DOWNSTREAM SIGNAL TRANSDUCTION%REACTOME DATABASE ID RELEASE 65%186763                                   |  0.0213012|  0.5913993|  1.942621|             6|
| HS-GAG BIOSYNTHESIS%REACTOME%R-HSA-2022928.1                                                            |  0.0213012|  0.5863204|  1.962934|             6|
| CELL-CELL COMMUNICATION%REACTOME%R-HSA-1500931.3                                                        |  0.0217038|  0.3838493|  1.662849|             7|
| CHONDROITIN SULFATE DERMATAN SULFATE METABOLISM%REACTOME DATABASE ID RELEASE 65%1793185                 |  0.0217038|  0.5015270|  1.841936|             7|
| N-GLYCAN TRIMMING IN THE ER AND CALNEXIN CALRETICULIN CYCLE%REACTOME%R-HSA-532668.2                     |  0.0217038|  0.5402788|  1.888971|             7|
| CHROMATIN ORGANIZATION%REACTOME%R-HSA-4839726.2                                                         |  0.0220338|  0.3022272|  1.476339|             8|
| CHROMATIN MODIFYING ENZYMES%REACTOME%R-HSA-3247509.4                                                    |  0.0220338|  0.3022272|  1.476339|             8|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                       |  0.0243093|  0.6333236|  1.900251|             9|
| INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING%REACTOME DATABASE ID RELEASE 65%6785807                      |  0.0243093|  0.4048856|  1.698653|            10|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                                    |  0.0243093|  0.5829550|  1.914883|             9|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                         |  0.0243093|  0.5166638|  1.842501|             9|
| SIGNALING BY TGF-BETA RECEPTOR COMPLEX%REACTOME DATABASE ID RELEASE 65%170834                           |  0.0243093|  0.4227119|  1.726517|            10|
| CLATHRIN-MEDIATED ENDOCYTOSIS%REACTOME%R-HSA-8856828.3                                                  |  0.0243093|  0.3553193|  1.589406|            10|
| INTRACELLULAR SIGNALING BY SECOND MESSENGERS%REACTOME DATABASE ID RELEASE 65%9006925                    |  0.0243830|  0.2992943|  1.461396|            11|
| TRANSCRIPTIONAL REGULATION BY THE AP-2 (TFAP2) FAMILY OF TRANSCRIPTION FACTORS%REACTOME%R-HSA-8864260.3 |  0.0243830|  0.5613401|  1.879303|            10|
| PARADOXICAL ACTIVATION OF RAF SIGNALING BY KINASE INACTIVE BRAF%REACTOME DATABASE ID RELEASE 65%6802955 |  0.0255765|  0.5317474|  1.844759|            11|
| SIGNALING BY MODERATE KINASE ACTIVITY BRAF MUTANTS%REACTOME DATABASE ID RELEASE 65%6802946              |  0.0255765|  0.5317474|  1.844759|            11|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                          |  0.0277383|  0.6638013|  1.898201|            12|
| ACTIVATION OF BH3-ONLY PROTEINS%REACTOME DATABASE ID RELEASE 65%114452                                  |  0.0308643|  0.5633874|  1.866249|            14|
| TP53 REGULATES TRANSCRIPTION OF CELL CYCLE GENES%REACTOME DATABASE ID RELEASE 65%6791312                |  0.0320367|  0.4637610|  1.720988|            15|
| LAMININ INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000157                                            |  0.0322724|  0.5939494|  1.831096|            15|
| VEGFA-VEGFR2 PATHWAY%REACTOME%R-HSA-4420097.3                                                           |  0.0323062|  0.3945066|  1.653474|            16|
| INTRINSIC PATHWAY FOR APOPTOSIS%REACTOME%R-HSA-109606.2                                                 |  0.0329625|  0.4877898|  1.769922|            16|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                           |  0.0329683|  0.5607128|  1.824955|            16|
| REGULATION OF FZD BY UBIQUITINATION%REACTOME%R-HSA-4641263.2                                            |  0.0341146|  0.6037169|  1.811418|            17|
| CALMODULIN INDUCED EVENTS%REACTOME DATABASE ID RELEASE 65%111933                                        |  0.0341146|  0.5727439|  1.786873|            17|
| CAM PATHWAY%REACTOME DATABASE ID RELEASE 65%111997                                                      |  0.0341146|  0.5727439|  1.786873|            17|
| L1CAM INTERACTIONS%REACTOME%R-HSA-373760.2                                                              |  0.0356714|  0.3913974|  1.644468|            19|
| ASPARAGINE N-LINKED GLYCOSYLATION%REACTOME%R-HSA-446203.4                                               |  0.0356714|  0.2836056|  1.402221|            20|
| DISASSEMBLY OF THE DESTRUCTION COMPLEX AND RECRUITMENT OF AXIN TO THE MEMBRANE%REACTOME%R-HSA-4641262.3 |  0.0378939|  0.5488056|  1.817946|            20|
| HEPARAN SULFATE HEPARIN (HS-GAG) METABOLISM%REACTOME%R-HSA-1638091.1                                    |  0.0395338|  0.4529951|  1.705219|            22|
| SIGNALING BY FGFR3 IN DISEASE%REACTOME DATABASE ID RELEASE 65%5655332                                   |  0.0430868|  0.6282799|  1.796624|            26|
| TGF-BETA RECEPTOR SIGNALING ACTIVATES SMADS%REACTOME DATABASE ID RELEASE 65%2173789                     |  0.0430868|  0.5199507|  1.765433|            26|
| DAG AND IP3 SIGNALING%REACTOME DATABASE ID RELEASE 65%1489509                                           |  0.0430868|  0.5378068|  1.781512|            26|
| CELL DEATH SIGNALLING VIA NRAGE, NRIF AND NADE%REACTOME DATABASE ID RELEASE 65%204998                   |  0.0430868|  0.3994354|  1.608117|            27|
| SIGNALING BY FGFR3 POINT MUTANTS IN CANCER%REACTOME%R-HSA-8853338.2                                     |  0.0430868|  0.6282799|  1.796624|            26|
| REGULATION OF PTEN GENE TRANSCRIPTION%REACTOME DATABASE ID RELEASE 65%8943724                           |  0.0448197|  0.4358255|  1.674393|            28|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                         |  0.0448197|  0.3887661|  1.611853|            29|
| ESR-MEDIATED SIGNALING%REACTOME DATABASE ID RELEASE 65%8939211                                          |  0.0459731|  0.3294227|  1.490957|            31|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                              |  0.0459731|  0.4581408|  1.682594|            29|
| INTRA-GOLGI AND RETROGRADE GOLGI-TO-ER TRAFFIC%REACTOME DATABASE ID RELEASE 65%6811442                  |  0.0466066|  0.3121420|  1.462149|            32|
| PRE-NOTCH TRANSCRIPTION AND TRANSLATION%REACTOME%R-HSA-1912408.3                                        |  0.0478051|  0.3876789|  1.602692|            32|
| CA-DEPENDENT EVENTS%REACTOME DATABASE ID RELEASE 65%111996                                              |  0.0489977|  0.5519660|  1.759067|            32|
| GROWTH HORMONE RECEPTOR SIGNALING%REACTOME%R-HSA-982772.1                                               |  0.0490638|  0.5877482|  1.763504|            33|
| RHO GTPASES ACTIVATE ROCKS%REACTOME DATABASE ID RELEASE 65%5627117                                      |  0.0490638|  0.5869809|  1.761202|            33|

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
| Kegg wnt signaling pathway                                | 0.0000000         | 0.0000000            |
| Kegg pathways in cancer                                   | 0.0324000         | 0.0000000            |
| Kegg mapk signaling pathway                               | 0.0081000         | 0.0000000            |
| Kegg focal adhesion                                       | 0.0432000         | 0.0040500            |
| Kegg basal cell carcinoma                                 | 0.0453600         | 0.0194400            |
| Kegg melanogenesis                                        | 0.1275750         | 0.0270000            |
| Kegg alzheimers disease                                   | 0.1472727         | 0.0671143            |
| Kegg arrhythmogenic right ventricular cardiomyopathy arvc | 0.1944000         | 0.0749250            |
| Kegg parkinsons disease                                   | 0.1566000         | 0.0774000            |
| Kegg ecm receptor interaction                             | 0.1365429         | 0.0793800            |

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

There are 37 up-regulated and 11 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                                        |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                                 |  0.0020785|  0.4056368|  1.976130|             0|
| KEGG\_ENDOCYTOSIS                                              |  0.0020785|  0.4037747|  1.882619|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY                                  |  0.0020785|  0.4299674|  1.929426|             0|
| KEGG\_AXON\_GUIDANCE                                           |  0.0020785|  0.4493756|  1.996689|             0|
| KEGG\_FOCAL\_ADHESION                                          |  0.0020785|  0.5544194|  2.600340|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                               |  0.0020785|  0.6108126|  2.464211|             0|
| KEGG\_NEUROTROPHIN\_SIGNALING\_PATHWAY                         |  0.0020785|  0.4347978|  1.923537|             0|
| KEGG\_REGULATION\_OF\_ACTIN\_CYTOSKELETON                      |  0.0020785|  0.3542808|  1.662650|             0|
| KEGG\_INSULIN\_SIGNALING\_PATHWAY                              |  0.0020785|  0.4063390|  1.805467|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                                     |  0.0020785|  0.3939500|  1.970658|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA                                   |  0.0020785|  0.5864775|  2.152822|             0|
| KEGG\_CHRONIC\_MYELOID\_LEUKEMIA                               |  0.0020785|  0.4665465|  1.903322|             0|
| KEGG\_SMALL\_CELL\_LUNG\_CANCER                                |  0.0020785|  0.4679651|  1.951182|             0|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY                             |  0.0025183|  0.5263141|  1.925664|             1|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY                            |  0.0025183|  0.4547765|  1.870971|             1|
| KEGG\_COLORECTAL\_CANCER                                       |  0.0025183|  0.5037259|  1.997648|             1|
| KEGG\_RENAL\_CELL\_CARCINOMA                                   |  0.0025183|  0.4860469|  1.964766|             1|
| KEGG\_PANCREATIC\_CANCER                                       |  0.0025183|  0.4877522|  1.976654|             1|
| KEGG\_BLADDER\_CANCER                                          |  0.0025183|  0.5453271|  1.948163|             1|
| KEGG\_ARRHYTHMOGENIC\_RIGHT\_VENTRICULAR\_CARDIOMYOPATHY\_ARVC |  0.0025183|  0.5006561|  1.979505|             1|
| KEGG\_DILATED\_CARDIOMYOPATHY                                  |  0.0025183|  0.4775598|  1.930459|             1|
| KEGG\_ERBB\_SIGNALING\_PATHWAY                                 |  0.0033367|  0.4317931|  1.797462|             2|
| KEGG\_MELANOGENESIS                                            |  0.0033367|  0.4276664|  1.777422|             2|
| KEGG\_ENDOMETRIAL\_CANCER                                      |  0.0033367|  0.4891885|  1.866184|             2|
| KEGG\_CHEMOKINE\_SIGNALING\_PATHWAY                            |  0.0042309|  0.3725686|  1.663794|             3|
| KEGG\_GNRH\_SIGNALING\_PATHWAY                                 |  0.0060322|  0.4131290|  1.717003|             5|
| KEGG\_NON\_SMALL\_CELL\_LUNG\_CANCER                           |  0.0060322|  0.4732919|  1.825857|             5|
| KEGG\_NOTCH\_SIGNALING\_PATHWAY                                |  0.0092461|  0.4997174|  1.834345|             9|
| KEGG\_HYPERTROPHIC\_CARDIOMYOPATHY\_HCM                        |  0.0102438|  0.4366952|  1.743585|            11|
| KEGG\_GAP\_JUNCTION                                            |  0.0142218|  0.3969463|  1.635315|            16|
| KEGG\_TYPE\_II\_DIABETES\_MELLITUS                             |  0.0249594|  0.4709110|  1.715529|            30|
| KEGG\_ACUTE\_MYELOID\_LEUKEMIA                                 |  0.0263453|  0.4398111|  1.692343|            33|
| KEGG\_GLYCOSAMINOGLYCAN\_BIOSYNTHESIS\_CHONDROITIN\_SULFATE    |  0.0282403|  0.5813376|  1.764570|            35|
| KEGG\_VASOPRESSIN\_REGULATED\_WATER\_REABSORPTION              |  0.0284813|  0.4731124|  1.690178|            37|
| KEGG\_B\_CELL\_RECEPTOR\_SIGNALING\_PATHWAY                    |  0.0292085|  0.3984707|  1.607557|            40|
| KEGG\_ADHERENS\_JUNCTION                                       |  0.0307605|  0.3871865|  1.579565|            43|
| KEGG\_PROSTATE\_CANCER                                         |  0.0358206|  0.3672851|  1.538065|            54|

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
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       | 0.000000          | 0.00000              |
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            | 0.000000          | 0.00000              |
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                | 0.000000          | 0.00990              |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  | 0.011880          | 0.01485              |
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens        | 0.007425          | 0.02376              |
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens | 0.009900          | 0.04950              |

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

There are 111 up-regulated and 6 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                                          |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| EGF/EGFR Signaling Pathway%WikiPathways\_20180810%WP437%Homo sapiens                                                             |  0.0015650|  0.3903432|  1.790745|             0|
| TGF-beta Receptor Signaling%WikiPathways\_20180810%WP560%Homo sapiens                                                            |  0.0015650|  0.5486288|  2.070489|             0|
| Integrated Breast Cancer Pathway%WikiPathways\_20180810%WP1984%Homo sapiens                                                      |  0.0015650|  0.4025444|  1.818882|             0|
| TGF-B Signaling in Thyroid Cells for Epithelial-Mesenchymal Transition%WikiPathways\_20180810%WP3859%Homo sapiens                |  0.0015650|  0.7395393|  2.142542|             0|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                                         |  0.0015650|  0.5363231|  2.498797|             0|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                                       |  0.0015650|  0.5407072|  2.140375|             0|
| TGF-beta Signaling Pathway%WikiPathways\_20180810%WP366%Homo sapiens                                                             |  0.0015650|  0.4403416|  1.956455|             0|
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                                          |  0.0015650|  0.5582691|  2.333432|             0|
| Arrhythmogenic Right Ventricular Cardiomyopathy%WikiPathways\_20180810%WP2118%Homo sapiens                                       |  0.0015650|  0.4991551|  1.954295|             0|
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens                          |  0.0015650|  0.5114115|  2.024409|             0|
| Leptin signaling pathway%WikiPathways\_20180810%WP2034%Homo sapiens                                                              |  0.0015650|  0.4841693|  1.967392|             0|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                                           |  0.0015650|  0.5048453|  2.201281|             0|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                                             |  0.0015650|  0.4872567|  2.088075|             0|
| Notch Signaling Pathway%WikiPathways\_20180810%WP61%Homo sapiens                                                                 |  0.0015650|  0.4922462|  1.917203|             0|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                                         |  0.0015650|  0.4580451|  2.072627|             0|
| miRNA targets in ECM and membrane receptors%WikiPathways\_20180810%WP2911%Homo sapiens                                           |  0.0015650|  0.7770215|  2.385923|             0|
| Integrin-mediated Cell Adhesion%WikiPathways\_20180810%WP185%Homo sapiens                                                        |  0.0015650|  0.4920657|  2.059595|             0|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                                              |  0.0015650|  0.4208167|  1.886583|             0|
| MAPK Signaling Pathway%WikiPathways\_20180810%WP382%Homo sapiens                                                                 |  0.0015650|  0.4004338|  1.916804|             0|
| Prolactin Signaling Pathway%WikiPathways\_20180810%WP2037%Homo sapiens                                                           |  0.0015650|  0.4912431|  2.000064|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                                     |  0.0015650|  0.6785853|  2.465279|             0|
| Endometrial cancer%WikiPathways\_20180810%WP4155%Homo sapiens                                                                    |  0.0015650|  0.5009958|  1.974217|             0|
| Association Between Physico-Chemical Features and Toxicity Associated Pathways%WikiPathways\_20180810%WP3680%Homo sapiens        |  0.0015650|  0.5023847|  1.974492|             0|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                                  |  0.0015650|  0.4641251|  1.985186|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                                    |  0.0015650|  0.5874748|  2.205979|             0|
| Ebola Virus Pathway on Host%WikiPathways\_20180810%WP4217%Homo sapiens                                                           |  0.0015650|  0.4229853|  1.844346|             0|
| Splicing factor NOVA regulated synaptic proteins%WikiPathways\_20180810%WP4148%Homo sapiens                                      |  0.0015650|  0.6007186|  2.135622|             0|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                                        |  0.0015650|  0.3913055|  1.915134|             0|
| IL-4 Signaling Pathway%WikiPathways\_20180810%WP395%Homo sapiens                                                                 |  0.0015650|  0.4950647|  1.881548|             0|
| IL-3 Signaling Pathway%WikiPathways\_20180810%WP286%Homo sapiens                                                                 |  0.0015650|  0.5810400|  2.065663|             0|
| VEGFA-VEGFR2 Signaling Pathway%WikiPathways\_20180810%WP3888%Homo sapiens                                                        |  0.0015650|  0.4380901|  2.107448|             0|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens                   |  0.0015650|  0.5432686|  2.244267|             0|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                                                            |  0.0015650|  0.3672727|  1.806463|             0|
| Breast cancer pathway%WikiPathways\_20180810%WP4262%Homo sapiens                                                                 |  0.0015650|  0.4285451|  1.910591|             0|
| MET in type 1 papillary renal cell carcinoma%WikiPathways\_20180810%WP4205%Homo sapiens                                          |  0.0015650|  0.5322077|  2.052050|             0|
| Ras Signaling%WikiPathways\_20180810%WP4223%Homo sapiens                                                                         |  0.0015650|  0.3782239|  1.728161|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens                             |  0.0015650|  0.4366771|  1.969274|             0|
| B Cell Receptor Signaling Pathway%WikiPathways\_20180810%WP23%Homo sapiens                                                       |  0.0021057|  0.4320285|  1.808303|             1|
| Angiopoietin Like Protein 8 Regulatory Pathway%WikiPathways\_20180810%WP3915%Homo sapiens                                        |  0.0021057|  0.4055949|  1.796202|             1|
| Wnt/beta-catenin Signaling Pathway in Leukemia%WikiPathways\_20180810%WP3658%Homo sapiens                                        |  0.0021057|  0.6818005|  2.118401|             1|
| miR-509-3p alteration of YAP1/ECM axis%WikiPathways\_20180810%WP3967%Homo sapiens                                                |  0.0021057|  0.7271337|  2.106601|             1|
| Corticotropin-releasing hormone signaling pathway%WikiPathways\_20180810%WP2355%Homo sapiens                                     |  0.0021057|  0.4538540|  1.870616|             1|
| DNA Damage Response (only ATM dependent)%WikiPathways\_20180810%WP710%Homo sapiens                                               |  0.0021057|  0.4218610|  1.799767|             1|
| Bladder Cancer%WikiPathways\_20180810%WP2828%Homo sapiens                                                                        |  0.0021057|  0.6401758|  2.109612|             1|
| Signaling Pathways in Glioblastoma%WikiPathways\_20180810%WP2261%Homo sapiens                                                    |  0.0021057|  0.4426197|  1.835760|             1|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                                                              |  0.0021057|  0.4987816|  1.928162|             1|
| RANKL/RANK (Receptor activator of NFKB (ligand)) Signaling Pathway%WikiPathways\_20180810%WP2018%Homo sapiens                    |  0.0021057|  0.5097552|  1.949425|             1|
| Apoptosis-related network due to altered Notch3 in ovarian cancer%WikiPathways\_20180810%WP2864%Homo sapiens                     |  0.0021057|  0.5095398|  1.913331|             1|
| Androgen receptor signaling pathway%WikiPathways\_20180810%WP138%Homo sapiens                                                    |  0.0021057|  0.4416392|  1.842018|             1|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP363%Homo sapiens                                                                  |  0.0021057|  0.5084747|  1.909332|             1|
| Pancreatic adenocarcinoma pathway%WikiPathways\_20180810%WP4263%Homo sapiens                                                     |  0.0021057|  0.4608888|  1.929101|             1|
| Chromosomal and microsatellite instability in colorectal cancer %WikiPathways\_20180810%WP4216%Homo sapiens                      |  0.0021057|  0.4620539|  1.872339|             1|
| Hepatitis C and Hepatocellular Carcinoma%WikiPathways\_20180810%WP3646%Homo sapiens                                              |  0.0028307|  0.5156964|  1.911342|             2|
| Insulin Signaling%WikiPathways\_20180810%WP481%Homo sapiens                                                                      |  0.0028307|  0.3701844|  1.692323|             2|
| MicroRNAs in cardiomyocyte hypertrophy%WikiPathways\_20180810%WP1544%Homo sapiens                                                |  0.0028307|  0.4283846|  1.776720|             2|
| Human Thyroid Stimulating Hormone (TSH) signaling pathway%WikiPathways\_20180810%WP2032%Homo sapiens                             |  0.0028307|  0.4648067|  1.826802|             2|
| Wnt Signaling Pathway and Pluripotency%WikiPathways\_20180810%WP399%Homo sapiens                                                 |  0.0036290|  0.4098185|  1.723313|             3|
| PDGF Pathway%WikiPathways\_20180810%WP2526%Homo sapiens                                                                          |  0.0036290|  0.5480423|  1.936998|             3|
| IL-6 signaling pathway%WikiPathways\_20180810%WP364%Homo sapiens                                                                 |  0.0036290|  0.5407779|  1.922526|             3|
| Estrogen signaling pathway%WikiPathways\_20180810%WP712%Homo sapiens                                                             |  0.0036404|  0.6564416|  1.987982|             3|
| Alpha 6 Beta 4 signaling pathway%WikiPathways\_20180810%WP244%Homo sapiens                                                       |  0.0044391|  0.5713602|  1.930198|             4|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                                               |  0.0052012|  0.5585430|  1.932829|             5|
| Interleukin-11 Signaling Pathway%WikiPathways\_20180810%WP2332%Homo sapiens                                                      |  0.0057404|  0.5155563|  1.872999|             6|
| T-Cell antigen Receptor (TCR) Signaling Pathway%WikiPathways\_20180810%WP69%Homo sapiens                                         |  0.0057404|  0.4262051|  1.731858|             6|
| Non-small cell lung cancer%WikiPathways\_20180810%WP4255%Homo sapiens                                                            |  0.0057404|  0.4373975|  1.746854|             6|
| IL17 signaling pathway%WikiPathways\_20180810%WP2112%Homo sapiens                                                                |  0.0064347|  0.5926364|  1.909433|             7|
| Oncostatin M Signaling Pathway%WikiPathways\_20180810%WP2374%Homo sapiens                                                        |  0.0064347|  0.4384284|  1.723129|             7|
| Regulation of Actin Cytoskeleton%WikiPathways\_20180810%WP51%Homo sapiens                                                        |  0.0067821|  0.3692273|  1.642033|             8|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                                     |  0.0067821|  0.3942640|  1.668600|             8|
| IL-5 Signaling Pathway%WikiPathways\_20180810%WP127%Homo sapiens                                                                 |  0.0075613|  0.5352106|  1.891646|             9|
| Kit receptor signaling pathway%WikiPathways\_20180810%WP304%Homo sapiens                                                         |  0.0081057|  0.4389567|  1.709651|            10|
| Nanoparticle-mediated activation of receptor signaling%WikiPathways\_20180810%WP2643%Homo sapiens                                |  0.0081063|  0.5658621|  1.849419|            10|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                                       |  0.0081063|  0.5589618|  1.841983|            10|
| PDGFR-beta pathway%WikiPathways\_20180810%WP3972%Homo sapiens                                                                    |  0.0094250|  0.5466798|  1.820490|            12|
| Chemokine signaling pathway%WikiPathways\_20180810%WP3929%Homo sapiens                                                           |  0.0100646|  0.3648770|  1.605027|            14|
| Myometrial Relaxation and Contraction Pathways%WikiPathways\_20180810%WP289%Homo sapiens                                         |  0.0100646|  0.3563647|  1.588788|            14|
| Lipid Metabolism Pathway%WikiPathways\_20180810%WP3965%Homo sapiens                                                              |  0.0111216|  0.5683300|  1.831119|            15|
| Hedgehog Signaling Pathway%WikiPathways\_20180810%WP4249%Homo sapiens                                                            |  0.0116001|  0.5100525|  1.802727|            16|
| TNF related weak inducer of apoptosis (TWEAK) Signaling Pathway%WikiPathways\_20180810%WP2036%Homo sapiens                       |  0.0127821|  0.4906609|  1.770717|            18|
| Notch Signaling Pathway%WikiPathways\_20180810%WP268%Homo sapiens                                                                |  0.0132772|  0.4822076|  1.751845|            19|
| Leptin Insulin Overlap%WikiPathways\_20180810%WP3935%Homo sapiens                                                                |  0.0137320|  0.6531517|  1.827734|            19|
| Spinal Cord Injury%WikiPathways\_20180810%WP2431%Homo sapiens                                                                    |  0.0145744|  0.3801514|  1.603360|            22|
| Photodynamic therapy-induced AP-1 survival signaling.%WikiPathways\_20180810%WP3611%Homo sapiens                                 |  0.0172769|  0.4655896|  1.725630|            26|
| Signaling of Hepatocyte Growth Factor Receptor%WikiPathways\_20180810%WP313%Homo sapiens                                         |  0.0172880|  0.5027303|  1.739690|            26|
| Viral Acute Myocarditis%WikiPathways\_20180810%WP4298%Homo sapiens                                                               |  0.0173323|  0.4020623|  1.620232|            27|
| Aryl Hydrocarbon Receptor Pathway%WikiPathways\_20180810%WP2873%Homo sapiens                                                     |  0.0175933|  0.5236745|  1.743880|            27|
| IL-1 signaling pathway%WikiPathways\_20180810%WP195%Homo sapiens                                                                 |  0.0182079|  0.4333018|  1.657049|            29|
| Aryl Hydrocarbon Receptor%WikiPathways\_20180810%WP2586%Homo sapiens                                                             |  0.0182079|  0.4709270|  1.699500|            29|
| Brain-Derived Neurotrophic Factor (BDNF) signaling pathway%WikiPathways\_20180810%WP2380%Homo sapiens                            |  0.0200708|  0.3383520|  1.516881|            34|
| Extracellular vesicle-mediated signaling in recipient cells%WikiPathways\_20180810%WP2870%Homo sapiens                           |  0.0240638|  0.5183525|  1.726157|            39|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                                              |  0.0257172|  0.3392901|  1.507479|            45|
| Thymic Stromal LymphoPoietin (TSLP) Signaling Pathway%WikiPathways\_20180810%WP2203%Homo sapiens                                 |  0.0270451|  0.4544407|  1.657251|            46|
| BMP Signaling Pathway in Eyelid Development%WikiPathways\_20180810%WP3927%Homo sapiens                                           |  0.0270451|  0.5936611|  1.747039|            44|
| miRs in Muscle Cell Differentiation%WikiPathways\_20180810%WP2012%Homo sapiens                                                   |  0.0271985|  0.5339054|  1.696730|            46|
| ErbB Signaling Pathway%WikiPathways\_20180810%WP673%Homo sapiens                                                                 |  0.0347143|  0.4194829|  1.594291|            62|
| AGE/RAGE pathway%WikiPathways\_20180810%WP2324%Homo sapiens                                                                      |  0.0378829|  0.3985061|  1.560234|            68|
| Simplified Interaction Map Between LOXL4 and Oxidative Stress Pathway%WikiPathways\_20180810%WP3670%Homo sapiens                 |  0.0379475|  0.5737966|  1.688581|            66|
| Serotonin Receptor 4/6/7 and NR3C Signaling%WikiPathways\_20180810%WP734%Homo sapiens                                            |  0.0387914|  0.5812404|  1.683930|            68|
| The effect of progerin on the involved genes in Hutchinson-Gilford Progeria Syndrome%WikiPathways\_20180810%WP4320%Homo sapiens  |  0.0394226|  0.4691903|  1.623625|            72|
| MAPK Cascade%WikiPathways\_20180810%WP422%Homo sapiens                                                                           |  0.0400117|  0.4896546|  1.630591|            74|
| White fat cell differentiation%WikiPathways\_20180810%WP4149%Homo sapiens                                                        |  0.0422985|  0.5005813|  1.636061|            80|
| Regulation of Wnt/B-catenin Signaling by Small Molecule Compounds%WikiPathways\_20180810%WP3664%Homo sapiens                     |  0.0422985|  0.5969873|  1.670567|            78|
| Microglia Pathogen Phagocytosis Pathway%WikiPathways\_20180810%WP3937%Homo sapiens                                               |  0.0422985|  0.5247597|  1.630465|            79|
| Calcium Regulation in the Cardiac Cell%WikiPathways\_20180810%WP536%Homo sapiens                                                 |  0.0445308|  0.3269138|  1.441622|            92|
| Preimplantation Embryo%WikiPathways\_20180810%WP3527%Homo sapiens                                                                |  0.0445308|  0.4265397|  1.561729|            87|
| TNF alpha Signaling Pathway%WikiPathways\_20180810%WP231%Homo sapiens                                                            |  0.0450534|  0.3570467|  1.494458|            91|
| GABA receptor Signaling%WikiPathways\_20180810%WP4159%Homo sapiens                                                               |  0.0451687|  0.5056628|  1.606976|            89|
| Inflammatory Response Pathway%WikiPathways\_20180810%WP453%Homo sapiens                                                          |  0.0452128|  0.5292745|  1.625191|            89|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                                            |  0.0452927|  0.3253178|  1.434584|            97|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                                        |  0.0473810|  0.2794731|  1.348367|           106|
| Fas Ligand (FasL) pathway and Stress induction of Heat Shock Proteins (HSP) regulation%WikiPathways\_20180810%WP314%Homo sapiens |  0.0491561|  0.4312218|  1.556210|           102|

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

| Pathway                                        | Corrected p-value | Corrected MF p-value |
|:-----------------------------------------------|:------------------|:---------------------|
| Extracellular matrix organization              | 0.0000000         | 0.0000000            |
| Class B 2 (Secretin family receptors)          | 0.0000000         | 0.0000000            |
| Neuronal System                                | 0.0000000         | 0.0206333            |
| G alpha (i) signalling events                  | 0.0000000         | 0.0309500            |
| Class A 1 (Rhodopsin-like receptors)           | 0.0619000         | 0.0495200            |
| Collagen formation                             | 0.1061143         | 0.1238000            |
| Platelet activation, signaling and aggregation | 0.0515833         | 0.1341167            |
| Collagen chain trimerization                   | 0.3156900         | 0.3156900            |
| ECM proteoglycans                              | 0.3172375         | 0.3249750            |
| Collagen biosynthesis and modifying enzymes    | 0.3207546         | 0.3320091            |

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

There are 25 up-regulated and 61 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                              |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                 |  0.0055795|  0.6411296|  2.020710|             0|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                |  0.0055795|  0.4229049|  1.742119|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                  |  0.0055795|  0.6223854|  2.090884|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                            |  0.0055795|  0.5076242|  2.116065|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                          |  0.0055795|  0.5986648|  2.094898|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                            |  0.0055795|  0.6009169|  1.940825|             0|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                      |  0.0055795|  0.6215197|  1.958904|             0|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                           |  0.0055795|  0.4408896|  1.799602|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090 |  0.0055795|  0.6139364|  1.998641|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                          |  0.0055795|  0.5305453|  1.948622|             0|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                  |  0.0067169|  0.4773152|  1.778265|             1|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                   |  0.0067169|  0.5061851|  1.832749|             1|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                        |  0.0067169|  0.4983576|  1.841735|             1|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                              |  0.0103774|  0.7521806|  1.922390|             2|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                         |  0.0137535|  0.4328159|  1.675652|             4|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                             |  0.0144072|  0.4911359|  1.732939|             4|
| COMPLEMENT CASCADE%REACTOME DATABASE ID RELEASE 65%166658                                            |  0.0190809|  0.6803816|  1.885419|             6|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                           |  0.0190809|  0.5639546|  1.828236|             6|
| DISEASES OF GLYCOSYLATION%REACTOME DATABASE ID RELEASE 65%3781865                                    |  0.0202920|  0.4987830|  1.720353|             7|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                    |  0.0233417|  0.6871840|  1.858908|             8|
| G ALPHA (Q) SIGNALLING EVENTS%REACTOME DATABASE ID RELEASE 65%416476                                 |  0.0334706|  0.4091554|  1.587097|            15|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                        |  0.0415905|  0.6020824|  1.778269|            17|
| CLASS B 2 (SECRETIN FAMILY RECEPTORS)%REACTOME DATABASE ID RELEASE 65%373080                         |  0.0440231|  0.5049330|  1.674479|            19|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                           |  0.0472782|  0.5007993|  1.666460|            21|
| TRAF6 MEDIATED INDUCTION OF NFKB AND MAP KINASES UPON TLR7 8 OR 9 ACTIVATION%REACTOME%R-HSA-975138.1 |  0.0479776|  0.4372266|  1.583070|            23|

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

| Pathway                                        | Corrected p-value | Corrected MF p-value |
|:-----------------------------------------------|:------------------|:---------------------|
| Kegg neuroactive ligand receptor interaction   | 0.0000000         | 0.0000000            |
| Kegg ecm receptor interaction                  | 0.0000000         | 0.0000000            |
| Kegg cell adhesion molecules cams              | 0.0000000         | 0.0000000            |
| Kegg cytokine cytokine receptor interaction    | 0.0216000         | 0.0364500            |
| Kegg calcium signaling pathway                 | 0.0162000         | 0.0388800            |
| Kegg focal adhesion                            | 0.0040500         | 0.0648000            |
| Kegg hematopoietic cell lineage                | 0.1512000         | 0.1984500            |
| Kegg taste transduction                        | 0.1393200         | 0.2013429            |
| Kegg systemic lupus erythematosus              | 0.2054842         | 0.2533091            |
| Kegg aldosterone regulated sodium reabsorption | 0.1825200         | 0.2557286            |

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

There are 4 up-regulated and 3 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                      |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_FOCAL\_ADHESION                        |  0.0132553|  0.4368927|  1.764752|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION             |  0.0132553|  0.5627785|  1.962808|             0|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES |  0.0223574|  0.6288321|  1.928557|             2|
| KEGG\_HEMATOPOIETIC\_CELL\_LINEAGE           |  0.0296202|  0.5695954|  1.824025|             4|

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

There are 11 up-regulated and 2 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| Hematopoietic Stem Cell Differentiation%WikiPathways\_20180810%WP2849%Homo sapiens                                        |  0.0110842|  0.6673470|  2.093985|             0|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                                       |  0.0110842|  0.4502829|  1.731594|             0|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                              |  0.0110842|  0.4983099|  1.833308|             0|
| Endothelin Pathways%WikiPathways\_20180810%WP2197%Homo sapiens                                                            |  0.0110842|  0.6915424|  1.957409|             0|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                                                     |  0.0110842|  0.3826395|  1.612594|             0|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                           |  0.0143534|  0.4600371|  1.706477|             1|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                              |  0.0201803|  0.5958494|  1.896321|             2|
| Complement and Coagulation Cascades%WikiPathways\_20180810%WP558%Homo sapiens                                             |  0.0276151|  0.6264701|  1.894223|             4|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                                 |  0.0276151|  0.3673466|  1.541883|             5|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                                                       |  0.0288667|  0.5285724|  1.782180|             5|
| Factors and pathways affecting insulin-like growth factor (IGF1)-Akt signaling%WikiPathways\_20180810%WP3850%Homo sapiens |  0.0467723|  0.6186584|  1.826011|            10|

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
| Collagen formation                                   | 0.0000000         | 0.0000000            |
| Collagen chain trimerization                         | 0.0000000         | 0.0000000            |
| Collagen biosynthesis and modifying enzymes          | 0.0000000         | 0.0000000            |
| ECM proteoglycans                                    | 0.0123800         | 0.0123800            |
| G alpha (i) signalling events                        | 0.0515833         | 0.1238000            |
| Signaling by NODAL                                   | 0.2033857         | 0.2033857            |
| Binding and Uptake of Ligands by Scavenger Receptors | 0.2166500         | 0.2243875            |
| Signaling by BMP                                     | 0.2407222         | 0.2407222            |
| Class A 1 (Rhodopsin-like receptors)                 | 0.2537900         | 0.2723600            |

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

There are 128 up-regulated and 117 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| LAMININ INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000157                                                                              |  0.0046751|  0.7271213|  2.196060|             0|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                                                      |  0.0046751|  0.7634609|  2.658627|             0|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                                                              |  0.0046751|  0.5364065|  1.944253|             0|
| REGULATION OF FZD BY UBIQUITINATION%REACTOME%R-HSA-4641263.2                                                                              |  0.0046751|  0.6978179|  2.054636|             0|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                                                         |  0.0046751|  0.6910985|  2.034852|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                                                       |  0.0046751|  0.7369268|  2.756805|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                                                 |  0.0046751|  0.5727275|  2.702206|             0|
| GLYCOSAMINOGLYCAN METABOLISM%REACTOME DATABASE ID RELEASE 65%1630316                                                                      |  0.0046751|  0.4614803|  1.962052|             0|
| SIGNALING BY NOTCH3%REACTOME DATABASE ID RELEASE 65%9012852                                                                               |  0.0046751|  0.5394003|  1.947690|             0|
| EPH-EPHRIN MEDIATED REPULSION OF CELLS%REACTOME%R-HSA-3928665.3                                                                           |  0.0046751|  0.5303114|  1.929447|             0|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                                                                |  0.0046751|  0.6000930|  1.978457|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                                                               |  0.0046751|  0.6796685|  2.663266|             0|
| EPH-EPHRIN SIGNALING%REACTOME%R-HSA-2682334.1                                                                                             |  0.0046751|  0.4717783|  1.942702|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                                                 |  0.0046751|  0.6656855|  2.383990|             0|
| SIGNALING BY PDGF%REACTOME DATABASE ID RELEASE 65%186797                                                                                  |  0.0046751|  0.5737726|  2.098727|             0|
| HS-GAG BIOSYNTHESIS%REACTOME%R-HSA-2022928.1                                                                                              |  0.0046751|  0.6647125|  2.170248|             0|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                                                        |  0.0046751|  0.5364088|  2.181969|             0|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                                                                |  0.0046751|  0.5966669|  2.146060|             0|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                                           |  0.0046751|  0.7033912|  2.449444|             0|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                                                |  0.0046751|  0.6060797|  2.247937|             0|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                                                         |  0.0046751|  0.6600929|  2.118562|             0|
| RHO GTPASE CYCLE%REACTOME DATABASE ID RELEASE 65%194840                                                                                   |  0.0046751|  0.5021043|  2.194329|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090                                      |  0.0046751|  0.7323880|  2.644539|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                                                               |  0.0046751|  0.5596483|  2.308530|             0|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                                                     |  0.0058748|  0.3726225|  1.732740|             1|
| SIGNALING BY NTRKS%REACTOME DATABASE ID RELEASE 65%166520                                                                                 |  0.0058748|  0.4218281|  1.762156|             1|
| HEPARAN SULFATE HEPARIN (HS-GAG) METABOLISM%REACTOME%R-HSA-1638091.1                                                                      |  0.0058748|  0.5241214|  1.928828|             1|
| SIGNALING BY VEGF%REACTOME DATABASE ID RELEASE 65%194138                                                                                  |  0.0058748|  0.4638137|  1.929282|             1|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                                                             |  0.0058748|  0.6276843|  1.998422|             1|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                                                             |  0.0058748|  0.4303519|  1.785298|             1|
| SIGNALING BY NTRK1 (TRKA)%REACTOME%R-HSA-187037.2                                                                                         |  0.0058748|  0.4692721|  1.871403|             1|
| DISEASES OF GLYCOSYLATION%REACTOME DATABASE ID RELEASE 65%3781865                                                                         |  0.0058748|  0.4907506|  1.893578|             1|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                                                                   |  0.0058748|  0.7379233|  2.032795|             1|
| CELL-EXTRACELLULAR MATRIX INTERACTIONS%REACTOME%R-HSA-446353.1                                                                            |  0.0058748|  0.7135915|  1.965767|             1|
| NRAGE SIGNALS DEATH THROUGH JNK%REACTOME DATABASE ID RELEASE 65%193648                                                                    |  0.0067913|  0.5086439|  1.886550|             2|
| N-GLYCAN TRIMMING IN THE ER AND CALNEXIN CALRETICULIN CYCLE%REACTOME%R-HSA-532668.2                                                       |  0.0067913|  0.5501566|  1.880378|             2|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                                                           |  0.0067913|  0.4428251|  1.792606|             2|
| MET PROMOTES CELL MOTILITY%REACTOME DATABASE ID RELEASE 65%8875878                                                                        |  0.0067913|  0.6098729|  1.957381|             2|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                                                              |  0.0067913|  0.4235370|  1.778252|             2|
| VEGFA-VEGFR2 PATHWAY%REACTOME%R-HSA-4420097.3                                                                                             |  0.0078036|  0.4313524|  1.765010|             3|
| CHONDROITIN SULFATE DERMATAN SULFATE METABOLISM%REACTOME DATABASE ID RELEASE 65%1793185                                                   |  0.0078036|  0.5070917|  1.823881|             3|
| SIGNALING BY EGFR%REACTOME DATABASE ID RELEASE 65%177929                                                                                  |  0.0078036|  0.5160705|  1.838951|             3|
| O-GLYCOSYLATION OF TSR DOMAIN-CONTAINING PROTEINS%REACTOME%R-HSA-5173214.1                                                                |  0.0078036|  0.5821632|  1.935568|             3|
| DEFECTIVE B3GALTL CAUSES PETERS-PLUS SYNDROME (PPS)%REACTOME%R-HSA-5083635.1                                                              |  0.0078036|  0.5731837|  1.889739|             3|
| NON-INTEGRIN MEMBRANE-ECM INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000171                                                            |  0.0078036|  0.5349328|  1.862815|             3|
| L1CAM INTERACTIONS%REACTOME%R-HSA-373760.2                                                                                                |  0.0089819|  0.4252722|  1.746959|             4|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                                                             |  0.0089942|  0.5940740|  1.925793|             4|
| RUNX2 REGULATES BONE DEVELOPMENT%REACTOME%R-HSA-8941326.1                                                                                 |  0.0089942|  0.5949921|  1.909622|             4|
| TP53 REGULATES TRANSCRIPTION OF CELL CYCLE GENES%REACTOME DATABASE ID RELEASE 65%6791312                                                  |  0.0089942|  0.4974348|  1.802997|             4|
| G ALPHA (12 13) SIGNALLING EVENTS%REACTOME%R-HSA-416482.3                                                                                 |  0.0099241|  0.4359478|  1.744848|             5|
| TOLL LIKE RECEPTOR TLR6:TLR2 CASCADE%REACTOME%R-HSA-168188.1                                                                              |  0.0107725|  0.4101515|  1.688933|             6|
| MYD88:MAL CASCADE INITIATED ON PLASMA MEMBRANE%REACTOME%R-HSA-166058.2                                                                    |  0.0107725|  0.4101515|  1.688933|             6|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                                                       |  0.0107725|  0.3986534|  1.669768|             6|
| TOLL LIKE RECEPTOR 2 (TLR2) CASCADE%REACTOME DATABASE ID RELEASE 65%181438                                                                |  0.0107725|  0.4101515|  1.688933|             6|
| TOLL LIKE RECEPTOR TLR1:TLR2 CASCADE%REACTOME DATABASE ID RELEASE 65%168179                                                               |  0.0107725|  0.4101515|  1.688933|             6|
| CLASS B 2 (SECRETIN FAMILY RECEPTORS)%REACTOME DATABASE ID RELEASE 65%373080                                                              |  0.0108042|  0.4816090|  1.780844|             6|
| DISEASES ASSOCIATED WITH O-GLYCOSYLATION OF PROTEINS%REACTOME%R-HSA-3906995.2                                                             |  0.0108916|  0.5061652|  1.779930|             6|
| INTRA-GOLGI AND RETROGRADE GOLGI-TO-ER TRAFFIC%REACTOME DATABASE ID RELEASE 65%6811442                                                    |  0.0116958|  0.3474726|  1.584855|             7|
| TRAF6 MEDIATED INDUCTION OF NFKB AND MAP KINASES UPON TLR7 8 OR 9 ACTIVATION%REACTOME%R-HSA-975138.1                                      |  0.0118122|  0.4133133|  1.681249|             7|
| SIGNALING BY WNT%REACTOME DATABASE ID RELEASE 65%195721                                                                                   |  0.0172980|  0.3056332|  1.480203|            12|
| SIGNALLING TO RAS%REACTOME DATABASE ID RELEASE 65%167044                                                                                  |  0.0173042|  0.6459158|  1.849818|            11|
| CLATHRIN-MEDIATED ENDOCYTOSIS%REACTOME%R-HSA-8856828.3                                                                                    |  0.0173042|  0.3596602|  1.571811|            12|
| CALNEXIN CALRETICULIN CYCLE%REACTOME%R-HSA-901042.2                                                                                       |  0.0173042|  0.5631319|  1.792901|            11|
| MET ACTIVATES PTK2 SIGNALING%REACTOME DATABASE ID RELEASE 65%8874081                                                                      |  0.0173042|  0.6485362|  1.857323|            11|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                                                                  |  0.0182930|  0.5650694|  1.813585|            12|
| TOLL LIKE RECEPTOR 7 8 (TLR7 8) CASCADE%REACTOME DATABASE ID RELEASE 65%168181                                                            |  0.0191960|  0.4047850|  1.649814|            14|
| MYD88 DEPENDENT CASCADE INITIATED ON ENDOSOME%REACTOME%R-HSA-975155.1                                                                     |  0.0191960|  0.4047850|  1.649814|            14|
| O-LINKED GLYCOSYLATION%REACTOME%R-HSA-5173105.3                                                                                           |  0.0200344|  0.4037059|  1.637278|            15|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                                                             |  0.0200344|  0.4638759|  1.696751|            15|
| TOLL LIKE RECEPTOR 9 (TLR9) CASCADE%REACTOME DATABASE ID RELEASE 65%168138                                                                |  0.0200344|  0.4011339|  1.645305|            15|
| ASPARAGINE N-LINKED GLYCOSYLATION%REACTOME%R-HSA-446203.4                                                                                 |  0.0200344|  0.3045354|  1.464234|            16|
| MYD88 CASCADE INITIATED ON PLASMA MEMBRANE%REACTOME%R-HSA-975871.1                                                                        |  0.0204726|  0.4014290|  1.628043|            16|
| TOLL LIKE RECEPTOR 10 (TLR10) CASCADE%REACTOME DATABASE ID RELEASE 65%168142                                                              |  0.0204726|  0.4014290|  1.628043|            16|
| TOLL LIKE RECEPTOR 5 (TLR5) CASCADE%REACTOME%R-HSA-168176.1                                                                               |  0.0204726|  0.4014290|  1.628043|            16|
| ACTIVATION OF MATRIX METALLOPROTEINASES%REACTOME%R-HSA-1592389.1                                                                          |  0.0204726|  0.6350531|  1.818709|            15|
| SIGNALING BY WNT IN CANCER%REACTOME DATABASE ID RELEASE 65%4791275                                                                        |  0.0210193|  0.5181554|  1.750161|            16|
| INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING%REACTOME DATABASE ID RELEASE 65%6785807                                                        |  0.0212632|  0.3982062|  1.633297|            17|
| NOTCH3 ACTIVATION AND TRANSMISSION OF SIGNAL TO THE NUCLEUS%REACTOME DATABASE ID RELEASE 65%9013507                                       |  0.0230685|  0.5749506|  1.794496|            18|
| TRIF(TICAM1)-MEDIATED TLR4 SIGNALING%REACTOME%R-HSA-937061.2                                                                              |  0.0235970|  0.3917665|  1.616023|            20|
| MYD88-INDEPENDENT TLR4 CASCADE%REACTOME%R-HSA-166166.2                                                                                    |  0.0235970|  0.3917665|  1.616023|            20|
| MAP KINASE ACTIVATION%REACTOME%R-HSA-450294.3                                                                                             |  0.0240976|  0.4285091|  1.653417|            20|
| SMOOTH MUSCLE CONTRACTION%REACTOME%R-HSA-445355.2                                                                                         |  0.0243711|  0.5225519|  1.737373|            20|
| TOLL LIKE RECEPTOR 3 (TLR3) CASCADE%REACTOME DATABASE ID RELEASE 65%168164                                                                |  0.0251328|  0.3905721|  1.608308|            22|
| MAPK TARGETS NUCLEAR EVENTS MEDIATED BY MAP KINASES%REACTOME DATABASE ID RELEASE 65%450282                                                |  0.0251328|  0.5219748|  1.735454|            21|
| CASPASE ACTIVATION VIA EXTRINSIC APOPTOTIC SIGNALLING PATHWAY%REACTOME%R-HSA-5357769.2                                                    |  0.0259692|  0.6251075|  1.790226|            22|
| PEPTIDE HORMONE METABOLISM%REACTOME DATABASE ID RELEASE 65%2980736                                                                        |  0.0259692|  0.4425518|  1.675515|            23|
| DEATH RECEPTOR SIGNALLING%REACTOME%R-HSA-73887.3                                                                                          |  0.0260123|  0.3487413|  1.521510|            24|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                                                                    |  0.0261562|  0.3752864|  1.573642|            24|
| SIGNALLING TO ERKS%REACTOME%R-HSA-187687.1                                                                                                |  0.0263562|  0.5403635|  1.751681|            23|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                                                |  0.0271854|  0.3126855|  1.440568|            26|
| INTRINSIC PATHWAY FOR APOPTOSIS%REACTOME%R-HSA-109606.2                                                                                   |  0.0290889|  0.4682451|  1.658333|            27|
| SIGNALING BY RAS MUTANTS%REACTOME%R-HSA-6802949.1                                                                                         |  0.0297692|  0.4587096|  1.656329|            28|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                                                            |  0.0321754|  0.6271626|  1.767251|            30|
| REGULATION OF PTEN GENE TRANSCRIPTION%REACTOME DATABASE ID RELEASE 65%8943724                                                             |  0.0321754|  0.4355678|  1.637082|            31|
| TGF-BETA RECEPTOR SIGNALING ACTIVATES SMADS%REACTOME DATABASE ID RELEASE 65%2173789                                                       |  0.0327911|  0.5140958|  1.709259|            31|
| TOLL LIKE RECEPTOR 4 (TLR4) CASCADE%REACTOME%R-HSA-166016.2                                                                               |  0.0327982|  0.3554246|  1.522914|            33|
| G-PROTEIN BETA:GAMMA SIGNALLING%REACTOME DATABASE ID RELEASE 65%397795                                                                    |  0.0336006|  0.5633308|  1.740422|            33|
| DOWNSTREAM SIGNAL TRANSDUCTION%REACTOME DATABASE ID RELEASE 65%186763                                                                     |  0.0336006|  0.5307263|  1.703361|            33|
| NOTCH4 INTRACELLULAR DOMAIN REGULATES TRANSCRIPTION%REACTOME DATABASE ID RELEASE 65%9013695                                               |  0.0336006|  0.5904751|  1.738579|            33|
| FCERI MEDIATED MAPK ACTIVATION%REACTOME%R-HSA-2871796.2                                                                                   |  0.0336006|  0.5523948|  1.724096|            33|
| INTERLEUKIN-17 SIGNALING%REACTOME DATABASE ID RELEASE 65%448424                                                                           |  0.0336006|  0.4152766|  1.618355|            34|
| TP53 REGULATES TRANSCRIPTION OF ADDITIONAL CELL CYCLE GENES WHOSE EXACT ROLE IN THE P53 PATHWAY REMAIN UNCERTAIN%REACTOME%R-HSA-6804115.1 |  0.0352214|  0.5888055|  1.733663|            35|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                                                           |  0.0363674|  0.4176251|  1.598121|            39|
| SIGNALING BY NOTCH1 PEST DOMAIN MUTANTS IN CANCER%REACTOME DATABASE ID RELEASE 65%2644602                                                 |  0.0363674|  0.4275683|  1.612304|            40|
| CONSTITUTIVE SIGNALING BY NOTCH1 HD+PEST DOMAIN MUTANTS%REACTOME DATABASE ID RELEASE 65%2894862                                           |  0.0363674|  0.4275683|  1.612304|            40|
| COPI-MEDIATED ANTEROGRADE TRANSPORT%REACTOME%R-HSA-6807878.1                                                                              |  0.0363674|  0.3899728|  1.569653|            41|
| SIGNALING BY NOTCH1 IN CANCER%REACTOME DATABASE ID RELEASE 65%2644603                                                                     |  0.0363674|  0.4275683|  1.612304|            40|
| SIGNALING BY NOTCH1 HD+PEST DOMAIN MUTANTS IN CANCER%REACTOME DATABASE ID RELEASE 65%2894858                                              |  0.0363674|  0.4275683|  1.612304|            40|
| CONSTITUTIVE SIGNALING BY NOTCH1 PEST DOMAIN MUTANTS%REACTOME%R-HSA-2644606.1                                                             |  0.0363674|  0.4275683|  1.612304|            40|
| NEURODEGENERATIVE DISEASES%REACTOME DATABASE ID RELEASE 65%8863678                                                                        |  0.0363674|  0.5811884|  1.732534|            39|
| DEREGULATED CDK5 TRIGGERS MULTIPLE NEURODEGENERATIVE PATHWAYS IN ALZHEIMER'S DISEASE MODELS%REACTOME%R-HSA-8862803.2                      |  0.0363674|  0.5811884|  1.732534|            39|
| CELL DEATH SIGNALLING VIA NRAGE, NRIF AND NADE%REACTOME DATABASE ID RELEASE 65%204998                                                     |  0.0375146|  0.4017010|  1.578812|            42|
| NOTCH3 INTRACELLULAR DOMAIN REGULATES TRANSCRIPTION%REACTOME%R-HSA-9013508.1                                                              |  0.0397003|  0.5615138|  1.714931|            44|
| PRE-NOTCH PROCESSING IN GOLGI%REACTOME%R-HSA-1912420.2                                                                                    |  0.0397003|  0.6137732|  1.729522|            44|
| SYNTHESIS OF PIPS AT THE PLASMA MEMBRANE%REACTOME DATABASE ID RELEASE 65%1660499                                                          |  0.0411732|  0.4358078|  1.603824|            48|
| P75 NTR RECEPTOR-MEDIATED SIGNALLING%REACTOME%R-HSA-193704.1                                                                              |  0.0421269|  0.3745126|  1.538446|            50|
| RUNX2 REGULATES OSTEOBLAST DIFFERENTIATION%REACTOME DATABASE ID RELEASE 65%8940973                                                        |  0.0421950|  0.5738155|  1.710555|            48|
| GAMMA CARBOXYLATION, HYPUSINE FORMATION AND ARYLSULFATASE ACTIVATION%REACTOME DATABASE ID RELEASE 65%163841                               |  0.0435337|  0.5140355|  1.666334|            50|
| G ALPHA (Q) SIGNALLING EVENTS%REACTOME DATABASE ID RELEASE 65%416476                                                                      |  0.0437914|  0.3331424|  1.458623|            55|
| ACTIVATION OF BH3-ONLY PROTEINS%REACTOME DATABASE ID RELEASE 65%114452                                                                    |  0.0440030|  0.5133241|  1.664028|            51|
| EPHRIN SIGNALING%REACTOME DATABASE ID RELEASE 65%3928664                                                                                  |  0.0451743|  0.5976574|  1.711612|            53|
| INTEGRIN SIGNALING%REACTOME DATABASE ID RELEASE 65%9006921                                                                                |  0.0456785|  0.5512139|  1.683474|            56|
| INTEGRIN ALPHAIIB BETA3 SIGNALING%REACTOME DATABASE ID RELEASE 65%354192                                                                  |  0.0456785|  0.5512139|  1.683474|            56|
| SIGNALING BY LIGAND-RESPONSIVE EGFR VARIANTS IN CANCER%REACTOME DATABASE ID RELEASE 65%5637815                                            |  0.0456785|  0.5804440|  1.689122|            55|
| CONSTITUTIVE SIGNALING BY LIGAND-RESPONSIVE EGFR CANCER VARIANTS%REACTOME DATABASE ID RELEASE 65%1236382                                  |  0.0456785|  0.5804440|  1.689122|            55|
| SIGNALING BY EGFR IN CANCER%REACTOME%R-HSA-1643713.1                                                                                      |  0.0456785|  0.5804440|  1.689122|            55|
| SIGNALING BY PTK6%REACTOME%R-HSA-8848021.2                                                                                                |  0.0464445|  0.4252362|  1.564919|            60|
| SIGNALING BY NON-RECEPTOR TYROSINE KINASES%REACTOME%R-HSA-9006927.2                                                                       |  0.0464445|  0.4252362|  1.564919|            60|

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
| Kegg pathways in cancer                     | 0.0054000         | 0.0000000            |
| Kegg focal adhesion                         | 0.0000000         | 0.0000000            |
| Kegg ecm receptor interaction               | 0.0000000         | 0.0000000            |
| Kegg wnt signaling pathway                  | 0.0162000         | 0.0081000            |
| Kegg mapk signaling pathway                 | 0.0439714         | 0.0129600            |
| Kegg axon guidance                          | 0.0162000         | 0.0135000            |
| Kegg cell adhesion molecules cams           | 0.0378000         | 0.0254571            |
| Kegg erbb signaling pathway                 | 0.2106000         | 0.0810000            |
| Kegg hypertrophic cardiomyopathy hcm        | 0.2400546         | 0.0824727            |
| Kegg cytokine cytokine receptor interaction | 0.0951750         | 0.0826200            |

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

There are 47 up-regulated and 15 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                                        |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                                 |  0.0019181|  0.3649497|  1.727791|             0|
| KEGG\_ENDOCYTOSIS                                              |  0.0019181|  0.3778583|  1.719172|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY                                  |  0.0019181|  0.4658496|  2.036442|             0|
| KEGG\_AXON\_GUIDANCE                                           |  0.0019181|  0.4539690|  1.972390|             0|
| KEGG\_FOCAL\_ADHESION                                          |  0.0019181|  0.5517971|  2.521566|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                               |  0.0019181|  0.6864207|  2.691792|             0|
| KEGG\_NEUROTROPHIN\_SIGNALING\_PATHWAY                         |  0.0019181|  0.5058686|  2.178317|             0|
| KEGG\_MELANOGENESIS                                            |  0.0019181|  0.4873528|  1.978194|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                                     |  0.0019181|  0.4080179|  1.980857|             0|
| KEGG\_COLORECTAL\_CANCER                                       |  0.0019181|  0.5189777|  1.998612|             0|
| KEGG\_RENAL\_CELL\_CARCINOMA                                   |  0.0019181|  0.4768893|  1.875558|             0|
| KEGG\_PANCREATIC\_CANCER                                       |  0.0019181|  0.5206546|  2.054588|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA                                   |  0.0019181|  0.6165568|  2.211266|             0|
| KEGG\_CHRONIC\_MYELOID\_LEUKEMIA                               |  0.0019181|  0.5024482|  1.998625|             0|
| KEGG\_SMALL\_CELL\_LUNG\_CANCER                                |  0.0019181|  0.5145797|  2.098850|             0|
| KEGG\_HYPERTROPHIC\_CARDIOMYOPATHY\_HCM                        |  0.0019181|  0.4956741|  1.922429|             0|
| KEGG\_ERBB\_SIGNALING\_PATHWAY                                 |  0.0027732|  0.4471755|  1.819063|             1|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY                             |  0.0027732|  0.5382415|  1.921788|             1|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY                            |  0.0027732|  0.4493235|  1.804956|             1|
| KEGG\_VASOPRESSIN\_REGULATED\_WATER\_REABSORPTION              |  0.0027732|  0.5789300|  2.021125|             1|
| KEGG\_DILATED\_CARDIOMYOPATHY                                  |  0.0027732|  0.4723732|  1.857797|             1|
| KEGG\_REGULATION\_OF\_ACTIN\_CYTOSKELETON                      |  0.0047378|  0.3536441|  1.617179|             3|
| KEGG\_INSULIN\_SIGNALING\_PATHWAY                              |  0.0047378|  0.3818684|  1.659129|             3|
| KEGG\_MTOR\_SIGNALING\_PATHWAY                                 |  0.0056728|  0.4843856|  1.785630|             4|
| KEGG\_GNRH\_SIGNALING\_PATHWAY                                 |  0.0086009|  0.4203573|  1.706255|             7|
| KEGG\_ACUTE\_MYELOID\_LEUKEMIA                                 |  0.0095751|  0.4792117|  1.791018|             8|
| KEGG\_CHEMOKINE\_SIGNALING\_PATHWAY                            |  0.0100211|  0.3639969|  1.586514|             9|
| KEGG\_GLYCOSAMINOGLYCAN\_BIOSYNTHESIS\_HEPARAN\_SULFATE        |  0.0112671|  0.5890263|  1.833633|            10|
| KEGG\_LYSOSOME                                                 |  0.0131692|  0.3674954|  1.573355|            13|
| KEGG\_ARRHYTHMOGENIC\_RIGHT\_VENTRICULAR\_CARDIOMYOPATHY\_ARVC |  0.0138665|  0.4456180|  1.712825|            14|
| KEGG\_TOLL\_LIKE\_RECEPTOR\_SIGNALING\_PATHWAY                 |  0.0143462|  0.4177962|  1.643151|            15|
| KEGG\_PRION\_DISEASES                                          |  0.0145998|  0.5582828|  1.773116|            15|
| KEGG\_PROGESTERONE\_MEDIATED\_OOCYTE\_MATURATION               |  0.0176451|  0.3950158|  1.595530|            21|
| KEGG\_VEGF\_SIGNALING\_PATHWAY                                 |  0.0219677|  0.4112727|  1.617495|            27|
| KEGG\_APOPTOSIS                                                |  0.0234220|  0.3916074|  1.573108|            31|
| KEGG\_T\_CELL\_RECEPTOR\_SIGNALING\_PATHWAY                    |  0.0234220|  0.3747279|  1.543508|            31|
| KEGG\_FC\_EPSILON\_RI\_SIGNALING\_PATHWAY                      |  0.0234220|  0.4195961|  1.615888|            31|
| KEGG\_ENDOMETRIAL\_CANCER                                      |  0.0239471|  0.4409380|  1.637600|            32|
| KEGG\_NON\_SMALL\_CELL\_LUNG\_CANCER                           |  0.0248811|  0.4341424|  1.630390|            35|
| KEGG\_PHOSPHATIDYLINOSITOL\_SIGNALING\_SYSTEM                  |  0.0256545|  0.3982205|  1.577354|            40|
| KEGG\_PROSTATE\_CANCER                                         |  0.0256545|  0.3724379|  1.528541|            40|
| KEGG\_INOSITOL\_PHOSPHATE\_METABOLISM                          |  0.0257748|  0.4388364|  1.621468|            40|
| KEGG\_VASCULAR\_SMOOTH\_MUSCLE\_CONTRACTION                    |  0.0258641|  0.3676556|  1.520980|            43|
| KEGG\_HEMATOPOIETIC\_CELL\_LINEAGE                             |  0.0369885|  0.4481770|  1.600213|            64|
| KEGG\_GAP\_JUNCTION                                            |  0.0372911|  0.3715653|  1.497363|            68|
| KEGG\_ADIPOCYTOKINE\_SIGNALING\_PATHWAY                        |  0.0392809|  0.3998956|  1.533167|            74|
| KEGG\_BLADDER\_CANCER                                          |  0.0429040|  0.4578896|  1.598556|            82|

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
| miRNA targets in ECM and membrane receptors%WikiPathways\_20180810%WP2911%Homo sapiens                         | 0.0074250         | 0.00990              |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  | 0.0148500         | 0.01485              |
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       | 0.0099000         | 0.01485              |
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                | 0.0118800         | 0.02376              |
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens | 0.0111375         | 0.04950              |

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

There are 108 up-regulated and 7 down-regulated Wiki pathways with fgsea

``` r
upPathwaysWiki %>% kable()
```

| pathway                                                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:--------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| EGF/EGFR Signaling Pathway%WikiPathways\_20180810%WP437%Homo sapiens                                                      |  0.0018661|  0.3989701|  1.783473|             0|
| TGF-beta Receptor Signaling%WikiPathways\_20180810%WP560%Homo sapiens                                                     |  0.0018661|  0.6176381|  2.262146|             0|
| Focal Adhesion%WikiPathways\_20180810%WP306%Homo sapiens                                                                  |  0.0018661|  0.5097105|  2.306513|             0|
| Primary Focal Segmental Glomerulosclerosis FSGS%WikiPathways\_20180810%WP2572%Homo sapiens                                |  0.0018661|  0.6442048|  2.479248|             0|
| TGF-beta Signaling Pathway%WikiPathways\_20180810%WP366%Homo sapiens                                                      |  0.0018661|  0.4500663|  1.944088|             0|
| Corticotropin-releasing hormone signaling pathway%WikiPathways\_20180810%WP2355%Homo sapiens                              |  0.0018661|  0.5046611|  2.021149|             0|
| Neural Crest Differentiation%WikiPathways\_20180810%WP2064%Homo sapiens                                                   |  0.0018661|  0.5158416|  2.101185|             0|
| Hair Follicle Development: Cytodifferentiation (Part 3 of 3)%WikiPathways\_20180810%WP2840%Homo sapiens                   |  0.0018661|  0.5036266|  1.938227|             0|
| Myometrial Relaxation and Contraction Pathways%WikiPathways\_20180810%WP289%Homo sapiens                                  |  0.0018661|  0.4348730|  1.889511|             0|
| Adipogenesis%WikiPathways\_20180810%WP236%Homo sapiens                                                                    |  0.0018661|  0.5203722|  2.204301|             0|
| ESC Pluripotency Pathways%WikiPathways\_20180810%WP3931%Homo sapiens                                                      |  0.0018661|  0.4622360|  1.926485|             0|
| miRNA targets in ECM and membrane receptors%WikiPathways\_20180810%WP2911%Homo sapiens                                    |  0.0018661|  0.8041657|  2.414899|             0|
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                                       |  0.0018661|  0.4296359|  1.855838|             0|
| Integrin-mediated Cell Adhesion%WikiPathways\_20180810%WP185%Homo sapiens                                                 |  0.0018661|  0.4787230|  1.955385|             0|
| MAPK Signaling Pathway%WikiPathways\_20180810%WP382%Homo sapiens                                                          |  0.0018661|  0.3575754|  1.662302|             0|
| Prolactin Signaling Pathway%WikiPathways\_20180810%WP2037%Homo sapiens                                                    |  0.0018661|  0.4870824|  1.930656|             0|
| Heart Development%WikiPathways\_20180810%WP1591%Homo sapiens                                                              |  0.0018661|  0.7275293|  2.566589|             0|
| Senescence and Autophagy in Cancer%WikiPathways\_20180810%WP615%Homo sapiens                                              |  0.0018661|  0.5388780|  2.220257|             0|
| MicroRNAs in cardiomyocyte hypertrophy%WikiPathways\_20180810%WP1544%Homo sapiens                                         |  0.0018661|  0.4495802|  1.813283|             0|
| Association Between Physico-Chemical Features and Toxicity Associated Pathways%WikiPathways\_20180810%WP3680%Homo sapiens |  0.0018661|  0.5103733|  1.953910|             0|
| Endochondral Ossification%WikiPathways\_20180810%WP474%Homo sapiens                                                       |  0.0018661|  0.6271042|  2.363919|             0|
| RANKL/RANK (Receptor activator of NFKB (ligand)) Signaling Pathway%WikiPathways\_20180810%WP2018%Homo sapiens             |  0.0018661|  0.5374950|  2.005611|             0|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                           |  0.0018661|  0.5332985|  2.216669|             0|
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                             |  0.0018661|  0.6085853|  2.219840|             0|
| Human Thyroid Stimulating Hormone (TSH) signaling pathway%WikiPathways\_20180810%WP2032%Homo sapiens                      |  0.0018661|  0.5014604|  1.919788|             0|
| Focal Adhesion-PI3K-Akt-mTOR-signaling pathway%WikiPathways\_20180810%WP3932%Homo sapiens                                 |  0.0018661|  0.4246874|  2.015931|             0|
| VEGFA-VEGFR2 Signaling Pathway%WikiPathways\_20180810%WP3888%Homo sapiens                                                 |  0.0018661|  0.4759047|  2.220924|             0|
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens            |  0.0018661|  0.5700028|  2.286284|             0|
| PI3K-Akt Signaling Pathway%WikiPathways\_20180810%WP4172%Homo sapiens                                                     |  0.0018661|  0.4132502|  1.971084|             0|
| Epithelial to mesenchymal transition in colorectal cancer%WikiPathways\_20180810%WP4239%Homo sapiens                      |  0.0018661|  0.4829492|  2.120591|             0|
| Angiopoietin Like Protein 8 Regulatory Pathway%WikiPathways\_20180810%WP3915%Homo sapiens                                 |  0.0023942|  0.3909708|  1.686176|             1|
| Hematopoietic Stem Cell Differentiation%WikiPathways\_20180810%WP2849%Homo sapiens                                        |  0.0023942|  0.5777307|  1.998813|             1|
| Wnt/beta-catenin Signaling Pathway in Leukemia%WikiPathways\_20180810%WP3658%Homo sapiens                                 |  0.0023942|  0.6734948|  2.043701|             1|
| miR-509-3p alteration of YAP1/ECM axis%WikiPathways\_20180810%WP3967%Homo sapiens                                         |  0.0023942|  0.6964823|  1.969846|             1|
| DNA Damage Response (only ATM dependent)%WikiPathways\_20180810%WP710%Homo sapiens                                        |  0.0023942|  0.4507485|  1.870473|             1|
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                                  |  0.0023942|  0.4089892|  1.797154|             1|
| Nanoparticle-mediated activation of receptor signaling%WikiPathways\_20180810%WP2643%Homo sapiens                         |  0.0023942|  0.6690084|  2.134760|             1|
| Aryl Hydrocarbon Receptor Pathway%WikiPathways\_20180810%WP2873%Homo sapiens                                              |  0.0023942|  0.6066411|  1.973806|             1|
| Canonical and Non-Canonical TGF-B signaling%WikiPathways\_20180810%WP3874%Homo sapiens                                    |  0.0023942|  0.7019356|  1.985269|             1|
| Apoptosis-related network due to altered Notch3 in ovarian cancer%WikiPathways\_20180810%WP2864%Homo sapiens              |  0.0023942|  0.5388780|  1.965580|             1|
| Wnt Signaling in Kidney Disease%WikiPathways\_20180810%WP4150%Homo sapiens                                                |  0.0023942|  0.6437870|  2.072201|             1|
| Interleukin-11 Signaling Pathway%WikiPathways\_20180810%WP2332%Homo sapiens                                               |  0.0023942|  0.5417170|  1.911078|             1|
| T-Cell antigen Receptor (TCR) Signaling Pathway%WikiPathways\_20180810%WP69%Homo sapiens                                  |  0.0023942|  0.4604169|  1.820539|             1|
| Breast cancer pathway%WikiPathways\_20180810%WP4262%Homo sapiens                                                          |  0.0023942|  0.4106147|  1.784110|             1|
| Chromosomal and microsatellite instability in colorectal cancer %WikiPathways\_20180810%WP4216%Homo sapiens               |  0.0023942|  0.4741848|  1.867889|             1|
| Ras Signaling%WikiPathways\_20180810%WP4223%Homo sapiens                                                                  |  0.0023942|  0.3830946|  1.703206|             1|
| Photodynamic therapy-induced AP-1 survival signaling.%WikiPathways\_20180810%WP3611%Homo sapiens                          |  0.0033864|  0.5379339|  1.936787|             2|
| Hypothesized Pathways in Pathogenesis of Cardiovascular Disease%WikiPathways\_20180810%WP3668%Homo sapiens                |  0.0034085|  0.6625199|  2.048994|             2|
| Regulation of toll-like receptor signaling pathway%WikiPathways\_20180810%WP1449%Homo sapiens                             |  0.0040636|  0.4109383|  1.714699|             3|
| MAPK Cascade%WikiPathways\_20180810%WP422%Homo sapiens                                                                    |  0.0040636|  0.5826328|  1.895690|             3|
| ErbB Signaling Pathway%WikiPathways\_20180810%WP673%Homo sapiens                                                          |  0.0040636|  0.5158083|  1.905505|             3|
| Pancreatic adenocarcinoma pathway%WikiPathways\_20180810%WP4263%Homo sapiens                                              |  0.0040636|  0.4176090|  1.705760|             3|
| MET in type 1 papillary renal cell carcinoma%WikiPathways\_20180810%WP4205%Homo sapiens                                   |  0.0040636|  0.4875766|  1.833300|             3|
| Simplified Interaction Map Between LOXL4 and Oxidative Stress Pathway%WikiPathways\_20180810%WP3670%Homo sapiens          |  0.0040900|  0.6747657|  1.943240|             3|
| Calcium Regulation in the Cardiac Cell%WikiPathways\_20180810%WP536%Homo sapiens                                          |  0.0045597|  0.3874143|  1.663334|             4|
| Vitamin D Receptor Pathway%WikiPathways\_20180810%WP2877%Homo sapiens                                                     |  0.0045597|  0.3887124|  1.668908|             4|
| Leptin signaling pathway%WikiPathways\_20180810%WP2034%Homo sapiens                                                       |  0.0045993|  0.4436301|  1.754162|             4|
| Splicing factor NOVA regulated synaptic proteins%WikiPathways\_20180810%WP4148%Homo sapiens                               |  0.0045993|  0.5451452|  1.886074|             4|
| IL-3 Signaling Pathway%WikiPathways\_20180810%WP286%Homo sapiens                                                          |  0.0045993|  0.5416330|  1.873923|             4|
| TGF-B Signaling in Thyroid Cells for Epithelial-Mesenchymal Transition%WikiPathways\_20180810%WP3859%Homo sapiens         |  0.0046042|  0.6823482|  1.929870|             4|
| Angiogenesis%WikiPathways\_20180810%WP1539%Homo sapiens                                                                   |  0.0046042|  0.5959847|  1.843219|             4|
| B Cell Receptor Signaling Pathway%WikiPathways\_20180810%WP23%Homo sapiens                                                |  0.0050496|  0.4121499|  1.683461|             5|
| Aryl Hydrocarbon Receptor%WikiPathways\_20180810%WP2586%Homo sapiens                                                      |  0.0051978|  0.5267824|  1.847505|             5|
| Insulin Signaling%WikiPathways\_20180810%WP481%Homo sapiens                                                               |  0.0056660|  0.3456344|  1.538057|             6|
| Differentiation Pathway%WikiPathways\_20180810%WP2848%Homo sapiens                                                        |  0.0058540|  0.5333665|  1.802313|             6|
| Chemokine signaling pathway%WikiPathways\_20180810%WP3929%Homo sapiens                                                    |  0.0062487|  0.3765263|  1.615195|             7|
| IL-6 signaling pathway%WikiPathways\_20180810%WP364%Homo sapiens                                                          |  0.0072538|  0.5356893|  1.853359|             8|
| Hedgehog Signaling Pathway%WikiPathways\_20180810%WP4249%Homo sapiens                                                     |  0.0072538|  0.5251672|  1.805721|             8|
| Signaling Pathways in Glioblastoma%WikiPathways\_20180810%WP2261%Homo sapiens                                             |  0.0076318|  0.4090615|  1.649860|             9|
| Alpha 6 Beta 4 signaling pathway%WikiPathways\_20180810%WP244%Homo sapiens                                                |  0.0078496|  0.5398652|  1.783160|             9|
| BMP Signaling Pathway in Eyelid Development%WikiPathways\_20180810%WP3927%Homo sapiens                                    |  0.0079073|  0.6462598|  1.861147|             9|
| Toll-like Receptor Signaling Pathway%WikiPathways\_20180810%WP75%Homo sapiens                                             |  0.0080464|  0.4319551|  1.691790|            10|
| Oncostatin M Signaling Pathway%WikiPathways\_20180810%WP2374%Homo sapiens                                                 |  0.0080464|  0.4495244|  1.720956|            10|
| Hepatitis C and Hepatocellular Carcinoma%WikiPathways\_20180810%WP3646%Homo sapiens                                       |  0.0088568|  0.4909576|  1.767653|            11|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                                       |  0.0090256|  0.3531775|  1.541385|            12|
| Spinal Cord Injury%WikiPathways\_20180810%WP2431%Homo sapiens                                                             |  0.0090464|  0.4004925|  1.644643|            12|
| Arrhythmogenic Right Ventricular Cardiomyopathy%WikiPathways\_20180810%WP2118%Homo sapiens                                |  0.0105183|  0.4427742|  1.687848|            14|
| IL-1 signaling pathway%WikiPathways\_20180810%WP195%Homo sapiens                                                          |  0.0111789|  0.4556654|  1.700271|            15|
| Wnt Signaling Pathway%WikiPathways\_20180810%WP363%Homo sapiens                                                           |  0.0118188|  0.4774436|  1.741496|            16|
| TNF related weak inducer of apoptosis (TWEAK) Signaling Pathway%WikiPathways\_20180810%WP2036%Homo sapiens                |  0.0136569|  0.4870180|  1.708045|            19|
| Regulation of Actin Cytoskeleton%WikiPathways\_20180810%WP51%Homo sapiens                                                 |  0.0136569|  0.3496077|  1.514128|            20|
| Kit receptor signaling pathway%WikiPathways\_20180810%WP304%Homo sapiens                                                  |  0.0139387|  0.4431447|  1.677332|            20|
| Structural Pathway of Interleukin 1 (IL-1)%WikiPathways\_20180810%WP2637%Homo sapiens                                     |  0.0145690|  0.4708771|  1.724623|            21|
| Signal Transduction of S1P Receptor%WikiPathways\_20180810%WP26%Homo sapiens                                              |  0.0151386|  0.5631428|  1.760016|            22|
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                                 |  0.0153681|  0.3040116|  1.422908|            25|
| Estrogen signaling pathway%WikiPathways\_20180810%WP712%Homo sapiens                                                      |  0.0156973|  0.6054016|  1.795044|            23|
| Ebola Virus Pathway on Host%WikiPathways\_20180810%WP4217%Homo sapiens                                                    |  0.0187479|  0.3590918|  1.521116|            30|
| IL-5 Signaling Pathway%WikiPathways\_20180810%WP127%Homo sapiens                                                          |  0.0216496|  0.4821209|  1.657712|            34|
| Hypertrophy Model%WikiPathways\_20180810%WP516%Homo sapiens                                                               |  0.0216496|  0.6215551|  1.757931|            33|
| Canonical and Non-canonical Notch signaling%WikiPathways\_20180810%WP3845%Homo sapiens                                    |  0.0237981|  0.5497633|  1.718200|            38|
| TNF alpha Signaling Pathway%WikiPathways\_20180810%WP231%Homo sapiens                                                     |  0.0237981|  0.3741306|  1.528169|            40|
| Notch Signaling Pathway%WikiPathways\_20180810%WP61%Homo sapiens                                                          |  0.0241276|  0.4203454|  1.591035|            40|
| Brain-Derived Neurotrophic Factor (BDNF) signaling pathway%WikiPathways\_20180810%WP2380%Homo sapiens                     |  0.0246255|  0.3362748|  1.467615|            43|
| Inflammatory Response Pathway%WikiPathways\_20180810%WP453%Homo sapiens                                                   |  0.0256410|  0.5757371|  1.728930|            42|
| Endometrial cancer%WikiPathways\_20180810%WP4155%Homo sapiens                                                             |  0.0288624|  0.4077715|  1.564159|            50|
| Lipid Metabolism Pathway%WikiPathways\_20180810%WP3965%Homo sapiens                                                       |  0.0295799|  0.5429885|  1.697027|            50|
| IL-4 Signaling Pathway%WikiPathways\_20180810%WP395%Homo sapiens                                                          |  0.0305311|  0.4267597|  1.576540|            53|
| Thymic Stromal LymphoPoietin (TSLP) Signaling Pathway%WikiPathways\_20180810%WP2203%Homo sapiens                          |  0.0317370|  0.4539192|  1.610629|            56|
| Extracellular vesicle-mediated signaling in recipient cells%WikiPathways\_20180810%WP2870%Homo sapiens                    |  0.0317370|  0.5064609|  1.647853|            56|
| RAC1/PAK1/p38/MMP2 Pathway%WikiPathways\_20180810%WP3303%Homo sapiens                                                     |  0.0347231|  0.3955655|  1.527949|            64|
| NOTCH1 regulation of human endothelial cell calcification%WikiPathways\_20180810%WP3413%Homo sapiens                      |  0.0351213|  0.6045801|  1.689952|            62|
| PDGF Pathway%WikiPathways\_20180810%WP2526%Homo sapiens                                                                   |  0.0354631|  0.4616787|  1.587424|            65|
| Signaling of Hepatocyte Growth Factor Receptor%WikiPathways\_20180810%WP313%Homo sapiens                                  |  0.0365279|  0.4700848|  1.588476|            68|
| Serotonin Receptor 2 and ELK-SRF/GATA4 signaling%WikiPathways\_20180810%WP732%Homo sapiens                                |  0.0401205|  0.5713155|  1.669647|            73|
| Non-small cell lung cancer%WikiPathways\_20180810%WP4255%Homo sapiens                                                     |  0.0452470|  0.3844523|  1.496266|            88|
| Imatinib and Chronic Myeloid Leukemia%WikiPathways\_20180810%WP3640%Homo sapiens                                          |  0.0452614|  0.5643316|  1.649237|            84|
| EPO Receptor Signaling%WikiPathways\_20180810%WP581%Homo sapiens                                                          |  0.0460100|  0.5277873|  1.632303|            87|
| Viral Acute Myocarditis%WikiPathways\_20180810%WP4298%Homo sapiens                                                        |  0.0492995|  0.3746350|  1.469384|           100|
