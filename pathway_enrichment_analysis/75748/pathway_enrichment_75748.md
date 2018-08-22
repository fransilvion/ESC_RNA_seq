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

| Pathway                                                                    | Corrected p-value | Corrected MF p-value |
|:---------------------------------------------------------------------------|:------------------|:---------------------|
| Extracellular matrix organization                                          | 0.0000000         | 0.0000000            |
| TCF dependent signaling in response to WNT                                 | 0.0000000         | 0.0619000            |
| Processing of DNA double-strand break ends                                 | 0.1664422         | 0.2174049            |
| Senescence-Associated Secretory Phenotype (SASP)                           | 0.1650667         | 0.2180568            |
| Transcriptional regulation by small RNAs                                   | 0.1604347         | 0.2210714            |
| Neurotransmitter release cycle                                             | 0.1641696         | 0.2228400            |
| NoRC negatively regulates rRNA expression                                  | 0.1702250         | 0.2231279            |
| B-WICH complex positively regulates rRNA expression                        | 0.1654635         | 0.2285538            |
| RUNX1 regulates transcription of genes involved in differentiation of HSCs | 0.2476000         | 0.2310933            |
| Chromosome Maintenance                                                     | 0.1745897         | 0.2329395            |

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

There are 89 up-regulated and 51 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:----------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| FACTORS INVOLVED IN MEGAKARYOCYTE DEVELOPMENT AND PLATELET PRODUCTION%REACTOME%R-HSA-983231.2             |  0.0081530|  0.3952818|  1.747861|             0|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                      |  0.0081530|  0.7081415|  2.533950|             0|
| NRAGE SIGNALS DEATH THROUGH JNK%REACTOME DATABASE ID RELEASE 65%193648                                    |  0.0081530|  0.5352238|  2.041896|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                       |  0.0081530|  0.6786716|  2.610652|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                 |  0.0081530|  0.4685203|  2.280334|             0|
| SIGNALING BY NOTCH3%REACTOME DATABASE ID RELEASE 65%9012852                                               |  0.0081530|  0.5592139|  2.080504|             0|
| CELLULAR SENESCENCE%REACTOME%R-HSA-2559583.2                                                              |  0.0081530|  0.3550414|  1.667183|             0|
| EPH-EPHRIN MEDIATED REPULSION OF CELLS%REACTOME%R-HSA-3928665.3                                           |  0.0081530|  0.5557962|  2.085483|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                               |  0.0081530|  0.6177198|  2.491723|             0|
| EPH-EPHRIN SIGNALING%REACTOME%R-HSA-2682334.1                                                             |  0.0081530|  0.4808375|  2.027737|             0|
| SIGNALING BY PDGF%REACTOME DATABASE ID RELEASE 65%186797                                                  |  0.0081530|  0.6334639|  2.385575|             0|
| RHO GTPASE CYCLE%REACTOME DATABASE ID RELEASE 65%194840                                                   |  0.0081530|  0.4909874|  2.191390|             0|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                              |  0.0081530|  0.4468272|  1.924674|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090      |  0.0081530|  0.6219693|  2.313980|             0|
| SIGNALING BY WNT%REACTOME DATABASE ID RELEASE 65%195721                                                   |  0.0081530|  0.3292186|  1.639591|             0|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                           |  0.0117967|  0.5294765|  2.083425|             1|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                              |  0.0117967|  0.5435469|  2.031360|             1|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                                |  0.0117967|  0.5877192|  1.986562|             1|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                |  0.0117967|  0.5225696|  1.993620|             1|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                               |  0.0117967|  0.4546404|  1.920331|             1|
| DOWNSTREAM SIGNAL TRANSDUCTION%REACTOME DATABASE ID RELEASE 65%186763                                     |  0.0143592|  0.5913993|  1.947170|             2|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                 |  0.0143592|  0.5624801|  2.070155|             2|
| HS-GAG BIOSYNTHESIS%REACTOME%R-HSA-2022928.1                                                              |  0.0143592|  0.5863204|  1.960843|             2|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                        |  0.0143592|  0.4564114|  1.902880|             2|
| SEMA4D IN SEMAPHORIN SIGNALING%REACTOME%R-HSA-400685.2                                                    |  0.0143592|  0.6456456|  2.040166|             2|
| MET PROMOTES CELL MOTILITY%REACTOME DATABASE ID RELEASE 65%8875878                                        |  0.0143592|  0.6095724|  2.007004|             2|
| NOTCH3 INTRACELLULAR DOMAIN REGULATES TRANSCRIPTION%REACTOME%R-HSA-9013508.1                              |  0.0154652|  0.6446422|  2.011912|             3|
| DOWNSTREAM SIGNALING OF ACTIVATED FGFR1%REACTOME%R-HSA-5654687.2                                          |  0.0154652|  0.6662548|  2.003806|             3|
| SEMA4D INDUCED CELL MIGRATION AND GROWTH-CONE COLLAPSE%REACTOME%R-HSA-416572.1                            |  0.0154652|  0.6602435|  1.985726|             3|
| ONCOGENIC MAPK SIGNALING%REACTOME%R-HSA-6802957.2                                                         |  0.0154652|  0.4635441|  1.869818|             3|
| SIGNALING BY NTRK1 (TRKA)%REACTOME%R-HSA-187037.2                                                         |  0.0154652|  0.4430739|  1.813553|             3|
| SIGNALING BY RAS MUTANTS%REACTOME%R-HSA-6802949.1                                                         |  0.0154652|  0.5381687|  2.002207|             3|
| NON-INTEGRIN MEMBRANE-ECM INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000171                            |  0.0154652|  0.5647645|  2.020902|             3|
| SIGNALING BY FGFR1%REACTOME DATABASE ID RELEASE 65%5654736                                                |  0.0165048|  0.5538595|  1.958506|             4|
| SIGNALING BY NTRKS%REACTOME DATABASE ID RELEASE 65%166520                                                 |  0.0165048|  0.4226017|  1.807815|             4|
| GLYCOSAMINOGLYCAN METABOLISM%REACTOME DATABASE ID RELEASE 65%1630316                                      |  0.0165048|  0.3918191|  1.702801|             4|
| PLASMA LIPOPROTEIN CLEARANCE%REACTOME DATABASE ID RELEASE 65%8964043                                      |  0.0165048|  0.5829550|  1.919367|             4|
| N-GLYCAN TRIMMING IN THE ER AND CALNEXIN CALRETICULIN CYCLE%REACTOME%R-HSA-532668.2                       |  0.0165048|  0.5402788|  1.896802|             4|
| PRE-NOTCH PROCESSING IN GOLGI%REACTOME%R-HSA-1912420.2                                                    |  0.0165048|  0.6905442|  1.989405|             4|
| O-LINKED GLYCOSYLATION%REACTOME%R-HSA-5173105.3                                                           |  0.0166500|  0.4268245|  1.774328|             5|
| PRE-NOTCH EXPRESSION AND PROCESSING%REACTOME DATABASE ID RELEASE 65%1912422                               |  0.0166500|  0.4081292|  1.738452|             5|
| G ALPHA (12 13) SIGNALLING EVENTS%REACTOME%R-HSA-416482.3                                                 |  0.0166500|  0.4293278|  1.762169|             5|
| INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING%REACTOME DATABASE ID RELEASE 65%6785807                        |  0.0166500|  0.4048856|  1.703868|             5|
| CHONDROITIN SULFATE DERMATAN SULFATE METABOLISM%REACTOME DATABASE ID RELEASE 65%1793185                   |  0.0166500|  0.5015270|  1.855231|             5|
| SIGNALING BY VEGF%REACTOME DATABASE ID RELEASE 65%194138                                                  |  0.0166500|  0.4062985|  1.730654|             5|
| LDL CLEARANCE%REACTOME%R-HSA-8964038.1                                                                    |  0.0168126|  0.6550377|  1.938471|             5|
| CELL-CELL COMMUNICATION%REACTOME%R-HSA-1500931.3                                                          |  0.0185632|  0.3838493|  1.664032|             6|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                           |  0.0187671|  0.5166638|  1.848783|             6|
| TRANSCRIPTIONAL REGULATION BY THE AP-2 (TFAP2) FAMILY OF TRANSCRIPTION FACTORS%REACTOME%R-HSA-8864260.3   |  0.0188120|  0.5613401|  1.877301|             6|
| PARADOXICAL ACTIVATION OF RAF SIGNALING BY KINASE INACTIVE BRAF%REACTOME DATABASE ID RELEASE 65%6802955   |  0.0199870|  0.5317474|  1.854208|             7|
| SIGNALING BY TGF-BETA RECEPTOR COMPLEX%REACTOME DATABASE ID RELEASE 65%170834                             |  0.0199870|  0.4227119|  1.725628|             7|
| SIGNALING BY MODERATE KINASE ACTIVITY BRAF MUTANTS%REACTOME DATABASE ID RELEASE 65%6802946                |  0.0199870|  0.5317474|  1.854208|             7|
| TCF DEPENDENT SIGNALING IN RESPONSE TO WNT%REACTOME%R-HSA-201681.1                                        |  0.0203980|  0.3211413|  1.543141|             8|
| CLATHRIN-MEDIATED ENDOCYTOSIS%REACTOME%R-HSA-8856828.3                                                    |  0.0205539|  0.3553193|  1.585872|             8|
| ACTIVATION OF BH3-ONLY PROTEINS%REACTOME DATABASE ID RELEASE 65%114452                                    |  0.0209757|  0.5633874|  1.869433|             8|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                         |  0.0209757|  0.6333236|  1.904763|             8|
| L1CAM INTERACTIONS%REACTOME%R-HSA-373760.2                                                                |  0.0235347|  0.3913974|  1.648328|            10|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                            |  0.0236603|  0.6638013|  1.912361|            10|
| CHROMATIN ORGANIZATION%REACTOME%R-HSA-4839726.2                                                           |  0.0236603|  0.3022272|  1.479312|            11|
| CHROMATIN MODIFYING ENZYMES%REACTOME%R-HSA-3247509.4                                                      |  0.0236603|  0.3022272|  1.479312|            11|
| VEGFA-VEGFR2 PATHWAY%REACTOME%R-HSA-4420097.3                                                             |  0.0244254|  0.3945066|  1.657060|            11|
| DISASSEMBLY OF THE DESTRUCTION COMPLEX AND RECRUITMENT OF AXIN TO THE MEMBRANE%REACTOME%R-HSA-4641262.3   |  0.0248731|  0.5488056|  1.821047|            11|
| INTRACELLULAR SIGNALING BY SECOND MESSENGERS%REACTOME DATABASE ID RELEASE 65%9006925                      |  0.0292896|  0.2992943|  1.463810|            15|
| INTRINSIC PATHWAY FOR APOPTOSIS%REACTOME%R-HSA-109606.2                                                   |  0.0292896|  0.4877898|  1.779119|            14|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                             |  0.0314797|  0.5607128|  1.829316|            15|
| DAG AND IP3 SIGNALING%REACTOME DATABASE ID RELEASE 65%1489509                                             |  0.0327454|  0.5378068|  1.784551|            16|
| HEPARAN SULFATE HEPARIN (HS-GAG) METABOLISM%REACTOME%R-HSA-1638091.1                                      |  0.0348140|  0.4529951|  1.717664|            18|
| TGF-BETA RECEPTOR SIGNALING ACTIVATES SMADS%REACTOME DATABASE ID RELEASE 65%2173789                       |  0.0353597|  0.5199507|  1.767828|            18|
| LAMININ INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000157                                              |  0.0370127|  0.5939494|  1.827885|            19|
| PRE-NOTCH TRANSCRIPTION AND TRANSLATION%REACTOME%R-HSA-1912408.3                                          |  0.0370127|  0.3876789|  1.603202|            21|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                           |  0.0370127|  0.3887661|  1.612117|            21|
| TP53 REGULATES TRANSCRIPTION OF CELL CYCLE GENES%REACTOME DATABASE ID RELEASE 65%6791312                  |  0.0370127|  0.4637610|  1.733182|            21|
| CALMODULIN INDUCED EVENTS%REACTOME DATABASE ID RELEASE 65%111933                                          |  0.0400752|  0.5727439|  1.787519|            23|
| CAM PATHWAY%REACTOME DATABASE ID RELEASE 65%111997                                                        |  0.0400752|  0.5727439|  1.787519|            23|
| REGULATION OF PTEN GENE TRANSCRIPTION%REACTOME DATABASE ID RELEASE 65%8943724                             |  0.0409665|  0.4358255|  1.683729|            25|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                                |  0.0409665|  0.4581408|  1.694738|            25|
| ASPARAGINE N-LINKED GLYCOSYLATION%REACTOME%R-HSA-446203.4                                                 |  0.0417972|  0.2836056|  1.403947|            28|
| INTRA-GOLGI AND RETROGRADE GOLGI-TO-ER TRAFFIC%REACTOME DATABASE ID RELEASE 65%6811442                    |  0.0417972|  0.3121420|  1.462114|            27|
| O-GLYCOSYLATION OF TSR DOMAIN-CONTAINING PROTEINS%REACTOME%R-HSA-5173214.1                                |  0.0428639|  0.5090161|  1.730650|            27|
| REGULATION OF FZD BY UBIQUITINATION%REACTOME%R-HSA-4641263.2                                              |  0.0434330|  0.6037169|  1.815719|            28|
| SIGNALING BY FGFR3 IN DISEASE%REACTOME DATABASE ID RELEASE 65%5655332                                     |  0.0438841|  0.6282799|  1.810026|            29|
| SIGNALING BY FGFR3 POINT MUTANTS IN CANCER%REACTOME%R-HSA-8853338.2                                       |  0.0438841|  0.6282799|  1.810026|            29|
| IRS-RELATED EVENTS TRIGGERED BY IGF1R%REACTOME%R-HSA-2428928.1                                            |  0.0452425|  0.4797510|  1.696451|            31|
| SIGNALING BY TYPE 1 INSULIN-LIKE GROWTH FACTOR 1 RECEPTOR (IGF1R)%REACTOME DATABASE ID RELEASE 65%2404192 |  0.0452425|  0.4797510|  1.696451|            31|
| IGF1R SIGNALING CASCADE%REACTOME%R-HSA-2428924.1                                                          |  0.0452425|  0.4797510|  1.696451|            31|
| CELL DEATH SIGNALLING VIA NRAGE, NRIF AND NADE%REACTOME DATABASE ID RELEASE 65%204998                     |  0.0458063|  0.3994354|  1.615989|            33|
| CA-DEPENDENT EVENTS%REACTOME DATABASE ID RELEASE 65%111996                                                |  0.0458242|  0.5519660|  1.766785|            32|
| HDACS DEACETYLATE HISTONES%REACTOME%R-HSA-3214815.2                                                       |  0.0491349|  0.3765636|  1.569977|            36|
| ESR-MEDIATED SIGNALING%REACTOME DATABASE ID RELEASE 65%8939211                                            |  0.0492621|  0.3294227|  1.488062|            37|

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
| Kegg wnt signaling pathway                                | 0.0081000         | 0.0000000            |
| Kegg pathways in cancer                                   | 0.0445500         | 0.0000000            |
| Kegg mapk signaling pathway                               | 0.0162000         | 0.0000000            |
| Kegg focal adhesion                                       | 0.0432000         | 0.0040500            |
| Kegg basal cell carcinoma                                 | 0.0648000         | 0.0129600            |
| Kegg melanogenesis                                        | 0.1296000         | 0.0297000            |
| Kegg parkinsons disease                                   | 0.1152000         | 0.0826200            |
| Kegg alzheimers disease                                   | 0.1708364         | 0.0833143            |
| Kegg ecm receptor interaction                             | 0.1296000         | 0.0870750            |
| Kegg arrhythmogenic right ventricular cardiomyopathy arvc | 0.2018769         | 0.0882000            |

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

There are 38 up-regulated and 11 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                                        |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_MAPK\_SIGNALING\_PATHWAY                                 |  0.0017512|  0.4056368|  1.974749|             0|
| KEGG\_ENDOCYTOSIS                                              |  0.0017512|  0.4037747|  1.886421|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY                                  |  0.0017512|  0.4299674|  1.932222|             0|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY                            |  0.0017512|  0.4547765|  1.864011|             0|
| KEGG\_AXON\_GUIDANCE                                           |  0.0017512|  0.4493756|  1.999892|             0|
| KEGG\_FOCAL\_ADHESION                                          |  0.0017512|  0.5544194|  2.607173|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                               |  0.0017512|  0.6108126|  2.448137|             0|
| KEGG\_NEUROTROPHIN\_SIGNALING\_PATHWAY                         |  0.0017512|  0.4347978|  1.922440|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                                     |  0.0017512|  0.3939500|  1.969469|             0|
| KEGG\_COLORECTAL\_CANCER                                       |  0.0017512|  0.5037259|  1.983489|             0|
| KEGG\_RENAL\_CELL\_CARCINOMA                                   |  0.0017512|  0.4860469|  1.954449|             0|
| KEGG\_PANCREATIC\_CANCER                                       |  0.0017512|  0.4877522|  1.967138|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA                                   |  0.0017512|  0.5864775|  2.144277|             0|
| KEGG\_CHRONIC\_MYELOID\_LEUKEMIA                               |  0.0017512|  0.4665465|  1.892074|             0|
| KEGG\_SMALL\_CELL\_LUNG\_CANCER                                |  0.0017512|  0.4679651|  1.944568|             0|
| KEGG\_ARRHYTHMOGENIC\_RIGHT\_VENTRICULAR\_CARDIOMYOPATHY\_ARVC |  0.0017512|  0.5006561|  1.961426|             0|
| KEGG\_DILATED\_CARDIOMYOPATHY                                  |  0.0017512|  0.4775598|  1.920321|             0|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY                             |  0.0028874|  0.5263141|  1.910763|             1|
| KEGG\_INSULIN\_SIGNALING\_PATHWAY                              |  0.0028874|  0.4063390|  1.808362|             1|
| KEGG\_ERBB\_SIGNALING\_PATHWAY                                 |  0.0038995|  0.4317931|  1.788431|             2|
| KEGG\_REGULATION\_OF\_ACTIN\_CYTOSKELETON                      |  0.0038995|  0.3542808|  1.667174|             2|
| KEGG\_ENDOMETRIAL\_CANCER                                      |  0.0047229|  0.4891885|  1.858982|             3|
| KEGG\_BLADDER\_CANCER                                          |  0.0047229|  0.5453271|  1.936910|             3|
| KEGG\_MELANOGENESIS                                            |  0.0055559|  0.4276664|  1.770169|             4|
| KEGG\_CHEMOKINE\_SIGNALING\_PATHWAY                            |  0.0062825|  0.3725686|  1.667353|             5|
| KEGG\_NON\_SMALL\_CELL\_LUNG\_CANCER                           |  0.0062825|  0.4732919|  1.814964|             5|
| KEGG\_GNRH\_SIGNALING\_PATHWAY                                 |  0.0087505|  0.4131290|  1.709997|             8|
| KEGG\_NOTCH\_SIGNALING\_PATHWAY                                |  0.0095743|  0.4997174|  1.827065|             9|
| KEGG\_HYPERTROPHIC\_CARDIOMYOPATHY\_HCM                        |  0.0138809|  0.4366952|  1.728986|            15|
| KEGG\_GLYCOSAMINOGLYCAN\_BIOSYNTHESIS\_CHONDROITIN\_SULFATE    |  0.0185765|  0.5813376|  1.765884|            21|
| KEGG\_ACUTE\_MYELOID\_LEUKEMIA                                 |  0.0209175|  0.4398111|  1.679516|            25|
| KEGG\_GAP\_JUNCTION                                            |  0.0210172|  0.3969463|  1.631465|            26|
| KEGG\_TYPE\_II\_DIABETES\_MELLITUS                             |  0.0218870|  0.4709110|  1.696237|            27|
| KEGG\_B\_CELL\_RECEPTOR\_SIGNALING\_PATHWAY                    |  0.0304359|  0.3984707|  1.597071|            40|
| KEGG\_VASOPRESSIN\_REGULATED\_WATER\_REABSORPTION              |  0.0311381|  0.4731124|  1.680416|            41|
| KEGG\_ADHERENS\_JUNCTION                                       |  0.0368189|  0.3871865|  1.570231|            51|
| KEGG\_PROSTATE\_CANCER                                         |  0.0371363|  0.3672851|  1.533986|            54|
| KEGG\_O\_GLYCAN\_BIOSYNTHESIS                                  |  0.0499109|  0.5008774|  1.615182|            76|

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
| Neuronal System                                      | 0.0000000         | 0.0412667            |
| G alpha (i) signalling events                        | 0.0464250         | 0.0495200            |
| Class A 1 (Rhodopsin-like receptors)                 | 0.0495200         | 0.0619000            |
| Collagen formation                                   | 0.0884286         | 0.1134833            |
| Platelet activation, signaling and aggregation       | 0.0722167         | 0.1326429            |
| Collagen chain trimerization                         | 0.3327125         | 0.3327125            |
| Binding and Uptake of Ligands by Scavenger Receptors | 0.3280700         | 0.3466400            |
| ECM proteoglycans                                    | 0.3438889         | 0.3507667            |

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

There are 22 up-regulated and 66 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                              |       padj|         ES|       NES|  nMoreExtreme|
|:-----------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                 |  0.0054835|  0.6411296|  2.002102|             0|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                |  0.0054835|  0.4229049|  1.735654|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                  |  0.0054835|  0.6223854|  2.078607|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                            |  0.0054835|  0.5076242|  2.111459|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                          |  0.0054835|  0.5986648|  2.101459|             0|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                           |  0.0054835|  0.4408896|  1.794462|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090 |  0.0054835|  0.6139364|  1.989166|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                          |  0.0068563|  0.5305453|  1.945435|             1|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                            |  0.0072743|  0.6009169|  1.927318|             1|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                   |  0.0086177|  0.5061851|  1.830426|             2|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                      |  0.0091009|  0.6215197|  1.940865|             2|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                  |  0.0101540|  0.4773152|  1.770368|             3|
| CLASS A 1 (RHODOPSIN-LIKE RECEPTORS)%REACTOME%R-HSA-373076.6                                         |  0.0101540|  0.4328159|  1.671655|             3|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                        |  0.0101540|  0.4983576|  1.835386|             3|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                           |  0.0128786|  0.5639546|  1.817642|             4|
| G ALPHA (Q) SIGNALLING EVENTS%REACTOME DATABASE ID RELEASE 65%416476                                 |  0.0137966|  0.4091554|  1.582657|             5|
| COMPLEMENT CASCADE%REACTOME DATABASE ID RELEASE 65%166658                                            |  0.0178956|  0.6803816|  1.855732|             6|
| PEPTIDE LIGAND-BINDING RECEPTORS%REACTOME%R-HSA-375276.4                                             |  0.0182613|  0.4911359|  1.736403|             7|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                              |  0.0198358|  0.7521806|  1.891347|             7|
| DISEASES OF GLYCOSYLATION%REACTOME DATABASE ID RELEASE 65%3781865                                    |  0.0202218|  0.4987830|  1.720278|             8|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                    |  0.0380949|  0.6871840|  1.833565|            16|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                        |  0.0421276|  0.6020824|  1.764698|            19|

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
| Kegg ecm receptor interaction                  | 0.0081000         | 0.0081000            |
| Kegg cell adhesion molecules cams              | 0.0162000         | 0.0162000            |
| Kegg cytokine cytokine receptor interaction    | 0.0254571         | 0.0526500            |
| Kegg calcium signaling pathway                 | 0.0259200         | 0.0648000            |
| Kegg focal adhesion                            | 0.0202500         | 0.0756000            |
| Kegg taste transduction                        | 0.1656000         | 0.2152286            |
| Kegg hematopoietic cell lineage                | 0.1701000         | 0.2430000            |
| Kegg complement and coagulation cascades       | 0.2052000         | 0.2606727            |
| Kegg aldosterone regulated sodium reabsorption | 0.1900800         | 0.2615143            |

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

There are 4 up-regulated and 4 down-regulated Reactome pathways with fgsea.

``` r
upPathwaysKEGG %>% kable()
```

| pathway                                      |       padj|         ES|       NES|  nMoreExtreme|
|:---------------------------------------------|----------:|----------:|---------:|-------------:|
| KEGG\_FOCAL\_ADHESION                        |  0.0131303|  0.4368927|  1.768565|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION             |  0.0131303|  0.5627785|  1.970268|             0|
| KEGG\_COMPLEMENT\_AND\_COAGULATION\_CASCADES |  0.0151487|  0.6288321|  1.931924|             1|
| KEGG\_HEMATOPOIETIC\_CELL\_LINEAGE           |  0.0180110|  0.5695954|  1.829783|             2|

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
| Collagen chain trimerization                         | 0.0154750         | 0.0154750            |
| Collagen biosynthesis and modifying enzymes          | 0.0206333         | 0.0206333            |
| ECM proteoglycans                                    | 0.0247600         | 0.0247600            |
| G alpha (i) signalling events                        | 0.0515833         | 0.1547500            |
| Signaling by BMP                                     | 0.2407222         | 0.2544778            |
| Signaling by NODAL                                   | 0.2829714         | 0.2630750            |
| Class A 1 (Rhodopsin-like receptors)                 | 0.2414100         | 0.2723600            |
| Binding and Uptake of Ligands by Scavenger Receptors | 0.2553375         | 0.3006571            |

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

There are 131 up-regulated and 120 down-regulated Reactome pathways with fgsea

``` r
upPathwaysReactome %>% kable()
```

| pathway                                                                                                                                   |       padj|         ES|       NES|  nMoreExtreme|
|:------------------------------------------------------------------------------------------------------------------------------------------|----------:|----------:|---------:|-------------:|
| LAMININ INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000157                                                                              |  0.0035573|  0.7271213|  2.206171|             0|
| COLLAGEN CHAIN TRIMERIZATION%REACTOME DATABASE ID RELEASE 65%8948216                                                                      |  0.0035573|  0.7634609|  2.666482|             0|
| NCAM SIGNALING FOR NEURITE OUT-GROWTH%REACTOME DATABASE ID RELEASE 65%375165                                                              |  0.0035573|  0.5364065|  1.957657|             0|
| PLATELET ACTIVATION, SIGNALING AND AGGREGATION%REACTOME%R-HSA-76002.2                                                                     |  0.0035573|  0.3726225|  1.735045|             0|
| REGULATION OF FZD BY UBIQUITINATION%REACTOME%R-HSA-4641263.2                                                                              |  0.0035573|  0.6978179|  2.062981|             0|
| SIGNAL TRANSDUCTION BY L1%REACTOME%R-HSA-445144.1                                                                                         |  0.0035573|  0.6910985|  2.043117|             0|
| NRAGE SIGNALS DEATH THROUGH JNK%REACTOME DATABASE ID RELEASE 65%193648                                                                    |  0.0035573|  0.5086439|  1.899381|             0|
| COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES%REACTOME DATABASE ID RELEASE 65%1650814                                                       |  0.0035573|  0.7369268|  2.777169|             0|
| EXTRACELLULAR MATRIX ORGANIZATION%REACTOME DATABASE ID RELEASE 65%1474244                                                                 |  0.0035573|  0.5727275|  2.709274|             0|
| GLYCOSAMINOGLYCAN METABOLISM%REACTOME DATABASE ID RELEASE 65%1630316                                                                      |  0.0035573|  0.4614803|  1.964547|             0|
| SIGNALING BY NOTCH3%REACTOME DATABASE ID RELEASE 65%9012852                                                                               |  0.0035573|  0.5394003|  1.956473|             0|
| EPH-EPHRIN MEDIATED REPULSION OF CELLS%REACTOME%R-HSA-3928665.3                                                                           |  0.0035573|  0.5303114|  1.940477|             0|
| COLLAGEN FORMATION%REACTOME%R-HSA-1474290.1                                                                                               |  0.0035573|  0.6796685|  2.676519|             0|
| HEPARAN SULFATE HEPARIN (HS-GAG) METABOLISM%REACTOME%R-HSA-1638091.1                                                                      |  0.0035573|  0.5241214|  1.935890|             0|
| EPH-EPHRIN SIGNALING%REACTOME%R-HSA-2682334.1                                                                                             |  0.0035573|  0.4717783|  1.947875|             0|
| ECM PROTEOGLYCANS%REACTOME DATABASE ID RELEASE 65%3000178                                                                                 |  0.0035573|  0.6656855|  2.394527|             0|
| SIGNALING BY PDGF%REACTOME DATABASE ID RELEASE 65%186797                                                                                  |  0.0035573|  0.5737726|  2.107241|             0|
| HS-GAG BIOSYNTHESIS%REACTOME%R-HSA-2022928.1                                                                                              |  0.0035573|  0.6647125|  2.180776|             0|
| SIGNALING BY VEGF%REACTOME DATABASE ID RELEASE 65%194138                                                                                  |  0.0035573|  0.4638137|  1.931879|             0|
| POST-TRANSLATIONAL PROTEIN PHOSPHORYLATION%REACTOME DATABASE ID RELEASE 65%8957275                                                        |  0.0035573|  0.5364088|  2.192457|             0|
| COLLAGEN DEGRADATION%REACTOME%R-HSA-1442490.3                                                                                             |  0.0035573|  0.6276843|  2.007113|             0|
| PLATELET DEGRANULATION%REACTOME DATABASE ID RELEASE 65%114608                                                                             |  0.0035573|  0.4303519|  1.789878|             0|
| NETRIN-1 SIGNALING%REACTOME%R-HSA-373752.2                                                                                                |  0.0035573|  0.5966669|  2.156177|             0|
| ELASTIC FIBRE FORMATION%REACTOME DATABASE ID RELEASE 65%1566948                                                                           |  0.0035573|  0.7033912|  2.456681|             0|
| INTEGRIN CELL SURFACE INTERACTIONS%REACTOME%R-HSA-216083.2                                                                                |  0.0035573|  0.6060797|  2.263226|             0|
| MOLECULES ASSOCIATED WITH ELASTIC FIBRES%REACTOME%R-HSA-2129379.1                                                                         |  0.0035573|  0.6600929|  2.129481|             0|
| RHO GTPASE CYCLE%REACTOME DATABASE ID RELEASE 65%194840                                                                                   |  0.0035573|  0.5021043|  2.188029|             0|
| SIGNALING BY NTRK1 (TRKA)%REACTOME%R-HSA-187037.2                                                                                         |  0.0035573|  0.4692721|  1.886616|             0|
| ASSEMBLY OF COLLAGEN FIBRILS AND OTHER MULTIMERIC STRUCTURES%REACTOME DATABASE ID RELEASE 65%2022090                                      |  0.0035573|  0.7323880|  2.656464|             0|
| REGULATION OF IGF ACTIVITY BY IGFBP%REACTOME%R-HSA-381426.2                                                                               |  0.0035573|  0.5596483|  2.311641|             0|
| SIGNALING BY NTRKS%REACTOME DATABASE ID RELEASE 65%166520                                                                                 |  0.0057681|  0.4218281|  1.763564|             1|
| G ALPHA (12 13) SIGNALLING EVENTS%REACTOME%R-HSA-416482.3                                                                                 |  0.0057681|  0.4359478|  1.755687|             1|
| DEGRADATION OF THE EXTRACELLULAR MATRIX%REACTOME DATABASE ID RELEASE 65%1474228                                                           |  0.0057681|  0.4428251|  1.799972|             1|
| SIGNALING BY TGF-BETA FAMILY MEMBERS%REACTOME DATABASE ID RELEASE 65%9006936                                                              |  0.0057681|  0.4235370|  1.785296|             1|
| VEGFA-VEGFR2 PATHWAY%REACTOME%R-HSA-4420097.3                                                                                             |  0.0070885|  0.4313524|  1.773915|             2|
| L1CAM INTERACTIONS%REACTOME%R-HSA-373760.2                                                                                                |  0.0070885|  0.4252722|  1.753663|             2|
| NCAM1 INTERACTIONS%REACTOME%R-HSA-419037.1                                                                                                |  0.0070885|  0.6000930|  1.985978|             2|
| DISEASES OF GLYCOSYLATION%REACTOME DATABASE ID RELEASE 65%3781865                                                                         |  0.0070885|  0.4907506|  1.909891|             2|
| SCAVENGING BY CLASS A RECEPTORS%REACTOME DATABASE ID RELEASE 65%3000480                                                                   |  0.0070885|  0.7379233|  2.035867|             2|
| O-GLYCOSYLATION OF TSR DOMAIN-CONTAINING PROTEINS%REACTOME%R-HSA-5173214.1                                                                |  0.0081714|  0.5821632|  1.936071|             3|
| MET PROMOTES CELL MOTILITY%REACTOME DATABASE ID RELEASE 65%8875878                                                                        |  0.0081714|  0.6098729|  1.967470|             3|
| TOLL LIKE RECEPTOR TLR6:TLR2 CASCADE%REACTOME%R-HSA-168188.1                                                                              |  0.0092925|  0.4101515|  1.693431|             4|
| MYD88:MAL CASCADE INITIATED ON PLASMA MEMBRANE%REACTOME%R-HSA-166058.2                                                                    |  0.0092925|  0.4101515|  1.693431|             4|
| TOLL LIKE RECEPTOR 2 (TLR2) CASCADE%REACTOME DATABASE ID RELEASE 65%181438                                                                |  0.0092925|  0.4101515|  1.693431|             4|
| TOLL LIKE RECEPTOR TLR1:TLR2 CASCADE%REACTOME DATABASE ID RELEASE 65%168179                                                               |  0.0092925|  0.4101515|  1.693431|             4|
| INTRA-GOLGI AND RETROGRADE GOLGI-TO-ER TRAFFIC%REACTOME DATABASE ID RELEASE 65%6811442                                                    |  0.0092925|  0.3474726|  1.585648|             4|
| BINDING AND UPTAKE OF LIGANDS BY SCAVENGER RECEPTORS%REACTOME%R-HSA-2173782.1                                                             |  0.0094085|  0.5940740|  1.932565|             4|
| RUNX2 REGULATES BONE DEVELOPMENT%REACTOME%R-HSA-8941326.1                                                                                 |  0.0103983|  0.5949921|  1.919464|             5|
| N-GLYCAN TRIMMING IN THE ER AND CALNEXIN CALRETICULIN CYCLE%REACTOME%R-HSA-532668.2                                                       |  0.0103983|  0.5501566|  1.888318|             5|
| TP53 REGULATES TRANSCRIPTION OF CELL CYCLE GENES%REACTOME DATABASE ID RELEASE 65%6791312                                                  |  0.0103983|  0.4974348|  1.815427|             5|
| TRAF6 MEDIATED INDUCTION OF NFKB AND MAP KINASES UPON TLR7 8 OR 9 ACTIVATION%REACTOME%R-HSA-975138.1                                      |  0.0103983|  0.4133133|  1.689330|             5|
| NON-INTEGRIN MEMBRANE-ECM INTERACTIONS%REACTOME DATABASE ID RELEASE 65%3000171                                                            |  0.0103983|  0.5349328|  1.868319|             5|
| CELL-EXTRACELLULAR MATRIX INTERACTIONS%REACTOME%R-HSA-446353.1                                                                            |  0.0103983|  0.7135915|  1.968738|             5|
| RESPONSE TO ELEVATED PLATELET CYTOSOLIC CA2+%REACTOME%R-HSA-76005.2                                                                       |  0.0114123|  0.3986534|  1.675281|             6|
| CHONDROITIN SULFATE DERMATAN SULFATE METABOLISM%REACTOME DATABASE ID RELEASE 65%1793185                                                   |  0.0114891|  0.5070917|  1.832478|             6|
| DEFECTIVE B3GALTL CAUSES PETERS-PLUS SYNDROME (PPS)%REACTOME%R-HSA-5083635.1                                                              |  0.0129070|  0.5731837|  1.896923|             7|
| TOLL LIKE RECEPTOR 9 (TLR9) CASCADE%REACTOME DATABASE ID RELEASE 65%168138                                                                |  0.0138446|  0.4011339|  1.650703|             8|
| CLASS B 2 (SECRETIN FAMILY RECEPTORS)%REACTOME DATABASE ID RELEASE 65%373080                                                              |  0.0138979|  0.4816090|  1.788984|             8|
| ASPARAGINE N-LINKED GLYCOSYLATION%REACTOME%R-HSA-446203.4                                                                                 |  0.0139286|  0.3045354|  1.468983|             9|
| SIGNALING BY WNT%REACTOME DATABASE ID RELEASE 65%195721                                                                                   |  0.0139286|  0.3056332|  1.482095|             9|
| TOLL LIKE RECEPTOR 7 8 (TLR7 8) CASCADE%REACTOME DATABASE ID RELEASE 65%168181                                                            |  0.0143260|  0.4047850|  1.658423|             9|
| INTERLEUKIN-4 AND INTERLEUKIN-13 SIGNALING%REACTOME DATABASE ID RELEASE 65%6785807                                                        |  0.0143260|  0.3982062|  1.638655|             9|
| MYD88 DEPENDENT CASCADE INITIATED ON ENDOSOME%REACTOME%R-HSA-975155.1                                                                     |  0.0143260|  0.4047850|  1.658423|             9|
| SIGNALING BY EGFR%REACTOME DATABASE ID RELEASE 65%177929                                                                                  |  0.0144350|  0.5160705|  1.847308|             9|
| O-LINKED GLYCOSYLATION%REACTOME%R-HSA-5173105.3                                                                                           |  0.0161779|  0.4037059|  1.645124|            11|
| TRIF(TICAM1)-MEDIATED TLR4 SIGNALING%REACTOME%R-HSA-937061.2                                                                              |  0.0161779|  0.3917665|  1.618201|            11|
| CLATHRIN-MEDIATED ENDOCYTOSIS%REACTOME%R-HSA-8856828.3                                                                                    |  0.0161779|  0.3596602|  1.567298|            11|
| MYD88-INDEPENDENT TLR4 CASCADE%REACTOME%R-HSA-166166.2                                                                                    |  0.0161779|  0.3917665|  1.618201|            11|
| TOLL LIKE RECEPTOR 3 (TLR3) CASCADE%REACTOME DATABASE ID RELEASE 65%168164                                                                |  0.0168290|  0.3905721|  1.612591|            12|
| SIGNALING BY BMP%REACTOME%R-HSA-201451.4                                                                                                  |  0.0171570|  0.5650694|  1.822933|            12|
| MET ACTIVATES PTK2 SIGNALING%REACTOME DATABASE ID RELEASE 65%8874081                                                                      |  0.0171575|  0.6485362|  1.864456|            12|
| SIGNALLING TO RAS%REACTOME DATABASE ID RELEASE 65%167044                                                                                  |  0.0181078|  0.6459158|  1.856922|            13|
| DISEASES ASSOCIATED WITH O-GLYCOSYLATION OF PROTEINS%REACTOME%R-HSA-3906995.2                                                             |  0.0190112|  0.5061652|  1.787307|            14|
| MYD88 CASCADE INITIATED ON PLASMA MEMBRANE%REACTOME%R-HSA-975871.1                                                                        |  0.0203439|  0.4014290|  1.635846|            16|
| TOLL LIKE RECEPTOR 10 (TLR10) CASCADE%REACTOME DATABASE ID RELEASE 65%168142                                                              |  0.0203439|  0.4014290|  1.635846|            16|
| TOLL LIKE RECEPTOR 5 (TLR5) CASCADE%REACTOME%R-HSA-168176.1                                                                               |  0.0203439|  0.4014290|  1.635846|            16|
| CALNEXIN CALRETICULIN CYCLE%REACTOME%R-HSA-901042.2                                                                                       |  0.0208536|  0.5631319|  1.800697|            16|
| NOTCH3 ACTIVATION AND TRANSMISSION OF SIGNAL TO THE NUCLEUS%REACTOME DATABASE ID RELEASE 65%9013507                                       |  0.0208536|  0.5749506|  1.796359|            16|
| MAP KINASE ACTIVATION%REACTOME%R-HSA-450294.3                                                                                             |  0.0213209|  0.4285091|  1.667661|            17|
| LIPOPROTEIN METABOLISM%REACTOME DATABASE ID RELEASE 65%174824                                                                             |  0.0237646|  0.4638759|  1.703634|            19|
| DEATH RECEPTOR SIGNALLING%REACTOME%R-HSA-73887.3                                                                                          |  0.0240983|  0.3487413|  1.518305|            20|
| CELL SURFACE INTERACTIONS AT THE VASCULAR WALL%REACTOME%R-HSA-202733.4                                                                    |  0.0248860|  0.3752864|  1.579677|            21|
| ACTIVATION OF MATRIX METALLOPROTEINASES%REACTOME%R-HSA-1592389.1                                                                          |  0.0255584|  0.6350531|  1.825693|            21|
| TOLL LIKE RECEPTOR 4 (TLR4) CASCADE%REACTOME%R-HSA-166016.2                                                                               |  0.0257435|  0.3554246|  1.521942|            22|
| PEPTIDE HORMONE METABOLISM%REACTOME DATABASE ID RELEASE 65%2980736                                                                        |  0.0257934|  0.4425518|  1.684304|            22|
| GPCR LIGAND BINDING%REACTOME DATABASE ID RELEASE 65%500792                                                                                |  0.0262608|  0.3126855|  1.442187|            24|
| REGULATION OF PTEN GENE TRANSCRIPTION%REACTOME DATABASE ID RELEASE 65%8943724                                                             |  0.0280424|  0.4355678|  1.647588|            25|
| SIGNALING BY WNT IN CANCER%REACTOME DATABASE ID RELEASE 65%4791275                                                                        |  0.0283213|  0.5181554|  1.750054|            25|
| INTERLEUKIN-17 SIGNALING%REACTOME DATABASE ID RELEASE 65%448424                                                                           |  0.0310076|  0.4152766|  1.629998|            29|
| SIGNALING BY RAS MUTANTS%REACTOME%R-HSA-6802949.1                                                                                         |  0.0311375|  0.4587096|  1.663798|            29|
| SMOOTH MUSCLE CONTRACTION%REACTOME%R-HSA-445355.2                                                                                         |  0.0313220|  0.5225519|  1.737825|            29|
| MAPK TARGETS NUCLEAR EVENTS MEDIATED BY MAP KINASES%REACTOME DATABASE ID RELEASE 65%450282                                                |  0.0313220|  0.5219748|  1.735905|            29|
| CASPASE ACTIVATION VIA EXTRINSIC APOPTOTIC SIGNALLING PATHWAY%REACTOME%R-HSA-5357769.2                                                    |  0.0318220|  0.6251075|  1.797101|            30|
| SIGNALLING TO ERKS%REACTOME%R-HSA-187687.1                                                                                                |  0.0318220|  0.5403635|  1.757841|            30|
| SIGNALING BY NOTCH1 PEST DOMAIN MUTANTS IN CANCER%REACTOME DATABASE ID RELEASE 65%2644602                                                 |  0.0320863|  0.4275683|  1.622810|            32|
| CONSTITUTIVE SIGNALING BY NOTCH1 HD+PEST DOMAIN MUTANTS%REACTOME DATABASE ID RELEASE 65%2894862                                           |  0.0320863|  0.4275683|  1.622810|            32|
| SIGNALING BY NOTCH1 IN CANCER%REACTOME DATABASE ID RELEASE 65%2644603                                                                     |  0.0320863|  0.4275683|  1.622810|            32|
| SIGNALING BY NOTCH1 HD+PEST DOMAIN MUTANTS IN CANCER%REACTOME DATABASE ID RELEASE 65%2894858                                              |  0.0320863|  0.4275683|  1.622810|            32|
| CONSTITUTIVE SIGNALING BY NOTCH1 PEST DOMAIN MUTANTS%REACTOME%R-HSA-2644606.1                                                             |  0.0320863|  0.4275683|  1.622810|            32|
| WNT LIGAND BIOGENESIS AND TRAFFICKING%REACTOME%R-HSA-3238698.1                                                                            |  0.0325327|  0.6271626|  1.763039|            32|
| CELL DEATH SIGNALLING VIA NRAGE, NRIF AND NADE%REACTOME DATABASE ID RELEASE 65%204998                                                     |  0.0325327|  0.4017010|  1.588386|            33|
| INTRINSIC PATHWAY FOR APOPTOSIS%REACTOME%R-HSA-109606.2                                                                                   |  0.0335993|  0.4682451|  1.664039|            34|
| G-PROTEIN BETA:GAMMA SIGNALLING%REACTOME DATABASE ID RELEASE 65%397795                                                                    |  0.0336199|  0.5633308|  1.747437|            34|
| COPI-MEDIATED ANTEROGRADE TRANSPORT%REACTOME%R-HSA-6807878.1                                                                              |  0.0336199|  0.3899728|  1.577987|            35|
| NOTCH4 INTRACELLULAR DOMAIN REGULATES TRANSCRIPTION%REACTOME DATABASE ID RELEASE 65%9013695                                               |  0.0336227|  0.5904751|  1.745640|            34|
| SEMAPHORIN INTERACTIONS%REACTOME%R-HSA-373755.1                                                                                           |  0.0346881|  0.4176251|  1.605683|            36|
| TGF-BETA RECEPTOR SIGNALING ACTIVATES SMADS%REACTOME DATABASE ID RELEASE 65%2173789                                                       |  0.0354579|  0.5140958|  1.709703|            37|
| TP53 REGULATES TRANSCRIPTION OF ADDITIONAL CELL CYCLE GENES WHOSE EXACT ROLE IN THE P53 PATHWAY REMAIN UNCERTAIN%REACTOME%R-HSA-6804115.1 |  0.0355881|  0.5888055|  1.740705|            38|
| NEURODEGENERATIVE DISEASES%REACTOME DATABASE ID RELEASE 65%8863678                                                                        |  0.0355881|  0.5811884|  1.736439|            38|
| DEREGULATED CDK5 TRIGGERS MULTIPLE NEURODEGENERATIVE PATHWAYS IN ALZHEIMER'S DISEASE MODELS%REACTOME%R-HSA-8862803.2                      |  0.0355881|  0.5811884|  1.736439|            38|
| FCERI MEDIATED MAPK ACTIVATION%REACTOME%R-HSA-2871796.2                                                                                   |  0.0360587|  0.5523948|  1.725886|            39|
| DOWNSTREAM SIGNAL TRANSDUCTION%REACTOME DATABASE ID RELEASE 65%186763                                                                     |  0.0365747|  0.5307263|  1.712141|            40|
| SYNTHESIS OF PIPS AT THE PLASMA MEMBRANE%REACTOME DATABASE ID RELEASE 65%1660499                                                          |  0.0369673|  0.4358078|  1.609695|            41|
| P75 NTR RECEPTOR-MEDIATED SIGNALLING%REACTOME%R-HSA-193704.1                                                                              |  0.0370032|  0.3745126|  1.544349|            42|
| PRE-NOTCH PROCESSING IN GOLGI%REACTOME%R-HSA-1912420.2                                                                                    |  0.0383262|  0.6137732|  1.725399|            42|
| G ALPHA (Q) SIGNALLING EVENTS%REACTOME DATABASE ID RELEASE 65%416476                                                                      |  0.0400039|  0.3331424|  1.454881|            47|
| RUNX2 REGULATES OSTEOBLAST DIFFERENTIATION%REACTOME DATABASE ID RELEASE 65%8940973                                                        |  0.0400039|  0.5738155|  1.714411|            45|
| NOTCH3 INTRACELLULAR DOMAIN REGULATES TRANSCRIPTION%REACTOME%R-HSA-9013508.1                                                              |  0.0414773|  0.5615138|  1.722435|            48|
| EPHRIN SIGNALING%REACTOME DATABASE ID RELEASE 65%3928664                                                                                  |  0.0423815|  0.5976574|  1.718186|            50|
| SIGNALING BY PTK6%REACTOME%R-HSA-8848021.2                                                                                                |  0.0459362|  0.4252362|  1.570648|            56|
| SIGNALING BY NON-RECEPTOR TYROSINE KINASES%REACTOME%R-HSA-9006927.2                                                                       |  0.0459362|  0.4252362|  1.570648|            56|
| SIGNALING BY LIGAND-RESPONSIVE EGFR VARIANTS IN CANCER%REACTOME DATABASE ID RELEASE 65%5637815                                            |  0.0466690|  0.5804440|  1.696929|            57|
| CONSTITUTIVE SIGNALING BY LIGAND-RESPONSIVE EGFR CANCER VARIANTS%REACTOME DATABASE ID RELEASE 65%1236382                                  |  0.0466690|  0.5804440|  1.696929|            57|
| SIGNALING BY EGFR IN CANCER%REACTOME%R-HSA-1643713.1                                                                                      |  0.0466690|  0.5804440|  1.696929|            57|
| GAMMA CARBOXYLATION, HYPUSINE FORMATION AND ARYLSULFATASE ACTIVATION%REACTOME DATABASE ID RELEASE 65%163841                               |  0.0478851|  0.5140355|  1.672194|            59|
| INTEGRIN SIGNALING%REACTOME DATABASE ID RELEASE 65%9006921                                                                                |  0.0478954|  0.5512139|  1.690840|            60|
| INTEGRIN ALPHAIIB BETA3 SIGNALING%REACTOME DATABASE ID RELEASE 65%354192                                                                  |  0.0478954|  0.5512139|  1.690840|            60|
| ACTIVATION OF BH3-ONLY PROTEINS%REACTOME DATABASE ID RELEASE 65%114452                                                                    |  0.0482793|  0.5133241|  1.669879|            61|
| SIGNALING BY ERBB4%REACTOME DATABASE ID RELEASE 65%1236394                                                                                |  0.0484902|  0.4664183|  1.600900|            62|
| INTRACELLULAR SIGNALING BY SECOND MESSENGERS%REACTOME DATABASE ID RELEASE 65%9006925                                                      |  0.0484902|  0.2818410|  1.343812|            67|
| SIGNALING BY NOTCH1%REACTOME%R-HSA-1980143.2                                                                                              |  0.0484902|  0.3836691|  1.520735|            63|

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
| Kegg pathways in cancer                     | 0.0194400         | 0.0000000            |
| Kegg focal adhesion                         | 0.0000000         | 0.0000000            |
| Kegg ecm receptor interaction               | 0.0000000         | 0.0000000            |
| Kegg wnt signaling pathway                  | 0.0162000         | 0.0040500            |
| Kegg cell adhesion molecules cams           | 0.0216000         | 0.0138857            |
| Kegg axon guidance                          | 0.0162000         | 0.0162000            |
| Kegg mapk signaling pathway                 | 0.0462857         | 0.0162000            |
| Kegg cytokine cytokine receptor interaction | 0.0749250         | 0.0738000            |
| Kegg hypertrophic cardiomyopathy hcm        | 0.2255538         | 0.0751091            |
| Kegg erbb signaling pathway                 | 0.2325857         | 0.0777600            |

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
| KEGG\_MAPK\_SIGNALING\_PATHWAY                                 |  0.0020088|  0.3649497|  1.725032|             0|
| KEGG\_WNT\_SIGNALING\_PATHWAY                                  |  0.0020088|  0.4658496|  2.035462|             0|
| KEGG\_AXON\_GUIDANCE                                           |  0.0020088|  0.4539690|  1.968068|             0|
| KEGG\_FOCAL\_ADHESION                                          |  0.0020088|  0.5517971|  2.520479|             0|
| KEGG\_ECM\_RECEPTOR\_INTERACTION                               |  0.0020088|  0.6864207|  2.688593|             0|
| KEGG\_NEUROTROPHIN\_SIGNALING\_PATHWAY                         |  0.0020088|  0.5058686|  2.179403|             0|
| KEGG\_VASOPRESSIN\_REGULATED\_WATER\_REABSORPTION              |  0.0020088|  0.5789300|  2.013662|             0|
| KEGG\_PATHWAYS\_IN\_CANCER                                     |  0.0020088|  0.4080179|  1.979692|             0|
| KEGG\_COLORECTAL\_CANCER                                       |  0.0020088|  0.5189777|  1.995777|             0|
| KEGG\_PANCREATIC\_CANCER                                       |  0.0020088|  0.5206546|  2.049767|             0|
| KEGG\_BASAL\_CELL\_CARCINOMA                                   |  0.0020088|  0.6165568|  2.210298|             0|
| KEGG\_CHRONIC\_MYELOID\_LEUKEMIA                               |  0.0020088|  0.5024482|  1.994959|             0|
| KEGG\_SMALL\_CELL\_LUNG\_CANCER                                |  0.0020088|  0.5145797|  2.090740|             0|
| KEGG\_HYPERTROPHIC\_CARDIOMYOPATHY\_HCM                        |  0.0020088|  0.4956741|  1.918201|             0|
| KEGG\_MELANOGENESIS                                            |  0.0035058|  0.4873528|  1.969683|             1|
| KEGG\_ENDOCYTOSIS                                              |  0.0043910|  0.3778583|  1.714274|             2|
| KEGG\_HEDGEHOG\_SIGNALING\_PATHWAY                             |  0.0043910|  0.5382415|  1.921601|             2|
| KEGG\_RENAL\_CELL\_CARCINOMA                                   |  0.0043910|  0.4768893|  1.872130|             2|
| KEGG\_DILATED\_CARDIOMYOPATHY                                  |  0.0043910|  0.4723732|  1.854401|             2|
| KEGG\_ERBB\_SIGNALING\_PATHWAY                                 |  0.0054853|  0.4471755|  1.812301|             3|
| KEGG\_TGF\_BETA\_SIGNALING\_PATHWAY                            |  0.0076035|  0.4493235|  1.800929|             5|
| KEGG\_REGULATION\_OF\_ACTIN\_CYTOSKELETON                      |  0.0083115|  0.3536441|  1.615089|             6|
| KEGG\_ACUTE\_MYELOID\_LEUKEMIA                                 |  0.0083682|  0.4792117|  1.788915|             6|
| KEGG\_MTOR\_SIGNALING\_PATHWAY                                 |  0.0100191|  0.4843856|  1.777463|             8|
| KEGG\_INSULIN\_SIGNALING\_PATHWAY                              |  0.0100191|  0.3818684|  1.655494|             8|
| KEGG\_GLYCOSAMINOGLYCAN\_BIOSYNTHESIS\_HEPARAN\_SULFATE        |  0.0107463|  0.5890263|  1.845614|             9|
| KEGG\_GNRH\_SIGNALING\_PATHWAY                                 |  0.0111960|  0.4203573|  1.698914|            10|
| KEGG\_ARRHYTHMOGENIC\_RIGHT\_VENTRICULAR\_CARDIOMYOPATHY\_ARVC |  0.0135432|  0.4456180|  1.707499|            13|
| KEGG\_LYSOSOME                                                 |  0.0168963|  0.3674954|  1.573178|            18|
| KEGG\_CHEMOKINE\_SIGNALING\_PATHWAY                            |  0.0172981|  0.3639969|  1.585819|            19|
| KEGG\_PRION\_DISEASES                                          |  0.0173955|  0.5582828|  1.786570|            19|
| KEGG\_TOLL\_LIKE\_RECEPTOR\_SIGNALING\_PATHWAY                 |  0.0200170|  0.4177962|  1.640147|            23|
| KEGG\_VEGF\_SIGNALING\_PATHWAY                                 |  0.0208111|  0.4112727|  1.614538|            27|
| KEGG\_NON\_SMALL\_CELL\_LUNG\_CANCER                           |  0.0208111|  0.4341424|  1.629891|            26|
| KEGG\_FC\_EPSILON\_RI\_SIGNALING\_PATHWAY                      |  0.0243439|  0.4195961|  1.613595|            33|
| KEGG\_ENDOMETRIAL\_CANCER                                      |  0.0243439|  0.4409380|  1.633848|            33|
| KEGG\_PROGESTERONE\_MEDIATED\_OOCYTE\_MATURATION               |  0.0246613|  0.3950158|  1.590286|            35|
| KEGG\_INOSITOL\_PHOSPHATE\_METABOLISM                          |  0.0281534|  0.4388364|  1.620132|            41|
| KEGG\_APOPTOSIS                                                |  0.0282507|  0.3916074|  1.569598|            43|
| KEGG\_PHOSPHATIDYLINOSITOL\_SIGNALING\_SYSTEM                  |  0.0291794|  0.3982205|  1.572390|            46|
| KEGG\_T\_CELL\_RECEPTOR\_SIGNALING\_PATHWAY                    |  0.0308263|  0.3747279|  1.540382|            50|
| KEGG\_VASCULAR\_SMOOTH\_MUSCLE\_CONTRACTION                    |  0.0319984|  0.3676556|  1.518502|            54|
| KEGG\_PROSTATE\_CANCER                                         |  0.0364540|  0.3724379|  1.523104|            64|
| KEGG\_HEMATOPOIETIC\_CELL\_LINEAGE                             |  0.0372850|  0.4481770|  1.600058|            65|
| KEGG\_BLADDER\_CANCER                                          |  0.0391522|  0.4578896|  1.592654|            70|
| KEGG\_ADIPOCYTOKINE\_SIGNALING\_PATHWAY                        |  0.0418527|  0.3998956|  1.528435|            79|
| KEGG\_GAP\_JUNCTION                                            |  0.0433245|  0.3715653|  1.492902|            84|
