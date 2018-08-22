Pathway enrichment analysis, GEO52158
================
German Novakovskiy
August 21, 2018

Pathways
--------

``` r
x <- scan("~/ESC_RNA_seq/pathway_enrichment_analysis/ReactomePathways.tsv", what = "", sep = "\n")
x <- strsplit(x, "[ \t]+")
max.col <- max(sapply(x, length))

## specify col.names as ?read.table suggests
cn <- paste("V", 1:max.col, sep = "")
reactome_pathways <- read.table("~/ESC_RNA_seq/pathway_enrichment_analysis/ReactomePathways.tsv", fill = TRUE, 
                 col.names = cn, sep = '\t', quote = "")

#reactome_pathways[1:5, 1:10] %>% kable()
```

Analysis of GSE52158
--------------------

Load DE data:

``` r
load("~/ESC_RNA_seq/analysis_of_public_data/GSE52158/DEgenes_0h_96h_52158.Rdata")
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

``` r
probes_to_genes <- genesym.probeid %>%
  filter(probe_id %in% rownames(DEgenes_0h_96h_52158))

topProbes <- DEgenes_0h_96h_52158 %>%
  rownames_to_column("probes") %>%
  filter(probes %in% probes_to_genes$probe_id)

probes_to_genes <- probes_to_genes %>% column_to_rownames('probe_id')
symbs <- probes_to_genes[topProbes$probes,]
topProbes$Symbol <- symbs

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

head(ermineInputProbeScores, 10)
```

    ##         absolute_logFC
    ## APOA2         4.717422
    ## PRSS2         4.631256
    ## SOX17         4.429966
    ## ERBB4         4.363262
    ## GATA6         4.139878
    ## MOBP          4.098615
    ## PRDM1         3.843296
    ## COLEC12       3.815424
    ## HHEX          3.791644
    ## SLC40A1       3.776331

#### Reactome pathways

``` r
enrichmentResultReactome <- precRecall(scores = ermineInputProbeScores,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/ReactomePathways.tsv",
                               minClassSize = 15,
                               maxClassSize = 300)

enrichmentResultReactome$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

| Name                                                                                                                        | ID            |  NumProbes|  NumGenes|   RawScore|    Pval|  CorrectedPvalue|  MFPvalue|  CorrectedMFPvalue|  Multifunctionality| Same as | GeneMembers                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|:----------------------------------------------------------------------------------------------------------------------------|:--------------|----------:|---------:|----------:|-------:|----------------:|---------:|------------------:|-------------------:|:--------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Plasma lipoprotein remodeling                                                                                               | R-HSA-8963899 |         16|        16|  0.0749959|  0.0000|        0.0000000|    0.0000|          0.0000000|              0.5970| NA      | ABCG1|APOA1|APOA2|APOA4|APOA5|APOE|FURIN|LCAT|LIPG|LMF1|LMF2|LPL|MBTPS1|MBTPS2|PCSK5|PCSK6|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Extracellular matrix organization                                                                                           | R-HSA-1474244 |        219|       219|  0.0527687|  0.0000|        0.0000000|    0.0001|          0.0432500|              0.4850| NA      | ACAN|ACTN1|ADAM10|ADAM15|ADAM17|ADAM19|ADAM8|ADAM9|ADAMTS1|ADAMTS14|ADAMTS16|ADAMTS18|ADAMTS2|ADAMTS3|ADAMTS4|ADAMTS5|ADAMTS8|ADAMTS9|AGRN|BCAN|BGN|BMP1|BMP2|BMP4|BMP7|BSG|CAPN1|CAPN10|CAPN12|CAPN13|CAPN15|CAPN2|CAPN5|CAPN6|CAPN7|CAPNS1|CAPNS2|CASK|CASP3|CD151|CD44|CD47|CDH1|COL11A1|COL11A2|COL12A1|COL13A1|COL16A1|COL18A1|COL1A1|COL1A2|COL20A1|COL21A1|COL22A1|COL23A1|COL26A1|COL27A1|COL2A1|COL4A1|COL4A2|COL4A5|COL4A6|COL5A1|COL5A2|COL5A3|COL6A1|COL6A2|COL6A3|COL7A1|COL8A1|COL9A1|COL9A2|COL9A3|COLGALT1|COLGALT2|CRTAP|CTRB2|CTSB|CTSD|CTSG|CTSK|CTSL|CTSV|DAG1|DCN|DDR2|DMD|DST|EFEMP1|EFEMP2|ELN|F11R|FBLN1|FBLN2|FBLN5|FBN1|FBN2|FBN3|FGF2|FMOD|FN1|FURIN|GDF5|HAPLN1|HSPG2|HTRA1|ICAM1|ICAM3|ITGA10|ITGA11|ITGA2|ITGA2B|ITGA3|ITGA5|ITGA6|ITGA7|ITGA9|ITGAE|ITGAM|ITGAV|ITGB1|ITGB3|ITGB4|ITGB5|ITGB7|ITGB8|JAM2|JAM3|KDR|KLK2|KLKB1|LAMA1|LAMA2|LAMA3|LAMA4|LAMA5|LAMB1|LAMB2|LAMB3|LAMC1|LAMC2|LAMC3|LOX|LOXL1|LOXL2|LOXL3|LOXL4|LRP4|LTBP1|LTBP2|LTBP3|LTBP4|MADCAM1|MATN3|MATN4|MFAP1|MFAP2|MFAP3|MFAP4|MFAP5|MMP11|MMP14|MMP15|MMP16|MMP17|MMP19|MMP2|MMP24|MMP25|MMP9|NCAM1|NCSTN|NID1|NRXN1|NTN4|P3H1|P3H2|P3H3|P4HA1|P4HA2|P4HB|PCOLCE|PDGFA|PDGFB|PECAM1|PHYKPL|PLEC|PLOD1|PLOD2|PLOD3|PPIB|PRKCA|PRSS2|PSEN1|PTPRS|PXDN|SCUBE3|SDC1|SDC2|SDC3|SDC4|SERPINE1|SERPINH1|SH3PXD2A|SPARC|SPP1|TGFB1|TGFB2|TGFB3|THBS1|TIMP1|TIMP2|TLL2|TMPRSS6|TNC|TNXB|TPSAB1|TTR|VCAN|                                                                                                                                                                                                                                                                                              |
| Plasma lipoprotein assembly, remodeling, and clearance                                                                      | R-HSA-174824  |         44|        44|  0.0442264|  0.0002|        0.0432500|    0.0002|          0.0576667|              0.7780| NA      | ABCA1|ABCG1|AMN|AP2A1|AP2A2|AP2B1|AP2M1|AP2S1|APOA1|APOA2|APOA4|APOA5|APOC1|APOE|BMP1|CES3|CLTA|CLTC|FURIN|HDLBP|LCAT|LDLR|LIPA|LIPG|LMF1|LMF2|LPL|LSR|MBTPS1|MBTPS2|MYLIP|NCEH1|NPC1|NPC2|P4HB|PCSK5|PCSK6|PCSK9|PRKACA|PRKACB|SCARB1|SOAT1|VLDLR|ZDHHC8|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| Asparagine N-linked glycosylation                                                                                           | R-HSA-446203  |        259|       259|  0.0538784|  0.0001|        0.0288333|    0.0013|          0.2811250|              0.2500| NA      | ACTR10|ACTR1A|ALG1|ALG10|ALG12|ALG13|ALG14|ALG2|ALG3|ALG5|ALG6|ALG8|ALG9|AMFR|ANK1|ANK2|ANK3|ANKRD28|ARCN1|ARF1|ARF3|ARF4|ARF5|ARFGAP1|ARFGAP2|ARFGAP3|ASGR1|ASGR2|B4GALT1|B4GALT2|B4GALT3|B4GALT4|B4GALT5|B4GALT6|BET1|BET1L|CALR|CANX|CAPZA1|CAPZA2|CAPZB|CD55|CD59|CHST10|CHST8|CMAS|CNIH1|CNIH2|CNIH3|COG1|COG2|COG3|COG4|COG5|COG6|COG7|COG8|COL7A1|COPA|COPB1|COPB2|COPE|COPG1|COPG2|COPZ1|COPZ2|CSNK1D|CTSA|CTSC|CTSZ|DAD1|DCTN1|DCTN2|DCTN3|DCTN4|DCTN5|DCTN6|DDOST|DERL1|DERL2|DHDDS|DOLK|DOLPP1|DPAGT1|DPM1|DPM2|DPM3|DYNC1H1|DYNC1I1|DYNC1I2|DYNC1LI1|DYNC1LI2|DYNLL1|DYNLL2|EDEM1|EDEM2|EDEM3|ENGASE|F8|FOLR1|FPGT|FUCA1|FUK|FUOM|FUT3|FUT8|GANAB|GBF1|GFPT1|GFPT2|GMDS|GMPPA|GNE|GNPNAT1|GOLGA2|GOLGB1|GORASP1|GOSR1|GOSR2|KDELR1|KDELR2|KDELR3|LMAN1|LMAN1L|LMAN2|LMAN2L|MAGT1|MAN1A1|MAN1A2|MAN1B1|MAN1C1|MAN2A1|MAN2A2|MANEA|MARCH6|MCFD2|MGAT1|MGAT2|MGAT3|MGAT4A|MGAT4B|MGAT4C|MGAT5|MIA3|MLEC|MOGS|MPI|MVD|NAGK|NANP|NANS|NAPA|NAPB|NAPG|NEU1|NEU3|NEU4|NGLY1|NPL|NSF|NUDT14|NUS1|OS9|PDIA3|PGM3|PMM1|PMM2|PPP6C|PPP6R1|PPP6R3|PREB|PRKCSH|PSMC1|RAB1A|RAB1B|RAD23B|RENBP|RFT1|RNF103|RNF139|RNF185|RPN1|RPN2|RPS27A|SAR1B|SCFD1|SEC13|SEC16A|SEC16B|SEC22A|SEC22B|SEC22C|SEC23A|SEC23IP|SEC24A|SEC24B|SEC24C|SEC24D|SEC31A|SEL1L|SLC17A5|SLC35A1|SLC35C1|SPTAN1|SPTB|SPTBN1|SPTBN2|SPTBN4|SRD5A3|ST3GAL1|ST3GAL2|ST3GAL3|ST3GAL4|ST3GAL5|ST3GAL6|ST6GAL1|ST6GAL2|ST6GALNAC2|ST6GALNAC3|ST6GALNAC4|ST6GALNAC5|ST6GALNAC6|ST8SIA2|ST8SIA3|ST8SIA4|ST8SIA5|STT3A|STX17|STX5|TBC1D20|TFG|TGFA|TMED10|TMED2|TMED3|TMED9|TMEM115|TRAPPC1|TRAPPC10|TRAPPC2|TRAPPC2L|TRAPPC3|TRAPPC5|TRAPPC6A|TRAPPC6B|TRAPPC9|TRIM13|TSTA3|TUSC3|UAP1|UBA52|UBB|UBC|UBXN1|UGGT1|UGGT2|USO1|VCP|YKT6| |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | R-HSA-381426  |         92|        92|  0.0276797|  0.0012|        0.1038000|    0.0013|          0.2249000|              0.0249| NA      | ADAM10|AFP|ANO8|APLP2|APOA1|APOA2|APOA5|APOE|APOL1|APP|BMP4|BPIFB2|C3|CALU|CDH2|CHGB|CHRDL1|CKAP4|CP|CSF1|CST3|CTSG|CYR61|DNAJC3|FAM20C|FBN1|FN1|FSTL1|FSTL3|FUCA2|GAS6|GOLM1|GPC3|HRC|HSP90B1|IGF1|IGFALS|IGFBP2|IGFBP3|IGFBP4|IGFBP5|IGFBP6|IGFBP7|IL6|KLK1|KLK2|KLK3|KTN1|LAMB1|LAMB2|LAMC1|LGALS1|LTBP1|MATN3|MBTPS1|MELTF|MEN1|MEPE|MFGE8|MGAT4A|MIA3|MMP2|MXRA8|NOTUM|NUCB1|P4HB|PAPPA|PAPPA2|PCSK9|PDIA6|PENK|PNPLA2|PRKCSH|PRSS23|QSOX1|RCN1|SCG2|SCG3|SDC2|SHISA5|SPARCL1|SPP1|STC2|TF|TGOLN2|TIMP1|TMEM132A|TNC|VCAN|VGF|VWA1|WFS1|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| PPARA activates gene expression                                                                                             | R-HSA-1989781 |         31|        31|  0.0416342|  0.0009|        0.0973125|    0.0016|          0.2306667|              0.2890| NA      | ABCA1|ACADM|ACOX1|ACSL1|ALAS1|ANGPTL4|ANKRD1|APOA1|APOA2|APOA5|CPT1A|CPT2|CYP4A11|FDFT1|FHL2|G0S2|GLIPR1|GRHL1|HMGCR|HMGCS1|ME1|NPAS2|PEX11A|PLIN2|PPARA|RGL1|SLC27A1|TIAM2|TNFRSF21|TRIB3|TXNRD1|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| Post-translational protein phosphorylation                                                                                  | R-HSA-8957275 |         81|        81|  0.0271749|  0.0016|        0.1258182|    0.0017|          0.2100714|              0.0340| NA      | ADAM10|AFP|ANO8|APLP2|APOA1|APOA2|APOA5|APOE|APOL1|APP|BMP4|BPIFB2|C3|CALU|CDH2|CHGB|CHRDL1|CKAP4|CP|CSF1|CST3|CYR61|DNAJC3|FAM20C|FBN1|FN1|FSTL1|FSTL3|FUCA2|GAS6|GOLM1|GPC3|HRC|HSP90B1|IGFBP3|IGFBP4|IGFBP5|IGFBP7|IL6|KTN1|LAMB1|LAMB2|LAMC1|LGALS1|LTBP1|MATN3|MBTPS1|MELTF|MEN1|MEPE|MFGE8|MGAT4A|MIA3|MXRA8|NOTUM|NUCB1|P4HB|PCSK9|PDIA6|PENK|PNPLA2|PRKCSH|PRSS23|QSOX1|RCN1|SCG2|SCG3|SDC2|SHISA5|SPARCL1|SPP1|STC2|TF|TGOLN2|TIMP1|TMEM132A|TNC|VCAN|VGF|VWA1|WFS1|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| Regulation of lipid metabolism by Peroxisome proliferator-activated receptor alpha (PPARalpha)                              | R-HSA-400206  |         48|        48|  0.0331364|  0.0005|        0.0720833|    0.0021|          0.2270625|              0.3860| NA      | ABCA1|ACADM|ACOX1|ACSL1|ALAS1|ANGPTL4|ANKRD1|APOA1|APOA2|APOA5|CHD9|CPT1A|CPT2|CREBBP|CYP4A11|FDFT1|FHL2|G0S2|GLIPR1|GRHL1|HDAC3|HELZ2|HMGCR|HMGCS1|ME1|MED1|NCOA1|NCOA2|NCOA6|NCOR1|NCOR2|NPAS2|PEX11A|PLIN2|PPARA|RGL1|RXRA|SIN3A|SIN3B|SLC27A1|SMARCD3|TBL1X|TBL1XR1|TGS1|TIAM2|TNFRSF21|TRIB3|TXNRD1|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Retinoid metabolism and transport                                                                                           | R-HSA-975634  |         33|        33|  0.0391143|  0.0003|        0.0519000|    0.0024|          0.2306667|              0.6900| NA      | AGRN|AKR1B10|AKR1C1|AKR1C3|APOA1|APOA2|APOA4|APOE|APOM|CLPS|GPC1|GPC2|GPC3|GPC4|GPC6|HSPG2|LDLR|LPL|LRAT|LRP1|LRP10|LRP12|LRP2|LRP8|PLB1|RBP1|RDH11|RETSAT|SDC1|SDC2|SDC3|SDC4|TTR|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Metabolism of fat-soluble vitamins                                                                                          | R-HSA-6806667 |         36|        36|  0.0363624|  0.0008|        0.0988571|    0.0026|          0.2249000|              0.7080| NA      | AGRN|AKR1B10|AKR1C1|AKR1C3|APOA1|APOA2|APOA4|APOE|APOM|CLPS|GPC1|GPC2|GPC3|GPC4|GPC6|HSPG2|LDLR|LPL|LRAT|LRP1|LRP10|LRP12|LRP2|LRP8|PLB1|RBP1|RDH11|RETSAT|SDC1|SDC2|SDC3|SDC4|TTR|UBIAD1|VKORC1|VKORC1L1|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |

``` r
enrichmentResultReactome$results %>% 
  dplyr::select(Name, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  head(10) %>% 
  kable(align = "l", col.names = c("Pathway", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| Pathway                                                                                                                     | Corrected p-value | Corrected MF p-value |
|:----------------------------------------------------------------------------------------------------------------------------|:------------------|:---------------------|
| Plasma lipoprotein remodeling                                                                                               | 0.0000000         | 0.0000000            |
| Extracellular matrix organization                                                                                           | 0.0000000         | 0.0432500            |
| Plasma lipoprotein assembly, remodeling, and clearance                                                                      | 0.0432500         | 0.0576667            |
| Post-translational protein phosphorylation                                                                                  | 0.1258182         | 0.2100714            |
| Metabolism of fat-soluble vitamins                                                                                          | 0.0988571         | 0.2249000            |
| Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) | 0.1038000         | 0.2249000            |
| Regulation of lipid metabolism by Peroxisome proliferator-activated receptor alpha (PPARalpha)                              | 0.0720833         | 0.2270625            |
| Retinoid metabolism and transport                                                                                           | 0.0519000         | 0.2306667            |
| PPARA activates gene expression                                                                                             | 0.0973125         | 0.2306667            |
| Activation of Matrix Metalloproteinases                                                                                     | 0.1791786         | 0.2437727            |
| \#\#\#\# KEGG pathways                                                                                                      |                   |                      |

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

enrichmentResultKEGG$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

| Name                                   | ID                                          |  NumProbes|  NumGenes|   RawScore|    Pval|  CorrectedPvalue|  MFPvalue|  CorrectedMFPvalue|  Multifunctionality| Same as | GeneMembers                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|:---------------------------------------|:--------------------------------------------|----------:|---------:|----------:|-------:|----------------:|---------:|------------------:|-------------------:|:--------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Kegg pathways in cancer                | KEGG\_PATHWAYS\_IN\_CANCER                  |        275|       275|  0.0802344|  0.3567|        0.8464970|    0.0000|          0.0000000|              0.9950| NA      | ABL1|AKT1|AKT2|AKT3|APC|APC2|APPL1|AR|ARAF|ARNT|ARNT2|AXIN1|AXIN2|BAD|BAX|BCL2|BCL2L1|BCR|BID|BIRC2|BIRC5|BMP2|BMP4|BRAF|BRCA2|CASP3|CASP8|CASP9|CBL|CBLB|CBLC|CCDC6|CCNA1|CCND1|CCNE1|CCNE2|CDC42|CDH1|CDK2|CDK4|CDK6|CDKN1A|CDKN1B|CDKN2B|CEBPA|CHUK|CKS1B|COL4A1|COL4A2|COL4A6|CREBBP|CRK|CRKL|CSF3R|CTBP1|CTBP2|CTNNA1|CTNNA2|CTNNA3|CTNNB1|CUL2|CYCS|DAPK1|DVL1|DVL2|DVL3|E2F1|E2F2|E2F3|EGF|EGFR|EGLN1|EGLN3|EP300|EPAS1|ERBB2|ETS1|FADD|FAS|FGF1|FGF11|FGF12|FGF13|FGF16|FGF17|FGF18|FGF19|FGF2|FGF22|FGF4|FGF6|FGF8|FGFR1|FGFR2|FGFR3|FH|FLT3LG|FN1|FOS|FOXO1|FZD1|FZD10|FZD2|FZD3|FZD4|FZD5|FZD6|FZD7|GLI1|GLI2|GLI3|GRB2|GSK3B|GSTP1|HDAC1|HDAC2|HGF|HIF1A|HRAS|HSP90AA1|HSP90AB1|HSP90B1|IGF1|IGF1R|IKBKB|IKBKG|IL6|ITGA2|ITGA2B|ITGA3|ITGA6|ITGAV|ITGB1|JAK1|JUN|JUP|KIT|KITLG|KLK3|KRAS|LAMA1|LAMA2|LAMA3|LAMA4|LAMA5|LAMB1|LAMB2|LAMB3|LAMC1|LAMC2|LAMC3|LEF1|MAP2K1|MAP2K2|MAPK1|MAPK10|MAPK3|MAPK8|MAPK9|MAX|MDM2|MECOM|MET|MITF|MLH1|MMP2|MMP9|MSH2|MSH3|MSH6|MTOR|MYC|NCOA4|NFKB1|NFKBIA|NKX3-1|NOS2|NRAS|PAX8|PDGFA|PDGFB|PDGFRA|PDGFRB|PGF|PIAS1|PIAS2|PIAS3|PIAS4|PIK3CA|PIK3CB|PIK3CD|PIK3R1|PIK3R2|PIK3R3|PIK3R5|PLCG1|PLCG2|PLD1|PML|PPARD|PRKCA|PRKCB|PRKCG|PTCH1|PTCH2|PTEN|PTK2|RAC1|RAC2|RAC3|RAD51|RAF1|RALA|RALB|RALBP1|RALGDS|RARA|RARB|RASSF1|RASSF5|RB1|RBX1|RELA|RET|RHOA|RUNX1T1|RXRA|RXRB|SHH|SKP2|SLC2A1|SMAD2|SMAD3|SMAD4|SMO|SOS1|SOS2|STAT1|STAT3|STAT5A|STAT5B|STK36|STK4|SUFU|TCF7|TCF7L1|TCF7L2|TFG|TGFA|TGFB1|TGFB2|TGFB3|TGFBR1|TGFBR2|TP53|TPM3|TPR|TRAF2|TRAF3|TRAF4|TRAF5|TRAF6|VEGFA|VEGFB|VHL|WNT10B|WNT16|WNT5A|WNT5B|WNT6|WNT7A|WNT7B|WNT8A|WNT9A|XIAP| |
| Kegg focal adhesion                    | KEGG\_FOCAL\_ADHESION                       |        175|       175|  0.0583071|  0.0404|        0.5839636|    0.0002|          0.0159000|              0.9840| NA      | ACTB|ACTN1|ACTN2|ACTN3|ACTN4|AKT1|AKT2|AKT3|ARHGAP35|ARHGAP5|BAD|BCAR1|BCL2|BIRC2|BRAF|CAPN2|CAV1|CAV2|CAV3|CCND1|CCND2|CCND3|CDC42|CHAD|COL11A1|COL11A2|COL1A1|COL1A2|COL2A1|COL4A1|COL4A2|COL4A6|COL5A1|COL5A2|COL5A3|COL6A1|COL6A2|COL6A3|CRK|CRKL|CTNNB1|DIAPH1|DOCK1|EGF|EGFR|ELK1|ERBB2|FLNA|FLNB|FLNC|FLT1|FN1|FYN|GRB2|GSK3B|HGF|HRAS|IGF1|IGF1R|ILK|ITGA10|ITGA11|ITGA2|ITGA2B|ITGA3|ITGA5|ITGA6|ITGA7|ITGA9|ITGAV|ITGB1|ITGB3|ITGB4|ITGB5|ITGB7|ITGB8|JUN|KDR|LAMA1|LAMA2|LAMA3|LAMA4|LAMA5|LAMB1|LAMB2|LAMB3|LAMC1|LAMC2|LAMC3|MAP2K1|MAPK1|MAPK10|MAPK3|MAPK8|MAPK9|MET|MYL10|MYL12A|MYL12B|MYL5|MYL7|MYL9|MYLK|MYLK2|MYLPF|PAK1|PAK2|PAK3|PAK4|PAK6|PARVA|PARVB|PARVG|PDGFA|PDGFB|PDGFC|PDGFD|PDGFRA|PDGFRB|PDPK1|PGF|PIK3CA|PIK3CB|PIK3CD|PIK3R1|PIK3R2|PIK3R3|PIK3R5|PIP5K1C|PPP1CA|PPP1CB|PPP1CC|PPP1R12A|PRKCA|PRKCB|PRKCG|PTEN|PTK2|PXN|RAC1|RAC2|RAC3|RAF1|RAP1A|RAP1B|RAPGEF1|RELN|RHOA|ROCK1|ROCK2|SHC1|SHC2|SHC3|SOS1|SOS2|SPP1|SRC|THBS1|THBS2|THBS3|THBS4|TLN1|TLN2|TNC|TNXB|VASP|VAV1|VAV2|VAV3|VCL|VEGFA|VEGFB|VWF|XIAP|ZYX|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| Kegg ppar signaling pathway            | KEGG\_PPAR\_SIGNALING\_PATHWAY              |         48|        48|  0.0444526|  0.0018|        0.1431000|    0.0011|          0.0583000|              0.2470| NA      | ACAA1|ACADM|ACOX1|ACOX2|ACOX3|ACSL1|ACSL3|ACSL4|ACSL5|ANGPTL4|APOA1|APOA2|APOA5|CPT1A|CPT1C|CPT2|CYP27A1|CYP4A11|DBI|EHHADH|FABP3|FABP5|FABP6|FABP7|FADS2|GK|ILK|LPL|ME1|NR1H3|PCK2|PDPK1|PLIN1|PLTP|PPARA|PPARD|RXRA|RXRB|SCD|SCD5|SCP2|SLC27A1|SLC27A2|SLC27A4|SLC27A5|SLC27A6|SORBS1|UBC|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| Kegg wnt signaling pathway             | KEGG\_WNT\_SIGNALING\_PATHWAY               |        129|       129|  0.0490600|  0.0130|        0.5167500|    0.0013|          0.0516750|              0.4140| NA      | APC|APC2|AXIN1|AXIN2|BTRC|CACYBP|CAMK2A|CAMK2B|CAMK2D|CAMK2G|CCND1|CCND2|CCND3|CER1|CHD8|CREBBP|CSNK1A1|CSNK1E|CSNK2A1|CSNK2A2|CSNK2B|CTBP1|CTBP2|CTNNB1|CTNNBIP1|CUL1|CXXC4|DAAM1|DAAM2|DKK1|DKK4|DVL1|DVL2|DVL3|EP300|FBXW11|FOSL1|FRAT1|FRAT2|FZD1|FZD10|FZD2|FZD3|FZD4|FZD5|FZD6|FZD7|GSK3B|JUN|LEF1|LRP5|LRP6|MAP3K7|MAPK10|MAPK8|MAPK9|MYC|NFAT5|NFATC1|NFATC2|NFATC3|NFATC4|NKD1|NKD2|NLK|PLCB1|PLCB2|PLCB3|PLCB4|PORCN|PPARD|PPP2CA|PPP2CB|PPP2R1A|PPP2R1B|PPP2R5A|PPP2R5B|PPP2R5C|PPP2R5D|PPP2R5E|PPP3CA|PPP3CB|PPP3CC|PPP3R1|PRICKLE1|PRICKLE2|PRKACA|PRKACB|PRKCA|PRKCB|PRKCG|PRKX|PSEN1|RAC1|RAC2|RAC3|RBX1|RHOA|ROCK1|ROCK2|RUVBL1|SENP2|SFRP1|SFRP2|SFRP5|SIAH1|SKP1|SMAD2|SMAD3|SMAD4|SOX17|TBL1X|TBL1XR1|TCF7|TCF7L1|TCF7L2|TP53|VANGL1|VANGL2|WIF1|WNT10B|WNT16|WNT5A|WNT5B|WNT6|WNT7A|WNT7B|WNT8A|WNT9A|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Kegg erbb signaling pathway            | KEGG\_ERBB\_SIGNALING\_PATHWAY              |         80|        80|  0.0300770|  0.0593|        0.6734786|    0.0020|          0.0636000|              0.9250| NA      | ABL1|ABL2|AKT1|AKT2|AKT3|ARAF|BAD|BRAF|CAMK2A|CAMK2B|CAMK2D|CAMK2G|CBL|CBLB|CBLC|CDKN1A|CDKN1B|CRK|CRKL|EGF|EGFR|EIF4EBP1|ELK1|ERBB2|ERBB3|ERBB4|GAB1|GRB2|GSK3B|HBEGF|HRAS|JUN|KRAS|MAP2K1|MAP2K2|MAP2K4|MAP2K7|MAPK1|MAPK10|MAPK3|MAPK8|MAPK9|MTOR|MYC|NCK1|NCK2|NRAS|NRG1|NRG2|NRG3|PAK1|PAK2|PAK3|PAK4|PAK6|PIK3CA|PIK3CB|PIK3CD|PIK3R1|PIK3R2|PIK3R3|PIK3R5|PLCG1|PLCG2|PRKCA|PRKCB|PRKCG|PTK2|RAF1|RPS6KB1|RPS6KB2|SHC1|SHC2|SHC3|SOS1|SOS2|SRC|STAT5A|STAT5B|TGFA|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Kegg neurotrophin signaling pathway    | KEGG\_NEUROTROPHIN\_SIGNALING\_PATHWAY      |        113|       113|  0.0357791|  0.1125|        0.8130682|    0.0024|          0.0636000|              0.9890| NA      | ABL1|AKT1|AKT2|AKT3|ARHGDIA|ARHGDIB|ATF4|BAD|BAX|BCL2|BDNF|BRAF|CALM1|CALM3|CAMK2A|CAMK2B|CAMK2D|CAMK2G|CAMK4|CDC42|CRK|CRKL|CSK|FOXO3|FRS2|GAB1|GRB2|GSK3B|HRAS|IKBKB|IRAK1|IRAK2|IRAK4|IRS1|IRS2|IRS4|JUN|KIDINS220|KRAS|MAGED1|MAP2K1|MAP2K2|MAP2K5|MAP2K7|MAP3K1|MAP3K3|MAP3K5|MAPK1|MAPK10|MAPK11|MAPK12|MAPK13|MAPK14|MAPK3|MAPK7|MAPK8|MAPK9|MAPKAPK2|NFKB1|NFKBIA|NFKBIB|NFKBIE|NGF|NRAS|NTF3|NTRK2|NTRK3|PDPK1|PIK3CA|PIK3CB|PIK3CD|PIK3R1|PIK3R2|PIK3R3|PIK3R5|PLCG1|PLCG2|PRDM4|PRKCD|PSEN1|PTPN11|RAC1|RAF1|RAP1A|RAP1B|RAPGEF1|RELA|RHOA|RIPK2|RPS6KA1|RPS6KA2|RPS6KA3|RPS6KA4|RPS6KA5|RPS6KA6|SH2B1|SH2B2|SH2B3|SHC1|SHC2|SHC3|SORT1|SOS1|SOS2|TP53|TRAF6|YWHAB|YWHAE|YWHAG|YWHAH|YWHAQ|YWHAZ|ZNF274|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Kegg axon guidance                     | KEGG\_AXON\_GUIDANCE                        |        113|       113|  0.0541916|  0.0003|        0.0477000|    0.0025|          0.0567857|              0.0376| NA      | ABL1|ABLIM1|ABLIM2|ARHGEF12|CDC42|CDK5|CFL1|CFL2|CXCL12|CXCR4|DPYSL2|EFNA1|EFNA2|EFNA3|EFNA4|EFNA5|EFNB1|EFNB2|EFNB3|EPHA1|EPHA2|EPHA4|EPHA7|EPHA8|EPHB1|EPHB2|EPHB3|EPHB4|EPHB6|FES|FYN|GNAI1|GNAI2|GNAI3|GSK3B|HRAS|ITGB1|KRAS|L1CAM|LIMK1|LIMK2|MAPK1|MAPK3|MET|NCK1|NCK2|NFAT5|NFATC1|NFATC2|NFATC3|NFATC4|NGEF|NRAS|NRP1|NTN1|NTN4|NTNG1|PAK1|PAK2|PAK3|PAK4|PAK6|PLXNA1|PLXNA2|PLXNA3|PLXNB1|PLXNB2|PLXNB3|PLXNC1|PPP3CA|PPP3CB|PPP3CC|PPP3R1|PTK2|RAC1|RAC2|RAC3|RASA1|RGS3|RHOA|RHOD|RND1|ROBO1|ROBO2|ROBO3|ROCK1|ROCK2|SEMA3A|SEMA3E|SEMA3F|SEMA3G|SEMA4A|SEMA4B|SEMA4C|SEMA4D|SEMA4F|SEMA4G|SEMA5A|SEMA5B|SEMA6A|SEMA6B|SEMA6C|SEMA6D|SEMA7A|SLIT1|SLIT2|SLIT3|SRGAP1|SRGAP3|UNC5A|UNC5B|UNC5C|UNC5D|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Kegg chemokine signaling pathway       | KEGG\_CHEMOKINE\_SIGNALING\_PATHWAY         |        129|       129|  0.0414801|  0.1466|        0.8965154|    0.0043|          0.0854625|              0.7420| NA      | ADCY1|ADCY2|ADCY3|ADCY4|ADCY5|ADCY6|ADCY7|ADCY9|AKT1|AKT2|AKT3|ARRB1|ARRB2|BCAR1|BRAF|CCL13|CCL2|CCL22|CCL26|CCL28|CCR5|CDC42|CHUK|CRK|CRKL|CSK|CX3CL1|CXCL1|CXCL12|CXCL14|CXCL3|CXCL5|CXCL6|CXCR1|CXCR2|CXCR3|CXCR4|CXCR5|DOCK2|ELMO1|FOXO3|GNAI1|GNAI2|GNAI3|GNB1|GNB2|GNB3|GNB4|GNB5|GNG11|GNG12|GNG2|GNG3|GNG4|GNG5|GNG7|GNGT2|GRB2|GRK1|GRK4|GRK5|GRK6|GSK3A|GSK3B|HRAS|IKBKB|IKBKG|JAK2|JAK3|KRAS|LYN|MAP2K1|MAPK1|MAPK3|NFKB1|NFKBIA|NFKBIB|NRAS|PAK1|PARD3|PF4|PIK3CA|PIK3CB|PIK3CD|PIK3R1|PIK3R2|PIK3R3|PIK3R5|PLCB1|PLCB2|PLCB3|PLCB4|PREX1|PRKACA|PRKACB|PRKCB|PRKCD|PRKCZ|PRKX|PTK2|PTK2B|PXN|RAC1|RAC2|RAF1|RAP1A|RAP1B|RASGRP2|RELA|RHOA|ROCK1|ROCK2|SHC1|SHC2|SHC3|SOS1|SOS2|STAT1|STAT2|STAT3|STAT5B|TIAM1|TIAM2|VAV1|VAV2|VAV3|WAS|WASL|XCL1|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Kegg t cell receptor signaling pathway | KEGG\_T\_CELL\_RECEPTOR\_SIGNALING\_PATHWAY |         85|        85|  0.0263703|  0.3490|        0.8670469|    0.0047|          0.0830333|              0.9350| NA      | AKT1|AKT2|AKT3|BCL10|CARD11|CBL|CBLB|CBLC|CD4|CD40LG|CD8B|CDC42|CDK4|CHUK|CSF2|CTLA4|DLG1|FOS|FYN|GRB2|GSK3B|HRAS|IKBKB|IKBKG|JUN|KRAS|LAT|LCK|MALT1|MAP2K1|MAP2K2|MAP2K7|MAP3K7|MAP3K8|MAPK1|MAPK11|MAPK12|MAPK13|MAPK14|MAPK3|MAPK9|NCK1|NCK2|NFAT5|NFATC1|NFATC2|NFATC3|NFATC4|NFKB1|NFKBIA|NFKBIB|NFKBIE|NRAS|PAK1|PAK2|PAK3|PAK4|PAK6|PDCD1|PDPK1|PIK3CA|PIK3CB|PIK3CD|PIK3R1|PIK3R2|PIK3R3|PIK3R5|PLCG1|PPP3CA|PPP3CB|PPP3CC|PPP3R1|PRKCQ|PTPN6|RAF1|RASGRP1|RELA|RHOA|SOS1|SOS2|TEC|VAV1|VAV2|VAV3|ZAP70|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| Kegg vegf signaling pathway            | KEGG\_VEGF\_SIGNALING\_PATHWAY              |         59|        59|  0.0225815|  0.0865|        0.7640833|    0.0083|          0.1319700|              0.9190| NA      | AKT1|AKT2|AKT3|BAD|CASP9|CDC42|HRAS|HSPB1|KDR|KRAS|MAP2K1|MAP2K2|MAPK1|MAPK11|MAPK12|MAPK13|MAPK14|MAPK3|MAPKAPK2|MAPKAPK3|NFAT5|NFATC1|NFATC2|NFATC3|NFATC4|NOS3|NRAS|PIK3CA|PIK3CB|PIK3CD|PIK3R1|PIK3R2|PIK3R3|PIK3R5|PLA2G12A|PLA2G2A|PLA2G3|PLA2G6|PLCG1|PLCG2|PPP3CA|PPP3CB|PPP3CC|PPP3R1|PRKCA|PRKCB|PRKCG|PTK2|PXN|RAC1|RAC2|RAC3|RAF1|SH2D2A|SHC2|SPHK1|SPHK2|SRC|VEGFA|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |

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
| Kegg pathways in cancer                | 0.8464970         | 0.0000000            |
| Kegg focal adhesion                    | 0.5839636         | 0.0159000            |
| Kegg wnt signaling pathway             | 0.5167500         | 0.0516750            |
| Kegg axon guidance                     | 0.0477000         | 0.0567857            |
| Kegg ppar signaling pathway            | 0.1431000         | 0.0583000            |
| Kegg erbb signaling pathway            | 0.6734786         | 0.0636000            |
| Kegg neurotrophin signaling pathway    | 0.8130682         | 0.0636000            |
| Kegg t cell receptor signaling pathway | 0.8670469         | 0.0830333            |
| Kegg chemokine signaling pathway       | 0.8965154         | 0.0854625            |
| Kegg dilated cardiomyopathy            | 0.9285600         | 0.1219000            |

#### Wikipathways

``` r
#converting symbols to EntrezID
geneSymbols <- rownames(ermineInputProbeScores)
geneEntrez <- bitr(geneSymbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
ermineInputGeneScoresRPA <- ermineInputProbeScores %>%
  rownames_to_column("SYMBOL") %>%
  filter(SYMBOL %in% geneEntrez$SYMBOL)

inputGSEA <- merge(ermineInputGeneScoresRPA, geneEntrez, sort=FALSE)

inputGSEA <- inputGSEA %>%
  dplyr::select(absolute_logFC, ENTREZID) %>%
  column_to_rownames('ENTREZID')
```

``` r
#entrezScores <- as.data.frame(test2)
enrichmentResultWiki <- precRecall(scores = inputGSEA,
                               scoreColumn = 1,
                               bigIsBetter = TRUE,
                               aspects = "B",
                               iterations = 10000,
                               geneSetDescription = NULL,
                               customGeneSets = "~/ESC_RNA_seq/pathway_enrichment_analysis/Wikipathways.gmt",
                               minClassSize = 15,
                               maxClassSize = 300)

enrichmentResultWiki$results %>% arrange(MFPvalue) %>% head(10) %>% kable()
```

| Name                                                 | ID                                                                                                             |  NumProbes|  NumGenes|   RawScore|    Pval|  CorrectedPvalue|  MFPvalue|  CorrectedMFPvalue|  Multifunctionality| Same as | GeneMembers                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
|:-----------------------------------------------------|:---------------------------------------------------------------------------------------------------------------|----------:|---------:|----------:|-------:|----------------:|---------:|------------------:|-------------------:|:--------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| <http://www.wikipathways.org/instance/WP2858_r94911> | Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                            |        122|       122|  0.0709612|  0.0000|        0.0000000|    0.0000|          0.0000000|             0.00668| NA      | 10013|1004|1006|10253|10257|10486|10634|11177|114815|114822|1487|1496|1499|1501|151613|152789|166|170302|1756|1896|1952|2263|2274|22871|2300|23133|23229|2534|25987|2627|27324|27352|2737|3148|3170|345557|3728|3801|3975|398|4061|4089|4204|440193|4609|4771|4772|4920|4926|4982|5010|5080|50937|51043|51222|51422|5150|51592|51701|51762|5236|5269|5292|5297|5362|5420|54207|5452|54521|5467|54806|54828|54880|54898|54903|56899|56913|56963|57154|57462|57476|57621|5783|6016|6238|6347|6386|6469|6498|652|657|6622|6657|6781|6785|6907|6929|7022|7226|7343|7593|7855|7903|79658|79731|79776|81576|8322|83439|8345|83933|8473|8495|8675|8820|8835|8848|9079|91653|9687|9924|9935|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| <http://www.wikipathways.org/instance/WP2853_r88152> | Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            |        129|       129|  0.0633946|  0.0000|        0.0000000|    0.0000|          0.0000000|             0.02000| NA      | 10128|10153|10451|10459|10973|1112|11169|11190|115825|121536|123169|131405|140597|1488|1499|154091|1601|1789|1844|1846|1847|1994|2005|2131|2146|2186|22823|22943|23007|2308|23168|23181|23286|23373|23528|23576|23682|2626|2627|26610|26986|27244|27245|27324|2736|2908|3087|31|3169|3170|3227|324|3251|351|3720|3915|3975|4087|4088|4089|4830|4838|4851|5015|5087|5090|51176|51366|51460|51701|54623|54799|54892|55809|55832|56946|57198|57669|57690|58499|5916|6001|63035|63978|6422|64321|6478|64864|655|657|6657|6671|6772|6875|6877|6925|6928|6932|7040|708|7080|7478|7547|7709|79035|79577|79668|79776|79923|79977|80155|80312|8320|83439|83595|83881|84295|8450|85363|8554|8928|9044|9241|9350|9425|9573|9646|9682|9760|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| <http://www.wikipathways.org/instance/WP2406_r89157> | Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  |         41|        41|  0.0636748|  0.0000|        0.0000000|    0.0001|          0.0096000|             0.47900| NA      | 132625|140885|1432|145873|1482|2247|22943|2626|2932|3084|3170|3479|3670|3791|3815|389421|4624|4684|4838|4851|4920|5080|5156|55897|6331|64321|649|652|6656|6657|6862|70|7040|7070|7137|7139|7852|79727|79923|83881|9241|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| <http://www.wikipathways.org/instance/WP430_r96008>  | Statin Pathway%WikiPathways\_20180810%WP430%Homo sapiens                                                       |         20|        20|  0.0674888|  0.0005|        0.0288000|    0.0007|          0.0504000|             0.45700| NA      | 116519|19|2222|3156|335|336|337|341|348|3931|3949|4023|4035|5360|64714|6646|6713|84532|8694|949|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| <http://www.wikipathways.org/instance/WP2882_r97489> | Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                      |        229|       229|  0.0783084|  0.0000|        0.0000000|    0.0032|          0.1843200|             0.65000| NA      | 10057|10062|1019|10252|10257|1027|1028|10486|10499|10560|10602|1066|10728|10769|10891|10998|11057|11140|11182|11214|11309|114134|115584|116519|121512|125206|1374|1376|154091|1544|1545|1548|1555|1576|1579|160728|1622|1728|177|1831|1836|1839|1956|1958|196|1962|1969|200010|201266|2033|2040|2099|211|218|2194|221079|2258|2289|22949|23054|2308|23491|23516|23657|23764|241|2495|2512|2539|25800|26509|2678|27063|27173|27190|2729|2730|283375|283848|2877|2878|2908|29103|2936|2938|2941|2944|2946|2947|2948|2949|2950|2953|29985|29988|30|3082|3084|3162|3280|329|3320|3326|3337|335|336|34|340024|347902|348932|3589|3725|3726|3727|376497|387509|388662|3895|405|4097|4199|4240|4257|4258|4259|4609|4616|4780|5052|51|51129|5142|51426|5155|5226|5272|5360|54566|5465|5467|55117|55224|5524|55334|55630|56262|56606|57007|5705|57181|57678|581|595|5997|6256|6319|6337|6342|6347|6397|64116|645|6513|6515|6518|6524|6526|6528|6530|6531|6533|6534|6535|6538|6591|6594|66035|6667|6714|6720|6774|6817|7039|7040|7042|7048|7049|7056|7128|7266|7295|7296|7421|7546|79056|7922|79571|80315|81031|8140|81706|8202|8204|8243|84002|8431|8454|84951|8507|8553|8648|8660|8714|873|874|8824|8856|8878|8884|89795|9049|91252|92086|92737|9588|9792|9817|983|9965| |
| <http://www.wikipathways.org/instance/WP3942_r94205> | PPAR signaling pathway%WikiPathways\_20180810%WP3942%Homo sapiens                                              |         45|        45|  0.0396740|  0.0009|        0.0432000|    0.0037|          0.1776000|             0.47700| NA      | 10062|10580|10998|10999|11001|116519|126129|1374|1376|1579|1593|1622|1962|2170|2171|2172|2173|2180|2181|2182|28965|30|335|336|34|3611|376497|4023|4199|51|5106|51129|5170|51703|5346|5360|5465|5467|6256|6257|6319|6342|8309|8310|9415|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| <http://www.wikipathways.org/instance/WP2878_r94794> | PPAR Alpha Pathway%WikiPathways\_20180810%WP2878%Homo sapiens                                                  |         20|        20|  0.0560545|  0.0024|        0.0691200|    0.0044|          0.1810286|             0.72800| NA      | 10062|1019|116519|1374|1376|1579|1622|1962|30|335|336|34|376497|4609|5360|5465|595|6256|6342|983|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| <http://www.wikipathways.org/instance/WP2857_r87780> | Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       |        131|       131|  0.0484928|  0.0020|        0.0640000|    0.0050|          0.1800000|             0.11600| NA      | 10124|10155|10413|10451|10637|10973|11169|11190|115825|121536|131405|1466|1488|154091|1789|2005|203228|2131|2132|2253|2260|22823|22943|2296|23007|2303|23181|23373|23499|23528|23576|2625|2627|26610|27244|27245|27324|29072|3064|31|3169|3170|3172|3251|3717|3720|4086|4087|4088|4089|4091|4780|4838|5080|5087|5090|51176|51366|51701|5308|54799|54892|5515|5566|55704|5573|55809|56946|57045|57198|57669|58499|5915|5916|595|6001|6169|63035|64321|652|655|657|659|6615|6657|6722|688|6899|6911|6925|6926|7003|708|7403|7546|7547|7855|79035|79668|79776|79923|79977|8030|80304|80312|8091|8312|8313|8320|8322|83439|83881|84159|84295|8450|8463|84667|85363|8553|8554|8646|8728|8928|90|9113|92|9241|93|9314|9573|9760|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| <http://www.wikipathways.org/instance/WP4258_r97136> | LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens |         78|        78|  0.0390898|  0.0015|        0.0617143|    0.0061|          0.1952000|             0.68800| NA      | 10023|1021|1024|10297|11197|11211|1452|1454|1457|1459|1460|147111|1487|1488|1499|1855|1856|1857|2146|22943|23401|2535|27121|2932|3190|3192|324|3725|4040|4041|407040|4609|467|4919|4920|51176|51384|51701|5176|5328|55506|56998|57680|59343|595|6259|6422|6423|6425|64321|64840|6885|6929|6932|6934|7474|7475|7476|7477|7480|7855|7976|80319|8061|81029|8312|8313|8321|8323|8324|83439|83999|85407|85409|8607|894|896|9350|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| <http://www.wikipathways.org/instance/WP428_r94896>  | Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                |         98|        98|  0.0411540|  0.0018|        0.0648000|    0.0083|          0.2390400|             0.77700| NA      | 10023|11197|11211|144165|1452|1454|1457|1459|1460|147111|1487|1488|1499|166336|1855|1856|1857|2239|22943|23002|23236|23401|23500|2535|27121|27130|2932|324|3725|387|4040|4041|4609|4772|4773|4775|4776|4919|4920|51176|51384|51701|5176|5328|5330|5331|5332|5530|5532|5533|5534|5578|5579|5582|5599|5601|56998|57216|57680|5879|59343|595|6259|6422|6423|6425|64321|64840|6885|6932|6934|7474|7475|7476|7477|7480|7855|7976|80319|8061|81029|815|816|817|818|81839|8312|8321|8323|8324|83439|83999|85407|85409|894|896|9350|9475|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |

``` r
enrichmentResultWiki$results %>% 
  dplyr::select(ID, CorrectedPvalue, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  head(10) %>% 
  kable(align = "l", col.names = c("ID", "Corrected p-value", 
                                   "Corrected MF p-value"))
```

| ID                                                                                                             | Corrected p-value | Corrected MF p-value |
|:---------------------------------------------------------------------------------------------------------------|:------------------|:---------------------|
| Endoderm Differentiation%WikiPathways\_20180810%WP2853%Homo sapiens                                            | 0.0000000         | 0.0000000            |
| Ectoderm Differentiation%WikiPathways\_20180810%WP2858%Homo sapiens                                            | 0.0000000         | 0.0000000            |
| Cardiac Progenitor Differentiation%WikiPathways\_20180810%WP2406%Homo sapiens                                  | 0.0000000         | 0.0096000            |
| Statin Pathway%WikiPathways\_20180810%WP430%Homo sapiens                                                       | 0.0288000         | 0.0504000            |
| PPAR signaling pathway%WikiPathways\_20180810%WP3942%Homo sapiens                                              | 0.0432000         | 0.1776000            |
| Mesodermal Commitment Pathway%WikiPathways\_20180810%WP2857%Homo sapiens                                       | 0.0640000         | 0.1800000            |
| PPAR Alpha Pathway%WikiPathways\_20180810%WP2878%Homo sapiens                                                  | 0.0691200         | 0.1810286            |
| Nuclear Receptors Meta-Pathway%WikiPathways\_20180810%WP2882%Homo sapiens                                      | 0.0000000         | 0.1843200            |
| LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways\_20180810%WP4258%Homo sapiens | 0.0617143         | 0.1952000            |
| Wnt Signaling Pathway%WikiPathways\_20180810%WP428%Homo sapiens                                                | 0.0648000         | 0.2390400            |

#### Trying to run everything with fgsea

``` r
# scores forfgsea (for absolute)
#ermineInputGeneScoresFGSEA <- DEgenes_0h_96h_75748 %>% 
#  rownames_to_column("gene") %>%
#  #mutate(absolute_logFC = abs(logFC)) %>% 
#  dplyr::select(gene, logFC) %>% 
#  na.omit() %>% 
#  as.data.frame() %>% 
#  arrange(desc(logFC)) %>% 
#  column_to_rownames("gene")


scoresFGSEA <- ermineInputProbeScores$absolute_logFC
names(scoresFGSEA) <- rownames(ermineInputProbeScores)

#for wikipathways
#ermineInputGeneScoresWiki <- ermineInputGeneScoresFGSEA %>%
#  rownames_to_column("SYMBOL") %>%
#  filter(SYMBOL %in% geneEntrez$SYMBOL)

#scoresWikiFGSEA <- merge(ermineInputGeneScoresWiki, geneEntrez, sort=FALSE)

#scoresWikiFGSEA <- scoresWikiFGSEA %>%
#  dplyr::select(logFC, ENTREZID) %>%
#  column_to_rownames('ENTREZID')

scoresWiki <- inputGSEA$absolute_logFC
names(scoresWiki) <- rownames(inputGSEA)
```

Reactome:

``` r
pathwaysReactome <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/ReactomePathways.tsv")

#ES – enrichment score, same as in Broad GSEA implementation;
#NES – enrichment score normalized to mean enrichment of random samples of the same size;
fgseaRes <- fgsea(pathwaysReactome, scoresFGSEA, minSize=15, maxSize=300, nperm=10000)
reactomeIdPathway <- reactome_pathways[,c(1,2)]
colnames(reactomeIdPathway) <- c("pathway", "Description") 
fgseaRes <- merge(fgseaRes, reactomeIdPathway, sort = FALSE)

#activated pathways
activPathwaysReactome <- fgseaRes %>% 
  arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05) %>% 
  dplyr::select(c("pathway", "Description", "padj", "ES", "NES", "nMoreExtreme"))

#inactivated pathways
inactivPathwaysReactome <- fgseaRes %>% 
  arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05) %>%
  dplyr::select(c("pathway", "Description", "padj", "ES", "NES", "nMoreExtreme"))
```

There are 15 activated and 56 inactivated pathways in Reactome for this data set.

All active pathways:

``` r
activPathwaysReactome %>% kable()
```

| pathway       | Description                                            |       padj|         ES|       NES|  nMoreExtreme|
|:--------------|:-------------------------------------------------------|----------:|----------:|---------:|-------------:|
| R-HSA-446203  | Asparagine N-linked glycosylation                      |  0.0160987|  0.5028584|  1.820588|             0|
| R-HSA-6807878 | COPI-mediated anterograde transport                    |  0.0160987|  0.5632953|  1.773994|             3|
| R-HSA-199977  | ER to Golgi Anterograde Transport                      |  0.0160987|  0.5217251|  1.748014|             3|
| R-HSA-1474244 | Extracellular matrix organization                      |  0.0160987|  0.4478595|  1.598895|             2|
| R-HSA-381070  | IRE1alpha activates chaperones                         |  0.0160987|  0.6601280|  1.952930|             0|
| R-HSA-174824  | Plasma lipoprotein assembly, remodeling, and clearance |  0.0160987|  0.6397869|  1.843138|             2|
| R-HSA-381119  | Unfolded Protein Response (UPR)                        |  0.0160987|  0.6524081|  2.025502|             0|
| R-HSA-381038  | XBP1(S) activates chaperone genes                      |  0.0160987|  0.6553692|  1.919402|             0|
| R-HSA-4791275 | Signaling by WNT in cancer                             |  0.0168293|  0.7041005|  1.897812|             3|
| R-HSA-416700  | Other semaphorin interactions                          |  0.0190234|  0.8076681|  1.906055|             4|
| R-HSA-114608  | Platelet degranulation                                 |  0.0190234|  0.5504052|  1.788574|             5|
| R-HSA-76005   | Response to elevated platelet cytosolic Ca2+           |  0.0191475|  0.5308110|  1.736947|             6|
| R-HSA-948021  | Transport to the Golgi and subsequent modification     |  0.0191475|  0.4763411|  1.634816|             6|
| R-HSA-4085001 | Sialic acid metabolism                                 |  0.0478185|  0.6632057|  1.753215|            21|
| R-HSA-1251985 | Nuclear signaling by ERBB4                             |  0.0488704|  0.6953696|  1.715482|            21|

KEGG:

``` r
pathwaysKEGG <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/KeggPathways.gmt")

fgseaRes <- fgsea(pathwaysKEGG, scoresFGSEA, minSize=15, maxSize=300, nperm=10000)

#activated pathways
activPathwaysKEGG <- fgseaRes %>% arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05)

#inactivated pathways
inactivPathwaysKEGG <- fgseaRes %>% arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05)
```

There are 0 activated and 0 inactivated pathways in Reactome for this data set.

All active pathways:

``` r
activPathwaysKEGG
```

    ## [1] pathway      pval         padj         ES           NES         
    ## [6] nMoreExtreme size         leadingEdge 
    ## <0 rows> (or 0-length row.names)

WikiPathways:

``` r
pathwaysWiki <- gmtPathways("~/ESC_RNA_seq/pathway_enrichment_analysis/Wikipathways.gmt")

fgseaRes <- fgsea(pathwaysWiki, scoresWiki, minSize=15, maxSize=300, nperm=10000)

#activated pathways
activPathwaysWiki <- fgseaRes %>% arrange(padj) %>% filter(NES > 0) %>% filter(padj <= 0.05)

#inactivated pathways
inactivPathwaysWiki <- fgseaRes %>% arrange(padj) %>% filter(NES < 0) %>% filter(padj <= 0.05)
```

There are 13 activated and 6 inactivated pathways in Reactome for this data set.

All active pathways:

``` r
activPathwaysWiki
```

    ##                                                                                                          pathway
    ## 1                                       Nuclear Receptors Meta-Pathway%WikiPathways_20180810%WP2882%Homo sapiens
    ## 2                                                         NRF2 pathway%WikiPathways_20180810%WP2884%Homo sapiens
    ## 3                                             Ectoderm Differentiation%WikiPathways_20180810%WP2858%Homo sapiens
    ## 4                                             Endoderm Differentiation%WikiPathways_20180810%WP2853%Homo sapiens
    ## 5               Photodynamic therapy-induced unfolded protein response%WikiPathways_20180810%WP3613%Homo sapiens
    ## 6                                                 Wnt Signaling Pathway%WikiPathways_20180810%WP428%Homo sapiens
    ## 7                                   Cardiac Progenitor Differentiation%WikiPathways_20180810%WP2406%Homo sapiens
    ## 8  LncRNA involvement in canonical Wnt signaling and colorectal cancer%WikiPathways_20180810%WP4258%Homo sapiens
    ## 9                                                    Heart Development%WikiPathways_20180810%WP1591%Homo sapiens
    ## 10                                                       Statin Pathway%WikiPathways_20180810%WP430%Homo sapiens
    ## 11                                Selenium Metabolism and Selenoproteins%WikiPathways_20180810%WP28%Homo sapiens
    ## 12                                        GPCRs, Class A Rhodopsin-like%WikiPathways_20180810%WP455%Homo sapiens
    ## 13                                                         Adipogenesis%WikiPathways_20180810%WP236%Homo sapiens
    ##            pval       padj        ES      NES nMoreExtreme size
    ## 1  0.0001186240 0.01690141 0.4791373 1.717150            0  229
    ## 2  0.0003994674 0.01690141 0.5316593 1.746136            2  106
    ## 3  0.0001299714 0.01690141 0.6453161 2.158787            0  122
    ## 4  0.0005149331 0.01690141 0.5105947 1.721432            3  129
    ## 5  0.0006417455 0.01690141 0.7293670 1.862191            3   24
    ## 6  0.0006721334 0.01690141 0.5237013 1.704111            4   98
    ## 7  0.0004517392 0.01690141 0.6500998 1.845202            2   41
    ## 8  0.0005585032 0.01690141 0.5569046 1.750251            3   78
    ## 9  0.0012334258 0.02537333 0.6597535 1.807221            7   34
    ## 10 0.0018050542 0.03061197 0.7230128 1.777805           10   20
    ## 11 0.0018069568 0.03061197 0.6119552 1.736935           11   41
    ## 12 0.0029264214 0.04630516 0.5149895 1.621557           20   79
    ## 13 0.0030548546 0.04630516 0.4763388 1.568483           22  108
    ##                                                                                                                                                                                                                                                                                                                                                   leadingEdge
    ## 1  336, 57007, 22949, 29103, 4258, 1958, 200010, 7056, 2877, 10257, 8507, 1028, 5272, 6720, 2878, 9792, 7040, 2194, 6515, 6347, 6319, 8824, 1728, 645, 11214, 23657, 1969, 1962, 10486, 3727, 4199, 51, 23054, 7546, 8878, 2953, 23764, 3084, 6714, 1027, 8204, 6526, 26509, 4097, 11182, 376497, 27063, 1622, 1836, 5465, 8553, 8660, 1839, 2539, 335, 91252
    ## 2                                                                                                                                                                                           22949, 4258, 1958, 200010, 2877, 10257, 2878, 7040, 6515, 8824, 1728, 645, 23657, 1969, 4199, 8878, 2953, 23764, 3084, 6526, 4097, 11182, 1839, 2539, 91252, 4780
    ## 3                                                                                      2627, 7903, 3975, 3170, 5362, 4920, 56899, 6781, 54898, 6469, 10257, 7855, 57462, 8495, 114815, 4061, 345557, 152789, 9687, 6347, 25987, 3728, 8322, 166, 10634, 10486, 8820, 1896, 6238, 57476, 56963, 5783, 10253, 4772, 10013, 8345, 81576, 8848, 6386, 91653, 3801
    ## 4                                                                                                                                                                                                                                                      64321, 2627, 3087, 22943, 3975, 3170, 6422, 83881, 1846, 9350, 8320, 5015, 2626, 7040, 3169, 655, 6772
    ## 5                                                                                                                                                                                                                                           1649, 57761, 2081, 4189, 23645, 7453, 467, 5611, 9695, 4780, 10130, 51726, 7184, 22926, 440, 3309, 811, 9451, 468
    ## 6                                                                                                                                                                                                                                             64321, 22943, 4920, 6422, 9350, 27121, 2535, 7855, 817, 81839, 64840, 27130, 2239, 5176, 324, 8324, 23002, 4772
    ## 7                                                                                                                                                                                                                                                                                                     64321, 7852, 22943, 3170, 4920, 83881, 2626, 3815, 7040
    ## 8                                                                                                                                                                                                                                                                       64321, 22943, 4920, 6422, 9350, 27121, 2535, 7855, 55506, 64840, 467, 5176, 324, 8324
    ## 9                                                                                                                                                                                                                                                                                        2627, 3170, 650, 2626, 2296, 6469, 7422, 4772, 659, 8553, 8928, 1499
    ## 10                                                                                                                                                                                                                                                                                                            336, 19, 6646, 341, 348, 3949, 64714, 335, 5360
    ## 11                                                                                                                                                                                                                                                           6414, 1735, 2877, 2878, 85465, 55829, 1491, 140606, 58515, 51714, 51091, 4780, 6301, 1390, 54952
    ## 12                                                                                                                                                                                                                                                                                7852, 57007, 887, 2863, 1909, 2861, 4986, 187, 624, 23432, 3354, 6915, 4886
    ## 13                                                                                                                                                                          83881, 2487, 650, 2626, 1592, 1649, 3977, 8788, 6095, 23175, 57761, 6720, 7040, 688, 6319, 3572, 80339, 1316, 6772, 28999, 8204, 6773, 10555, 9021, 5465, 8660, 5106, 1647, 10135
