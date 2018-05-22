library(data.table)
suppressMessages(suppressWarnings(library(hgu133a.db)))
suppressMessages(suppressWarnings(library(annotate)))
suppressMessages(suppressWarnings(library(ermineR)))
suppressMessages(suppressWarnings(library(ggplot2)))
library(dplyr)

#loading data 
load("analysis_of_public_data/GSE109658/upRegulated109658.Rdata") #data.frame
load("analysis_of_public_data/GSE109658/downRegulated109658.Rdata") #data.frame
load("analysis_of_public_data/GSE75748/upRegulated75748.Rdata") #data.frame
load("analysis_of_public_data/GSE75748/downRegulated75748.Rdata") #data.frame
load("analysis_of_public_data/GSE52158/upRegulated52158.Rdata") #character vector
load("analysis_of_public_data/GSE52158/downRegulated52158.Rdata") #character vector

upGenes109658 <- upRegulated109658$gene #1619
downGenes109658 <- downRegulated109658$gene #990
upGenes75748 <- upRegulated75748$gene #1984
downGenes75748 <- downRegulated75748$gene #2152
upGenes52158 <- upRegulated52158
downGenes52158 <- downRegulated52158

common_upRegulated_1 <- intersect(upGenes109658, upGenes75748) #735
common_downRegulated_1 <- intersect(downGenes109658, downGenes75748) #433

common_upRegulated_2 <- intersect(upGenes52158, upGenes75748) #258
common_downRegulated_2 <- intersect(downGenes52158, downGenes75748) #123

common_upRegulated_3 <- intersect(upGenes109658, upGenes52158) #295
common_downRegulated_3 <- intersect(downGenes109658, downGenes52158) #151

#GSC and FOXA2 are not here
all_upRegulated <- intersect(upGenes109658, intersect(upGenes75748, upGenes52158)) #208
#POU5F1, DPPA4, CDH1 are not here
all_downRegulated <- intersect(downGenes109658, intersect(downGenes75748, downGenes52158)) #84

#take genes that are present in at least two intersection sets:
confident_up <- c()
confident_up <- c( confident_up , common_upRegulated_1[ ! common_upRegulated_1 %chin% confident_up ] ) #735
confident_up <- c( confident_up , common_upRegulated_2[ ! common_upRegulated_2 %chin% confident_up ] ) #785
confident_up <- c( confident_up , common_upRegulated_3[ ! common_upRegulated_3 %chin% confident_up ] ) #872

confident_down <- c()
confident_down <- c( confident_down , common_downRegulated_1[ ! common_downRegulated_1 %chin% confident_down ] ) #433
confident_down <- c( confident_down , common_downRegulated_2[ ! common_downRegulated_2 %chin% confident_down ] ) #472
confident_down <- c( confident_down , common_downRegulated_3[ ! common_downRegulated_3 %chin% confident_down ] ) #539

#saving all up and down regulated
save(all_upRegulated, file="allUpRegulated.Rdata")
save(all_downRegulated, file="allDownRegulated.Rdata")

#saving common up and down regulated in 109658, 75748 #1
save(common_upRegulated_1, file="common_upRegulated_109658_75748.Rdata")
save(common_downRegulated_1, file="common_downRegulated_109658_75748.Rdata")

#saving common up and down regulated in 52158, 75748 #2
save(common_upRegulated_2, file="common_upRegulated_52158_75748.Rdata")
save(common_downRegulated_2, file="common_downRegulated_52158_75748.Rdata")

#saving common up and down regulated in 109658, 52158 #3
save(common_upRegulated_3, file="common_upRegulated_52158_109658.Rdata")
save(common_downRegulated_3, file="common_downRegulated_52158_109658.Rdata")

#saving confident up and down
save(confident_up, file="in_at_least_two_sets.Rdata")
save(confident_down, file="in_at_least_two_sets.Rdata")

#==============================================
#converting golden markers files into probes
x <- hgu133aSYMBOL
# Get the probe identifiers - gene symbol mappings
mapped_probes <- mappedkeys(x)
# Convert to a dataframe
genesym.probeid <- as.data.frame(x[mapped_probes])
#head(genesym.probeid)

up_vector <- scan("~/ESC_RNA_seq/FILES_FOR_LINCS/Golden_list_Up.grp", character(), quote = "") 
down_vector <- scan("~/ESC_RNA_seq/FILES_FOR_LINCS/Golden_list_Down.grp", character(), quote = "") 

probe_up <- genesym.probeid %>%
  filter(symbol %in% up_vector) 

probe_up <- probe_up$probe_id

probe_down <- genesym.probeid %>%
  filter(symbol %in% down_vector) 

probe_down <- probe_down$probe_id

fileConn2<-file("~/ESC_RNA_seq/FILES_FOR_LINCS/Golden_list_Down_Probes.grp")
writeLines(probe_down, fileConn2)
close(fileConn2)

fileConn2<-file("~/ESC_RNA_seq/FILES_FOR_LINCS/Golden_list_Up_Probes.grp")
writeLines(probe_up, fileConn2)
close(fileConn2)

#==============================================
#best of the best DE markers
up_list <- c("SNAI1", "CER1", "SOX17", "CXCR4", "FOXA2", "CDH2", "GATA4", "GATA6", "KIT", "KRT19", "ITGA5", "CD55")
down_list <- c("FGF8", "WNT3", "T", "FGF4", "MIXL1", "LHX1", "SOX7", "CDX1", "SOX1", "ZIC1", "POU5F1",
               "SOX2", "FOXF1", "MEOX1", "FLK1", "BMP4", "CDX2")

probe_up_2 <- genesym.probeid %>%
  filter(symbol %in% up_list) 

probe_up_2 <- probe_up_2$probe_id

probe_down_2 <- genesym.probeid %>%
  filter(symbol %in% down_list) 

probe_down_2 <- probe_down_2$probe_id

fileConn2<-file("~/ESC_RNA_seq/FILES_FOR_LINCS/Best_list_Down_Probes.grp")
writeLines(probe_down_2, fileConn2)
close(fileConn2)

fileConn2<-file("~/ESC_RNA_seq/FILES_FOR_LINCS/Best_list_Up_Probes.grp")
writeLines(probe_up_2, fileConn2)
close(fileConn2)


#==============================================
#enrichment analysis for compounds (do it for every file in Compounds folder)
#no significant results for KI_8751, periplocymarin, testosterone, WYE_354, hydroxyfasudil
compound <- read.delim("Compounds/hydroxyfasudil_A549.gct", header = F)

compound <- compound[-c(1:11),]
compound <- compound[,-c(4:8)]
colnames(compound) <- c("id", "symbol", "function", "space", "cell_1", "cell_2")

compound$cell_1 <- as.numeric(levels(compound$cell_1))[compound$cell_1]
compound$cell_2 <- as.numeric(levels(compound$cell_2))[compound$cell_2]

compound$average <- (compound$cell_1 + compound$cell_2)/2
#compound$average <- compound$cell_1

sorted_for_enrichment <- compound[with(compound, order(-average)), ]
sorted_for_enrichment <- sorted_for_enrichment[,c("symbol", "average")]
rownames(sorted_for_enrichment) <-  c()
sorted_for_enrichment <- sorted_for_enrichment %>% column_to_rownames("symbol")

enrichmentResult <- precRecall(scores = sorted_for_enrichment, 
                               scoreColumn = 1, # column 1 is the scores 
                               bigIsBetter = TRUE, # larger logFC should be ranked higher
                               annotation = "analysis_of_public_data/GSE109658/Generic_human", # ask ermineJ to use the Generic_human annotation file (will automatically download)
                               aspects = "B", # look at only biological processes 
                               iterations = 10000, # 10K sampling iterations so that results are stable
                               geneSetDescription = "analysis_of_public_data/GSE109658/GO.xml") # use the GO XML file in current directory

x <- enrichmentResult$results %>% arrange(MFPvalue) %>% head(10)
Enrichment <- enrichmentResult$results
Enrichment$Name <- as.factor(Enrichment$Name)

Enrichment %>% 
  dplyr::select(Name, NumGenes, CorrectedMFPvalue) %>% 
  arrange(CorrectedMFPvalue) %>% 
  filter(CorrectedMFPvalue <= 5e-2) %>% 
  head(25) %>% 
  ggplot(aes(x = fct_reorder(Name, CorrectedMFPvalue), 
             y = NumGenes, fill = CorrectedMFPvalue)) +
  geom_col() +
  labs(title = "Biological Processes - fasudil_A549.gct", 
       x = "", y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() +
  theme_light() 
