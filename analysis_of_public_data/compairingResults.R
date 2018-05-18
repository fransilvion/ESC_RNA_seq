library(data.table)

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

save(all_upRegulated, file="allUpRegulated.Rdata")
save(all_downRegulated, file="allDownRegulated.Rdata")
