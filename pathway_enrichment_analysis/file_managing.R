# CUSTOM ANNOTATION FILE

## read in the file to determine the max number of columns
x <- scan("~/ESC_RNA_seq/pathway_enrichment_analysis/Generic_human_original", what = "", sep = "\n")
x <- strsplit(x, "[ \t]+")
max.col <- max(sapply(x, length))

## specify col.names as ?read.table suggests
cn <- paste("V", 1:max.col, sep = "")
z1 <- read.table("~/ESC_RNA_seq/pathway_enrichment_analysis/Generic_human_original", fill = TRUE, col.names = cn)

## creating a new annotation with just two columns: ProbeName and GeneSymbols

new_annotation <- z1[,1:2]
new_annotation <- new_annotation[-1,]
colnames(new_annotation) <- c("ProbeName", "GeneSymbols")
write.table(new_annotation, file = "~/ESC_RNA_seq/pathway_enrichment_analysis/Generic_human", 
            row.names = FALSE, quote = FALSE, sep = "\t")

#################################################################################################
# ANALYZING REACTOME PATHWAYS FILE AND CHANGING IT A LITTLE BIT

## read in the file to determine the max number of columns
x <- scan("~/ESC_RNA_seq/pathway_enrichment_analysis/ReactomePathways.gmt", what = "", sep = "\n")
x <- strsplit(x, '\t')
max.col <- max(sapply(x, length))

## specify col.names as ?read.table suggests
cn <- paste("V", 1:max.col, sep = "")
z1 <- read.table("~/ESC_RNA_seq/pathway_enrichment_analysis/ReactomePathways.gmt", fill = TRUE, col.names = cn,
                 sep = "\t", quote="")

## reformatting the file according to specification:
## One one line, the fields are:
##  
##  1. A unique gene set identifier
##  2. A description (can be blank, but cannot be ommitted)
##  3. The remaining fields are interpreted as gene symbols (keyed to the second column of your annotation file).

## deleting Reactome pathway column
z1 <- z1[,-3]
## changing the order of columns
z1 <- z1[,c(2,1,3:ncol(z1))]

## turning blank cells into NAs for futher easy removing
z1 <- apply(z1, 2, function(x) gsub("^$|^ $", NA, x))

## I hate for loops but here it's just for simplicity (don't blame me too much)
## probably doing it with lapply should dramaticaly increase the speed
for (i in 1:nrow(z1)){
  
  one_row <- z1[i,]
  one_row <- one_row[!is.na(one_row)]
  one_row <- as.data.frame(t(one_row))
  write.table(one_row, file = "~/ESC_RNA_seq/pathway_enrichment_analysis/ReactomePathways.tsv", 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
}

