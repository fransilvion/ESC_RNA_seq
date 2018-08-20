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
            row.names = FALSE, quote = FALSE)

# ANALYZING REACTOME PATHWAYS FILE AND CHANGING IT A LITTLE BIT

## read in the file to determine the max number of columns
x <- scan("~/ESC_RNA_seq/pathway_enrichment_analysis/ReactomePathways.gmt", what = "", sep = "\n")
x <- strsplit(x, '\t')
max.col <- max(sapply(x, length))

## specify col.names as ?read.table suggests
cn <- paste("V", 1:max.col, sep = "")
z1 <- read.table("~/ESC_RNA_seq/pathway_enrichment_analysis/ReactomePathways.gmt", fill = TRUE, col.names = cn,
                 sep = "\t")

## reformatting the file according to specification:
## One one line, the fields are:
##  
##  1. A unique gene set identifier
##  2. A description (can be blank, but cannot be ommitted)
##  3. The remaining fields are interpreted as gene symbols (keyed to the second column of your annotation file).

