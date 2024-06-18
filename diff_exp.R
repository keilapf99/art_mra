#
#3.differential_expression_DEseq2.R
#This script perform 

#Libraries --- ---

#BiocManager::install("DESeq2")

pacman::p_load("dplyr", 
               "DESeq2", 
               "pheatmap", 
               "ggplot2", 
               "ggrepel", 
               "biomaRt")

#Define functions --- ---

#Function to translate gene names
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Connecting to the Ensembl database through biomaRt


# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = ensembl)
}

#Get data --- ---

counts <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.txt") %>% as.data.frame()
counts <- counts[ -c(1:4),] #Delete alignment stats
dim(counts)

#
rownames(counts) <- counts$feature
dim(counts)
#[1] 60603   1142

metadata <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
dim(metadata)
#[1] 1141   41

table(metadata$dicho_NIA_reagan, useNA = "ifany")

#  0    1 <NA> 
#  307  573  261  

#Filter to obtain only the ones with NIA-Reagan dicho

metadata <- metadata %>% filter(!is.na(dicho_NIA_reagan))

#Only samples with metadata
counts <- counts %>% dplyr::select(all_of(metadata$specimenID))
dim(counts)
#[1] 60603   880

#Differential expression --- ---

#Experimental design

coldata <- as.data.frame(colnames(counts))
colnames(coldata) <- "specimenID"

coldata <- coldata %>% left_join(metadata, by = "specimenID")
coldata <- coldata %>% dplyr::select("specimenID", "msex", "sequencingBatch",  "cogdx", "ceradsc", "dicho_NIA_reagan")
coldata$dicho_NIA_reagan <- as.factor(coldata$dicho_NIA_reagan)  #Convert to factor

#DESeqData object

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ sequencingBatch + msex + dicho_NIA_reagan)

#Pre-filtering

# A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 treated samples.
smallestGroupSize <- 3
#Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples.
#The count of 10 is a reasonable choice for bulk RNA-seq.
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds)
#[1] 30657   880

#Specify conditions to compare
#dds$condition <- factor(dds$condition, levels = c("untreated","treated"))

dds$dicho_NIA_reagan <- factor(dds$dicho_NIA_reagan, levels = c("0", "1"))

#Differential expression analysis --- ---

dds <- DESeq(dds)  #Little slow
#Results
res <- results(dds, alpha=0.05) #alpha is the FDR limit

#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
#[1] 4914

#How many adjusted p-values were less than 0.05?
sum(res$padj < 0.05, na.rm=TRUE)
#[1] 2984

#Order results by p-value
res.df <- res[order(res$padj),] %>% as.data.frame()

#Remove version gene
rownames(res.df) <- str_remove(rownames(res.df), "\\..*$")

# add a column of NAs
res.df$diffexpressed <- "NO"
# if log2Foldchange > 0.1 and pvalue < 0.05, set as "UP" 
res.df$diffexpressed[res.df$log2FoldChange > 1 & res.df$padj < 0.05] <- "UP"
# if log2Foldchange < -0.1 and pvalue < 0.05, set as "DOWN"
res.df$diffexpressed[res.df$log2FoldChange < -1 & res.df$padj < 0.05] <- "DOWN"

#Add gene names in SYMBOL --- ---

#Create dictionary
symbol <- convert_ens_to_symbol(rownames(res.df))
symbol$external_gene_name <- ifelse(symbol$external_gene_name == "", symbol$ensembl_gene_id, symbol$external_gene_name)
#merge

# Convertir rownames en columna
res.df <- res.df %>% mutate(ensembl_gene_id = rownames(.), .before = 1) 
res.df <- res.df %>% left_join(symbol, by = "ensembl_gene_id")

#Extract DEGs --- ---

DEGS <- res.df %>% filter(diffexpressed != "NO")
rownames(DEGS) <- DEGS$ensembl_gene_id

#Create matrix of DEG counts --- ---

DEG_mat <- BiocGenerics::counts(dds, normalized = T) %>% as.data.frame()
rownames(DEG_mat)  <- str_remove(rownames(DEG_mat) , "\\..*$")

DEG_mat <- DEG_mat %>% filter(rownames(DEG_mat) %in% DEGS$ensembl_gene_id)
dim(DEG_mat)

# Create correlation matrix

DEG_mat.z <- t(apply(DEG_mat, 1, scale))
#Add sample names
colnames(DEG_mat.z) <- colnames(DEG_mat) 

#Separe by condition

condition <- data.frame(specimenID = colnames(DEG_mat.z))
condition <- condition %>% left_join(coldata, by = "specimenID")
condition <- condition %>% dplyr::select(specimenID, dicho_NIA_reagan)

#Heat

ComplexHeatmap::Heatmap(DEG_mat.z, cluster_rows = T, cluster_columns = T, 
                        column_labels = colnames(DEG_mat.z),column_title = "Samples", row_title = "DEGs",
                       # column_split = list(specimenID = condition$specimenID, dicho_NIA_reagan = condition$dicho_NIA_reagan),
                        name = "Z-score", row_labels = DEGS[rownames(DEG_mat.z),]$external_gene_name)

#Vulcano plot --- ---

#Create labels for Vplot
res.df <- res.df %>% mutate(delabel = ifelse(diffexpressed != "NO", external_gene_name, NA))

#
vplot <- ggplot(data=res.df, aes(x=log2FoldChange, y=-log10(pvalue),  col=diffexpressed, label=delabel)) +
  geom_point() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") + #log2FoldChange threshold is 0.6
  geom_hline(yintercept=-log10(0.05), col="red") + #p-value threshold is 0.05
  scale_color_manual(values=c("#4E8098", "gray", "#A31621")) +
  theme_minimal() +
  geom_text_repel() +
  labs(title = "Differential expression", 
        sub = "Using dichotomic NIA-Reagan Criteria")
  
#Pheatmap 

ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:10]

pheatmap.df <- as.data.frame(colData(dds)[,c("dicho_NIA_reagan","msex")])

pheatmap(assay(ntd)[select,], cluster_rows = TRUE, show_rownames = FALSE,
         cluster_cols = TRUE , annotation_col = pheatmap.df)

#Save data --- ---

#Save DEGs data
vroom::vroom_write()
