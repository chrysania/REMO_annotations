library(dplyr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(hdf5r)
options(scipen = 999)

# read REMO
remo <- read.table('../../data/REMOv1_GRCh38.bed.gz', col.names = c("chr", "start", "end", "Module"))

# add GC content
get_GC <- function(x, genome) {
    sequences <- Biostrings::getSeq(x = genome, x)
    gc <- Biostrings::letterFrequency(sequences, letters = 'CG') / width(sequences) * 100
    colnames(gc) <- 'GC.percent'
    return(gc)
}

remo_gr <- makeGRangesFromDataFrame(remo, keep.extra.columns = TRUE)
gc_content <- get_GC(remo_gr, BSgenome.Hsapiens.UCSC.hg38)
remo_gr$GC <- as.numeric(gc_content[, 1])

remo <- as.data.frame(remo_gr)

# count modules and total bases
remo <- remo %>% group_by(Module) %>% summarize(CREs = n(), Bases = sum(width), Chromosome = seqnames, GC_mean = mean(GC)) %>% unique

remo <- as.data.frame(remo)
remo$GC_mean <- round(remo$GC_mean, 2)
rownames(remo) <- remo$Module



# read all REMO modules with bio annotations
remo$Bioterms <- NA
all_chromosome <- paste0("chr", c(1:22, "X"))
for (i in all_chromosome) {
    message(i)
    annot <- readr::read_tsv(paste0("remo_annotation_whole_body_", i, ".tsv"), col_names = c('Module', 'Terms'))
    annot <- annot[annot$Module %in% rownames(remo), ]
    remo[annot$Module, 'Bioterms'] <- annot$Terms
}

rownames(remo) <- NULL

# sort
remo <- remo[gtools::mixedorder(remo$Module, decreasing = TRUE), ]

write.table(x = remo, file = "../../data/whole_body/REMOv1_GRCh38_metadata.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# build module x term sparse matrix
# get unique terms
bioterms_list <- strsplit(remo$Bioterms, " \\| ")
bioterms_flattened <- unlist(bioterms_list)
unique_terms <- unique(bioterms_flattened)

# build triplet matrix
row_indices <- vector(mode = 'numeric')
col_indices <- vector(mode = 'numeric')

for (i in 1:nrow(remo)) {
  current_bioterms <- bioterms_list[[i]]
  for (term in current_bioterms) {
    col_idx <- fastmatch::fmatch(term, unique_terms)
    row_indices <- c(row_indices, i)
    col_indices <- c(col_indices, col_idx)
  }
}

membership <- sparseMatrix(
    i = row_indices,
    j = col_indices,
    x = 1,
    dims = c(nrow(remo), length(unique_terms)),
    dimnames = list(remo$Module, unique_terms)
)

# write as mtx file
writeMM(membership, "../../data/whole_body/REMOv1_GRCh38_bioterms/counts.mtx")
writeLines(rownames(membership), "../../data/whole_body/REMOv1_GRCh38_bioterms/rownames.txt")
writeLines(colnames(membership), "../../data/whole_body/REMOv1_GRCh38_bioterms/colnames.txt")

# write as h5 file
h5file <- H5File$new("../../data/whole_body/REMOv1_GRCh38_bioterms.h5", mode = "w")
h5file.grp <- h5file$create_group("annotations")

h5file.grp[["i"]] <- as.integer(row_indices)
h5file.grp[["j"]] <- as.integer(col_indices)
h5file.grp[["x"]] <- as.integer(membership@x)
h5file.grp[["dims"]] <- as.integer(dim(membership))
h5file.grp[["row_names"]] <- rownames(membership)
h5file.grp[["col_names"]] <- colnames(membership)
h5file$close()

# write as rds file
saveRDS(membership, "../../data/whole_body/REMOv1_GRCh38_bioterms.rds")
