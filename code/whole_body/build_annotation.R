library(peakcluster)
library(edgeR)
library(Matrix)
source("../utilities.R")

load_atac <- function(
    chromosome,
    whole_body = "~/scratch/MetaPeak/scATAC_atlas/data/zhang/cpeak_cre/chr/"
) {
    # adult whole body
    scatac <- t(read.table(paste0(whole_body, chromosome, ".csv.gz"), sep = ","))
    colnames(scatac) <- readLines(paste0(whole_body, chromosome, "_colnames.txt"))
    rownames(scatac) <- readLines("~/scratch/MetaPeak/scATAC_atlas/data/zhang/cpeak_cre/matrix/scatac_rownames.txt")

    return(scatac)
}

normalize_atac <- function(m, df_thresh = 0.7) { # modified normalize_atac with no ccre min accessibility filter
    # filter low coverage celltypes
    df <- rowSums(m > 0) / ncol(m)
    ct_use <- df > df_thresh
    m <- m[ct_use, , drop = FALSE]

#    # remove CREs with less than minimum accessibility
#    m <- m[, colSums(m) > cs_thresh, drop = FALSE]
    
    # normalize
    m_norm <- peakcluster:::normalize(m, l2 = FALSE, depthnorm = TRUE)

    return(m_norm)
}

top_frequency_experiments <- function(
    scatac_use, #encode_normalized,
    remo,
    all_ccre,
    n = 5
) {
    # get CCREs for REMO module
    idents <- all_ccre[all_ccre$REMO == remo, "id"]

    experiments <- vector(mode = "character")
    for (i in idents) {
        # get ccre column only
        all_targets <- scatac_use[, i]
        # get top 10 experiments 
        top_experiments <- order(all_targets, decreasing = TRUE)[1:10]
        celltype_names <- names(all_targets[top_experiments])
        
        experiments <- c(experiments, celltype_names)
    }
    # frequency of each experiment in vector
    exp_freq <- sort(table(experiments), decreasing = TRUE)
    return(head(names(exp_freq), n = n))
}

# load CRE information
ccre <- read.table("../../data/REMOv1_GRCh38.bed.gz")
ccre$peak <- paste(ccre$V1, ccre$V2, ccre$V3, sep="-")
ccre_all <- read.table("~/scratch/MetaPeak/ENCODE/combined_cre.tsv", sep = "\t")
ccre_all$peak <- paste(ccre_all$V1, ccre_all$V2, ccre_all$V3, sep="-")
pk <- ccre_all$V4
names(pk) <- ccre_all$peak
ccre$id <- pk[ccre$peak]
ccre$REMO <- ccre$V4

# get REMO annotations for all chromosomes
all_chromosome <- paste0("chr", c(1:22, "X"))
for (chr in all_chromosome) {
    gc()
    message(chr)

    # load and normalize brain pseudobulk matrix
    scatac_use <- load_atac(chromosome = chr)
    scatac_use <- normalize_atac(m = scatac_use, df_thresh = 0.7)

    # list of modules in the current chromosome
    chromosome_modules <- unique(ccre[ccre$V1 == chr, "REMO"])
    annot_vec <- vector(mode = "character", length = length(chromosome_modules))

    # annotate each REMO in the chromosome
    for (i in seq_along(chromosome_modules)) {
        message(i)
        annot <- top_frequency_experiments(
            scatac_use = scatac_use,
            remo = chromosome_modules[i],
            all_ccre = ccre,
            n = 5
        )
        annot_vec[i] <- paste0(annot, collapse = " | ")
    }
    df <- data.frame("REMO" = chromosome_modules, "Annotation" = annot_vec)
    write.table(df, paste0("./remo_annotation_whole_body_", chr, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
}