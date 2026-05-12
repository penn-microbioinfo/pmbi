renv::load(file.path(Sys.getenv("HOME"), "dev/pmbi/snakemake/pipelines/kallisto/"))

library(tidyverse)
library(tximport)
library(argparse)


parser <- ArgumentParser(description='Process kallisto output into a count matrix.')
parser$add_argument('-t', '--t2g', help='Path to the reference transcripts_to_genes file.')
parser$add_argument('-p', '--project', help='Project name.')
parser$add_argument('-o', '--outdir', help='Where to output the count matrix.')
parser$add_argument('kallisto_outdir', help='Path to kallisto output directory.')
args = parser$parse_args()

t2g=args$t2g
rundir=args$kallisto_outdir

# import transcript to gene mapping
tx2gene.df <- read.table(t2g, sep="\t")

# import kallisto files
sample.names <- dir(rundir)
# sample.names = sample.names[sample.names!="kallisto_logs"]
file.paths <- file.path(rundir, sample.names, "abundance.h5")
names(file.paths) <- sample.names
txi.res <- tximport(file.paths, type="kallisto", countsFromAbundance="lengthScaledTPM", tx2gene = tx2gene.df)

##### 1) EXPORT FULL TXI.RES object #####
saveRDS(txi.res, file=file.path(args$outdir, glue::glue("{args$project}_kallisto_txi.rds")))

##### 2) CREATE ENSEMBL GENE ID TO GENE NAME MAPPING
ensembl.gene.id.to.name.df <- distinct(tx2gene.df, V2, V3)
colnames(ensembl.gene.id.to.name.df) <- c("ensembl.gene.id", "gene.name")

##### 3) EXPORT lengthScaledTPM COUNTS #####

## exclude genes with 0 counts across all samples
counts.df <- as.data.frame(txi.res$counts[rowSums(txi.res$counts) > 0, ])

# %% Wrtite out the count matrix without summing across gene symbols
#    This preserves each ENSEMBL id as its own row
out.name <- file.path(args$outdir, glue::glue("{args$project}_ensembl_counts_matrix.tsv"))
write_delim(x = as_tibble(counts.df, rownames = "gene.id"),
            file = out.name,
            delim="\t",
            quote='none')

## add gene names
counts_named_by_geneid <- merge(counts.df, ensembl.gene.id.to.name.df, by.x="row.names", by.y="ensembl.gene.id")

## If you want to see which genes are duplicated, and their per-transcript counts
#duplicate_gene_names <- counts_named_by_geneid %>%
#    group_by(gene.name) %>%
#    summarise(n = n()) %>%
#    filter(n!=1)
#print(counts_named_by_geneid %>% filter(gene.name %in% duplicate_gene_names$gene.name))

# Sum across duplicate gene names, but distinct ensembl gene ids
counts_named_by_geneid_duplicate_summed <- counts_named_by_geneid %>%
    group_by(gene.name) %>%
    summarise_at(vars(!Row.names), sum)

## make gene.name the first column
counts.df.out <- counts_named_by_geneid_duplicate_summed %>% 
    select(gene.name, everything()) %>% 
    as_tibble

## write output
# out.name <- file.path(args$outdir, sprintf("%s_geneSymbolSummed__counts_matrix_n%ix%i.tsv", args$project, nrow(counts.df.out), ncol(counts.df.out)-1))
out.name <- file.path(args$outdir, glue::glue("{args$project}_geneSymbolSummed_counts_matrix.tsv"))
write_delim(x = counts.df.out,
            file = out.name,
            delim="\t",
            quote='none')

# 4) EXPORT TPM (stored in abundance)
## exclude genes with zero counts
# tpm.df <- as.data.frame(txi.res$abundance[rowSums(txi.res$abundance) > 0, ])

## add gene names
# tpm.df.w.gene.name <- merge(tpm.df, ensembl.gene.id.to.name.df, by.x="row.names", by.y="ensembl.gene.id")

# ## deal with duplicate gene names
# tpm.df.out <- tpm.df.w.gene.name %>%
#     group_by(gene.name) %>%
#     summarise_at(vars(!Row.names), sum)

## write output
# out.name <- file.path(args$outdir, sprintf("%s_kallisto_tpm_matrix_n%ix%i.tsv", args$project, nrow(tpm.df.out), ncol(tpm.df.out)-1))
# out.name <- file.path(args$outdir, sprintf("%s_kallisto_tpm_matrix_n%ix%i.tsv", args$project, nrow(tpm.df.out), ncol(tpm.df.out)-1))
# write.table(tpm.df.out, out.name, sep="\t", quote=F, row.names=F)

