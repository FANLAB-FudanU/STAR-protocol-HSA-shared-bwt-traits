library(tidyverse)
library(data.table)
library(argparse)

getwd()

### get arguments
parser <- ArgumentParser()
parser$add_argument("--fm_stats_per_locus_file", type = "character", required = TRUE, help='Output file of scripts/bin/summarize_fm_for_loci_file.sh script')
parser$add_argument("--out_prefix", type = "character", required = TRUE)
args <- parser$parse_args()

fm_stats_per_locus_file <- args$fm_stats_per_locus_file
out_prefix <- args$out_prefix
cat("fm_stats_per_locus_file is", fm_stats_per_locus_file, "\n")
cat("out_prefix is", out_prefix, "\n")

cat('Reading', fm_stats_per_locus_file, '\n')
fm_stats_per_locus <- fread(fm_stats_per_locus_file)
dim(fm_stats_per_locus)
head(fm_stats_per_locus)

# get num cs counts
num_cs_counts <- fm_stats_per_locus %>% 
    count(num_CSs, name='num_loci')
head(num_cs_counts)

out_num_cs_counts_file <- paste0(out_prefix, '.num_cs_counts.tsv')
write_tsv(num_cs_counts, out_num_cs_counts_file)

# get cs size counts
head(fm_stats_per_locus)

cs_sizes <- data.frame('cs_size' = as.numeric(unlist(strsplit(fm_stats_per_locus$CS_sizes, ','))))
dim(cs_sizes)
head(cs_sizes)
summary(cs_sizes$cs_size)

cs_size_counts <- cs_sizes %>% 
    count(cs_size, name='num_CSs')
dim(cs_size_counts)
head(cs_size_counts)

out_cs_size_counts_file <- paste0(out_prefix, '.cs_size_counts.tsv')
write_tsv(cs_size_counts, out_cs_size_counts_file)

# get summary of num cs and cs size
cat('Getting summary\n')
num_cs_smmry <- data.frame('num_cs' = unclass(summary(fm_stats_per_locus$num_CSs)), check.names = FALSE)
cs_size_smmry <- data.frame('cs_size' = unclass(summary(cs_sizes$cs_size)), check.names = FALSE)

smmry <- cbind(num_cs_smmry, cs_size_smmry)
smmry$num_cs <- round(smmry$num_cs, 1)
smmry$cs_size <- round(smmry$cs_size, 1)
smmry <- smmry %>% 
    mutate(stats = rownames(smmry), .before=1)
smmry

out_smmry_file <- paste0(out_prefix, '.summary.tsv')
write_tsv(smmry, out_smmry_file)

# Done
cat('Done. out_num_cs_counts_file is', out_num_cs_counts_file, '. out_cs_size_counts_file is', out_cs_size_counts_file, '\n')
cat('out_smmry_file is', out_smmry_file, '\n')



