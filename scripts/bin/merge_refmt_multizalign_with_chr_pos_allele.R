library(tidyverse)
library(data.table)
library(argparse)

### get arguments
parser <- ArgumentParser()
parser$add_argument("--snp_chr_pos_file", type = "character", required = TRUE, help='1-based position')
parser$add_argument("--refmtted_multizalign_file", type = "character", required = TRUE, help='0-based position')
parser$add_argument("--out_file", type = "character", required = TRUE)
args <- parser$parse_args()

snp_chr_pos_file <- args$snp_chr_pos_file
refmtted_multizalign_file <- args$refmtted_multizalign_file
out_file <- args$out_file
cat("snp_chr_pos_file is", snp_chr_pos_file, "\n")
cat("refmtted_multizalign_file is", refmtted_multizalign_file, "\n")
cat("out_file is", out_file, "\n")

cat('Reading file:', snp_chr_pos_file, refmtted_multizalign_file, '\n')
snp_chr_pos <- fread(snp_chr_pos_file)
head(snp_chr_pos)

refmtted_multizalign <- fread(refmtted_multizalign_file)
head(refmtted_multizalign)
refmtted_multizalign <- refmtted_multizalign %>% 
    mutate(human.pos = human.pos + 1) # convert 0-based to 1-based

# join
head(snp_chr_pos)
cat('Merging\n')
merged_primates_allele <- snp_chr_pos %>% 
    left_join(refmtted_multizalign, by=c('chrom' = 'human.chrom', 'start'='human.pos'))

cat('shape of merged data:\n')
print(dim(merged_primates_allele))

# write
write_tsv(merged_primates_allele, out_file)
cat('Done. out_file is', out_file, '\n')


