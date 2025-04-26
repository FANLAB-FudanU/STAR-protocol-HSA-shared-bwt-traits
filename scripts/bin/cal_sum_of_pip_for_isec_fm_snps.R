library(data.table)
library(argparse)

getwd()

### get arguments
parser <- ArgumentParser()
parser$add_argument("--isec_fm_snps_file", type = "character", required = TRUE, help='Must contain the following columns: chrom, locus_start, locus_stop, cs-HT, cs_specific_prob-HT, cs-BMR, cs_specific_prob-BMR')
parser$add_argument("--out_dir", type = "character", required = TRUE)
args <- parser$parse_args()

isec_fm_snps_file <- args$isec_fm_snps_file
out_dir <- args$out_dir
cat("isec_fm_snps_file is", isec_fm_snps_file, "\n")
cat("out_dir is", out_dir, "\n")

if (! dir.exists(out_dir)) {
    dir.create(out_dir)
}

cat('Reading', isec_fm_snps_file, '\n')
isec_fm_snps <- fread(isec_fm_snps_file)
dim(isec_fm_snps)

head(isec_fm_snps)

cat('Calculating sum of pips\n')
pip_sums_HT <- isec_fm_snps[, .(pip_sum = sum(`cs_specific_prob-HT`)), by=.(chrom, locus_start, locus_stop, `cs-HT`)]
dim(pip_sums_HT)
head(pip_sums_HT)

pip_sums_BMR <- isec_fm_snps[, .(pip_sum = sum(`cs_specific_prob-BMR`)), by=.(chrom, locus_start, locus_stop, `cs-BMR`)]
dim(pip_sums_BMR)
head(pip_sums_BMR)

# write into file
out_pip_sums_HT_file <- paste0(out_dir, '/PIP_sums.HT.tsv')
out_pip_sums_BMR_file <- paste0(out_dir, '/PIP_sums.BMR.tsv')

fwrite(pip_sums_HT, out_pip_sums_HT_file, sep = '\t', quote = F, row.names = F)
fwrite(pip_sums_BMR, out_pip_sums_BMR_file, sep = '\t', quote = F, row.names = F)

## get cs with pip sum greater than 0.5 or 0.9
cs_with_pip_gt05_HT <- pip_sums_HT[pip_sum >= 0.5]
cs_with_pip_gt09_HT <- pip_sums_HT[pip_sum >= 0.9]

cs_with_pip_gt05_BMR <- pip_sums_BMR[pip_sum >= 0.5]
cs_with_pip_gt09_BMR <- pip_sums_BMR[pip_sum >= 0.9]

out_cs_with_pip_gt05_HT_file <- paste0(out_dir, '/cs_with_pip_gt05.HT.tsv')
out_cs_with_pip_gt09_HT_file <- paste0(out_dir, '/cs_with_pip_gt09.HT.tsv')

out_cs_with_pip_gt05_BMR_file <- paste0(out_dir, '/cs_with_pip_gt05.BMR.tsv')
out_cs_with_pip_gt09_BMR_file <- paste0(out_dir, '/cs_with_pip_gt09.BMR.tsv')

fwrite(cs_with_pip_gt05_HT, out_cs_with_pip_gt05_HT_file, sep = '\t', quote = F, row.names = F)
fwrite(cs_with_pip_gt09_HT, out_cs_with_pip_gt09_HT_file, sep = '\t', quote = F, row.names = F)

fwrite(cs_with_pip_gt05_BMR, out_cs_with_pip_gt05_BMR_file, sep = '\t', quote = F, row.names = F)
fwrite(cs_with_pip_gt09_BMR, out_cs_with_pip_gt09_BMR_file, sep = '\t', quote = F, row.names = F)

cat('Done\n')
cat('out dir is', out_dir, '\n')
cat('out_pip_sums_HT_file is', out_pip_sums_HT_file, '\n')
cat('out_pip_sums_BMR_file is', out_pip_sums_BMR_file, '\n')


