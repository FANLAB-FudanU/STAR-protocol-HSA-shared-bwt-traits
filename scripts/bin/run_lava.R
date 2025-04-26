library(tidyverse)
library(data.table)
library(argparse)
library(LAVA)
library(parallel)

getwd()

# Rscript $script --input_info_file $input_info_file --sample_overlap_file $sample_overlap_file --ref_prefix $ref_prefix --loci_file $loci_file --threads $threads --out_prefix $out_prefix

### get arguments
parser <- ArgumentParser()
parser$add_argument("--input_info_file", type = "character", required = TRUE)
parser$add_argument("--sample_overlap_file", type = "character", required = TRUE)
parser$add_argument("--ref_bfile", type = "character", required = TRUE)
parser$add_argument("--loci_file", type = "character", required = TRUE)
parser$add_argument("--threads", type = "character", required = TRUE)
parser$add_argument("--out_prefix", type = "character", required = TRUE)
args <- parser$parse_args()

input_info_file <- args$input_info_file
sample_overlap_file <- args$sample_overlap_file
ref_bfile <- args$ref_bfile
loci_file <- args$loci_file
threads <- args$threads
out_prefix <- args$out_prefix

cat("input_info_file is", input_info_file, "\n")
cat("sample_overlap_file is", sample_overlap_file, "\n")
cat("ref_bfile is", ref_bfile, "\n")
cat("loci_file is", loci_file, "\n")
cat("threads is", threads, "\n")
cat("out_prefix is", out_prefix, "\n")

### Read in summary statistics and related info
cat('Reading summary statistics and reference file\n')
input = process.input(input.info.file=input_info_file, sample.overlap.file=sample_overlap_file, ref.prefix=ref_bfile)

### Read loci
loci <- read.loci(loci_file)
num_loci <- nrow(loci)
cat('Read', num_loci, 'loci\n')

### cal univ threshold
univ_pthresh <- 0.05 / num_loci
cat(paste0('Set P threshold of univariate test to 0.05/', num_loci, ' (', univ_pthresh, ')', '\n'))

### parallel for each locus
cat(paste("Starting LAVA analysis for", num_loci,"loci\n"))

u = b = NULL
start_time <- Sys.time() 
out <- mclapply(seq(1, num_loci), mc.cores=threads, function(idx){
    cat(paste0("Processing ", idx, " / ", num_loci, " locus\n"))

	# process locus
	locus <- process.locus(loci[idx, ], input)

	# in some cases the locus object cannot be created due to e.g too few SNPs or negative variances in all analysed phenotypes, hence this check
	if (!is.null(locus)) {
		# extract general locus info for output
		loc.info <- data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)

		# run univ & bivar analysis functions & store results
		ub <- run.univ.bivar(locus, univ.thresh=univ_pthresh)
		u <- cbind(loc.info, ub$univ)
		if (!is.null(ub$bivar)) {
			b <- cbind(loc.info, ub$bivar)
		}
	}
	return(list(univ=u, bivar=b))
})
end_time <- Sys.time()
end_time - start_time

### merge results of all loci
cat('Mering univ and bivar tests\n')
univ_df <- do.call(rbind, lapply(out, "[[","univ"))
bivar_df <- do.call(rbind, lapply(out, "[[","bivar"))

# add pthresh
num_bivar_tests <- nrow(bivar_df)
bivar_pthresh <- 0.05 / num_bivar_tests
cat(paste0('Using P threshold of bivariate test of 0.05/', num_bivar_tests, ' (', bivar_pthresh, ')', '\n'))

bivar_df['bivar_pthreshold'] <- bivar_pthresh

sig_bivar_df <- bivar_df %>% 
	filter(p <= bivar_pthreshold)
dim(sig_bivar_df)

out_univ_file <- paste0(out_prefix, '.univ.tsv')
out_bivar_file <- paste0(out_prefix, '.bivar.tsv')
out_sig_bivar_file <- paste0(out_prefix, '.sig_bivar.tsv')

write.table(univ_df, out_univ_file, row.names=F, quote=F, sep='\t')
write.table(bivar_df, out_bivar_file, row.names=F, quote=F, sep='\t')
write.table(sig_bivar_df, out_sig_bivar_file, row.names=F, quote=F, sep='\t')

cat('Done. Univariate tests are written into', out_univ_file, '\n')
cat('Bivariate tests are written into', out_bivar_file, '\n')
cat('Significant bivariate tests are written into', out_sig_bivar_file, '\n')


