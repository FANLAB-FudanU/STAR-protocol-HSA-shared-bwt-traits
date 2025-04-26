library(tidyverse)
library(data.table)
library(argparse)

getwd()

### get arguments
parser <- ArgumentParser()
parser$add_argument("--rg_summary", type = "character", required = TRUE, help='rg_summary: outputs by script: /public/home/fan_lab/wangjie/PublicDB/GWAS/UKBB_YangLab/04.genetic_correlation/scripts/bin/cal_pair_wise_rg.v2.sh')
parser$add_argument("--out_file", type = "character", required = TRUE)
args <- parser$parse_args()

rg_summary <- args$rg_summary
out_file <- args$out_file
cat("rg_summary is", rg_summary, "\n")
cat("out_file is", out_file, "\n")

cat('Reading', rg_summary, '\n')
scor <- read.table(rg_summary, header=T, sep='\t')              # read in
scor <- scor[,c("p1","p2","gcov_int")]             # retain key headers
scor$p1 <- as.character(scor$p1)
scor$p2 <- as.character(scor$p2)

head(scor)

### 这一步目的是将p1和p2中的文件名后缀去掉，由于我的文件不包括这些信息，因此不需要运行这一步
# scor$p1 = gsub("_munge.sumstats.gz","",scor$p1)   # assuming the munged files have format [phenotype]_munge.sumstats.gz
# scor$p2 = gsub("_munge.sumstats.gz","",scor$p2)   # (adapt as necessary)

phen <- unique(c(scor$p1, scor$p2))                  # obtain list of all phenotypes (assuming all combinations have been analysed)
n <- length(phen)

mat <- matrix(NA, n, n)                    # create matrix
rownames(mat) = colnames(mat) = phen    # set col/rownames

cat('Creating sample overlap file\n')
for (idx in seq(1, nrow(scor))) {
    row <- scor[idx, ]
    p1 <- row$p1
    p2 <- row$p2
    gcov_int <- row$gcov_int

    mat[p1, p2] <- gcov_int
}

diag(mat) <- 1

## 将下三角的数组设置为上三角的数值
mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

mat_std <- round(cov2cor(mat),5) # standardise

write.table(mat_std, out_file, quote=F)   # save

cat('Done. out file is', out_file, '\n')

