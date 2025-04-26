#!/usr/bin/env Rscript

## Main: perform fine-mapping using susie_rss function 
# inputs: beta, se, LD matrix, sample size

# copied and modified from UKBB_FinnGen_BBJ_finemapping_MasahiroKanai/finemapping-pipeline/R/run_susieR.R

options(stringsAsFactors = F)

library(argparse)
library(data.table)
library(dplyr)
library(stringr)
library(susieR)

load_R <- function(path, num_snps, dominant = FALSE, triangular_ld_matrix = FALSE) {
    if (tools::file_ext(path) == "ld") {
        R <- fread(path, header = FALSE)
    } else if (tools::file_ext(path) %in% c("gz", "bgz")) {
        R <- fread(cmd = paste("zcat", path), header = FALSE)
    } else {
        stop("LD matrix file extension is not supported")
    }

    R <- as.matrix(R)
    if (num_snps != dim(R)[1]) {
        stop("LD matrix size doesn't match the number of SNPs")
    }
    if (dominant) {
        warning("--dominant is specified. LD matrix is squared.")
        R <- R**2
    }

    if (triangular_ld_matrix) {
        warning("--triangular-ld-matrix is specified. LD matrix is re-densified.")
        R[lower.tri(R)] <- t(R)[lower.tri(R)]
    }

    colnames(R) <- seq(1, num_snps)
    rownames(R) <- seq(1, num_snps)

    return(R)
}

compute_yty <- function(beta, se, p, R, n, k) {
    # beta and se should be standarized
    beta_s <- beta * sqrt(2 * p * (1 - p))
    se_s <- se * sqrt(2 * p * (1 - p))

    # Y'Y =  Bj^2 (Xj'Xj) + Var(Bj)(Xj'Xj)(N - k)
    XjtXj <- (n - 1) * diag(R)
    yty <- beta_s**2 * XjtXj + se_s**2 * XjtXj * (n - k)

    return(median(yty))
}

# cf. https://github.com/stephenslab/susieR/blob/master/R/summary.susie.R
summarize.susie.cs <- function(object, orig_vars, R, ..., low_purity_threshold = 0.5) {
    if (is.null(object$sets)) {
        stop("Cannot summarize SuSiE object because credible set information is not available")
    }

    # variable_prob 即每个SNP的PIP
    variables <- data.frame(cbind(1:length(object$pip), object$pip, -1, NA, NA, NA))
    colnames(variables) <- c("variable", "variable_prob", "cs", "cs_specific_prob", "low_purity", "lead_r2")
    rownames(variables) <- NULL
    added_vars <- c()
    if (object$null_index > 0) variables <- variables[-object$null_index, ]
    if (!is.null(object$sets$cs)) {
        cs <- data.frame(matrix(NA, length(object$sets$cs), 5))
        colnames(cs) <- c("cs", "cs_log10bf", "cs_avg_r2", "cs_min_r2", "variable")
        for (i in 1:length(object$sets$cs)) {
            if (any(object$sets$cs[[i]] %in% added_vars)) {
                print(
                    sprintf("Skipping cs %d as there is an overlap between variants in this cs and previous credible sets", i)
                )
                print("Removed cs variants:")
                print(orig_vars[object$sets$cs[[i]], ], max = length(object$sets$cs[[i]]))
                next
            } else {
                added_vars <- append(added_vars, object$sets$cs[[i]])
            }
            in_cs_idx <- which(variables$variable %in% object$sets$cs[[i]])
            variables$cs[in_cs_idx] <- object$sets$cs_index[[i]]
            variables[in_cs_idx, "cs_specific_prob"] <- object$alpha[object$sets$cs_index[[i]], object$sets$cs[[i]]]
            variables$low_purity[in_cs_idx] <- object$sets$purity$min.abs.corr[i] < low_purity_threshold
            lead_pip_idx <- in_cs_idx[which.max(variables$variable_prob[in_cs_idx])]
            variables$lead_r2[in_cs_idx] <- R[lead_pip_idx, in_cs_idx]^2

            cs$cs[i] <- object$sets$cs_index[[i]]
            cs$cs_log10bf[i] <- log10(exp(object$lbf[cs$cs[i]]))
            cs$cs_avg_r2[i] <- object$sets$purity$mean.abs.corr[i]^2
            cs$cs_min_r2[i] <- object$sets$purity$min.abs.corr[i]^2
            cs$low_purity <- object$sets$purity$min.abs.corr[i] < low_purity_threshold
            cs$variable[i] <- paste(object$sets$cs[[i]], collapse = ",")
        }
        variables <- variables[order(variables$variable_prob, decreasing = T), ]
    } else {
        cs <- NULL
    }
    return(list(vars = variables, cs = na.omit(cs)))
}

susie_ss_wrapper <- function(df, R, n, L, var_y = 1.0, prior_weights = NULL, estimate_residual_variance = FALSE, min_abs_corr = 0.0, low_purity_threshold = 0.5) {
    beta <- df$beta
    se <- df$se
    fitted_bhat <- susie_rss(
        bhat = beta,
        shat = se,
        R = R,
        n = n,
        var_y = var_y,
        L = L,
        prior_weights = prior_weights,
        scaled_prior_variance = 0.1,
        estimate_residual_variance = estimate_residual_variance,
        estimate_prior_variance = TRUE,
        standardize = TRUE,
        check_input = FALSE,
        min_abs_corr = min_abs_corr
    )

    cs_summary <- summarize.susie.cs(fitted_bhat, df, R, low_purity_threshold = low_purity_threshold)
    variables <-
        cs_summary$vars %>%
        rename(prob = variable_prob) %>% # prob即原始的PIP
        arrange(variable) %>%
        mutate(
            mean = susie_get_posterior_mean(fitted_bhat),
            sd = susie_get_posterior_sd(fitted_bhat)
        )
    cs <- cs_summary$cs

    return(list(
        susie_obj = fitted_bhat,
        variables = variables,
        cs = cs
    ))
}

parser <- ArgumentParser()

parser$add_argument("--z", "-z", help="", type = "character", required = TRUE)
parser$add_argument("--ld", type = "character")
parser$add_argument("--out", type = "character")
parser$add_argument("--snp", type = "character")
parser$add_argument("--cred", type = "character")
parser$add_argument("--log", type = "character")
parser$add_argument("--susie-obj", type = "character")
parser$add_argument("--pheno", type = "character")
parser$add_argument("--n-samples", "-n", type = "integer", required = TRUE)
parser$add_argument("--var-y", type = "double", default = 1.0)
parser$add_argument("--L", type = "integer", default = 10)
parser$add_argument("--yty", type = "double")
parser$add_argument("--in-sample-ld", action = 'store_true', help='If sepcified, then in LD is in-sample LD and estimate_residual_variance is TRUE in susie_rss')
parser$add_argument("--min-cs-corr", default = 0.5, type = "double")
parser$add_argument("--low-purity-threshold", default = 0.5, type = "double")
parser$add_argument("--compute-yty", action = "store_true")
parser$add_argument("--n-covariates", "-k", type = "integer")
parser$add_argument("--prior-weights", type = "character")
parser$add_argument("--write-alpha", action = "store_true")
parser$add_argument("--write-single-effect", action = "store_true")
parser$add_argument("--save-susie-obj", action = "store_true")
parser$add_argument("--dominant", action = "store_true")
parser$add_argument("--triangular-ld-matrix", action = "store_true", help = "Use triangular LD matrix")

# run in command line
args <- parser$parse_args()

if (is.null(args$out)) {
    args$out <- "tmp"
}
if (is.null(args$snp)) {
    args$snp <- paste0(args$out, ".susie.snp")
    args$snp_in_cs <- paste0(args$out, ".susie.snp_in_cs.tsv")
}
if (is.null(args$cred)) {
    args$cred <- paste0(args$out, ".susie.cred")
}
if (is.null(args$log)) {
    args$log <- paste0(args$out, ".susie.log")
}
if (is.null(args$susie_obj)) {
    args$susie_obj <- paste0(args$out, ".susie.rds")
}

print(args)

if (is.null(args$out) & any(sapply(list(args$snp, args$cred, args$log), is.null))) {
    stop("Either --out or all of --snp, --cred, and --log should be specified.")
}
if ((args$var_y != 1.0 | !is.null(args$yty)) & args$compute_yty) {
    warning("--compute-yty will override the specified --var-y/--yty.")
}
if (args$compute_yty & is.null(args$n_covariates)) {
    stop("--n-covariates is required for --compute-yty.")
}

print("Analysis started")

df <- fread(args$z)
n <- args$n_samples
L <- args$L

if (!is.null(args$prior_weights)) {
	prior_weights <- fread(args$prior_weights)
	df$v <- str_c(df$chromosome, df$position, df$allele1, df$allele2, sep = ":")
	prior_weights <- prior_weights$mean_binpred[match(df$v, prior_weights$SNP)]
	if (is.null(prior_weights)) {
		stop("'min_binpred' is missing (--prior-weights).")
	}
	if (any(is.na(prior_weights))) {
		stop("All SNPs should have prior weights (--prior-weights).")
	}
} else {
	prior_weights <- NULL
}

num_snps <- nrow(df)
if (!is.null(args$ld)) {
	R <- load_R(args$ld, num_snps, dominant = args$dominant, triangular_ld_matrix = args$triangular_ld_matrix)
} else {
	warning("No LD file is provided.")
	# require(Matrix)
	# R <- Matirx::.sparseDiagonal(nrow(df))
	R <- diag(1, nrow(df), nrow(df))
}

if (args$compute_yty) {
	yty <- compute_yty(df$beta, df$se, df$maf, R, n, args$n_covariates)

	# we scale var_y accordingly, instead of tweaking susie's code
	# https://github.com/stephenslab/susieR/blob/master/R/susie_ss.R#L334
	# susie_ss(XtX = XtX, Xty = Xty, yty = var_y * (n-1), n = n, ...)
	var_y <- yty / (n - 1)
} else if (!is.null(args$yty)) {
	var_y <- args$yty / (n - 1)
} else if (!is.null(args$pheno)) {
	pheno <- fread(args$pheno)
	var_y <- var(pheno[, 0])
} else if (!is.null(args$var_y)) {
	var_y <- args$var_y
} else {
	var_y <- NULL
}

# set estimate_residual_variance
if (args$in_sample_ld) {
    estimate_residual_variance = TRUE
} else {
    estimate_residual_variance = FALSE
}

cat('Running susie\n')
res <- susie_ss_wrapper(
	df, R, n, L, var_y, prior_weights,
    estimate_residual_variance = estimate_residual_variance, 
	min_abs_corr = args$min_cs_corr, low_purity_threshold = args$low_purity_threshold
)

susie_obj <- res$susie_obj
variables <- cbind(df, res$variables[c("mean", "sd", "prob", "cs", "cs_specific_prob", "low_purity", "lead_r2")])
cs <- res$cs
if (!is.null(cs)) {
	cs <- cs %>%
		mutate(cs_size = unlist(lapply(str_split(variable, ","), length))) %>%
		select(-variable)
}

cat('Saving RDS to', args$susie_obj, '\n')
if (args$save_susie_obj) {
    saveRDS(susie_obj, file = args$susie_obj)
}

# get snps in cs
snps_in_cs <- variables %>% 
    filter(cs != -1)

if(nrow(snps_in_cs) == 0){
    cat('Warning! No SNPs in credible set!\n')
} else {
    write.table(snps_in_cs, args$snp_in_cs, sep = "\t", row.names = F, quote = F)
    write.table(cs, args$cred, sep = "\t", row.names = F, quote = F)
}

write.table(variables, args$snp, sep = "\t", row.names = F, quote = F)

# write snp pos
snps_pos_in_cs_without_chr <- snps_in_cs %>% 
    select(chromosome, position) %>% 
    mutate(end = position)

snps_pos_in_cs_with_chr <- snps_pos_in_cs_without_chr %>% 
    mutate(chromosome = paste0('chr', chromosome))

out_snps_pos_in_cs_with_chr_file <- paste0(args$out, '.susie.snp_in_cs.pos.with_chr.tsv')
out_snps_pos_in_cs_without_chr_file <- paste0(args$out, '.susie.snp_in_cs.pos.without_chr.tsv')

write.table(snps_pos_in_cs_without_chr, out_snps_pos_in_cs_without_chr_file, sep = "\t", row.names = F, quote = F, col.names = F)
write.table(snps_pos_in_cs_with_chr, out_snps_pos_in_cs_with_chr_file, sep = "\t", row.names = F, quote = F, col.names = F)

# done
cat('out file is', args$snp, args$snp_in_cs, args$cred, '\n')
cat('All done\n')


