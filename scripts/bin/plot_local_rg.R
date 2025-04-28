library(tidyverse)
library(data.table)
library(argparse)
library(RColorBrewer)

get_file_type <- function(filename){
    x <- unlist(strsplit(a, '\\.'))
    x[length(x)]
}

#' Plot heatmap across the whole genome
#' 
#' @param data A data.frame with four columns (chrom, start, end, value). No `chr` string in chrom column.
#' @param col A vector of color names used to generate colors.
#' @returns A plot.
#' @export 
#' @examples
#' genome_heatmap(df)
genome_heatmap <- function(
    data, 
    col = c("darkgreen", "yellow", "red"), 
    band = 3,
    width = 5,
    wh = 9, 
    ht = 6, 
    dpi = 300, 
    legend.cex = 1,
    main.cex = 1.5,
    main.font = 2,
    main = NULL,
    axis.cex = 1,
    axis.lwd = 1.5, # width of axis
    xticks.pos = 1,
    lower_limit = -1,
    upper_limit = 1,
    legend_length = 10,
    legend.y.intersp = 1,
    legend.x.intersp = 1, 
    save = FALSE, 
    use_RColorBrewer = FALSE, 
    palette = "Blues", 
    use_point_unit = FALSE, 
    axis.size = 7, # in point unit 
    main.size = 8, # in point unit 
    legend.size = 7, 
    filename = NULL
){
    n_regions <- nrow(data)
    data <- as.matrix(data)
    cat('Plotting heatmap of', n_regions, 'regions across the genome\n')
    ### order regions by chromosome and position
    order_index <- order(data[, 1], data[, 2])
    data <- data[order_index, ]

    ### get chrom and pos
    chrom <- data[, 1]
    start_pos <- data[, 2]
    stop_pos <- data[, 3]
    value <- data[, 4]

    ### set main
    if(is.null(main))   main <- "The heatmap of genomic regions"

    ### get chrom info
    max.chrom <- max(chrom)
    chrom.num <- unique(chrom)
    chorm.maxlen <- max(stop_pos)
    num_chroms <- length(chrom.num)

    ### get tick unit
    tick_unit <- ifelse(chorm.maxlen < 1e3, 1, ifelse(chorm.maxlen < 1e6, 1e3, 1e6))
    unit_label <- ifelse(tick_unit == 1, "bp", ifelse(tick_unit == 1e3, "Kb", "Mb"))

    ### get legend breaks
    min_value <- min(value)
    max_value <- max(value)

    if (is.null(lower_limit)) { lower_limit <- min_value }
    if (is.null(upper_limit)) { upper_limit <- max_value }

    # set number of units of legend
    if(upper_limit >1 & upper_limit <= legend_length) { legend_length <- upper_limit }

    legend_unit <- (upper_limit - lower_limit) / legend_length
    legend_breaks <- seq(lower_limit, upper_limit, by=legend_unit)

    ### get cols
    if (use_RColorBrewer) {
        cols <- brewer.pal(legend_length, palette)
    } else {
        cols <- colorRampPalette(col)(legend_length)
    }
    cols <- c('grey95', cols)

    # bin value
    value_bins <- cut(value, breaks=legend_breaks)

    # bin legend breaks
    legend_breaks_bins <- cut(legend_breaks, breaks=legend_breaks)

    # find cols of each value
    value_breaks <- legend_breaks[match(value_bins, legend_breaks_bins)]
    value_cols <- cols[match(value_breaks, legend_breaks)]

    ### get start pos, stop pos, max pos and colors of each chrom
    chrom_start_pos <- list()
    chrom_stop_pos <- list()
    chrom_region_cols <- list()
    chrom_max_pos <- c()

    for (chrom_idx in seq(1, num_chroms)) {
        chrom_start_pos[[chrom_idx]] <- start_pos[chrom == chrom.num[chrom_idx]]
        chrom_stop_pos[[chrom_idx]] <- stop_pos[chrom == chrom.num[chrom_idx]]
        chrom_region_cols[[chrom_idx]] <- value_cols[chrom == chrom.num[chrom_idx]]

        max_pos <- max(chrom_stop_pos[[chrom_idx]])
        chrom_max_pos <- c(chrom_max_pos, max_pos)

    }

    ### set canvas
    cat('Setting canvas\n')
    if (save) {
        file_type <- get_file_type(filename)
        if(file_type=="jpg"){ jpeg(filename, width=wh*dpi, height=ht*dpi, res=dpi, quality=100) }
        if(file_type=="pdf"){ jpeg(filename, width=wh, height=ht) }
        if(file_type=="tiff"){ jpeg(filename, width=wh*dpi, height=ht*dpi, res=dpi) }
        if(file_type=="png"){ jpeg(filename, width=wh*dpi, height=ht*dpi, res=dpi, bg=NA) }
    }

    if(use_point_unit){
        main_size <- main.size / 12
    } else {
        main_size <- main.cex
    }
    plot(NULL, xlim=c(0, chorm.maxlen + chorm.maxlen/10), ylim=c(0, length(chrom.num) * band + band), 
        main=main, cex.main=main_size, font.main=main.font, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")

    ### plot each chrom
    for (chrom_idx in seq(1, num_chroms)) {
        ybottom <- -width/5 - band * (chrom_idx - num_chroms - 1)
        ytop <- width/5 - band * (chrom_idx - num_chroms - 1)
        rect(xleft=0, ybottom=ybottom, xright=chrom_max_pos[chrom_idx], ytop=ytop, 
            col="grey95", border="grey80")

        region_start_pos <- chrom_start_pos[[chrom_idx]]
        region_stop_pos <- chrom_stop_pos[[chrom_idx]]
        region_cols <- chrom_region_cols[[chrom_idx]]

        # plot all regions
        rect(xleft=region_start_pos, ybottom=ybottom, xright=region_stop_pos, ytop=ytop, 
            col=region_cols, border='NA')
    }

    ### add chrom labels 
    cat('Adding chrom labels canvas\n')
    chrom_labels <- rev(chrom.num)
    if(use_point_unit){
        axis_size <- axis.size / 12
    } else {
        axis_size <- axis.cex
    }
    mtext(text=paste("Chr", chrom_labels, sep=""), 
        at=seq(band, num_chroms * band, band), side=2, las=2, font=1, 
        cex=axis_size, line=0.2, xpd=TRUE)

    ### add ticks
    xticks <- seq(0, chorm.maxlen / tick_unit, length=10)

    # round xticks
    if(round(xticks[2]) <= 10){
        xticks=seq(0, chorm.maxlen / tick_unit, round(xticks[2], 1))
    }else{
        xticks=seq(0, chorm.maxlen / tick_unit, round(xticks[2]))    
    }

    if((chorm.maxlen/tick_unit - max(xticks)) > 0.5*xticks[2]){
        xticks=c(xticks, round(chorm.maxlen / tick_unit))
    }

    # add xticks
    if(use_point_unit){
        axis_size <- axis.size / 12
    } else {
        axis_size <- axis.cex
    }
    axis(3, at=xticks*tick_unit, labels=paste(xticks, unit_label, sep=""), 
        font=1, lwd=axis.lwd, tck=0.01, padj=1.2, 
        cex.axis=axis_size, 
        mgp=c(3, xticks.pos, 0)
    )

    ### add legend
    cat('Adding legend\n')
    legend_labels <- legend_breaks
    legend_cols <- cols
    legend_x_pos <- chorm.maxlen + chorm.maxlen/50
    legend_y_pos <- -width/2.5 + band

    if(use_point_unit){
        legend_size <- legend.size / 12
    } else {
        legend_size <- legend.cex
    }
    legend(x=legend_x_pos, y=legend_y_pos, legend=legend_labels, title='', col=legend_cols, 
        pch=15, pt.cex=legend_size*3, cex=legend_size, 
        y.intersp=legend.y.intersp, x.intersp=legend.x.intersp, 
        xjust=0, yjust=0, xpd=TRUE, bty='n'
    )

    if (save) {
        dev.off()
        cat('Heatmap is saved to', filename, '\n')
    }
}

getwd()

### get arguments
parser <- ArgumentParser()
parser$add_argument("--sig_local_rg_file", help="File with significant local rg. Columns should be: locus, chr, start, stop, n.snps, rho, r2, p", type = "character", required = TRUE)
parser$add_argument("--title", help="Title of output file", type = "character", required = FALSE, default='Distribution of local rg across the genome')
parser$add_argument("--unit_mm", help="Whether to use millimetre (mm) as unit of output figure size.", type = "logical", required = FALSE, default=TRUE)
parser$add_argument("--width", help="Width of output figure in unit mm (or inches if unit_mm is FALSE)", type = "numeric", required = FALSE, default=210)
parser$add_argument("--height", help="Height of output figure in unit mm (or inches if unit_mm is FALSE)", type = "numeric", required = FALSE, default=180)
parser$add_argument("--out_file", help="Output PDF file", type = "character", required = TRUE)
args <- parser$parse_args()

sig_local_rg_file <- args$sig_local_rg_file
title <- args$title
unit_mm <- args$unit_mm
width <- args$width
height <- args$height
out_file <- args$out_file

cat("sig_local_rg_file is", sig_local_rg_file, "\n")
cat("title is", title, "\n")
cat("unit_mm is", unit_mm, "\n")
cat("width is", width, "\n")
cat("height is", height, "\n")
cat("out_file is", out_file, "\n")

# tmp
tmp <- function(){

    sig_local_rg_file <- '02.local_rg/03.lava_outputs/02.HT_BMR_local_rg.sig_bivar.tsv'
    title <- 'Heatmap'
    unit_mm <- TRUE
    width <- 210
    height <- 180
    out_file <- 'tmp/a.pdf'

}

cat('Reading', sig_local_rg_file, '\n')
sig_local_rg <- fread(sig_local_rg_file)
dim(sig_local_rg)
head(sig_local_rg)

## plot using genome_heatmap
cat('Plotting local rg using genome_heatmap\n')
local_rg <- sig_local_rg[, c(2,3,4,9)]

if(unit_mm){
    width <- width / 25.4
    height <- height / 25.4
}

pdf(out_file, width = width, height = height)
genome_heatmap(
    local_rg, use_RColorBrewer = TRUE, palette = 'Blues', legend_length=5, lower_limit = 0, upper_limit = 1, 
    use_point_unit = TRUE, axis.size = 7, main.size = 8, legend.size = 7, 
    axis.lwd = 1, main = title
)

dev.off()

cat('Done. out_file is', out_file, '\n')





