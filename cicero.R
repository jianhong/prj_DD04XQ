# conda activate R4.12
library(cicero)
library(monocle3)

# You can substitute the data path below to your scATAC-seq data.
data_folder <- "outs/filtered_peak_bc_matrix"

# Create a folder to save results
output_folder <- "cicero_output"
dir.create(output_folder)

# Read in matrix data using the Matrix package
indata <- Matrix::readMM(paste0(data_folder, "/matrix.mtx")) 
# Binarize the matrix
indata@x[indata@x > 0] <- 1

# Format cell info
cellinfo <- read.table(paste0(data_folder, "/barcodes.tsv"))
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

# Format peak info
peakinfo <- read.table(paste0(data_folder, "/peaks.bed"))
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# Make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))

input_cds <- monocle3::detect_genes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

# Visualize peak_count_per_cell
pdf(file.path(output_folder, "hist.plot.pdf"))
hist(Matrix::colSums(exprs(input_cds)))
dev.off()

# Filter cells by peak_count
# Please set an appropriate threshold values according to your data 
max_count <-  15000
min_count <- 2000
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) >= min_count] 
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) <= max_count] 

# Data preprocessing
set.seed(2023)

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP
detach("package:monocle3",unload = TRUE)
library(monocle)
make_cicero_cds <- function (cds, reduced_coordinates, k = 50, summary_stats = NULL, 
          size_factor_normalize = TRUE, silent = FALSE) {
    #assertthat::assert_that(is(cds, "CellDataSet"))
    assertthat::assert_that(is.data.frame(reduced_coordinates) | 
                                is.matrix(reduced_coordinates))
    assertthat::assert_that(assertthat::are_equal(nrow(reduced_coordinates), 
                                                  nrow(pData(cds))))
    assertthat::assert_that(setequal(row.names(reduced_coordinates), 
                                     colnames(cds)))
    assertthat::assert_that(assertthat::is.count(k) & k > 1)
    assertthat::assert_that(is.character(summary_stats) | is.null(summary_stats))
    if (!is.null(summary_stats)) {
        assertthat::assert_that(all(summary_stats %in% names(pData(cds))), 
                                msg = paste("One of your summary_stats is missing", 
                                            "from your pData table. Either add a", "column with the name in", 
                                            "summary_stats, or remove the name", "from the summary_stats parameter.", 
                                            collapse = " "))
        assertthat::assert_that(sum(vapply(summary_stats, function(x) {
            !(is(pData(cds)[, x], "numeric") | is(pData(cds)[, 
                                                             x], "integer"))
        }, 1)) == 0, msg = paste("All columns in summary_stats must be", 
                                 "of class numeric or integer.", collapse = " "))
    }
    assertthat::assert_that(is.logical(size_factor_normalize))
    assertthat::assert_that(is.logical(silent))
    reduced_coordinates <- as.data.frame(reduced_coordinates)
    reduced_coordinates <- reduced_coordinates[colnames(cds), 
    ]
    nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
    row.names(nn_map) <- row.names(reduced_coordinates)
    nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
    good_choices <- seq_len(nrow(nn_map))
    choice <- sample(seq_len(length(good_choices)), size = 1, 
                     replace = FALSE)
    chosen <- good_choices[choice]
    good_choices <- good_choices[good_choices != good_choices[choice]]
    it <- 0
    k2 <- k * 2
    get_shared <- function(other, this_choice) {
        k2 - length(union(cell_sample[other, ], this_choice))
    }
    while (length(good_choices) > 0 & it < 5000) {
        it <- it + 1
        choice <- sample(seq_len(length(good_choices)), size = 1, 
                         replace = FALSE)
        new_chosen <- c(chosen, good_choices[choice])
        good_choices <- good_choices[good_choices != good_choices[choice]]
        cell_sample <- nn_map[new_chosen, ]
        others <- seq_len(nrow(cell_sample) - 1)
        this_choice <- cell_sample[nrow(cell_sample), ]
        shared <- sapply(others, get_shared, this_choice = this_choice)
        if (max(shared) < 0.9 * k) {
            chosen <- new_chosen
        }
    }
    cell_sample <- nn_map[chosen, ]
    if (!silent) {
        combs <- combn(nrow(cell_sample), 2)
        shared <- apply(combs, 2, function(x) {
            k2 - length(unique(as.vector(cell_sample[x, ])))
        })
        message(paste0("Overlap QC metrics:\nCells per bin: ", 
                       k, "\nMaximum shared cells bin-bin: ", max(shared), 
                       "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ", 
                       median(shared)))
        if (mean(shared)/k > 0.1) 
            warning("On average, more than 10% of cells are shared between paired bins.")
    }
    exprs_old <- exprs(cds)
    mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in% 
                       cell_sample[x, , drop = FALSE])
    mask <- Matrix::Matrix(mask)
    new_exprs <- exprs_old %*% mask
    new_exprs <- Matrix::t(new_exprs)
    new_exprs <- as.matrix(new_exprs)
    pdata <- pData(cds)
    new_pcols <- "agg_cell"
    if (!is.null(summary_stats)) {
        new_pcols <- c(new_pcols, paste0("mean_", summary_stats))
    }
    new_pdata <- plyr::adply(cell_sample, 1, function(x) {
        sub <- pdata[x, ]
        df_l <- list()
        df_l["temp"] <- 1
        for (att in summary_stats) {
            df_l[paste0("mean_", att)] <- mean(sub[, att])
        }
        data.frame(df_l)
    })
    new_pdata$agg_cell <- paste("agg", chosen, sep = "")
    new_pdata <- new_pdata[, new_pcols, drop = FALSE]
    row.names(new_pdata) <- new_pdata$agg_cell
    row.names(new_exprs) <- new_pdata$agg_cell
    new_exprs <- as.matrix(t(new_exprs))
    fdf <- fData(cds)
    new_pdata$temp <- NULL
    fd <- new("AnnotatedDataFrame", data = as.data.frame(fdf))
    pd <- new("AnnotatedDataFrame", data = as.data.frame(new_pdata))
    cicero_cds <- suppressWarnings(newCellDataSet(new_exprs, 
                                                  phenoData = pd, featureData = fd, expressionFamily = negbinomial.size(), 
                                                  lowerDetectionLimit = 0))
    cicero_cds <- monocle::detectGenes(cicero_cds, min_expr = 0.1)
    cicero_cds <- BiocGenerics::estimateSizeFactors(cicero_cds)
    if (any(!c("chr", "bp1", "bp2") %in% names(fData(cicero_cds)))) {
        fData(cicero_cds)$chr <- NULL
        fData(cicero_cds)$bp1 <- NULL
        fData(cicero_cds)$bp2 <- NULL
        fData(cicero_cds) <- cbind(fData(cicero_cds), df_for_coords(row.names(fData(cicero_cds))))
    }
    if (size_factor_normalize) {
        Biobase::exprs(cicero_cds) <- t(t(Biobase::exprs(cicero_cds))/Biobase::pData(cicero_cds)$Size_Factor)
    }
    cicero_cds
}

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Save cds object if you want
saveRDS(cicero_cds, paste0(output_folder, "/cicero_cds.Rds"))


# !!Please make sure that the reference genome information below matches the reference genome of your scATAC-seq data.

# If your scATAC-seq uses mm10 reference genome, you can read chromosome length file with the following command.
download.file(url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes",
              destfile = "./hg38_chromosome_length.txt")
chromosome_length <- read.table("./hg38_chromosome_length.txt")

# For mm9 genome, you can use the following command.
#data("mouse.mm9.genome")
#chromosome_length <- mouse.mm9.genome

# For hg19 genome, you can use the following command.
#data("human.hg19.genome")
#chromosome_length <- mhuman.hg19.genome

# Run the main function
conns <- run_cicero(cicero_cds, chromosome_length) # Takes a few minutes to run

# Check results
head(conns)

# Save results if you want
saveRDS(conns, paste0(output_folder, "/cicero_connections.Rds"))

all_peaks <- row.names(input_cds)
write.csv(x = all_peaks, file = paste0(output_folder, "/all_peaks.csv"))
write.csv(x = conns, file = paste0(output_folder, "/cicero_connections.csv"))

