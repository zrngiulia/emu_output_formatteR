# Libraries
require(ggplot2)
require(vegan)
require(paletteer)
require(dplyr)

# Function to create counts tables from emu standard output
emu_counts_taxid <- function(barcodeslist, path) {
  
  # Suffix for the path to the non-collapsed taxonomy tables files
  pt_suff1 <- "_rel-abundance.tsv"
  
  # Check format of inputs
  
  # Barcodeslist is a vector of characters containing barcode numbers
  # (e.g. `c("01", "02")`)
  if(!is.character(barcodeslist)) stop("barcodeslist wrong format")
  if(length(grep("barcode", barcodeslist)) != 0) stop("barcodeslist wrong format")
  
  # Path is valid
  if(!dir.exists(path)) stop("path does not exist")
  if(!file.exists(paste0(path, "/barcode", barcodeslist[1], pt_suff1)))
    stop("files not found in path")
  
  # Pre-processing
  barcodeslist <- paste0("barcode", barcodeslist)
  path <- paste0(path, "/")
  n <- length(barcodeslist)
  
  # Read abundance tables in a list
  all_barcodes <- list()
  for (i in 1:n)
    all_barcodes[[i]] <- read.table(header = T, sep = "\t",
                                    file = paste0(path, "/", barcodeslist[i], pt_suff1),
                                    stringsAsFactors=FALSE)
  # Create a table of counts
  taxid_vector <- c()
  for (i in 1:n)
    taxid_vector <- unique(c(taxid_vector, all_barcodes[[i]]$tax_id))
  taxid_counts <- data.frame("tax_id"=taxid_vector)
  colnames(taxid_counts) <- "tax_id"
  
  for (i in 1:n) taxid_counts <-
    merge(taxid_counts, all_barcodes[[i]][,c("tax_id", "estimated.counts")],
          by="tax_id", all.x = TRUE)
  colnames(taxid_counts) <- c("tax_id", barcodeslist)
  
  # Substitute NA counts with 0
  taxid_counts[is.na(taxid_counts)] <- 0
  
  # Create taxid_counts_info table
  taxid_counts_info <- all_barcodes[[1]][1,c(1,3:10)]
  taxid_counts_info <- taxid_counts_info[-1,]
  for (i in 1:n) taxid_counts_info <-
    merge(taxid_counts_info, all_barcodes[[i]][,c(1, 3:10)], all.x = TRUE,
          all.y = TRUE)
  # Add all barcodes counts
  for (i in 1:n) taxid_counts_info <-
    merge(taxid_counts_info, all_barcodes[[i]][,c(1, 14)], by="tax_id",
          all.x = TRUE)
  # Change colnames
  colnames(taxid_counts_info) <- c(colnames(taxid_counts_info)[1:9], barcodeslist)
  # Substitute NA counts with 0
  taxid_counts_info[is.na(taxid_counts_info)] <- 0
  
  # Create a taxa table only, with taxid as a first column
  taxid_info_only <- taxid_counts_info %>% select(!barcodeslist)
  
  # Create a taxid counts table with rounded values
  nums <- vapply(taxid_counts, is.numeric, FUN.VALUE = logical(1))
  counts_rounded <- taxid_counts
  counts_rounded[,nums] <- round(counts_rounded[,nums], digits = 0)
  
  return(list(taxid_counts, taxid_counts_info, taxid_info_only, counts_rounded))
}