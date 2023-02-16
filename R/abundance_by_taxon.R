# Libraries
require(ggplot2)
require(vegan)
require(paletteer)
require(dplyr)

#################################################################################
# Function to create abundance tables from collapsed emu output

emu_abund_tables <- function(barcodeslist, path,
                             taxon=c('species', 'genus', 'family', 'order', 'class', 
                                     'phylum', 'clade', 'superkingdom')) {
  
  # Suffix for the path to the collapsed taxonomy tables files
  pt_suff <- paste0("_rel-abundance-", taxon, ".tsv")
  
  # Check format of inputs
  
  # Barcodeslist is a vector of characters containing barcode numbers
  # (e.g. `c("01", "02")`)
  if(!is.character(barcodeslist)) stop("barcodeslist wrong format")
  if(length(grep("barcode", barcodeslist)) != 0) stop("barcodeslist wrong format")
  
  # Path is valid
  if(!dir.exists(path)) stop("path does not exist")
  if(!file.exists(paste0(path, "/", taxon, "/barcode", barcodeslist[1], pt_suff))) 
    stop("files not found in path")
  
  # Taxon is one of the accepted values
  if (!taxon %in% c('species', 'genus', 'family', 'order', 'class', 
                    'phylum', 'clade', 'superkingdom'))
    stop("Taxon is not among accepted values ('species', 'genus', 'family', 'order', 'class', 
                                     'phylum', 'clade', 'superkingdom')")
  
  # Pre-processing
  barcodeslist <- paste0("barcode", barcodeslist)
  path <- paste0(path, "/")
  n <- length(barcodeslist)
  
  # Read collapsed abundance tables in a list
  all_barcodes <- list()
  for (i in 1:n)
    all_barcodes[[i]] <- read.table(header = T, sep = "\t",
                                    file = paste0(path, "/", taxon, "/", barcodeslist[i], pt_suff),
                                    stringsAsFactors=FALSE)
  
  # Create table of abundances
  taxa_vector <- c()
  for (i in 1:n) 
    taxa_vector <- unique(c(taxa_vector, all_barcodes[[i]][, taxon]))
  # Note: the column name is inside the variable "taxon". $variable does not
  # work as intended so the column is accessed by [, variable]
  
  taxa_abundance <- data.frame("taxa"=taxa_vector)
  colnames(taxa_abundance) <- taxon
  
  for (i in 1:n) taxa_abundance <-
    merge(taxa_abundance, all_barcodes[[i]][,c(taxon, "abundance")],
          by=taxon, all.x = TRUE)
  colnames(taxa_abundance) <- c(taxon, barcodeslist)
  
  taxa_abundance <- subset(taxa_abundance, !is.na(taxa_abundance[,taxon]))
  taxa_abundance <- subset(taxa_abundance, taxa_abundance[,taxon] != "")
  
  # Substitute NA abundances with 0
  taxa_abundance[is.na(taxa_abundance)] <- 0
  
  # Transpose table
  taxa_abundance_t <- cbind(Sample = colnames(taxa_abundance)[-1], as.data.frame(t(taxa_abundance[, -1])))
  colnames(taxa_abundance_t) <- c("Sample", as.character(taxa_abundance[,1]))
  
  return(list(taxa_abundance, taxa_abundance_t))
}
