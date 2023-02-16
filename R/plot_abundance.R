# Libraries
require(ggplot2)
require(vegan)
require(paletteer)
require(dplyr)
#################################################################################
# Function to create stacked barplots for taxa relative abundances
plot_abundance <- function(abundance_table, metadata, 
                           taxon=c('species', 'genus', 'family', 'order', 'class', 
                                   'phylum', 'clade', 'superkingdom'),
                           variable = c("Sample", "Name", "Type", "Sampling_date",
                                        "Batch", "Storage_time", "Pore_size")){
  
  # Check inputs
  # Check abundance table format
  if(!is.data.frame(abundance_table)) stop("abundance_table not a data.frame")
  # Check if taxon is a character string
  if (!is.character(taxon)) stop("taxon not a character string")
  if(!taxon %in% c('species', 'genus', 'family', 'order', 'class', 
                   'phylum', 'clade', 'superkingdom'))
    stop("taxon not among accepted values ('species', 'genus', 'family', 'order', 'class', 
                                                    'phylum', 'clade', 'superkingdom')")
  # Metadata must be in the correct format
  if(!is.data.frame(metadata)) stop("metadata is not a data.frame")
  if(!all(colnames(metadata) == c( "Sample", "Name", "Type", "Sampling_date",
                                   "Batch", "Storage_time", "Pore_size")))
    stop("Incorrect colnames in metadata")
  # Check if the variable is among the ones in the metadata
  if(!variable %in% colnames(metadata)) stop("variable not in metadata")
  if(!is.character(variable)) stop("variable is not a character string")
  
  
  # Transform the transposed table of abundance into the long format
  taxa_abundance_long <- reshape2::melt(abundance_table, id.vars = "Sample",
                                        variable.name = taxon)
  # Add metadata to the long table
  taxa_abundance_long <- merge(metadata, taxa_abundance_long, by.y = "Sample")
  
  # Remove factors
  for (i in 1:length(colnames(taxa_abundance_long))) {
    if(is.factor(taxa_abundance_long[,i])) taxa_abundance_long[,i] <- as.character(taxa_abundance_long[,i])
  }
  
  # Filter value != 0
  taxa_abundance_long <- subset(taxa_abundance_long, taxa_abundance_long$value != 0)
  
  # Plot
  pl <- ggplot(taxa_abundance_long, aes(x = Sample, y = value, fill = taxa_abundance_long[,taxon])) + 
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_text(angle = 90),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.text = element_text(size=10),
          legend.key.size = unit(0.2, "cm"))+
    facet_wrap(~taxa_abundance_long[, variable], scales="free_x", nrow=1)+
    labs(x="Sample",y="Relative Abundance (%)", fill = taxon) +
    scale_fill_manual(values = paletteer_c("grDevices::Spectral",
                                           n = length(unique(taxa_abundance_long[,taxon]))))
  
  return(pl)
}