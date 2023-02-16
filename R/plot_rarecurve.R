# Libraries
require(ggplot2)
require(vegan)
require(paletteer)
require(dplyr)

#################################################################################
# Function to plot rarefaction curves
plot_rarecurve <- function(reads_table, controls_barcodes){
  
  # Check if reads_table is a data.frame
  
  # Check if controls_barcodes is a vector of characters
  
  
  # Round counts to integers
  #reads_table <- round(reads_table, digits = 0)
  # Does not work for a dataframe
  
  nums <- vapply(reads_table, is.numeric, FUN.VALUE = logical(1))
  reads_table[,nums] <- round(reads_table[,nums], digits = 0)
  
  # Exclude controls
  controls <- paste0("barcode", controls_barcodes)
  reads_table <- subset(reads_table, select = !(names(reads_table) %in% controls))
  
  # Transpose
  reads_table_t <- cbind(Sample = colnames(reads_table)[-1], as.data.frame(t(reads_table[, -1])))
  colnames(reads_table_t) <- c("Sample", as.character(reads_table[,1]))
  
  #count the number of species
  S <- specnumber(reads_table_t[,-1])
  
  # Find the minimum number of reads:
  # sum all reads per row (thus reads per sample)
  # then find the minimum sum (reads number of the sample with less reads)
  raremax<-min(rowSums(reads_table_t[,-1]))
  #Rarefaction of the samples to the minimum number of reads in the experiment
  Srare <- rarefy(reads_table_t[,-1], raremax)
  
  #plot with vegan
  col <- paletteer_c("viridis::viridis", n = length(reads_table_t$Sample))
  lty <- c("solid", "dashed", "dotdash")
  lwd <- c(1, 2)
  pars <- expand.grid(col = col, lty = lty, lwd = lwd, 
                      stringsAsFactors = FALSE)
  head(pars)
  
  out <- with(pars[1:length(reads_table_t$Sample), ],
              rarecurve(reads_table_t[,-1], step = 20, sample = raremax, col = col,
                        lty = lty, label = TRUE))
  
  return(out)
  
}