library(preprocessCore)

args <- commandArgs(trailingOnly = TRUE)

data = read.table("K:/2024_NC/Calculate_volume/KR_MCFS/NHEK/3.TAD_volume_sort", sep = "\t")

d = as.matrix(as.numeric(data[,'TCI']))
normalize.quantiles.use.target(d, target = rnorm(50000), copy=FALSE, subset=NULL)
d <- as.character(d)
data[,'TCI_Qnor'] = d

data <- as.vector(data)

write.table( data, "K:/2024_NC/Calculate_volume/KR_MCFS/NHEK/8.TAD_volume_sort_normalized",quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)



