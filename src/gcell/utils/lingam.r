# Template R script for lingam
library(pcalg)

dataset <- read.csv(file='{FOLDER}{FILE}', header=FALSE, sep=",");
estDAG <- lingam(dataset, verbose = {VERBOSE})
write.csv(as.matrix(estDAG$Bpruned),row.names = FALSE, file = '{FOLDER}{OUTPUT}');
