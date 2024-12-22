library(pcalg)

dataset <- read.csv(file='data.csv', header=FALSE, sep=",");
estDAG <- lingam(dataset, verbose = FALSE)
write.csv(as.matrix(estDAG$Bpruned),row.names = FALSE, file = 'result.csv');
