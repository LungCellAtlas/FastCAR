###############################################################################
# FastCAR package
# Marijn Berg
# m.berg@umcg.nl
###############################################################################
# FastCAR removes ambient RNA from 10X data by looking at the profile of the
# empty droplets and removing per gene the highest value found in the ambient
# RNA if it's found to contaminate more than a certain threshold of droplets
###############################################################################
# This script contains functions to profile the Ambient RNA to suggest
# good settings to use to find genes to correct for

###############################################################################
library(Matrix)
library(Seurat)
library(qlcMatrix)
###############################################################################

# describe the number of genes identified in the background
# and the number of genes failing the contaminiation chance threshold
#
describe.ambient.RNA.sequence = function(fullCellMatrix, start, stop, by, contaminationChanceCutoff){
  genesInBackground  = vector(mode = "numeric", length = length(seq(start, stop, by)))
  genesContaminating = vector(mode = "numeric", length = length(seq(start, stop, by)))
  nEmptyDroplets     = vector(mode = "numeric", length = length(seq(start, stop, by)))

  ambientDescriptions = data.frame(nEmptyDroplets, genesInBackground, genesContaminating)
  rownames(ambientDescriptions) = seq(start, stop, by)
  for(emptyCutoff in seq(start, stop, by)){
    nEmpty = table((Matrix::colSums(fullCellMatrix) < emptyCutoff) &(Matrix::colSums(fullCellMatrix) > 0))[2]

    occurences = rowSums(fullCellMatrix[,Matrix::colSums(fullCellMatrix) < emptyCutoff] !=0)

    #probability of a background read of a gene ending up in a cell
    probabiltyCellContaminationPerGene = occurences / nEmpty
    nFailingThreshold = sum(probabiltyCellContaminationPerGene > contaminationChanceCutoff)

    nGenes = sum(occurences != 0)
    ambientDescriptions[as.character(emptyCutoff), c(1,2,3)] = c(nEmpty ,nGenes, nFailingThreshold)
  }
  return(ambientDescriptions)
}



plot.ambient.profile = function(ambientProfile){
  par(mfrow = c(3,1))
  plot(as.numeric(rownames(ambientProfile)), ambientProfile[,1],
       main = "Total number of empty droplets at cutoffs",
       xlab = "empty droplet UMI cutoff",
       ylab = "Number of empty droplets")

  plot(as.numeric(rownames(ambientProfile)), ambientProfile[,2],
       main = "Number of genes in ambient RNA",
       xlab = "empty droplet UMI cutoff",
       ylab = "Genes in empty droplets")

  plot(as.numeric(rownames(ambientProfile)), ambientProfile[,3],
       main = "number of genes to correct",
       xlab = "empty droplet UMI cutoff",
       ylab = "Genes identified as contamination")
}



