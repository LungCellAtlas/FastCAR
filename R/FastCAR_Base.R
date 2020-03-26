###############################################################################
# FastCAR package
# Marijn Berg
# m.berg@umcg.nl
###############################################################################
# FastCAR removes ambient RNA from 10X data by looking at the profile of the
# empty droplets and removing per gene the highest value found in the ambient
# RNA if it's found to contaminate more than a certain threshold of droplets
###############################################################################
# This file contains base functions
###############################################################################
library(Matrix)
library(Seurat)
library(qlcMatrix)
###############################################################################

# had to make this function to efficiently modify sparse matrices on a per gene basis
# A dgCMatrix object has the following elements that matter for this function
# i: a sequence of the row locations of each non-zero element
# x: the non-zero values in the matrix
# p: the where in 'i' and 'x' a new column starts
# Dimnames: The names of the rows and columns
remove.background = function(geneCellMatrix, ambientRNAprofile){
  for(gene in names(ambientProfile[ambientProfile > 0])){
    # Determine the locations where the gene is not zero, therefore referenced in i
    iLocs = which(geneCellMatrix@Dimnames[[1]] == gene)
    # Determine the location of the actual values,
    xLocs = which(geneCellMatrix@i == iLocs-1) # -1 because of 0 and 1 based counting systems
    # Remove the contaminating RNA
    geneCellMatrix@x[xLocs] = geneCellMatrix@x[xLocs] - ambientProfile[gene]
  }
  # correct for instances where the expression was corrected to below zero
  geneCellMatrix@x[geneCellMatrix@x < 0] = 0
  # remove the zeroes and return the corrected matrix
  return(drop0(geneCellMatrix, is.Csparse = TRUE))
}
##############

determine.background.to.remove = function(fullCellMatrix, cellMatrix, emptyDropletCutoff, contaminationChanceCutoff){

  # determines the highest expression value found for every gene in the droplets that we're sure don't contain cells
  backGroundMax   = as.vector(rowMax(fullCellMatrix[,Matrix::colSums(fullCellMatrix) < emptyDropletCutoff]))
  names(backGroundMax) = rownames(fullCellMatrix)
  nCell = ncol(cellMatrix)

  # droplets that are empty but not unused barcodes, unused barcodes have zero reads assigned to them.
  nEmpty = table((Matrix::colSums(fullCellMatrix) < emptyDropletCutoff) &(Matrix::colSums(fullCellMatrix) > 0))[2]
  # rowSum on a logical statement returns the number of TRUE occurences
  occurences = rowSums(fullCellMatrix[,Matrix::colSums(fullCellMatrix) < emptyDropletCutoff] !=0)

  #probability of a background read of a gene ending up in a cell
  probabiltyCellContaminationPerGene = occurences / nEmpty

  # if the probablity of a gene contaminating a cell is too low we set the value to zero so it doesn't get corrected
  backGroundMax[probabiltyCellContaminationPerGene < contaminationChanceCutoff] = 0
  return(backGroundMax)
}

##############

read.cell.matrix = function(cellFolderLocation){
  cellMatrix = Read10X(cellFolderLocation)
  return(cellMatrix)
}

##############

read.full.matrix = function(fullFolderLocation){
  fullMatrix = Read10X(fullFolderLocation)
  return(fullMatrix)
}

##############

write.corrected.matrix = function(correctedMatrix, folderLocation, correctedSignal){
  dir.create(folderLocation)

  writeMM(obj = correctedMatrix, file=paste(folderLocation, "/matrix.mtx", sep = ""))

  # save genes and cells names
  write(x = rownames(correctedMatrix), file = paste(folderLocation, "/genes.tsv", sep = ""))
  write(x = colnames(correctedMatrix), file = paste(folderLocation, "/barcodes.tsv", sep = ""))

  correctedSignal = correctedSignal[correctedSignal > 0]
  write.table(list(correctedSignal), file = paste(folderLocation, "/genesCorrectedFor.csv", sep = ""), row.names = TRUE, col.names = FALSE)

}

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












