###############################################################################
# FastCAR package
# Marijn Berg
# m.berg@umcg.nl
###############################################################################
# FastCAR removes ambient RNA from 10X data by looking at the profile of the
# empty droplets and removing per gene the highest value found in the ambient
# RNA if it's found to contaminate more than a certain threshold of droplets
###############################################################################

# This script runs a simple correction with standard settings
library(Matrix)
library(Seurat)
library(qlcMatrix)
source("R/FastCAR_Base.R")

emptyDropletCutoff        = 100  # Droplets with fewer than this number of UMIs are considered empty
contaminationChanceCutoff = 0.05 # This is the maximum fraction of droplets containing the contamination

cellExpressionFolder  = c("/home/marijn/Analysis_Drive/R_Projects/scRNA_Background_correction/Cellranger_output/4951STDY7487591_ARMS054/filtered_feature_bc_matrix/")
fullMatrixFolder      = c("/home/marijn/Analysis_Drive/R_Projects/scRNA_Background_correction/Cellranger_output/4951STDY7487591_ARMS054/raw_feature_bc_matrix/")

# This folder will contain the corrected cell matrix
correctedMatrixFolder = c("/home/marijn/Analysis_Drive/R_Projects/scRNA_Background_correction/Cellranger_output/4951STDY7487591_ARMS054/corrected_feature_bc_matrix")

# these are literally wrappers for the Seurat functions but it's good practice in case the function changes later
cellMatrix     = read.cell.matrix(cellExpressionFolder)
fullMatrix     = read.full.matrix(fullMatrixFolder)

###############################################################################
# This is an optional function that will show the effects of different cutoff values
# start, stop, and by all refer to number of UMI, it determines the background ate steps of by, from start to stop
ambProfile = describe.ambient.RNA.sequence(fullCellMatrix = fullMatrix, start = 10, stop = 500, by = 10, contaminationChanceCutoff = 0.05)
plot.ambient.profile(ambProfile)
###############################################################################

ambientProfile = determine.background.to.remove(fullMatrix, cellMatrix, emptyDropletCutoff, contaminationChanceCutoff)

cellMatrix     = remove.background(cellMatrix, ambientProfile)

write.corrected.matrix(cellMatrix, correctedMatrixFolder, ambientProfile)
