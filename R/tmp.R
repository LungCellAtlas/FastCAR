cutoffEmpty = 100   # set this to NA to determine dynamicallly

###############################################################################

# remove this or things will keep getting added to it
rm(biopsySeurat, backGroundOverviews)

for(folder in dataFolders[1]){
  gc()
  donor  = str_sub(folder, -7,-1)
  sample = str_sub(folder, 19)
  sample = str_sub(sample, 1, -9)

  allExpression  = Read10X(paste(folder, "/raw_feature_bc_matrix", sep = ""))
  colnames(allExpression) = paste(datasetInfo[sample, "prefix"], colnames(allExpression), sep = "")
  # change barcode ids to match r object

  if(is.na(cutoffEmpty)){
    ##########################################
    # determine the number of genes in the empty droplets at various cutoffs

    foundNumberOfGenes = vector(mode = "numeric", length = length(seq(kneePointStart, kneePointStop, kneePointSteps)))
    names(foundNumberOfGenes) = as.character(seq(kneePointStart, kneePointStop, kneePointSteps))

    for(cutoffValue in seq(kneePointStart, kneePointStop, kneePointSteps)){
      nExpressedCells = Matrix::rowSums(allExpression[, Matrix::colSums(allExpression) < cutoffValue])
      nExpressedCells[nExpressedCells > 0] = 1
      foundNumberOfGenes[as.character(cutoffValue)]= sum(nExpressedCells)
    }
    write.table(foundNumberOfGenes, file = paste("Foundgenes_cutoff_", donor, ".csv"))
    ##########################################
    # Find the point where increase in the number of genes decreases the most

    cutoffEmpty = which(diff(foundNumberOfGenes) == max(diff(foundNumberOfGenes)[find_peaks(diff(foundNumberOfGenes))]))

    # cutoffEmpty = as.numeric(names(which(diff(foundNumberOfGenes) == min(diff(foundNumberOfGenes)))))[1] # have to add this at the end, multiple changes can be the same value

    ##########################################
    png(paste("Empty_droplets_content_", donor, ".png", sep = ""), height = 400, width = 800)
    plot(as.numeric(names(foundNumberOfGenes)), foundNumberOfGenes,   main = paste(donor, " unique genes per cutoff"), xlab = "cutoff <")
    segments(cutoffEmpty, 0 ,cutoffEmpty, 25000)
    dev.off()
  }

  # read this in as a full matrix, haven't figgered out how to do this on a sparse matrix
  cellExpression = allExpression[,cellIDs[cellIDs %in% colnames(allExpression)]]

  if(!exists("biopsySeurat")){
    biopsySeurat <<- CreateSeuratObject(cellExpression, project = "Bronchial_Biopsy")
    biopsySeurat[["orig.ident"]] = donor
    biopsySeurat[["percent.mt"]] <- PercentageFeatureSet(biopsySeurat, pattern = "^MT-")
    # biopsySeurat <- subset(biopsySeurat, subset = percent.mt < quantile(biopsySeurat@meta.data$percent.mt, c(.9)) )
  }else{
    tmpbiopsySeurat = CreateSeuratObject(cellExpression, project = "Bronchial_Biopsy")
    tmpbiopsySeurat[["orig.ident"]] = donor

    tmpbiopsySeurat[["percent.mt"]] <- PercentageFeatureSet(tmpbiopsySeurat, pattern = "^MT-")
    # tmpbiopsySeurat <- subset(tmpbiopsySeurat, subset = percent.mt < quantile(tmpbiopsySeurat@meta.data$percent.mt, c(.9), na.rm = TRUE))

    biopsySeurat = merge(biopsySeurat, tmpbiopsySeurat)
  }


  backgroundTotal  = Matrix::rowSums(allExpression[,Matrix::colSums(allExpression) < cutoffEmpty])
  backGroundMax   = as.vector(rowMax(allExpression[names(backgroundTotal),Matrix::colSums(allExpression) < cutoffEmpty]))

  nCell = ncol(cellExpression)
  # droplets that are empty but not unused barcodes, unused barcodes have zero reads assigned to them.
  nEmpty = table((Matrix::colSums(allExpression) < cutoffEmpty) &(Matrix::colSums(allExpression) > 0))[2]
  occurences = rowSums(allExpression[,Matrix::colSums(allExpression) < cutoffEmpty] !=0)

  #probability of a background read of a gene ending up in a cell
  pCell = occurences / nEmpty

  write.csv(pCell, file = paste("pCell_", donor, ".csv", sep = ""))
  write.csv(backGroundMax, file = paste("bgMax_", donor, ".csv", sep = ""))



  # set the background value for those genes that are unlikely to be a significant contaminant to 0 (zero)
  backGroundMax[pCell < acceptableFractionContaminated] = 0

  # remove the background
  cellExpression = apply(cellExpression , 2, '-', backGroundMax)
  cellExpression[cellExpression < 0] = 0
  # cellExpression = as.sparse(cellExpression)

  if(!exists("biopsySeuratAC")){
    biopsySeuratAC <<- CreateSeuratObject(cellExpression, project = "AC_Bronchial_Biopsy")
    biopsySeuratAC[["orig.ident"]] = donor
    biopsySeuratAC[["percent.mt"]] <- PercentageFeatureSet(biopsySeuratAC, pattern = "^MT-")
    # biopsySeuratAC <- subset(biopsySeuratAC, subset = percent.mt < quantile(biopsySeuratAC@meta.data$percent.mt, c(.9)) )
  }else{
    tmpbiopsySeuratAC = CreateSeuratObject(cellExpression, project = "AC_Bronchial_Biopsy")
    tmpbiopsySeuratAC[["orig.ident"]] = donor

    tmpbiopsySeuratAC[["percent.mt"]] <- PercentageFeatureSet(tmpbiopsySeuratAC, pattern = "^MT-")
    # tmpbiopsySeuratAC <- subset(tmpbiopsySeuratAC, subset = percent.mt < quantile(tmpbiopsySeuratAC@meta.data$percent.mt, c(.9), na.rm = TRUE))

    biopsySeuratAC = merge(biopsySeuratAC, tmpbiopsySeuratAC)
  }
}
