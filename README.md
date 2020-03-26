# FastCAR

FastCAR is an R package to remove ambient RNA from cells in droplet based single cell RNA sequencing data.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

FastCAR can be install from git with the following command.

```
devtools::install_git("https://git.web.rug.nl/P278949/FastCAR")
```

Running FastCAR is quite simple.
First load the library and dependencies.

```
library(Matrix)
library(Seurat)
library(qlcMatrix)
library(FastCAR)

```



```

cellExpressionFolder  = c("Cellranger_output/sample1/filtered_feature_bc_matrix/")
fullMatrixFolder      = c("Cellranger_output/sample1/raw_feature_bc_matrix/")

```



```
# This folder will contain the corrected cell matrix
correctedMatrixFolder = c("Cellranger_output/sample1/corrected_feature_bc_matrix")


cellMatrix     = read.cell.matrix(cellExpressionFolder)
fullMatrix     = read.full.matrix(fullMatrixFolder)

```

The following functions give an idea of the effect that different settings have on the ambient RNA profile

```
ambProfile = describe.ambient.RNA.sequence(fullCellMatrix = fullMatrix, start = 10, stop = 500, by = 10, contaminationChanceCutoff = 0.05)
plot.ambient.profile(ambProfile)

``` 
![picture](https://git.web.rug.nl/P278949/FastCAR/src/branch/master/Images/Example_profile.png)



Set 

```

emptyDropletCutoff        = 100 
contaminationChanceCutoff = 0.05

```

```
ambientProfile = determine.background.to.remove(fullMatrix, cellMatrix, emptyDropletCutoff, contaminationChanceCutoff)
cellMatrix     = remove.background(cellMatrix, ambientProfile)


```

Finally write the corrected cell/gene matrix to a file, this matrix can be used in Seurat the same way as any other cell/gene matrix.

```

write.corrected.matrix(cellMatrix, correctedMatrixFolder, ambientProfile)

```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests


## Authors

* **Marijn Berg** - *Initial work* 

## License

This project is licensed under the GPL-3 License - see the [LICENSE.md](LICENSE.md) file for details

