#http://spatial.libd.org/spatialLIBD/
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
  #  install.packages("BiocManager")
  #}
#BiocManager::install("LieberInstitute/spatialLIBD")
## Load the package
library("spatialLIBD")

## Download the spot-level data
spe <- fetch_data(type = "spe")

## This is a SpatialExperiment object
spe
#> class: SpatialExperiment
#> dim: 33538 47681
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#>   ENSG00000268674
#> rowData names(9): source type ... gene_search is_top_hvg
#> colnames(47681): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
#>   TTGTTTCCATACAACT-1 TTGTTTGTGTAAATTC-1
#> colData names(66): sample_id Cluster ... spatialLIBD ManualAnnotation
#> reducedDimNames(6): PCA TSNE_perplexity50 ... TSNE_perplexity80
#>   UMAP_neighbors15
#> mainExpName: NULL
#> altExpNames(0):
#> spatialData names(3) : in_tissue array_row array_col
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor

## Note the memory size
lobstr::obj_size(spe)
#> 2.04 GB

## Remake the logo image with histology information
vis_clus(
  spe = spe,
  clustervar = "spatialLIBD",
  sampleid = "151673",
  colors = libd_layer_colors,
  ... = " DLPFC Human Brain Layers\nMade with research.libd.org/spatialLIBD/"
)
