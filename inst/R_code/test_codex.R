library(Giotto)

# 1. set working directory
results_folder = 'C:/Document/Serieux/Travail/Data_analysis_and_papers/spatial_data/'

# 2. set giotto python path
# set python path to your preferred python version path
# set python path to NULL if you want to automatically install (only the 1st time) and use the giotto miniconda environment
python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}



# download data to working directory
# use method = 'wget' if wget is available. This should be much faster.
# if you run into authentication issues with wget, then add " extra = '--no-check-certificate' "
getSpatialDataset(dataset = 'codex_spleen', directory = results_folder)



expr_path = paste0(results_folder, "codex_BALBc_3_expression.txt.gz")
loc_path = paste0(results_folder, "codex_BALBc_3_coord.txt")
meta_path = paste0(results_folder, "codex_BALBc_3_annotation.txt")


instrs = createGiottoInstructions(show_plot = FALSE,
                                  save_plot = TRUE,
                                  save_dir = results_folder,
                                  python_path = python_path)


# expression info
codex_expression = readExprMatrix(expr_path, transpose = F)
# cell coordinate info
codex_locations = data.table::fread(loc_path)
# metadata
codex_metadata = data.table::fread(meta_path)



## stitch x.y tile coordinates to global coordinates
xtilespan = 1344;
ytilespan = 1008;
# TODO: expand the documentation and input format of stitchTileCoordinates. Probably not enough information for new users.
stitch_file = stitchTileCoordinates(location_file = codex_metadata, Xtilespan = xtilespan, Ytilespan = ytilespan);
codex_locations = stitch_file[,.(Xcoord, Ycoord)]

# create Giotto object
codex_test <- createGiottoObject(raw_exprs = codex_expression,
                                 spatial_locs = codex_locations,
                                 instructions = instrs,
                                 cell_metadata = codex_metadata)

# subset Giotto object
cell_meta = pDataDT(codex_test)
cell_IDs_to_keep = cell_meta[Imaging_phenotype_cell_type != "dirt" & Imaging_phenotype_cell_type != "noid" & Imaging_phenotype_cell_type != "capsule",]$cell_ID
codex_test = subsetGiotto(codex_test, cell_ids = cell_IDs_to_keep)

## filter
codex_test <- filterGiotto(gobject = codex_test,
                           expression_threshold = 1,
                           gene_det_in_min_cells = 10,
                           min_det_genes_per_cell = 2,
                           expression_values = c('raw'),
                           verbose = T)

codex_test <- normalizeGiotto(gobject = codex_test, scalefactor = 6000, verbose = T,
                              log_norm = FALSE,library_size_norm = FALSE,
                              scale_genes = FALSE, scale_cells = TRUE)

## add gene & cell statistics
codex_test <- addStatistics(gobject = codex_test,expression_values = "normalized")

## adjust expression matrix for technical or known variables
codex_test <- adjustGiottoMatrix(gobject = codex_test,
                                 expression_values = c('normalized'),
                                 batch_columns = NULL,
                                 covariate_columns = NULL,
                                 return_gobject = TRUE,
                                 update_slot = c('custom'))

## visualize
spatPlot(gobject = codex_test,point_size = 0.1,
         coord_fix_ratio = NULL,point_shape = 'no_border',
         save_param = list(save_name = '2_a_spatPlot'))



spatPlot(gobject = codex_test, point_size = 0.2,
         coord_fix_ratio = 1, cell_color = 'sample_Xtile_Ytile',
         legend_symbol_size = 3,legend_text = 5,
         save_param = list(save_name = '2_b_spatPlot'))
