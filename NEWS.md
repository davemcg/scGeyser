# scGeyser 0.03

Fix bug in dot plot size. Print mgif YAML document on loading screen. 
Fix display bug where the headers for the loading screen were getting covered by the data tables.
Now takes mgif var field "active_gene_column" to show display name for genes (to help avoid situations where the var_name is "ENSG[digit]"). 

# scGeyser 0.02

Add GUI interface for "preloaded" mgif files (`scGeyser(mgif_dir = '/folder/where/mgif/files/are'))

# scGeyser 0.01

"scGeyser", a visualization tool for SingleCellExperiment 
(SCE). Also supports Seurat, hdf5-backed SCE, and the "mgif" format. 
