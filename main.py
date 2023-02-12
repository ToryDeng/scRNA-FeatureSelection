from experiments.run_experiments import run_cell_clustering, run_cell_classification, run_computation_time


fs_methods = ['giniclust3', 'deviance', 'seurat_v3', 'm3drop', 'feast', 'scran', 'triku', 'var', 'sc3', 'scmap']


run_cell_clustering(
    fs_methods=fs_methods,
    use_saved_genes=True,
    all_genes=False
)

# run_computation_time(fs_methods=fs_methods)

# 'giniclust3', 'deviance', 'seurat_v3', 'm3drop', 'feast', 'scran', 'triku', 'var', 'sc3', 'cv', 'scmap',
# 'pagest1w', 'pagest2w'

