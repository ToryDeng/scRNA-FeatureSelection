from experiments.run_experiments import run_cell_clustering


run_cell_clustering(
    fs_methods=['pagest1w', 'pagest2w'],
    use_saved_genes=False
)
# 'cv2', 'deviance', 'seurat_v3', 'm3drop', 'feast', 'scran', 'triku', 'pagest1w', 'pagest2w'
