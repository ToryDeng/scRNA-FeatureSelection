from experiments.run_experiments import run_cell_clustering


run_cell_clustering(
    fs_methods=['sct', 'giniclust3'],
    use_saved_genes=True,
    all_genes=True
)
# 'cv2', 'deviance', 'seurat_v3', 'm3drop', 'feast', 'scran', 'triku', 'sct', 'giniclust3', 'pagest1w', 'pagest2w'
# TODO: run SCT and giniclust3


