from utils.evaluation import evaluate_feature_selection_methods

# evaluate_feature_selection_methods(measurements=['population_demixing', 'marker_discovery',
#                                                  'classification', 'clustering',
#                                                  'batch_correction', 'computation_time'],
#                                    methods=['var', 'seurat', 'deviance', 'm3drop', 'scmap',
#                                             'rf', 'lgb', 'xgb', 'nsc', 'fisher_score', 'scGeneFit',
#                                             'cellranger', 'feast', 'rf+fisher_score', 'seurat+deviance'])
evaluate_feature_selection_methods(measurements=['classification'],
                                   methods=['var', 'cv2', 'seurat', 'deviance', 'm3drop', 'scmap', 'rf', 'scGeneFit',
                                            'lgb', 'xgb', 'nsc', 'fisher_score', 'cellranger', 'feast',
                                            'rf+fisher_score', 'seurat+deviance']
                                   )
