class ClassificationConfig:
    def __init__(self):
        self.methods_on = {'rf': 'raw',
                           'lgb': 'raw',
                           'xgb': 'raw',
                           'var': 'raw',
                           'cv2': 'raw',
                           'nsc': 'raw',
                           'seurat': 'raw',
                           'deviance': 'raw',
                           'm3drop': 'raw',
                           'cellassign': 'raw'}


class ClusteringConfig:
    def __init__(self):
        self.methods_on = {
            'var': 'raw',
            'cv2': 'raw',
            'seurat': 'raw',
            'deviance': 'raw',
            'm3drop': 'raw',
        }
