import scanpy as sc
import os
from differential_expression_analysis import Differential_Expression_Analysis as dea

class Cell_Selection:

    def __init__(self, adata, tsv, output_path, full_sample_name):
        self.adata = adata
        self.tsv = tsv
        self.output_path = output_path
        self.full_sample_name = full_sample_name


    def splittsv(self):
        pass

    # This method might be removed, but I think it is good to check for the 3 important markers?
    # Q: ask maurits 3 important markers or all of them?
    # If one or even all of them are missing --> remove the cluster!
    # end with boolean value? so like true or false --> if true, then continue
    # if false, remove cluster end loop.
    def markers_present(self):
        neurons = ['MAP2', 'DCX', 'NEUROG2']
        astrocytes = ['VIM', 'S100B', 'SOX9']
        neuronscheck = dea.checkpoint(self.adata, neurons)
        astrocytescheck = dea.checkpoint(self.adata, astrocytes)
        if len(neurons) == len(neuronscheck):
            check1 = True
        if len(astrocytes) == (astrocytescheck):
            check2 = True
            
        if check1 and check2:
            return True
        else:
            return False


    # If markers present is true --> check for positive log fold
    # if all of them are pos --> also true return
    def pos_log_fold(self):
        return True


    # check if both methods above are true
    # removes clusters
    # provides a list with clusters that are good to go
    # return list?
    def check_markers(self):
        pass


    # do the filtering
    def filter_clusters(self, list):
        pass


    # Q: should you do the integration on pca h5ad file again but now with clustered results?
    # ask maurits
    def perform_analysis(self, list):
        pass
# use final anndata after dea instead of pca only!
# redo pca after integration
# 20 pc in knn plot
# filteren op blA blC matrix raw counts (cell selectie)
# na cell selectie weer pca knn leiden en plots
#  nog een umap 
# wacht met umap totdat annotatie gedaan 

    