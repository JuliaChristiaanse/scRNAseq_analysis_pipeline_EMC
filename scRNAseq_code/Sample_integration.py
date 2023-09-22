import scanpy as sc
import os
import pandas as pd
import harmonypy as hm
import anndata
import numpy as np
from matplotlib import pyplot as plt
from Differential_expression_analysis import Differential_Expression_Analysis as dea

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=120, dpi_save=150, facecolor='white', color_map='Blues')

class Sample_integration:
    def __init__(self, co_sample, mono_sample, full_name, output_dir, sample_output, markerpath):
        self.co_sample = co_sample
        self.mono_sample = mono_sample
        self.output_dir = output_dir
        self.full_name = full_name
        self.concat_anndata = self.concatinate_adatas()
        self.sample_output = sample_output
        self.markerpath = markerpath
        self.run()


    # Directories are created
    def makedirs(self):
        os.makedirs(os.path.join(self.sample_output))
        os.chdir(os.path.join(self.sample_output))
        os.makedirs('DEA')
        os.makedirs('UMAPS')
        os.makedirs('AnnData_storage')


    # The mono and co samples are concatinated in one anndata object and returned
    def concatinate_adatas(self):
        self.co_adata = sc.read_h5ad(os.path.join(self.output_dir, self.co_sample,'AnnData_storage', self.co_sample+'.h5ad'))
        self.mono_adata = sc.read_h5ad(os.path.join(self.output_dir, self.mono_sample,'AnnData_storage', self.mono_sample+'.h5ad'))
        self.co_adata.obs['batch'] = self.co_sample
        self.mono_adata.obs['batch'] = self.mono_sample
        return anndata.concat([self.co_adata, self.mono_adata], index_unique='-', axis=0, join='inner', 
                              merge=None, uns_merge=None, label=None, keys=[self.co_sample, self.mono_sample],
                                fill_value=None, pairwise=None)
    

    # harmony is ran according to the vignettes & github issues found
    # plots are created
    # TO DO: Create plots for solo samples in for loop and add to axis (try in Notebook)
    def run_harmony(self, adataobj, umaptitle):
        sc.pp.pca(adataobj, n_comps=50)
        data_mat = adataobj.obsm['X_pca']
        meta_data = adataobj.obs
        vars_use = ['batch']
        title = self.co_sample+' '+self.mono_sample+' '+umaptitle
        ho = hm.run_harmony(data_mat, meta_data, vars_use,
                             theta=None, lamb=None, sigma=0.1,
                             nclust=None, tau=0, block_size=0.5,
                             max_iter_harmony=10, max_iter_kmeans=20,
                             epsilon_cluster=1e-5, epsilon_harmony=1e-4,
                             plot_convergence=True, verbose=True, reference_values=None,
                             cluster_prior=None, random_state=0)
        adjusted_pcs = pd.DataFrame(ho.Z_corr).T
        adataobj.obsm['X_pca']=adjusted_pcs.values
        sc.pp.neighbors(adataobj, n_pcs=20)
        sc.tl.leiden(adataobj, resolution=0.5)
        sc.tl.umap(adataobj)
        fig, axs = plt.subplots(2, 2, figsize=(10,8),constrained_layout=True)
        # change to lower case
        cocult = adataobj[adataobj.obs['batch'] == self.co_sample]
        sample = adataobj[adataobj.obs['batch'] == self.mono_sample]
        sc.pl.umap(adataobj, color="batch", title=title, ax=axs[0,0], show=False)
        sc.pl.umap(adataobj, color="leiden", title="Leidenalg UMAP", ax=axs[0,1], show=False)
        sc.pl.umap(sample, color="leiden", title=f"{self.mono_sample} sample only", ax=axs[1,0], show=False)
        sc.pl.umap(cocult, color="leiden", title=f"{self.co_sample} sample only", ax=axs[1,1], show=False)
        plt.savefig(os.path.join(self.sample_output, 'UMAPS', title+'.png'))


    # This function performs the cell selection
    def cell_selection(self, adata_DE, adataobj, genes, dp):
        dotdf = dp.dot_size_df
        pos = dotdf[(dotdf[genes[0]] > .2) & (dotdf[genes[1]] > .2) & (dotdf[genes[2]] > .2)]
        good_clusters = list(pos.index)
        self.adata_subset = adataobj[adataobj.obs['leiden'].isin(good_clusters)]
        self.run_harmony(self.adata_subset, "UMAPS after cell selection")
        self.path = os.path.join(self.sample_output, 'AnnData_storage', f'{self.full_name}_subset_anndata.h5ad')
        self.adata_subset.write(self.path)


    # This function calls all other functions & runs them
    def run(self):
        self.makedirs()
        self.run_harmony(self.concat_anndata, "UMAPS before cell selection")
        self.adata_DE = self.concat_anndata.raw.to_adata()
        deado = dea(self.adata_DE, self.sample_output, self.full_name, self.markerpath)
        deado.perform_dea()
        deado.basic_dea_plots()
        if self.mono_sample == 'BL_N':
            genes = ['DCX', 'MAP2', 'RBFOX3']
        else:
            genes = ['VIM', 'S100B', 'SOX9']
        # dp = sc.pl.rank_genes_groups_dotplot(deado.adata, deado.cluster, color_map='Blues',
        #                                       var_names=genes, return_fig=True, show=False)
        dp = deado.rank_genes_dotplot_overview(genes)
        self.cell_selection(self.adata_DE, self.concat_anndata, genes, dp)
        # We don't do anything with the integrated anndata object so don't save?
        # self.path_full_anndata = os.path.join(self.sample_output, 'AnnData_storage', f'{self.full_name}_concat_anndata.h5ad')
        # self.concat_anndata.write(os.path.join(self.path_full_anndata))
        
