import scanpy as sc
import os
import pandas as pd
import harmonypy as hm
import anndata
import numpy as np
from matplotlib import pyplot as plt
from differential_expression_analysis import Differential_Expression_Analysis as dea

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=120, dpi_save=150, facecolor='white', color_map='Blues')

class Sample_integration:

    def __init__(self, co_sample, mono_sample, output_dir, markerpath):
        self.co_sample = co_sample
        self.mono_sample = mono_sample
        self.output_dir = output_dir
        self.full_name = self.co_sample+'_'+self.mono_sample
        self.concat_anndata = self.concatinate_adatas()
        self.sample_output = os.path.join(output_dir, 'integration_dir_'+self.co_sample+self.mono_sample)
        self.markerpath = markerpath
        self.run()


    # Directories are created
    def makedirs(self):
        os.makedirs(os.path.join(self.sample_output))
        os.chdir(os.path.join(self.sample_output))
        os.makedirs('DEA')
        os.makedirs('UMAPS')
        os.makedirs('CSVtjes')


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
        sc.pp.neighbors(adataobj, n_pcs=20) # add 20 components
        sc.tl.leiden(adataobj, resolution=0.5) # Wanna keep resolution consistant throughout the pipeline? Found 1 a bit too correlated. 0.5 too sparse?
        sc.tl.umap(adataobj)
        fig, axs = plt.subplots(2, 2, figsize=(10,8),constrained_layout=True)
        Cocult = adataobj[adataobj.obs['batch'] == self.co_sample]
        Sample = adataobj[adataobj.obs['batch'] == self.mono_sample]
        sc.pl.umap(adataobj, color="batch", title=title, ax=axs[0,0], show=False)
        sc.pl.umap(adataobj, color="leiden", title="Leidenalg UMAP", ax=axs[0,1], show=False)
        sc.pl.umap(Sample, color="leiden", title=f"{self.mono_sample} sample only", ax=axs[1,0], show=False)
        sc.pl.umap(Cocult, color="leiden", title=f"{self.co_sample} sample only", ax=axs[1,1], show=False)
        plt.savefig(os.path.join(self.sample_output, 'UMAPS', title+'.png'))


    # cell selection is performed
    # EDIT: the creation of the .txt files and the .tsv files will be removed later!
    def perform_selection(self):
        if self.mono_sample == 'BL_A':
            marker_genes = ['VIM', 'FABP7', 'S100B']
        else:
            marker_genes = ['DCX', 'MAP2', 'NEUROG2']
        goede_anndatas = []
        with open(os.path.join(self.sample_output,'CSVtjes','clusters_goed_fout.txt'), 'w') as f:
            f.write('GOEDFOUT\tlen_pos\tlen_totaal\tpercentage\tcluster_nummer\n')
            
            for cluster in np.unique(self.adata_DE.obs['leiden']):
                subcluster = self.adata_DE[self.adata_DE.obs['leiden']==cluster]
                justmarkers = subcluster[:, marker_genes]
                df = justmarkers.to_df()
                df.to_csv(os.path.join(self.sample_output, 'CSVtjes', f'{cluster}{self.mono_sample}_BL_C.tsv'), sep='\t')
                pos = df[(df[marker_genes[0]] > 0) & (df[marker_genes[1]] > 0) & (df[marker_genes[2]] > 0)]
                pos.to_csv(os.path.join(self.sample_output, 'CSVtjes', f'{cluster}POS_{self.mono_sample}_BL_C.tsv'), sep='\t')
                som = int((len(pos)/len(df)*100))
                if som >= 20:
                    goede_anndatas.append(subcluster)
                    f.write('GOED'+'\t'+str(len(pos))+'\t\t'+str(len(df))+'\t'+str(som)+'\t'+cluster+'\n')
                else:
                    f.write('FOUT'+'\t'+str(len(pos))+'\t\t'+str(len(df))+'\t'+str(som)+'\t'+cluster+'\n')
            print("LENGTE**** "+str(len(goede_anndatas)))
        return goede_anndatas


    # This function checks whether the cells of the clusters with a high enough expression level
    # off all 3 marker genes are present in the solo samples
    def check_barcodes(self, goede_anndatas):
        cluster_barcodes, check = [], []
        for ga in goede_anndatas:
            a = list(ga.obs.index)
            for x in a:
                cluster_barcodes.append(x)
        adata_DE_barcodes = list(self.adata_DE.obs.index)
        for b in adata_DE_barcodes:
            for d in cluster_barcodes:
                if b == d:
                    check.append(d)
        return check


    # the selection found in check_barcodes() is applied to the two mono-samples
    def apply_selection(self, check):
        mask = self.adata_DE.obs_names.isin(check)
        filtered_adata = self.adata_DE[mask]
        BL_mono_barcodes = [barcode + f"-{self.mono_sample}" for barcode in list(self.mono_adata.obs_names)]
        BL_co_barcodes = [barcode + f"-{self.co_sample}" for barcode in list(self.co_adata.obs_names)]
        self.mono_adata.obs['edited_barcodes'] = BL_mono_barcodes
        self.co_adata.obs['edited_barcodes'] = BL_co_barcodes
        self.filtered_mono_adata = self.mono_adata[self.mono_adata.obs['edited_barcodes'].isin(list(filtered_adata[filtered_adata.obs['batch']==f'{self.mono_sample}'].obs_names))]
        self.filtered_co_adata = self.co_adata[self.co_adata.obs['edited_barcodes'].isin(list(filtered_adata[filtered_adata.obs['batch']==f'{self.co_sample}'].obs_names))]


    # This function concats the sub-selected solo anndata objects of the mono and co-culture
    # it performs a pca
    # run_harmony is called again with a new plot title
    def re_run_integration(self):
        adata = anndata.concat([self.filtered_co_adata, self.filtered_mono_adata], index_unique=None, axis=0,
                                join='inner', merge=None, uns_merge=None, label=None, keys=None, 
                                fill_value=None, pairwise=None)
        sc.pp.pca(adata, n_comps=50)
        self.run_harmony(adata, 'UMAPS_after_cell_selection')


    # This function calls all other functions & runs them
    def run(self):
        self.makedirs()
        self.run_harmony(self.concat_anndata, "UMAPS_before_cell_selection")
        self.adata_DE = self.concat_anndata.raw.to_adata()
        deado = dea(self.adata_DE, self.sample_output, self.full_name, self.markerpath)
        deado.perform_dea()
        deado.basic_dea_plots()
        good_anndatas = self.perform_selection()
        checked_barcodes = self.check_barcodes(good_anndatas)
        # print(checked_barcodes)
        # print(len(checked_barcodes))
        self.apply_selection(checked_barcodes)
        self.re_run_integration()
        
