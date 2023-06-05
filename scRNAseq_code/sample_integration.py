import scanpy as sc
import os
import pandas as pd
import harmonypy as hm
import anndata
from matplotlib import pyplot as plt
from individual_sample_analysis import Sample_Analysis

class Data_integration:

    def __init__(self, co_sample, mono_sample, output_dir):
        self.co_sample = co_sample
        self.mono_sample = mono_sample
        self.output_dir = output_dir
        self.full_name = self.co_sample+'_'+self.mono_sample
        self.concat_anndata = self.concatinate_adatas()
        self.sample_output = os.path.join(output_dir, 'integration_dir_'+self.co_sample+self.mono_sample)
        self.run()


    def makedirs(self):
        os.makedirs(os.path.join(self.sample_output))
        os.chdir(os.path.join(self.sample_output))
        os.makedirs('DEA')


    def concatinate_adatas(self):
        co_adata = sc.read_h5ad(os.path.join(self.output_dir, self.co_sample,'AnnData_storage', 'PCA_'+self.co_sample+'.h5ad'))
        mono_adata = sc.read_h5ad(os.path.join(self.output_dir, self.mono_sample,'AnnData_storage', 'PCA_'+self.mono_sample+'.h5ad'))
        co_adata.obs['batch'] = self.co_sample
        mono_adata.obs['batch'] = self.mono_sample
        return anndata.concat([co_adata, mono_adata], index_unique='_', axis=0, join='inner', 
                              merge=None, uns_merge=None, label=None, keys=None,
                                fill_value=None, pairwise=None)
    

    def run_harmony(self):
        data_mat = self.concat_anndata.obsm['X_pca']
        meta_data = self.concat_anndata.obs
        vars_use = ['batch']
        title = self.co_sample+'_'+self.mono_sample+'harmonypy_UMAP'
        ho = hm.run_harmony(data_mat, meta_data, vars_use,
                             theta=None, lamb=None, sigma=0.1,
                             nclust=None, tau=0, block_size=0.5,
                             max_iter_harmony=10, max_iter_kmeans=20,
                             epsilon_cluster=1e-5, epsilon_harmony=1e-4,
                             plot_convergence=True, verbose=True, reference_values=None,
                             cluster_prior=None, random_state=0)
        adjusted_pcs = pd.DataFrame(ho.Z_corr).T
        self.concat_anndata.obsm['X_pca']=adjusted_pcs.values
        sc.pp.neighbors(self.concat_anndata)
        sc.tl.leiden(self.concat_anndata, resolution=0.5) # Wanna keep resolution consistant throughout the pipeline? Found 1 a bit too correlated. 0.5 too sparse?
        sc.tl.umap(self.concat_anndata)
        sc.pl.umap(self.concat_anndata, color=['batch', 'leiden'], palette='tab20', color_map='magma', title=title, show=False)
        plt.savefig(os.path.join(self.sample_output, title+'.png'))


    # def DEA_analysis(self, adata, sample_output, full_name):
    #     sc.tl.rank_genes_groups(self.DE_adata, 'leiden', method='wilcoxon', corr_method='bonferroni', key='wilcoxon', pts=True)
    #     sc.tl.filter_rank_genes_groups(self.DE_adata, groupby='leiden', min_in_group_fraction=0.1, min_fold_change=1)
    #     sc.pl.rank_genes_groups(self.DE_adata, sharey=False)
    #     plt.savefig(os.path.join(sample_output, 'DEA', 'DEA_wilcoxon'+full_name+'.png'))
    #     cluster = range(len(self.DE_adata.obs['leiden'].unique()))
    #     cluster = [str(x) for x in cluster]
    #     rank_genes_df = sc.get.rank_genes_groups_df(self.DE_adata, group=cluster, key='rank_genes_groups', pval_cutoff=0.5, log2fc_min=0, log2fc_max=None, gene_symbols=None)
    #     rank_genes_df.to_csv(os.path.join(sample_output, full_name+'rank_genes.tsv', sep='\t', encoding='utf-8'))
        


    def run(self):
        self.makedirs()
        self.run_harmony()
        self.adata_DE = self.concat_anndata.raw.to_adata()
        self.cluster = Sample_Analysis.calculate_DE_genes(self, self.adata_DE, self.sample_output, self.full_name)
        # Sample_Analysis.perform_DEA(self, self.adata_DE, self.sample_output)