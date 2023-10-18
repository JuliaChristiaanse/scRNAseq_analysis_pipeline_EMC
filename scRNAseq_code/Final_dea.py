import scanpy as sc
import pandas as pd
import cellrank as cr
from matplotlib import pyplot as plt
from bioinfokit import visuz
from itertools import chain
import os
from Differential_expression_analysis import Differential_Expression_Analysis

class Final_Dea:

    def __init__(self, subset_adata_path, output_dir, markerpath):
        self.adata = sc.read(subset_adata_path)
        self.output_dir = output_dir
        self.markerpath = markerpath
        samples = self.adata.obs['batch'].unique()
        self.sample_name = samples[0]+'_'+samples[1]
        self.sample_output = os.path.join(self.output_dir, f'{self.sample_name}_dea_cellrank_output')
        self.run()

    
    def makedirs(self):
        os.makedirs(self.sample_output)
        os.chdir(self.sample_output)
        os.makedirs('dea_output')
        os.makedirs('cellrank_output')
        os.makedirs('anndata_storage')


    def color_picker(self, marker_group):
        if marker_group == 'astrocytes':
            color='Blues'
        elif marker_group == 'neurons':
            color='Blues'
        else:
            color='Blues'
        return color


    def dea(self):
        self.adata_DE = self.adata.raw.to_adata()
        sc.tl.rank_genes_groups(self.adata_DE, 'batch', groups=['BL_C'],
                        method='wilcoxon', corr_method='bonferroni', pts=True)
        self.rank_genes_df = sc.get.rank_genes_groups_df(self.adata_DE, group=['BL_C'], key='rank_genes_groups', pval_cutoff=None, log2fc_min=None, log2fc_max=None, gene_symbols=None)
        self.rank_genes_df.to_csv(os.path.join(self.sample_output, 'dea_output', f'{self.sample_name}_rank_genes_df.tsv'), sep='\t', encoding='utf-8')
        self.marker_dict = Differential_Expression_Analysis.read_markergenes(self, self.adata_DE, self.markerpath)
        os.chdir(os.path.join(self.sample_output, 'dea_output'))
        for set in self.marker_dict.items():
            self.marker_group, self.markers = set[0], set[1]
            os.makedirs(self.marker_group)
            with plt.rc_context({'figure.figsize': (3, 3)}):
                sc.pl.umap(self.adata_DE, color=self.markers, color_map=self.color_picker(self.marker_group), s=50, frameon=False,
                            ncols=4, vmax='p99', show=False)
                plt.savefig(os.path.join(self.sample_output, 'dea_output', self.marker_group,
                                          f'Features_overview_{self.sample_name}_{self.marker_group}.png'))

    def volcano_plot(self):
        os.chdir(os.path.join(self.sample_output, 'dea_output'))
        os.makedirs('volcano')
        os.chdir(os.path.join(self.sample_output, 'dea_output', 'volcano'))
        # access genes to show:
        genes_to_show = self.marker_dict.values()
        genes_to_show = tuple(chain(*genes_to_show))
        print(genes_to_show)
        visuz.GeneExpression.volcano(df=self.rank_genes_df, lfc="logfoldchanges", pv="pvals_adj", geneid="names", genenames=genes_to_show,
                                     gstyle=2, xlm=(-10,11,2), ylm=(-10, 300,50), figname=f'{self.sample_name}_volcano', figtype='png', axtickfontsize=10,
                                     lfc_thr=(1, 1), pv_thr=(0.05, 0.05), sign_line=False, color=('blue', 'grey', '#9ADCFF'), show=False)
        

    def run(self):
        self.makedirs()
        self.dea()
        self.volcano_plot()