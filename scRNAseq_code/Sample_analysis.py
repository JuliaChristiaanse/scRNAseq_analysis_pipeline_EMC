"""
Author:         	Julia Christiaanse
Date:               29-11-2023
Script:             Class Sample_Analysis
Python version:     3.10.9
Imports & settings:
"""
import os
from matplotlib import pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
from Differential_expression_analysis import Differential_Expression_Analysis as dea
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=150, facecolor='white')


class Sample_Analysis:
    """
    Class Sample_Analysis performs the first 3 steps of the pipeline, which include preprocessing, (including normalization, logtransformation, filtering & doublet detection),
    clustering, and cluster-based differential expression analysis (DEA).
    """ 


    def __init__(self, sample_dir, sample_name, output_dir, markerpath):
        """
        __init__ of the class Sample_Analysis.
        

        Args:
            sample_dir (String): Location of the 10X Genomics output folder containing all samples
            sample_name (String): Name of the individual sample
            output_dir (String): Path to the master output directory
            markerpath (String): Path to markergenes.txt
        """        
        self.sample_dir = sample_dir
        self.sample_name = sample_name
        self.output_dir = output_dir
        self.full_path = os.path.join(self.sample_dir, self.sample_name, 'filtered_feature_bc_matrix')
        self.sample_output = os.path.join(self.output_dir, self.sample_name)
        self.adata = self.create_AnnData()
        self.markerpath = markerpath
        self.run()


    def create_folders(self):
        """
        Generate necessary sub-folders within sample_output.
        """        
        os.makedirs(self.sample_output)
        os.chdir(self.sample_output)
        os.makedirs('QC')
        os.makedirs('PCA')
        os.makedirs('DEA')
        os.makedirs('Clusters')
        os.makedirs('AnnData_storage')
        os.makedirs('AnnData_test')

    
    def create_AnnData(self):
        """
        Creates an AnnData object of the sample, 
        using self.full_path to find the 10X Genomics output of the sample.

        Returns:
            AnnData : An AnnData object of the sample
        """        
        print(f'Creation of {self.sample_name} AnnData object pending...')
        return sc.read_10x_mtx(path=self.full_path, var_names='gene_symbols', cache=True)


    def calculate_qc(self):
        """
        Performs preprocessing on the AnnData object
        Filter cells from the AnnData object and calculate quality control based on amount of mitochondrial RNA. 
        """        
        sc.pp.filter_cells(self.adata, min_genes=700)   # Keep cells with at least 700 genes expressed
        sc.pp.filter_genes(self.adata, min_cells=3)     # Keep genes that are expressed in atleast 3 cells
        self.adata.uns[self.sample_name] = self.sample_name
        self.adata.var['mt'] = self.adata.var_names.str.startswith('MT-')  # indicate mitochondrial RNA
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        

    def create_QC_plots(self):
        """
        Create plots depicting the quality control of the AnnData object.
        Calculates and prints the correlation coefficient to the terminal.
        Figures are saved in the correct sub-folders.
        """        
        # Recent Scanpy bug caused error in the creation of a violin plot, commented out untill it's fixed.
        # sc.pl.violin(self.adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
        #              scale='width', color='#9ADCFF', show_legend=False, multi_panel=True, show=False)
        # plt.savefig(os.path.join(self.sample_output, 'QC', 'QC_nfeatn_count_percMT_'+self.sample_name+'.png'))

        sc.pl.scatter(self.adata, 'pct_counts_mt', 'total_counts', color='pct_counts_mt',
                      title="title", color_map='Blues', show=False)
        plt.savefig(os.path.join(self.sample_output, 'QC', 'Pct_MT_RNA_in_counts'+self.sample_name+'.png'))

        # Compute correlation coefficient, which is only printed to the user, not stored.
        x = self.adata.obs['n_genes_by_counts']
        y = self.adata.obs['total_counts']
        corr_coef = np.corrcoef(x, y)
        sc.pl.scatter(self.adata, 'n_genes_by_counts', 'total_counts', color='n_genes_by_counts',
                      title="", color_map='Blues', show=False)
        print(self.sample_name, f'R-squared: {round(corr_coef[0][1], 2):.2f}')
        plt.savefig(os.path.join(self.sample_output, 'QC', 'ngenes_by_counts'+self.sample_name+'.png'))


    def detect_doublets(self):
        """
        Detects doublet using the Scrublet package, using default settings.
        Stores score distribution figure to the correct sub-folder.
        """        
        # Note from 15-11-2023: with default setting normalize variance = True, 
        # normalize the data such that each gene has a variance of 1.
        # This normalization DOES NOT impact adata.X.
        # It's adviced to keep this setting as is.
        sc.external.pp.scrublet(self.adata, adata_sim=None, batch_key=None, sim_doublet_ratio=2.0, expected_doublet_rate=0.06, stdev_doublet_rate=0.02, synthetic_doublet_umi_subsampling=1.0, knn_dist_metric='euclidean',
                         normalize_variance=True, log_transform=False, mean_center=True, n_prin_comps=30, use_approx_neighbors=True, get_doublet_neighbor_parents=False, n_neighbors=None, threshold=None,
                           verbose=True, copy=False, random_state=0) # auto-set threshold. Do NOT manually generate doublets.
        sc.external.pl.scrublet_score_distribution(self.adata, show=False)
        plt.savefig(os.path.join(self.sample_output, 'QC', 'doublet_score_distribution_'+self.sample_name+'.png'))
        # rename 'predicted_doublet' column 'doublet_info'
        self.adata.obs['doublet_info'] = self.adata.obs["predicted_doublet"].astype(str)
        # Strore detected doublets in different AnnData object
        self.doublets_found = self.adata[self.adata.obs['doublet_info']=='True',:]
        print(self.doublets_found.var_names)
        # Remove the detected doublets from adata.
        self.adata = self.adata[self.adata.obs['doublet_info'] == 'False',:]


    def normalization_HVG(self):
        """
        Calls the self.detect_doublets() function to detect and remove potential doublets from adata.
        Stores rawcounts in an extra layer within adata.
        Normalizes, logtransforms, finds highly variable genes (HVG), subsets adata to HVG, regresses out and scales adata.
        Figure is stored to the correct sub-folders.
        """        
        self.detect_doublets()
        self.adata.layers['rawcounts'] = self.adata.X
        self.adata.raw = self.adata
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
        sc.pp.highly_variable_genes(self.adata, flavor='cell_ranger', subset=False)
        sc.pl.highly_variable_genes(self.adata, show=False)
        plt.savefig(os.path.join(self.sample_output, 'QC', 'Highly_variable_genes_'+self.sample_name+'.png'))
        self.adata = self.adata[:, self.adata.var.highly_variable]
        sc.pp.regress_out(self.adata, ['total_counts', 'pct_counts_mt']) # regress out sequencing depth and % MT-RNA
        sc.pp.scale(self.adata)


    def run_PCA(self):
        """
        Performs a Principal Component Analysis (PCA) with 50 components.
        Stores figures in correct sub-directory.
        """        
        sc.tl.pca(self.adata, n_comps=50, random_state=0)   
        sc.pl.pca(self.adata, annotate_var_explained=True, na_color='#9ADCFF', title="", show=False)
        plt.savefig(os.path.join(self.sample_output, 'PCA', 'PCA_Scores_'+self.sample_name+'.png'))
        sc.pl.pca_loadings(self.adata, components=[1,2], show=False)
        plt.savefig(os.path.join(self.sample_output, 'PCA', 'PCA_loadings_plot_'+self.sample_name+'.png'))
        sc.pl.pca_variance_ratio(self.adata, n_pcs=19, show=False)
        plt.savefig(os.path.join(self.sample_output, 'PCA', 'PCA_Variance_elbow_'+self.sample_name+'.png'))


    def unsupervised_clustering(self):
        """
        Performs a K-nearest neighbors (KNN) analysis, UMAP, and Louvain clustering algorithm on adata.
        Stores the UMAP in the correct sub-folder.
        Merges adata and doublets_found together to create doublets_included.
        Repeats KNN, UMAP, and Louvain clustering on doublets_included,
        then stores the UMAP in the correct sub-folder, to show the user which doublets were removed from adata.
        """        
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20, random_state=0)
        sc.tl.umap(self.adata, random_state=0)
        sc.tl.louvain(self.adata, resolution=0.5, random_state=0)
        title=f'Unsupervised Louvain Cluster {self.sample_name}'
        sc.pl.umap(self.adata, color=['louvain'], title="", legend_loc='on data', legend_fontsize=8, show=False)
        plt.savefig(os.path.join(self.sample_output, 'Clusters', 'Unsupervised_UMAP_'+self.sample_name+'.png'))
       
        # Create doublets_included, a concatinated AnnData object of doublets_found and adata.
        self.adata.obs['doublet?'] = 'No doublet'
        self.doublets_found.obs['doublet?'] = 'Doublet'
        self.doublets_included = ad.concat([self.adata, self.doublets_found], keys=['No doublet', 'Doublet'])
        # Re-run preprocessing steps on doublets_included
        sc.pp.pca(self.doublets_included, random_state=0)
        sc.pp.neighbors(self.doublets_included, n_neighbors=10, n_pcs=20, random_state=0)
        sc.tl.umap(self.doublets_included, random_state=0)
        sc.tl.louvain(self.doublets_included, resolution=0.5, random_state=0)   
        title=f'Doublets detected in dataset {self.sample_name}'
        # Create a UMAP to show the doublets vs not doublets and store in folder.
        sc.pl.umap(self.doublets_included, color=['doublet_score', 'doublet_info'], title=title, legend_loc='right margin', legend_fontsize=8, show=False)
        plt.savefig(os.path.join(self.sample_output, 'Clusters', 'Doublet_umap_'+self.sample_name+'.png'))


    def run(self):
        """
        Runs all class functions in chronological order. 
        Created for aesthetic purposes, to avoid running every function in the __init__.
        Accesses the raw counts of adata, then calls the external class Differential_Expression_Analysis to perform the DEA on the Louvain clusters.
        Adata is written to disk, in the correct sub-folder AnnData_storage.
        """        
        self.create_folders()
        self.calculate_qc()
        self.create_QC_plots() 
        self.normalization_HVG()
        df_obs = sc.get.obs_df(self.adata, ['n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'])
        df_vars = sc.get.var_df(self.adata, ['gene_ids', 'feature_types', 'n_cells', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'])
        df_obs.to_csv(self.sample_output+'/adata_obs', sep='\t', encoding='utf-8')
        df_vars.to_csv(self.sample_output+'/adata_vars', sep='\t', encoding='utf-8')
        self.run_PCA()
        self.unsupervised_clustering()
        self.adata = self.adata[self.adata.obs['doublet_info'] == 'False',:]
        self.adata_DE = self.adata.raw.to_adata()
        deado = dea(self.adata_DE, self.sample_output, self.sample_name, self.markerpath)
        deado.perform_dea()
        deado.basic_dea_plots()
        adata_storage = os.path.join(self.sample_output, 'AnnData_storage')
        self.adata.write(os.path.join(adata_storage, self.sample_name+'.h5ad'))
