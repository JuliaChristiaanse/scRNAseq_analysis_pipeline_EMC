import os
from matplotlib import pyplot as plt
import scanpy as sc
import math
import multiprocessing as mp

# Set scanpy settings, turned figure settings off for now
# To do: set 'magma' as standard color map for the pipeline.
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
# sc.settings.set_figure_params(dpi=150, facecolor='white')

class Sample_Analysis:

    # Init that sets variables, calls create_AnnData() method and runs all code below with run() method
    def __init__(self, sample_dir, sample_name, output_dir):
        self.sample_dir = sample_dir
        self.sample_name = sample_name
        self.output_dir = output_dir
        self.full_path = os.path.join(self.sample_dir, self.sample_name, 'filtered_feature_bc_matrix')
        self.sample_output = os.path.join(self.output_dir, self.sample_name)
        self.adata = self.create_AnnData()
        self.run()


    # Create master folder/sub folder structure
    def create_folders(self):
        os.makedirs(self.sample_output)
        os.chdir(self.sample_output)
        os.makedirs('QC')
        os.makedirs('PCA')
        os.makedirs('DEA')
        os.makedirs('Clusters')
        os.makedirs('AnnData_storage')

    
    # Create AnnData oject, store in cache for faster future access
    def create_AnnData(self):
        print(f'Creation of {self.sample_name} AnnData object pending...')
        return sc.read_10x_mtx(path=self.full_path, var_names='gene_symbols', cache=True)


    # filter_cells and calculate qc based on mitochondrial RNA 
    def calculate_qc(self):
        sc.pp.filter_cells(self.adata, min_genes=700)   # Keep cells with at least 700 genes expressed
        sc.pp.filter_genes(self.adata, min_cells=3)     # Keep genes that are expressed in atleast 3 cells
        self.adata.uns[self.sample_name] = self.sample_name
        self.adata.var['mt'] = self.adata.var_names.str.startswith('MT-')    # Filter out MT- RNA
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        

    # Save QC plots to correct folders    
    def create_QC_plots(self):    
        sc.pl.violin(self.adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                     scale='width', color='#9ADCFF', multi_panel=True, show=False)
        plt.savefig(os.path.join(self.sample_output, 'QC', 'QC_nfeatn_count_percMT_'+self.sample_name+'.png'))
        sc.pl.scatter(self.adata, 'pct_counts_mt', 'total_counts', color='pct_counts_mt',
                      title='Percentage of MT RNA in counts '+self.sample_name, color_map='Blues', show=False)
        plt.savefig(os.path.join(self.sample_output, 'QC', 'Pct_MT_RNA_in_counts'+self.sample_name+'.png'))
        sc.pl.scatter(self.adata, 'n_genes_by_counts', 'total_counts', color='n_genes_by_counts',
                      title='Amount of genes by counts '+self.sample_name, color_map='Blues', show=False)
        plt.savefig(os.path.join(self.sample_output, 'QC', 'ngenes_by_counts'+self.sample_name+'.png'))

    
    # Subset, Regress out bad stuff and normalize, find variablefeatures and do some other stuff like SCTransform
    def normalization_HVG(self):
        self.adata = self.adata[self.adata.obs.pct_counts_mt < 20, :] # EDIT 24-4-2023: slice away MT- counts that are too high
        sc.pp.normalize_total(self.adata, target_sum=1e4) # EDIT: removed this: still included highly expressed data for now
        sc.pp.log1p(self.adata)
        # sc.pp.highly_variable_genes(self.adata, flavor='cell_ranger', n_top_genes=2000, subset=False)
        sc.pp.highly_variable_genes(self.adata, flavor='cell_ranger', subset=False) # EDIT: added this line for testing
        sc.pl.highly_variable_genes(self.adata, show=False)
        plt.savefig(os.path.join(self.sample_output, 'QC', 'Highly_variable_genes_'+self.sample_name+'.png'))
        self.adata.raw = self.adata
        self.adata = self.adata[:, self.adata.var.highly_variable] # Actually do the slicing
        sc.pp.regress_out(self.adata, ['total_counts', 'pct_counts_mt']) # regress out sequencing depth and % MT-RNA
        sc.pp.scale(self.adata)


    # Run a PCA and plot output
    def run_PCA(self):
        sc.tl.pca(self.adata, n_comps=50)       # Create 50 components
        sc.pl.pca(self.adata, annotate_var_explained=True, na_color='#9ADCFF', title=f'PCA scoringsplot {self.sample_name}', show=False)
        plt.savefig(os.path.join(self.sample_output, 'PCA', 'PCA_Scores_'+self.sample_name+'.png'))
        sc.pl.pca_loadings(self.adata, components=[1,2], show=False) # Only show first 2.
        plt.savefig(os.path.join(self.sample_output, 'PCA', 'PCA_loadings_plot_'+self.sample_name+'.png'))
        sc.pl.pca_variance_ratio(self.adata, n_pcs=19, show=False)
        plt.savefig(os.path.join(self.sample_output, 'PCA', 'PCA_Variance_elbow_'+self.sample_name+'.png'))


    # Perform unsupervised clustering on the un-annotated sample
    def unsupervised_clustering(self):
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20) # calculate neighbors with 20 npc's
        sc.tl.umap(self.adata)
        sc.tl.leiden(self.adata, resolution=0.5)      # resolution default for scanpy = 1. resolution used in seurat = 0.5.
        title=f'Unsupervised Leiden Cluster {self.sample_name}'
        sc.pl.umap(self.adata, color=['leiden'], title=title, legend_loc='on data', legend_fontsize=8, show=False)
        plt.savefig(os.path.join(self.sample_output, 'Clusters', 'Unsupervised_UMAP_'+self.sample_name+'.png'))


    # Calculate differentially expressed genes
    def calculate_DE_genes(self):
        self.adata_DE = self.adata.raw.to_adata()
        sc.tl.rank_genes_groups(self.adata_DE, 'leiden', method='wilcoxon', corr_method='bonferroni', key='wilcoxon', pts=True, )
        sc.tl.filter_rank_genes_groups(self.adata_DE, groupby='leiden', min_in_group_fraction=0.1, min_fold_change=1)
        sc.pl.rank_genes_groups(self.adata_DE, sharey=False, show=False)
        plt.savefig(os.path.join(self.sample_output, 'DEA', 'DEA_wilcoxon'+self.sample_name+'.png'))
        self.cluster = range(len(self.adata_DE.obs['leiden'].unique()))
        self.cluster = [str(x) for x in self.cluster]
        rank_genes_df = sc.get.rank_genes_groups_df(self.adata_DE, group=self.cluster, key='rank_genes_groups', pval_cutoff=0.5, log2fc_min=0, log2fc_max=None, gene_symbols=None)
        rank_genes_df.to_csv(os.path.join(self.sample_output, 'DEA', 'rank_genes_df.tsv'), sep='\t', encoding='utf-8')


    # Simple function that stores the markergenes, might be deleted later
    def store_markergenes(self):
        markergenes = {
            'astrocyte_interest' : ["GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1", "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1", "APOE"], 
            'astrocyte_mature' : ["CD44", "FABP7", "VIM", "SOX9", "TOP2A", "S100B", "GJA", "SLC1A3", "IGFBP7", "ALDH1L1", "APOE"],
            'neuron_interest' : ["TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1"],
            'neuron_mature' : ["NEUROG2", "DCX", "MAP2", "RBFOX3", "SYN1", "SNAP25", "SYT1", "APOE"],
            'schema_psych_interest' : ["SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4", "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11"],
            'sloan_2017_interest' : ["AQP4", "ALDH1L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B", "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1", "HOPX", "FAM107A", "AGT"],
            'interneuron_interest' : ["SST", "PVALB", "GAD1"]
            }
        self.perform_DEA(markergenes)
        
    # Perform differential expression analysis
    # Look into running this in parallel
    def perform_DEA(self, markergenes):
        os.chdir(os.path.join(self.sample_output, 'DEA'))
        for set in markergenes.items():
            self.name, self.markers = set[0], self.checklist(set[1])
            os.makedirs(self.name)
            dirname = 'individual_features'
            parent = os.path.join(self.sample_output, 'DEA', self.name)
            path = os.path.join(parent, dirname)
            os.makedirs(path)
            self.features_plots(path)
            self.big_features_plot()
            self.dotplots()
            self.violinplots()
            self.heatmap()
            self.tracksplot()


    def checklist(self, markergenes_list):
        genes_passed = []
        for gene in markergenes_list:
            match = self.adata_DE[:, self.adata_DE.var_names.str.match('^'+gene+'$')]
            if match.n_vars > 0:
                genes_passed.append(gene)
            else:
                pass
        return genes_passed


    def features_plots(self, path):
        for gene in self.markers:
            sc.pl.scatter(self.adata_DE, color=gene, legend_loc='none', color_map='Blues', size=5, basis='umap', show=False)
            plt.savefig(os.path.join(path, gene+'.png'))


    # Create a big features plot --> overview of all umaps 
    def big_features_plot(self):
        with plt.rc_context({'figure.figsize': (3, 3)}):
            sc.pl.umap(self.adata_DE, color=self.markers, color_map='Blues', s=50, frameon=False, ncols=4, vmax='p99', show=False)
            plt.savefig(os.path.join(self.sample_output, 'DEA', self.name, 'Big_Features_overview_'+self.sample_name+'.png'))


    # Create Dotplot
    def dotplots(self):
        sc.pl.rank_genes_groups_dotplot(self.adata_DE, self.cluster, color_map='Blues', var_names=self.markers, show=False)
        plt.savefig(os.path.join(self.sample_output, 'DEA', self.name, 'Dotplot_Dendogram'+self.sample_name+'.png'))


    # Create Custom big violin plot with scanpy function & matplotlib    
    def violinplots(self):
        fig = plt.figure(figsize=(16, 16), constrained_layout=True)
        fig.suptitle('Violin plots of marker gene set', fontsize=18)
        gridspace = fig.add_gridspec(math.ceil(len(self.markers)/3), 3)
        for n, ticker in enumerate(self.markers):
            ax = fig.add_subplot(gridspace[n])
            sc.pl.violin(self.adata_DE, keys=ticker, groupby='leiden', palette = 'Blues', show=False, ax=ax)
            ax.set_title(ticker.upper())
        plt.savefig(os.path.join(self.sample_output, 'DEA', self.name, 'Violin_Plots'+self.sample_name+'.png'))


    # Create heatmap
    def heatmap(self):
        sc.pl.heatmap(self.adata_DE, var_names=self.markers, groupby='leiden', cmap='Blues', show=False)
        plt.savefig(os.path.join(self.sample_output, 'DEA', self.name, 'Heatmap'+self.sample_name+'.png'))


    # Create tracksplot
    def tracksplot(self):
        sc.pl.tracksplot(self.adata_DE, self.markers, 'leiden', log=False, color_map='Blues', show=False)
        plt.savefig(os.path.join(self.sample_output, 'DEA', self.name, 'Tracksplot'+self.sample_name+'.png'))


    # run all functions at once and write AnnData
    def run(self):
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
        self.calculate_DE_genes()
        self.store_markergenes()
        self.adata.write(os.path.join(self.sample_output, 'AnnData_storage', self.sample_name+'.h5ad'))
