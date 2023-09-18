import scanpy as sc
from matplotlib import pyplot as plt
import math
import os

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=120, dpi_save=150, facecolor='white', color_map='Blues')

class Differential_Expression_Analysis:

    def __init__(self, adata, output_path, full_sample_name, markerpath):
        self.adata = adata
        self.output_path = output_path
        self.full_sample_name = full_sample_name
        self.markerpath = markerpath


    # perform dea with wilcoxon & bonferroni correction
    # write results, grouped by clusters found, to csv for downstream analysis
    def perform_dea(self):
        sc.tl.rank_genes_groups(self.adata, 'leiden', method='wilcoxon', corr_method='bonferroni', key='wilcoxon', pts=True, )
        sc.tl.filter_rank_genes_groups(self.adata, groupby='leiden', min_in_group_fraction=0.1, min_fold_change=1)
        sc.pl.rank_genes_groups(self.adata, sharey=False, show=False)
        plt.savefig(os.path.join(self.output_path, 'DEA', self.full_sample_name+'_DEA_wilcoxon.png'))
        self.cluster = range(len(self.adata.obs['leiden'].unique()))
        self.cluster = [str(x) for x in self.cluster]
        rank_genes_df = sc.get.rank_genes_groups_df(self.adata, group=self.cluster, key='rank_genes_groups', pval_cutoff=0.5, log2fc_min=0, log2fc_max=None, gene_symbols=None)
        rank_genes_df.to_csv(os.path.join(self.output_path, 'DEA', self.full_sample_name+'_rank_genes_df.tsv'), sep='\t', encoding='utf-8')
        # return self.cluster # remove this line?


    # read .txt file containing markergenes and put them in a dictionary for downstream
    def read_markergenes(self, adata, markerpath):
        marker_dict = {}
        with open(markerpath, "r") as filestream:
            for line in filestream:
                marker = line.strip().split(",")
                marker_dict[marker[0]] = Differential_Expression_Analysis.checkpoint(self, adata, list(marker[1:]))
        return marker_dict


    # Checks of all markergenes are present in the list, if not, removes specific geens
    def checkpoint(self, adata, markergenes_list):
        genes_passed = []
        for gene in markergenes_list:
            match = adata[:, adata.var_names.str.match('^'+gene+'$')]
            if match.n_vars > 0:
                genes_passed.append(gene)
            else:
                pass
        return genes_passed
    

    # Individual features plots --> umap
    def features_plots(self, path, markers, adata):
        for gene in markers:
            sc.pl.scatter(adata, color=gene, legend_loc='none', color_map='Blues', size=5, basis='umap', show=False)
            plt.savefig(os.path.join(path, gene+'.png'))


    # Create a big features plot --> overview of all umaps created in features_plots() 
    def big_features_plot(self, path, markers, name, adata):
        with plt.rc_context({'figure.figsize': (3, 3)}):
            sc.pl.umap(adata, color=markers, color_map='Blues', s=50, frameon=False, ncols=4, vmax='p99', show=False)
            plt.savefig(os.path.join(path, 'DEA', name, 'Big_Features_overview_'+self.full_sample_name+'.png'))


    # Create Dotplot
    def dotplots(self, path, markers, name, adata):
        sc.pl.rank_genes_groups_dotplot(adata, self.cluster, color_map='Blues',
                                              var_names=markers, show=False)
        plt.savefig(os.path.join(path, 'DEA', name, 'Dotplot_Dendogram'+self.full_sample_name+'.png'))
    

    # Create Custom big violin plot with scanpy function & matplotlib    
    def violinplots(self, path, markers, name, adata):
        fig = plt.figure(figsize=(16, 16), constrained_layout=True)
        fig.suptitle('Violin plots of marker gene set', fontsize=18)
        gridspace = fig.add_gridspec(math.ceil(len(markers)/3), 3)
        for n, ticker in enumerate(markers):
            ax = fig.add_subplot(gridspace[n])
            sc.pl.violin(adata, keys=ticker, groupby='leiden', palette='Blues', show=False, ax=ax)
            ax.set_title(ticker.upper())
        plt.savefig(os.path.join(path, 'DEA', name, 'Violin_Plots'+self.full_sample_name+'.png'))

    
    # Create heatmap
    def heatmap(self, path, markers, name, adata):
        sc.pl.heatmap(adata, var_names=markers, groupby='leiden', cmap='Blues', show=False)
        plt.savefig(os.path.join(path, 'DEA', name, 'Heatmap'+self.full_sample_name+'.png'))


    # Create tracksplot
    def tracksplot(self, path, markers, name, adata):
        sc.pl.tracksplot(adata, markers, 'leiden', log=False, color_map='Blues', show=False)
        plt.savefig(os.path.join(path, 'DEA', name, 'Tracksplot'+self.full_sample_name+'.png'))

    
    def rank_genes_dotplot_overview(self, genes):
        dp = sc.pl.rank_genes_groups_dotplot(self.adata, self.cluster, color_map='Blues',
                                              var_names=genes, return_fig=True, show=False)
        return dp

    # Create dictionairy with marker genes and name of set
    # Seperate those in name and a list of markers
    # Create directories for the names'
    # Create plots & store away accordingly
    def basic_dea_plots(self):
        # put self.adata self.adata.raw.to_adata() here.
        marker_dict = Differential_Expression_Analysis.read_markergenes(self, self.adata, self.markerpath)
        os.chdir(os.path.join(self.output_path, 'DEA'))
        for set in marker_dict.items():
            self.name, self.markers = set[0], set[1]
            os.makedirs(self.name)
            dirname = 'individual_features'
            parent = os.path.join(self.output_path, 'DEA', self.name)
            path = os.path.join(parent, dirname)
            os.makedirs(path)
            Differential_Expression_Analysis.features_plots(self, path, self.markers, self.adata)
            Differential_Expression_Analysis.big_features_plot(self, self.output_path, self.markers, self.name, self.adata)
            Differential_Expression_Analysis.dotplots(self, self.output_path, self.markers, self.name, self.adata)
            Differential_Expression_Analysis.violinplots(self, self.output_path, self.markers, self.name, self.adata)
            Differential_Expression_Analysis.heatmap(self, self.output_path, self.markers, self.name, self.adata)
            Differential_Expression_Analysis.tracksplot(self, self.output_path, self.markers, self.name, self.adata)