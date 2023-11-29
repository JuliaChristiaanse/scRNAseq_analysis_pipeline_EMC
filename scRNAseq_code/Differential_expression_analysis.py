"""
Author:         	Julia Christiaanse
Date:               29-11-2023
Script:             Class Differential_Expression_Analysis
Python version:     3.10.9
Imports & settings:
"""
import scanpy as sc
from matplotlib import pyplot as plt
import math
import os
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=120, dpi_save=150, facecolor='white', color_map='Blues')


class Differential_Expression_Analysis:
    """
    In this class a differntial expression analysis (DEA) between clusters is performed.
    """


    def __init__(self, adata, output_path, full_sample_name, markerpath):
        """
        The __init__ of the class Differential_Expression_Analysis.

        Args:
            adata (AnnData Object): The AnnData object of the sample
            output_path (String): Location of the output directory
            full_sample_name (String): The full name of the sample
            markerpath (String): Path to markergenes.txt
        """
        self.adata = adata
        self.output_path = output_path
        self.full_sample_name = full_sample_name
        self.markerpath = markerpath


    def perform_dea(self):
        """
        This function performs a DEA on the Louvain clusters defined during Sample_Analysis.py.
        A Pandas DataFrame of the DEA results is obtained and turned into a .csv file.
        The output .csv is written to disk.
        """
        sc.tl.rank_genes_groups(self.adata, 'louvain', method='wilcoxon', corr_method='bonferroni', key='wilcoxon', pts=True, )
        sc.tl.filter_rank_genes_groups(self.adata, groupby='louvain', min_in_group_fraction=0.1, min_fold_change=1)
        sc.pl.rank_genes_groups(self.adata, sharey=False, show=False)
        plt.savefig(os.path.join(self.output_path, 'DEA', self.full_sample_name+'_DEA_wilcoxon.png'))
        # Obtain a list of all (unique) clusters present within the AnnData object.
        self.cluster = range(len(self.adata.obs['louvain'].unique()))
        self.cluster = [str(x) for x in self.cluster]
        rank_genes_df = sc.get.rank_genes_groups_df(self.adata, group=self.cluster, key='rank_genes_groups', pval_cutoff=0.5, log2fc_min=0, log2fc_max=None, gene_symbols=None)
        rank_genes_df.to_csv(os.path.join(self.output_path, 'DEA', self.full_sample_name+'_rank_genes_df.tsv'), sep='\t', encoding='utf-8')


    def read_markergenes(self, adata, markerpath):
        """
        This function opens the markergenes.txt file where the marker genes are stored and runs it by self.checkpoint(), 
        before storing its contents in a dictionary and returning it.

        Args:
            adata (AnnData Object): The AnnData object of the sample
            markerpath (String): Path to markergenes.txt

        Returns:
            Dict : marker_dict, dictionary of all marker genes found in the .txt file containing the marker genes.
            The key is the category, the values are a list of genes ( both are strings).
            Three possible categories are: Neurons, Astrocytes and Interesting_Genes.
        """
        marker_dict = {}
        with open(markerpath, "r") as filestream:
            for line in filestream:
                marker = line.strip().split(",")
                # The category becomes the key, the values are a list of genes present within the AnnData object
                marker_dict[marker[0]] = Differential_Expression_Analysis.checkpoint(self, adata, list(marker[1:]))
        return marker_dict


    def checkpoint(self, adata, markergenes_list):
        """This function "checks" if a marker gene of the .txt file is in fact present within the genes of the AnnData object.
        If not, this gene is removed.

        Args:
            adata (AnnData object): The AnnData object of the sample
            markergenes_list (String): Path to .txt file containing marker genes

        Returns:
            List: genes_passed, genes that pass the checkpoint, a.k.a are present in the AnnData object
        """
        genes_passed = []
        for gene in markergenes_list:
            match = adata[:, adata.var_names.str.match('^'+gene+'$')]
            if match.n_vars > 0:
                genes_passed.append(gene)
            else:
                pass
        return genes_passed
    

    def features_plots(self, path, markers, adata):
        """
        Create features plots of the sample and every marker gene.
        Args:
            path (String): Path to output directory
            markers (List): List of marker genes to plot
            adata (AnnData): The AnnData object of the sample
        """
        for gene in markers:
            sc.pl.scatter(adata, color=gene, legend_loc='none', color_map='Blues', size=5, basis='umap', show=False)
            plt.savefig(os.path.join(path, gene+'.png'))


    def big_features_plot(self, path, markers, name, adata):
        """
        Creates a big overview of all features plots.

        Args:
            path (String): Path to output directory
            markers (List): List of marker genes to plot
            name (String): Name of the sample
            adata (AnnData): The AnnData object of the sample
        """
        with plt.rc_context({'figure.figsize': (3, 3)}):
            sc.pl.umap(adata, color=markers, color_map='Blues', s=50, frameon=False, ncols=4, vmax='p99', show=False)
            plt.savefig(os.path.join(path, 'DEA', name, 'Big_Features_overview_'+self.full_sample_name+'.png'))


    def dotplots(self, path, markers, name, adata):
        """
        Creates a dotplot of the marker genes to visualize their expression levels per cluster.

        Args:
            path (String): Path to output directory
            markers (List): List of marker genes to plot
            name (String): Name of the sample
            adata (AnnData): The AnnData object of the sample
        """
        sc.pl.rank_genes_groups_dotplot(adata, self.cluster, color_map='Blues',
                                              var_names=markers, show=False)
        plt.savefig(os.path.join(path, 'DEA', name, 'Dotplot_Dendogram'+self.full_sample_name+'.png'))
    

    def violinplots(self, path, markers, name, adata):
        """
        Creates a custom made violin plot with Matplotlib.

        Args:
            path (String): Path to output directory
            markers (List): List of marker genes to plot
            name (String): Name of the sample
            adata (AnnData): The AnnData object of the sample
        """
        fig = plt.figure(figsize=(16, 16), constrained_layout=True)
        fig.suptitle('Violin plots of marker gene set', fontsize=18)
        # Define the gridspace so plots line out nicely.
        gridspace = fig.add_gridspec(math.ceil(len(markers)/3), 3)
        # For every marker, make a figure and place it in the grid.
        for n, ticker in enumerate(markers):
            ax = fig.add_subplot(gridspace[n])
            sc.pl.violin(adata, keys=ticker, groupby='louvain', palette='Blues', show=False, ax=ax)
            ax.set_title(ticker.upper())
        plt.savefig(os.path.join(path, 'DEA', name, 'Violin_Plots'+self.full_sample_name+'.png'))

    
    def heatmap(self, path, markers, name, adata):
        """
        Creates a heatmap.

        Args:
            path (String): Path to output directory
            markers (List): List of marker genes to plot
            name (String): Name of the sample
            adata (AnnData): The AnnData object of the sample
        """
        sc.pl.heatmap(adata, var_names=markers, groupby='louvain', cmap='Blues', show=False)
        plt.savefig(os.path.join(path, 'DEA', name, 'Heatmap'+self.full_sample_name+'.png'))


    def tracksplot(self, path, markers, name, adata):
        """Creates a tracksplot.

        Args:
            path (String): Path to output directory
            markers (List): List of marker genes to plot
            name (String): Name of the sample
            adata (AnnData Object): The AnnData object of the sample
        """
        sc.pl.tracksplot(adata, markers, 'louvain', log=False, color_map='Blues', show=False)
        plt.savefig(os.path.join(path, 'DEA', name, 'Tracksplot'+self.full_sample_name+'.png'))

    
    def rank_genes_dotplot_overview(self, genes):
        """
        This dotplot is created for the cell selection step.

        Args:
            genes (List): Small list of marker panels.

        Returns:
            Dotplot Figure: dp, the dotplot figure
        """
        dp = sc.pl.rank_genes_groups_dotplot(self.adata, self.cluster, color_map='Blues',
                                              var_names=genes, return_fig=True, show=False)
        return dp


    def basic_dea_plots(self):
        """
        This function calls the other functions to:
            I: Create the marker dict
            II: Loop over every key & value category pair in the dict
            III: Creates plot per category of the dict
        """
        marker_dict = Differential_Expression_Analysis.read_markergenes(self, self.adata, self.markerpath)
        os.chdir(os.path.join(self.output_path, 'DEA'))
        for set in marker_dict.items():
            self.name, self.markers = set[0], set[1]
            # Make directories for every marker gene to store the individual features plots in
            os.makedirs(self.name)
            dirname = 'individual_features'
            parent = os.path.join(self.output_path, 'DEA', self.name)
            path = os.path.join(parent, dirname)
            # Go back to a lower level in the folder structure and make the rest of the plots
            os.makedirs(path)
            Differential_Expression_Analysis.features_plots(self, path, self.markers, self.adata)
            Differential_Expression_Analysis.big_features_plot(self, self.output_path, self.markers, self.name, self.adata)
            Differential_Expression_Analysis.dotplots(self, self.output_path, self.markers, self.name, self.adata)
            Differential_Expression_Analysis.violinplots(self, self.output_path, self.markers, self.name, self.adata)
            Differential_Expression_Analysis.heatmap(self, self.output_path, self.markers, self.name, self.adata)
            Differential_Expression_Analysis.tracksplot(self, self.output_path, self.markers, self.name, self.adata)