"""
Author:         	Julia Christiaanse
Date:               29-11-2023
Script:             Class Final_Dea
Python version:     3.10.9
Imports:
"""
import scanpy as sc
import pandas as pd
import cellrank as cr
from matplotlib import pyplot as plt
import os
import rpy2.robjects
from rpy2.robjects import r
import rpy2.ipython.html
import rpy2	
import gseapy
from gseapy import Msigdb
import numpy as np
import networkx as nx
from gseapy.scipalette import SciPalette
r('library(Seurat)')
r('library(ggplot2)')
r('library(dplyr)')
r('library(EnhancedVolcano)')
from Differential_expression_analysis import Differential_Expression_Analysis


class Final_Dea:
    """
    Class Final DEA performs step 7B of the pipeline.
    Performs a final differential expression analysis (DEA), gene set enrichment (GSEA) analysis,
    and over-representation analysis (ORA) on the annotated AnnData object.
    During this final DEA cell cultures are compared.
    """


    def __init__(self, sample_name, annotated_adata_loc, final_dea_output_loc, output_dir, markerpath):
        """
        The __init__ of the class Final_Dea.

        Args:
            sample_name (String): Name of the sample
            annotated_adata_loc (String): Path to annotated AnnData object
            final_dea_output_loc (String): Path to subfolder to store output files
            output_dir (String): Path to the master output directory
            markerpath (String): Path to markergenes.txt
        """
        self.sample_name = sample_name
        self.sample_output = final_dea_output_loc
        self.output_dir = output_dir
        self.markerpath = markerpath
        self.makedirs(self.sample_output)
        self.adata = sc.read(os.path.join(annotated_adata_loc, f'{self.sample_name}_annotated.h5ad'))
        # Access the Molecular Signatures Database
        msig = Msigdb()
        # Collect the 'c5.go.bp' collection from the human '2023.2.Hs' version as a dictionary.
        gmt = msig.get_gmt(category='c5.go.bp', dbver='2023.2.Hs')
        # Run the DEA
        adata, marker_dict = self.dea(sample_name=self.sample_name, adata=self.adata,
                  markerpath=self.markerpath, sample_output=self.sample_output)
        # Plot results
        self.volcano_plot(sample_name=self.sample_name,
                          marker_dict=marker_dict, sample_output=self.sample_output)
        # Run GSEA and over representation analysis (ORA)
        self.gsea(adata=adata, gmt=gmt, sample_output=self.sample_output)
        self.ora(adata=adata, gmt=gmt, sample_output=self.sample_output)

    
    def makedirs(self, sample_output):
        """
        Generate necessary sub-folders within sample_output.

        Args:
            sample_output (String): Path to sample-specific sub-directory to store results in
        """
        os.makedirs(sample_output)
        os.chdir(sample_output)
        os.makedirs('dea_output')
        os.makedirs('AnnData_storage')
        os.makedirs('GSEA')
        os.makedirs('ORA')


    def dea(self, sample_name, adata, markerpath, sample_output):
        """
        Here DEA is performed which compares cell cultures.
        e.g., BL_A with BL_C & BL_N with BL_C.
        BL_C is taken as the main comparison party.
        A Pandas DataFrame containing the results is obtained.
        Additionally, the marker genes are loaded from the .txt file with the read_markergenes() function of the Differential_Expression_Analysis class.
        Overview features plots with all highlighted marker genes are created per category, and written to disk.
        The Pandas DataFrame is turned into a .csv file and written to disk.

        Args:
            sample_name (String): Name of the sample
            adata (AnnData): The AnnData object of the sample
            markerpath (String): Path to markergenes.txt
            sample_output (String): Path to sample-specific sub-directory to store results in

        Returns:
            AnnData: adata, the AnnData object of the sample
            Dict: marker_dict, the marker genes dictionary
        """
        sc.tl.rank_genes_groups(adata, 'Sample', groups=['BL_C'],
                        method='wilcoxon', corr_method='bonferroni', pts=True)
        rank_genes_df = sc.get.rank_genes_groups_df(adata, group=['BL_C'], key='rank_genes_groups', pval_cutoff=None, log2fc_min=None, log2fc_max=None, gene_symbols=None)
        rank_genes_df.to_csv(os.path.join(sample_output, 'dea_output', f'{sample_name}_rank_genes_df.tsv'), sep='\t', encoding='utf-8')
        marker_dict = Differential_Expression_Analysis.read_markergenes(self, adata, markerpath)
        os.chdir(os.path.join(sample_output, 'dea_output'))
        for set in marker_dict.items():
            marker_group, markers = set[0], set[1]
            os.makedirs(marker_group)
            with plt.rc_context({'figure.figsize': (3, 3)}):
                sc.pl.umap(adata, color=markers, color_map='Blues', s=50, frameon=False,
                            ncols=4, vmax='p99', show=False)
                plt.savefig(os.path.join(sample_output, 'dea_output', marker_group,
                                          f'Features_overview_{sample_name}_{marker_group}.png'))
        rank_genes_df.to_csv(os.path.join(sample_output, 'dea_output', f'{sample_name}_df.csv'))
        return adata, marker_dict
    

    def volcano_plot(self, sample_name, marker_dict, sample_output):
        """
        Calls R code to create a volcano plot. 
        This volcano plot is created to visualize the results of the DEA that compares the co-culture genes to genes of a the current mono-culture.
        It stores the figure in the correct folder.

        Args:
            sample_name (String): Name of the sample
            marker_dict (Dict): The marker genes dictionary
            sample_output (String): Path to sample-specific sub-directory to store results in
        """
        os.chdir(os.path.join(sample_output, 'dea_output'))
        os.makedirs('volcano')
        os.chdir(os.path.join(sample_output, 'dea_output', 'volcano'))
        genes_to_show = []
        for i in marker_dict.values():
            for j in i:
                genes_to_show.append(j)
        rpy2.robjects.globalenv['path_to_csv'] = os.path.join(sample_output, 'dea_output', f'{sample_name}_df.csv')
        rpy2.robjects.globalenv['genes_to_show'] = genes_to_show
        rpy2.robjects.globalenv['output_path'] = os.path.join(sample_output, 'dea_output', 'volcano', f'{sample_name}_volcano.png')
        r(
            '''
            data <- read.csv(path_to_csv)
            genes_to_show <- c(genes_to_show)
            genes_to_show
            data <- data %>% 
                mutate(significant = ifelse(abs(logfoldchanges) > 1 & pvals_adj < 0.05, "Significant", "Not significant"))
            keyvals <- ifelse(data$pvals_adj > 0.5, 'grey30', ifelse(data$logfoldchanges < -1, '#9ADCFF', ifelse(data$logfoldchanges > 1, '#FFBC9A', 'grey30')))
            keyvals[is.na(keyvals)] <- 'grey30'
            names(keyvals)[keyvals == '#FFBC9A'] <- 'Up-regulated genes'
            names(keyvals)[keyvals == 'grey30'] <- 'Not Significant'
            names(keyvals)[keyvals == '#9ADCFF'] <- 'Down-regulated genes'
            
            EnhancedVolcano(data,
            lab = data$names,
            selectLab = genes_to_show,
            x = 'logfoldchanges',
            y = 'pvals_adj',
            xlim =c(-27, 27),
            ylim =c(0, 350),
            xlab = bquote(~Log[2]~ 'fold change'),
            pCutoff = 0.05,
            FCcutoff = 1,
            cutoffLineType = 'solid',
            cutoffLineWidth = 0.3,
            colCustom = keyvals,
            pointSize = 1.0,
            labSize = 4.0,
            boxedLabels=TRUE,
            colAlpha = 1,
            drawConnectors = TRUE,
            typeConnectors = 'closed',
            widthConnectors = 0.5,
            lengthConnectors = 5,
            colConnectors = 'black',
            arrowheads = FALSE,
            legendLabels=c("Not significant, no L2FC", "Not significant, with L2FC", "p-value <0.05, no L2FC", "p-value <0.05, with L2FC"),
            legendPosition = 'right',
            legendLabSize = 10,
            legendIconSize = 5.0,
            gridlines.minor = FALSE,
            gridlines.major= FALSE)
            ggplot2::ggsave(file=paste0(output_path), width = 30, height = 20, units = "cm")
            '''
        )
        

    def gsea(self, adata, gmt, sample_output):
        """
        Performs GSEA with gseapy using the gene set dictionary.

        Args:
            adata (AnnData): The AnnData object of the sample
            gmt (Dict): Dictionary containing the 'c5.go.bp' gene set collection
            sample_output (String): Path to sample-specific sub-directory to store results in
        """
        # Obtain a Pandas Dataframe created as a result of the culture-based DEA.
        df = sc.get.rank_genes_groups_df(adata, group='BL_C', key='rank_genes_groups', pval_cutoff=None, log2fc_min=None, log2fc_max=None, gene_symbols=None)
        df = df.sort_values('logfoldchanges', ascending=False)
        # Create a smaller DataFrame with only gene names and log2foldchanges.
        ranking = df[['names', 'logfoldchanges']]
        # Run GSEA.
        pre_res = gseapy.prerank(rnk = ranking, gene_sets = gmt, seed = 0)
        # Obtain Pandas DataFrame with GSEA results.
        pre_res.res2d.rename(columns={pre_res.res2d.columns[8]: "gene_perc" }, inplace = True)
        # Manipulate 'gene_perc' column to perform a filtering step.
        pre_res.res2d['gene_perc'] = pre_res.res2d['gene_perc'].str.replace('%', '')
        pre_res.res2d['gene_perc'] = pre_res.res2d[['gene_perc']].astype('float')
        # Filter Pandas DataFrame to only contain gene sets that have at least 10 percent genes present from the sample in the gene sets.
        filter_gene_perc = pre_res.res2d[pre_res.res2d.gene_perc>10]
        # Sort the new Pandas Dataframe based on q-value.
        q_val_sort = filter_gene_perc.sort_values('FDR q-val')
        # Create and sort new DataFrames with positive and negative NES scores.
        pos = q_val_sort[q_val_sort.NES>0]
        neg = q_val_sort[q_val_sort.NES<0]
        pos = pos.sort_values(['FDR q-val','NOM p-val', 'NES'], ascending=[True, False, False])
        neg = neg.sort_values(['FDR q-val', 'NOM p-val','NES'], ascending=[True, True, True])
        # Subset DataFrames to the top 15 positive and negative NES scores.
        pos_NES_15 = pos.head(15).Term
        neg_NES_15 = neg.head(15).Term
        # Plot top 5 positive and negative pathways to show, and save figures.
        pos_plot = pre_res.plot(terms=pos_NES_15[0:5], ofname=os.path.join(sample_output, 'GSEA', 'positive_NES_pathways.png'))
        neg_plot = pre_res.plot(terms=neg_NES_15[0:5], ofname=os.path.join(sample_output, 'GSEA', 'negative_NES_pathways.png'))
        # Combine the top 15 positive and negative terms into a smaller DataFrame.
        combined_terms = list(pos_NES_15) + list(neg_NES_15)
        smaller_df = pre_res.res2d[pre_res.res2d['Term'].isin(combined_terms)]
        sci = SciPalette()
        NbDr = sci.create_colormap()
        # Create a dotplot from the smaller DataFrame and store it to disk.
        gseapy.dotplot(df=smaller_df,column='NES',y='Term', cmap=NbDr,figsize=(3,12), top_term=30,
                        show_ring=True, ofname=os.path.join(sample_output, 'GSEA', 'pos_neg_NES_dotplot.png'))
        # The networkplot appears to be bugged. This is commented out for now.
        #self.networkplot(df=smaller_df, sample_output=sample_output, sort_on='NES', cutoff=2, name='NES')        

        
    def networkplot(self, df, sort_on, cutoff, sample_output, name):
        """
        Creates a network plot with the top 15 positive and negative NES-score pathways.

        Args:
            df (Pandas DataFrame): DataFrame containing GSEA pathways
            sort_on (String): Column in df to map the node colors
            cutoff (Int): nodes with sort_on value < cut-off are shown. Work only for ("Adjusted P-value", "P-value", "NOM p-val", "FDR q-val")
            sample_output (String): Path to sample-specific sub-directory to store results in
            name (String): Name of column to plot
        """
        df = df.sort_values(sort_on)
        nodes, edges = gseapy.enrichment_map(df, column=sort_on, top_term=30, cutoff=cutoff)
        G = nx.from_pandas_edgelist(edges,
                            source='src_idx',
                            target='targ_idx',
                            edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])
        fig, ax = plt.subplots(figsize=(40, 20))
        pos=nx.layout.spiral_layout(G)
        nx.draw_networkx_nodes(G,
                            pos=pos,
                            cmap=plt.cm.RdYlBu,
                            node_color=list(nodes.NES),
                            node_size=list(nodes.Hits_ratio *1000))
        nx.draw_networkx_labels(G,
                                pos=pos,
                                labels=nodes.Term.to_dict())
        edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
        nx.draw_networkx_edges(G,
                            pos=pos,
                            width=list(map(lambda x: x*10, edge_weight)),
                            edge_color='#CDDBD4')
        plt.savefig(os.path.join(sample_output, 'GSEA', f'network_plot_{name}.png'))



    def ora(self, adata, gmt, sample_output):
        """
        Performs ORA on significant genes using enricher of gseapy.
        Selections are based on up-regulated and down-regulated genes and stored in Pandas DataFrames
        Enricher is applied to these selections.
        Plots are stored in corresponding output figures.

        Args:
            adata (AnnData): The AnnData object of the sample
            gmt (dict): Dictionary containing the 'c5.go.bp' gene set collection
            sample_output (String): Path to sample-specific sub-directory to store results in
        """

        df = sc.get.rank_genes_groups_df(adata, group='BL_C', key='rank_genes_groups', pval_cutoff=None, log2fc_min=None, log2fc_max=None, gene_symbols=None)
        sign_genes = df[df.pvals_adj<0.05]
        genes_up_reg = sign_genes[sign_genes.logfoldchanges>0]
        genes_down_reg = sign_genes[sign_genes.logfoldchanges<0]
        enr_up = gseapy.enrichr(genes_up_reg.names, gene_sets=gmt, outdir=None)
        enr_pos = enr_up.res2d.sort_values(['Adjusted P-value','Combined Score'], ascending=[True, True])
        up = gseapy.dotplot(enr_pos, figsize=(3,12), title="Up", top_term=30, cmap = plt.cm.autumn_r, ofname=os.path.join(sample_output, 'ORA', 'enr_up_dotplot.png'))
        enr_dw = gseapy.enrichr(genes_down_reg.names, gene_sets=gmt, outdir=None)
        enr_neg = enr_dw.res2d.sort_values(['Adjusted P-value','Combined Score'], ascending=[True, True])
        down = gseapy.dotplot(enr_neg, figsize=(3,12), title="Down", top_term=30, cmap = plt.cm.winter_r, size=5, ofname=os.path.join(sample_output, 'ORA', 'enr_down_dotplot.png'))
        enr_pos['UP_OR_DOWN'] = 'UP'
        enr_neg['UP_OR_DOWN'] = 'DOWN'
        enr_res = pd.concat([enr_pos, enr_neg])
        enr_res = enr_res.sort_values(['Adjusted P-value', 'Combined Score'], ascending=[True, False])        
        sci = SciPalette()
        NbDr = sci.create_colormap()
        gseapy.dotplot(enr_res,figsize=(3,5), x='UP_OR_DOWN', x_order = ["UP","DOWN"], top_term=30,
                        title="GO_BP", cmap = NbDr.reversed(), size=3, show_ring=True,
                          ofname=os.path.join(sample_output, 'ORA', 'up_and_down_pathways_dotplot.png'))
        
