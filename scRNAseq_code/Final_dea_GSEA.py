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

    def __init__(self, sample_name, annotated_adata_loc, final_dea_output_loc, output_dir, markerpath):
        # EDIT: DO we need raw counts for dea? Or can they be normalized? In this case they are normalized.
        self.sample_name = sample_name
        self.sample_output = final_dea_output_loc
        self.output_dir = output_dir
        self.markerpath = markerpath
        self.makedirs(self.sample_output)
        self.adata = sc.read(os.path.join(annotated_adata_loc, f'{self.sample_name}_annotated.h5ad'))
        msig = Msigdb()
        gmt = msig.get_gmt(category='c5.go.bp', dbver='2023.2.Hs')
        adata, marker_dict = self.dea(sample_name=self.sample_name, adata=self.adata,
                  markerpath=self.markerpath, sample_output=self.sample_output)
        self.volcano_plot(sample_name=self.sample_name,
                          marker_dict=marker_dict, sample_output=self.sample_output)
        self.gsea(adata=adata, gmt=gmt, sample_output=self.sample_output)
        self.ora(adata=adata, gmt=gmt, sample_output=self.sample_output)

    
    def makedirs(self, sample_output):
        os.makedirs(sample_output)
        os.chdir(sample_output)
        os.makedirs('dea_output')
        os.makedirs('AnnData_storage')
        os.makedirs('GSEA')
        os.makedirs('ORA')


    def dea(self, sample_name, adata, markerpath, sample_output):
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
        df = sc.get.rank_genes_groups_df(adata, group='BL_C', key='rank_genes_groups', pval_cutoff=0.05, log2fc_min=None, log2fc_max=None, gene_symbols=None)
        df['Rank'] = -np.log10(df.pvals_adj)*df.logfoldchanges
        df = df.sort_values('Rank', ascending=False)
        ranking = df[['names', 'Rank']]
        pre_res = gseapy.prerank(rnk = ranking, gene_sets = gmt, seed = 0)
        out_df = pre_res.res2d.sort_values('NES')
        neg_terms_to_graph = out_df.iloc[0:5].Term
        pos_terms_to_graph = out_df.tail(5).Term
        pre_res.plot(terms=pos_terms_to_graph, ofname=os.path.join(sample_output, 'GSEA', 'positive_NES_pathways.png'))
        pre_res.plot(terms=neg_terms_to_graph, ofname=os.path.join(sample_output, 'GSEA', 'negative_NES_pathways.png'))
        pre_res.res2d = pre_res.res2d.sort_values('FDR q-val')
        nodes, edges = gseapy.enrichment_map(pre_res.res2d)
        G = nx.from_pandas_edgelist(edges,
                            source='src_idx',
                            target='targ_idx',
                            edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])
        fig, ax = plt.subplots(figsize=(20, 8))
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
        plt.savefig(os.path.join(sample_output, 'GSEA', 'network_plot.png'))
        


    def ora(self, adata, gmt, sample_output):
        df = sc.get.rank_genes_groups_df(adata, group='BL_C', key='rank_genes_groups', pval_cutoff=0.05, log2fc_min=None, log2fc_max=None, gene_symbols=None)
        genes_up_reg = df[df.logfoldchanges>0]
        genes_down_reg = df[df.logfoldchanges<0]
        enr_up = gseapy.enrichr(genes_up_reg.names, gene_sets=gmt, outdir=None)
        enr_dw = gseapy.enrichr(genes_down_reg.names, gene_sets=gmt, outdir=None)
        enr_up.res2d['UP_OR_DOWN'] = 'UP'
        enr_dw.res2d['UP_OR_DOWN'] = 'DOWN'
        enr_res = pd.concat([enr_up.res2d, enr_dw.res2d])
        enr_res = enr_res.sort_values('Adjusted P-value')
        sci = SciPalette()
        NbDr = sci.create_colormap()
        gseapy.dotplot(enr_res,figsize=(3,5), x='UP_OR_DOWN', x_order = ["UP","DOWN"],
                        title="GO_BP", cmap = NbDr.reversed(), size=3, show_ring=True,
                          ofname=os.path.join(sample_output, 'ORA', 'up_and_down_pathways_dotplot.png'))
        
