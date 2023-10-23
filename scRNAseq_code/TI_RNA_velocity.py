import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import cellrank as cr
import scvelo as scv
from matplotlib import pyplot as plt
import os
import rpy2.robjects
from rpy2.robjects import r
import rpy2.ipython.html
import rpy2
import anndata2ri
anndata2ri.activate()
r('library(Seurat)')
r('library(scuttle)')
r('library(scran)')
r('library(reticulate)')
r('use_virtualenv("C:/Users/julia/.virtualenvs/project-VjJne3mB")')
r('library(SeuratDisk)')
r('library(SeuratData)')
r('library(patchwork)')
r('library(sceasy)')
r('library(monocle3)')
r('library(SeuratWrappers)')

sc.settings.verbosity = 3
sc.logging.print_header()

class TI_RNA_Velocity:

    def __init__(self, sample_name, annotated_adata_loc, tirv_output_loc, path_to_loomp_C, mono_loomp_files):
        self.sample_name = sample_name
        self.sample_output = tirv_output_loc
        #self.makedirs(self.sample_output)
        self.adata = self.open_and_process_data(sample_name=self.sample_name,
                                                annotated_adata_loc=annotated_adata_loc,
                                                path_to_loomp_C=path_to_loomp_C,
                                                mono_loomp_files=mono_loomp_files)
        self.gene_set = self.define_marker_gene_set(self.sample_name)
        #self.perform_RNA_velocity(self.adata, self.gene_set)
        root_cell =  self.find_root_cell(adata=self.adata, gene_set=self.gene_set)
        #self.transform_an_to_seur(annotated_adata_loc=annotated_adata_loc, sample_name=self.sample_name)
        self.perform_TI(root_cell=root_cell, gene_set=self.gene_set)


    def makedirs(self, tirv_output_loc):
        os.makedirs(tirv_output_loc)
        os.chdir(tirv_output_loc)
        os.mkdir('Figures_output')
        os.mkdir('AnnData_storage')
        os.mkdir('RDS_storage')

    
    def define_marker_gene_set(self, sample_name):
        if sample_name == 'BL_A_BL_C':
            gene_set = ['VIM', 'S100B', 'SOX9']
        else:
            gene_set = ['RBFOX3', 'NEUROG2', 'MAP2']
        return gene_set


    def open_and_process_data(self, sample_name, annotated_adata_loc, path_to_loomp_C, mono_loomp_files):
        if sample_name == 'BL_A_BL_C':
            mono_loomp = mono_loomp_files[0]
            replace = 'BL_A'
        else:
            mono_loomp = mono_loomp_files[1]
            replace = 'BL_N'
        loomp_mono = sc.read_loom(mono_loomp)
        loomp_BLC = sc.read_loom(path_to_loomp_C)
        adata = sc.read_h5ad(os.path.join(annotated_adata_loc, f'{sample_name}_annotated.h5ad'))
        adata.raw = adata
        loomp_mono.var_names_make_unique()
        loomp_BLC.var_names_make_unique()
        loomp_mono.obs_names = loomp_mono.obs_names.drop_duplicates()
        loomp_BLC.obs_names = loomp_BLC.obs_names.drop_duplicates()
        loomp_combined = ad.concat([loomp_mono, loomp_BLC])
        split_barcodes= [barcode.split('_') for barcode in adata.obs_names]
        new_barcodes = []
        for set in split_barcodes:
            if set[1] == '1':
                new_barcode = f'{set[0]}-{replace}'
            else:
                new_barcode = f'{set[0]}-BL_C'
            new_barcodes.append(new_barcode)
        adata.obs_names = new_barcodes
        split_barcodes = [barcode.split(':') for barcode in loomp_combined.obs_names]
        new_barcodes = []
        for set in split_barcodes:
            set[1] = set[1][:-1]
            new_barcode = f'{set[1]}-1-{set[0]}'
            new_barcodes.append(new_barcode)
        loomp_combined.obs_names = new_barcodes
        common_cells = adata.obs_names.intersection(loomp_combined.obs_names)
        common_genes = adata.var_names.intersection(loomp_combined.var_names)
        new_loomp_combined = loomp_combined[common_cells, common_genes]
        adata = adata[common_cells, common_genes]
        adata.layers['spliced'] = new_loomp_combined.layers['spliced']
        adata.layers['unspliced'] = new_loomp_combined.layers['unspliced']
        del loomp_mono
        del loomp_BLC
        return adata


    def perform_RNA_velocity(self, adata, gene_set):
        adata.obs['annotated_clusters'] = adata.obs['kriegstein.seurat.custom.clusters.mean']
        scv.pl.proportions(adata, groupby='annotated_clusters', show=False)
        plt.savefig(os.path.join(self.sample_output, 'Figures_output', 'splicing_proportions.png'))
        adata.raw = adata
        sc.pp.highly_variable_genes(adata)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
        scv.pp.moments(adata, n_pcs=20, n_neighbors=20, use_highly_variable=True)
        scv.tl.recover_dynamics(adata, n_jobs=2)
        scv.tl.velocity(adata, mode="dynamical", filter_genes=True, use_highly_variable=True)
        vk = cr.kernels.VelocityKernel(adata)
        vk.compute_transition_matrix()
        ck = cr.kernels.ConnectivityKernel(adata)
        ck.compute_transition_matrix()
        vk.plot_projection(color='annotated_clusters', show=False, save=os.path.join(self.sample_output, 'Figures_output', 'RNA_splicing_projection.png'))

    
    def transform_an_to_seur(self, annotated_adata_loc, sample_name):
        inpath=os.path.join(annotated_adata_loc, f'{sample_name}_annotated.h5ad')
        self.outpath=os.path.join(self.sample_output, 'RDS_storage', 'seuratobj.rds')
        rpy2.robjects.globalenv['inpath'] = inpath # inherrit root cell to R 
        rpy2.robjects.globalenv['outpath'] = self.outpath # inherrit root cell to R 
        r(
            '''
            antoseur_Converter <- function(inpath, outpath) {
                sceasy::convertFormat(
                inpath,
                from = "anndata",
                to = "seurat",
                outFile = outpath
                )
            }
            '''
        )
        r(
            '''
            initialize_and_run <- function(inpath, outpath) {
            # Attempt to run the function once
            tryCatch({
                antoseur_Converter(inpath, outpath)
            }, error = function(e) {
                cat("First attempt failed. Retrying...\n")
            })

            # Retry the function with the same parameters
            antoseur_Converter(inpath, outpath)
            }

            initialize_and_run(inpath, outpath)
            '''
        )
    
    def perform_TI(self, root_cell, gene_set):
        rpy2.robjects.globalenv['root_cell'] = root_cell
        rpy2.robjects.globalenv['path_to_rds'] = os.path.join(self.sample_output, 'RDS_storage', 'seuratobj.rds')
        rpy2.robjects.globalenv['figure_output'] = os.path.join(self.sample_output, 'Figures_output')
        rpy2.robjects.globalenv['gene_set'] = gene_set
        r(
            '''
            adata <- readRDS(file=path_to_rds)
            cds <- SeuratWrappers::as.cell_data_set(adata)
            cds <- cluster_cells(cds)
            p1 <- plot_cells(cds, show_trajectory_graph=FALSE, color_cells_by='partition')
            ggplot2::ggsave(file.path(figure_output, paste0("sample_partition.png")), plot = p1)
            cds <- learn_graph(cds, use_partition=FALSE)
            cds <- order_cells(cds, root_cells=root_cell)
            p2 <- plot_cells(cds, color_cells_by="pseudotime", label_branch_points=FALSE, label_leaves=FALSE)
            ggplot2::ggsave(file.path(figure_output, paste0("cds_pseudotime.png")), plot = p2)
            rowData(cds)$gene_name <- rownames(cds)
            rowData(cds)$gene_short_name <- rowData(cds)$gene_name
            p3 <- plot_cells(cds, genes=c(gene_set), label_cell_groups=FALSE, show_trajectory_graph=FALSE, min_expr=3)
            ggplot2::ggsave(file.path(figure_output, paste0("cds_genes_pseudotime.png")), plot = p3)          
            '''
        )

    
    def find_root_cell(self, adata, gene_set):
        adata_raw = adata.raw.to_adata()
        adata_raw = adata_raw[:, gene_set]
        df = adata_raw.to_df()
        df['sum'] = df[gene_set[0]]+df[gene_set[1]]+df[gene_set[2]]
        sorted_data = df.sort_values(by='sum')
        top_100_lowest_sums = sorted_data.head(100)
        root_cell_barcode = list(top_100_lowest_sums.index)[0]
        if root_cell_barcode.endswith("-BL_C"):
            root_cell_barcode=root_cell_barcode[:-5]+"_2"
        else:
            root_cell_barcode=root_cell_barcode[:-5]+"_1"
        return root_cell_barcode