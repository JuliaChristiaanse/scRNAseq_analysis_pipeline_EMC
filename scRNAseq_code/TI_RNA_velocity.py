"""
Author:         	Julia Christiaanse
Date:               29-11-2023
Script:             Class TI_RNA_Velocity
Python version:     3.10.9
Imports & settings:
"""
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
    """
    Class TI_RNA_Velocity performs step 7A of the pipeline.
    Performs trajectory inference (TI) and RNA velocity analysis.
    """


    def __init__(self, sample_name, annotated_adata_loc, tirv_output_loc, path_to_loomp_C, mono_loomp_files):
        """
        __init__ of class TI_RNA_Velocity. 

        Args:
            sample_name (String): Name of the sample
            annotated_adata_loc (String): Path to annotated AnnData object
            tirv_output_loc (String): Path to subfolder to store output files
            path_to_loomp_C (String): Path to co-culture .loom file
            mono_loomp_files (List): List of paths to mono-culture .loom files
        """        
        self.sample_name = sample_name
        self.sample_output = tirv_output_loc
        self.makedirs(self.sample_output)
        self.adata = self.open_and_process_data(sample_name=self.sample_name,
                                                annotated_adata_loc=annotated_adata_loc,
                                                path_to_loomp_C=path_to_loomp_C,
                                                mono_loomp_files=mono_loomp_files)
        self.gene_set = self.define_marker_gene_set(self.sample_name)
        self.perform_RNA_velocity(self.adata)
        root_cell =  self.find_root_cell(adata=self.adata, gene_set=self.gene_set)
        self.transform_an_to_seur(annotated_adata_loc=annotated_adata_loc, sample_name=self.sample_name)
        self.perform_TI(root_cell=root_cell, gene_set=self.gene_set)


    def makedirs(self, tirv_output_loc):
        """
        Generate necessary sub-folders within tirv_output_loc.

        Args:
            tirv_output_loc (String): Path to subfolder to store output files
        """        
        os.makedirs(tirv_output_loc)
        os.chdir(tirv_output_loc)
        os.mkdir('Figures_output')
        os.mkdir('AnnData_storage')
        os.mkdir('RDS_storage')

    
    def define_marker_gene_set(self, sample_name):
        """
        Defines a marker gene set based on the mono-culture currently being processed.

        Args:
            sample_name (String): Name of the sample

        Returns:
            List: Gene_set, list of marker genes that correspond to the name of the sample
        """        
        if sample_name == 'BL_A_BL_C':
            gene_set = ['VIM', 'S100B', 'SOX9']
        else:
            gene_set = ['RBFOX3', 'NEUROG2', 'MAP2']
        return gene_set


    def open_and_process_data(self, sample_name, annotated_adata_loc, path_to_loomp_C, mono_loomp_files):
        """
        Loads paths to .loom files and annotated AnnData object into AnnData objects.
        Processess and filter the .loom files so they become compatible with the annotated AnnData object.
        Adds the splicing information to the cells of the annotated AnnData object for downstream analysis.

        Args:
            sample_name (String): Name of the sample
            annotated_adata_loc (String): Path to annotated AnnData object
            path_to_loomp_C (String): Path to co-culture .loom file
            mono_loomp_files (List): List of paths to mono-culture .loom files

        Returns:
            AnnData: The AnnData object
        """
        # Select mono .loom file corresponding to the sample name.
        if sample_name == 'BL_A_BL_C':
            mono_loomp = mono_loomp_files[0]
            replace = 'BL_A'
        else:
            mono_loomp = mono_loomp_files[1]
            replace = 'BL_N'
        # Load files into AnnData objects.
        loomp_mono = sc.read_loom(mono_loomp)
        loomp_BLC = sc.read_loom(path_to_loomp_C)
        adata = sc.read_h5ad(os.path.join(annotated_adata_loc, f'{sample_name}_annotated.h5ad'))
        adata.raw = adata
        # Preprocess .loom files AnnData objects.
        loomp_mono.var_names_make_unique()
        loomp_BLC.var_names_make_unique()
        loomp_mono.obs_names = loomp_mono.obs_names.drop_duplicates()
        loomp_BLC.obs_names = loomp_BLC.obs_names.drop_duplicates()
        # Combine both .loom file AnnData objects into one object.
        loomp_combined = ad.concat([loomp_mono, loomp_BLC])
        split_barcodes= [barcode.split('_') for barcode in adata.obs_names]
        new_barcodes = []
        # Replace numbers at the end of the .loom AnnData barcodes with sample names.
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
        # Find common cells and genes between the annotated AnnData object and the combined AnnData object of the .loom files.
        common_cells = adata.obs_names.intersection(loomp_combined.obs_names)
        common_genes = adata.var_names.intersection(loomp_combined.var_names)
        # Subset the combined AnnData object to the common cells and genes.
        new_loomp_combined = loomp_combined[common_cells, common_genes]
        # Subset the annotated AnnData object to the common cells and genes.
        adata = adata[common_cells, common_genes]
        # Inherrit splicing information to annotated AnnData object from combined AnnData object.
        adata.layers['spliced'] = new_loomp_combined.layers['spliced']
        adata.layers['unspliced'] = new_loomp_combined.layers['unspliced']
        # delete the old .loom file AnnData objects.
        del loomp_mono
        del loomp_BLC
        return adata


    def perform_RNA_velocity(self, adata):
        """
        Performs RNA velocity analysis.
        Preprocesses the AnnData object and computes a velocity kernel and transition matrix.
        Plots the final projection and stores the figure.

        Args:
            adata (AnnData): The AnnData object
        """
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
        vk.plot_projection(color='annotated_clusters', show=False, legend_loc='right margin', save=os.path.join(self.sample_output, 'Figures_output', 'RNA_splicing_projection.png'))

    
    def transform_an_to_seur(self, annotated_adata_loc, sample_name):
        """
        Transforms annotated AnnData object to a Seurat object to prepare it for TI analysis with Monocle3.

        Args:
            annotated_adata_loc (String): Path to annotated AnnData object
            sample_name (String): Name of the sample
        """        
        inpath=os.path.join(annotated_adata_loc, f'{sample_name}_annotated.h5ad')
        self.outpath=os.path.join(self.sample_output, 'RDS_storage', 'seuratobj.rds')
        # Inherrit Python variables to R.
        rpy2.robjects.globalenv['inpath'] = inpath
        rpy2.robjects.globalenv['outpath'] = self.outpath
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
        """
        Performs TI with Monocle3.

        Args:
            root_cell (String): Barcode of the root cell
            gene_set (List):  list of strings, that are the genes of the marker gene panel
        """ 
        # Inherrit Python variables to R.   
        rpy2.robjects.globalenv['root_cell'] = root_cell
        rpy2.robjects.globalenv['path_to_rds'] = os.path.join(self.sample_output, 'RDS_storage', 'seuratobj.rds')
        rpy2.robjects.globalenv['figure_output'] = os.path.join(self.sample_output, 'Figures_output')
        rpy2.robjects.globalenv['gene_set'] = gene_set
        r(
            '''
            # Read in RDS object.
            adata <- readRDS(file=path_to_rds)

            # Follow Monocle3 vignette to perform TI, create & save plots.
            cds <- SeuratWrappers::as.cell_data_set(adata)
            cds <- cluster_cells(cds)
            p1 <- plot_cells(cds, show_trajectory_graph=FALSE, color_cells_by='partition')
            ggplot2::ggsave(file.path(figure_output, paste0("sample_partition.png")), plot = p1)
            cds <- learn_graph(cds, use_partition=FALSE)
            cds <- order_cells(cds, root_cells=root_cell)
            p2 <- plot_cells(cds, color_cells_by="pseudotime", label_branch_points=FALSE,label_leaves=FALSE)
            ggplot2::ggsave(file.path(figure_output, paste0("cds_pseudotime.png")), plot = p2)
            rowData(cds)$gene_name <- rownames(cds)
            rowData(cds)$gene_short_name <- rowData(cds)$gene_name
            p3 <- plot_cells(cds, genes=c(gene_set), label_cell_groups=FALSE, show_trajectory_graph=FALSE, min_expr=3)
            ggplot2::ggsave(file.path(figure_output, paste0("cds_genes_pseudotime.png")), plot = p3)          
            '''
        )

    
    def find_root_cell(self, adata, gene_set):
        """
        Computes a root cell based on the lowest expression value of genes in the gene_set.

        Args:
            adata (AnnData): The AnnData object
            gene_set (List): list of strings, that are the genes of the marker gene panel

        Returns:
            String: root_cell_barcode, barcode of the root cell
        """        
        adata_raw = adata.raw.to_adata()
        adata_raw = adata_raw[:, gene_set]
        df = adata_raw.to_df()
        # Compute 'sum' column, that is the sum of the expression value of all three genes of the marker panel.
        df['sum'] = df[gene_set[0]]+df[gene_set[1]]+df[gene_set[2]]
        # Sort the column and access the top 100 genes with the lowest expression value.
        sorted_data = df.sort_values(by='sum')
        top_100_lowest_sums = sorted_data.head(100)
        # Pick the top cell as root cell.
        root_cell_barcode = list(top_100_lowest_sums.index)[0]
        # Transform the end of the barcode to either '_2', or '_1' depending on the sample.
        # This is to ensure compatibility with barcodes from Seurat object, which end on numbers instead of sample names.
        if root_cell_barcode.endswith("-BL_C"):
            root_cell_barcode=root_cell_barcode[:-5]+"_2"
        else:
            root_cell_barcode=root_cell_barcode[:-5]+"_1"
        return root_cell_barcode