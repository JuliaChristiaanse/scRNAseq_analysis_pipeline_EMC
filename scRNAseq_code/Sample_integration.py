"""
Author:         	Julia Christiaanse
Date:               29-11-2023
Script:             Class Sample_integration
Python version:     3.10.9
Imports & settings:
"""
import scanpy as sc
import os
from matplotlib import pyplot as plt
import rpy2.robjects
from rpy2.robjects import r
import rpy2.ipython.html
import rpy2
import anndata2ri
anndata2ri.activate()
from Differential_expression_analysis import Differential_Expression_Analysis as dea
r('library(Seurat)')
r('library(SingleR)')
r('library(scuttle)')
r('library(scran)')
r('library(reticulate)')
r('use_virtualenv("C:/Users/julia/.virtualenvs/project-VjJne3mB")')
r('library(SeuratDisk)')
r('library(SeuratData)')
r('library(patchwork)')
r('library(sceasy)')
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=120, dpi_save=150, facecolor='white', color_map='Blues')


class Sample_integration:
    """
    Class Sample_integration performs steps 4 & 5 of the pipeline. 
    These steps include integration of two samples, cell selection based on marker genes, and re-integration of the remaining cells.
    """    


    def __init__(self, co_sample, mono_sample, full_name, output_dir, sample_output, markerpath):
        """
        __init__ of the class Sample_integration.

        Args:
            co_sample (String): Name of the co-culture sample
            mono_sample (String): Name of the mono-culture sample
            full_name (String): Combined name of co-culture and mono-culture sample
            output_dir (string): Path to the master output directory
            sample_output (String): Path to the sample-specific output sub-folder within output_dir
            markerpath (String): Path to markergenes.txt
        """        
        self.co_sample = co_sample
        self.mono_sample = mono_sample
        self.output_dir = output_dir
        self.full_name = full_name
        self.sample_output = sample_output
        self.markerpath = markerpath
        self.run()


    def makedirs(self):
        """
        Generate necessary sub-folders within the master folder, sample_output.
        """        
        os.makedirs(os.path.join(self.sample_output))
        os.chdir(os.path.join(self.sample_output))
        os.makedirs('DEA')
        os.makedirs('UMAPS')
        os.makedirs('AnnData_storage')
        os.makedirs('h5ad_RDS_storage')

    
    def anndata_to_rds(self):
        """
        Transforms the AnnData object to a Seurat (RDS) object to ensure data compatibility with the R ecosystem.
        Full paths to both the co-culture as the current mono-culture samples are provided.
        The 'Sample' column is added for sample recognition, and raw counts are accessed.
        R variables of the full paths to both samples are defined.
        Using Rpy2, R code is called and the AnnData object is transformed to a Seurat object.
        The object is stored on disk, in the 'head_RDS_storage' sub-folder.
        """        
        self.co_adata = sc.read_h5ad(os.path.join(self.output_dir, self.co_sample,'AnnData_storage', self.co_sample+'.h5ad'))
        self.mono_adata = sc.read_h5ad(os.path.join(self.output_dir, self.mono_sample,'AnnData_storage', self.mono_sample+'.h5ad'))
        
        self.co_adata.obs['Sample'] = self.co_sample
        self.mono_adata.obs['Sample'] = self.mono_sample

        self.co_adata = self.co_adata.raw.to_adata()
        self.mono_adata = self.mono_adata.raw.to_adata()

        self.co_adata.write(os.path.join(self.sample_output, 'h5ad_RDS_storage', f'{self.co_sample}_raw.h5ad'))
        self.mono_adata.write(os.path.join(self.sample_output, 'h5ad_RDS_storage', f'{self.mono_sample}_raw.h5ad'))
        
        mono_inpath = os.path.join(self.sample_output, 'h5ad_RDS_storage', f'{self.mono_sample}_raw.h5ad')
        mono_outpath = os.path.join(self.sample_output, 'h5ad_RDS_storage', f'{self.mono_sample}_raw.rds')
        
        co_inpath = os.path.join(self.sample_output, 'h5ad_RDS_storage', f'{self.co_sample}_raw.h5ad')
        co_outpath = os.path.join(self.sample_output, 'h5ad_RDS_storage', f'{self.co_sample}_raw.rds')
        
        rpy2.robjects.globalenv['mono_inpath'] = mono_inpath
        rpy2.robjects.globalenv['mono_outpath'] = mono_outpath
        rpy2.robjects.globalenv['co_inpath'] = co_inpath
        rpy2.robjects.globalenv['co_outpath'] = co_outpath

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
            # Attempt to run the function once, use tryCatch in case of error.
            tryCatch({
                antoseur_Converter(inpath, outpath)
            }, error = function(e) {
                cat("First attempt failed. Retrying...\n")
            })
            antoseur_Converter(inpath, outpath)
            }

            initialize_and_run(mono_inpath, mono_outpath)
            initialize_and_run(co_inpath, co_outpath)
            '''
        )


    def run_cca(self):
        """
        Runs first round of Seurat Canonical Correlation Analysis (CCA) on co-culture and current mono-culture sample.
        R variables 'mono_outpath', 'co_outpath' are created, containing both paths of the previously transformed Seurat object of both samples.
        New paths to store the first integrated embeddings as RDS and H5ad files are then created.
        R variables 'combined_rds_path', 'combined_h5ad_path' are made containing both paths respectively.
        R code is ran with Rpy2. See further documentation below.
        """        
        mono_outpath = os.path.join(self.sample_output, 'h5ad_RDS_storage', f'{self.mono_sample}_raw.rds')
        co_outpath = os.path.join(self.sample_output, 'h5ad_RDS_storage', f'{self.co_sample}_raw.rds')
        rpy2.robjects.globalenv['mono_outpath'] = mono_outpath
        rpy2.robjects.globalenv['co_outpath'] = co_outpath
        combined_rds_path = os.path.join(self.sample_output, 'h5ad_RDS_storage', 'combined_cca1.rds')
        combined_h5ad_path = os.path.join(self.sample_output, 'h5ad_RDS_storage', 'combined_cca1.h5ad')
        rpy2.robjects.globalenv['combined_rds_path'] = combined_rds_path
        rpy2.robjects.globalenv['combined_h5ad_path'] = combined_h5ad_path
        r(
            '''
            # Read in both samples as Seurat objects
            sample_files <- list(mono_outpath, co_outpath)
            samples.list <- lapply(X=sample_files, FUN=function(x) {
            readRDS(file=x)
            })
            # Preprocess both samples
            samples.list <- lapply(X=samples.list, FUN=function(x) {
            x <- NormalizeData(x)
            x <- FindVariableFeatures(x, selection.method="vst", nfeatures=3000)
            })
            # Find similar genes between samples
            features <- SelectIntegrationFeatures(object.list = samples.list)

            # Integrate both samples into 'combined'
            anchors <- FindIntegrationAnchors(object.list = samples.list, anchor.features = features)
            combined <- IntegrateData(anchorset = anchors)
            DefaultAssay(combined) <- "integrated"

            # Run scale, PCA, UMAP, KNN, and Louvain clusters on 'combined'
            combined <- ScaleData(combined, verbose = FALSE)
            combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
            combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
            combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
            combined <- FindClusters(combined, resolution = 0.5, algorithm=1)

            # Write Seurat object to 'h5ad_RDS_storage'
            saveRDS(combined, file = file.path(combined_rds_path))

            # Convert Seurat object to AnnData object and write to 'h5ad_RDS_storage' 
            sceasy::convertFormat(
            combined,
            from = "seurat",
            to = "anndata",
            outFile = combined_h5ad_path
            )
            '''
        )

    
    def perform_cell_selection(self):
        """
        Performs cell selection after the first round of CCA has been completed.
        This is done to select only clusters containing cell types from the co-culture corresponding to the mono-culture.
        Starts off with loading the h5ad file into an AnnData object, called adata.
        Performs a differential expression analysis (DEA) on adata.
        Selects marker panels, which correspond to the current mono-culture.
        Performs cell selection based on expression levels of the genes of that marker panel.
        Clusters expressing at least 20% of all three genes of the marker panel are kept.
        A new AnnData object of these clusters is then created, transformed to a Seurat object and stored to 'h5ad_RDS_storage'.
        """        
        adata = sc.read_h5ad(os.path.join(self.sample_output, 'h5ad_RDS_storage', 'combined_cca1.h5ad'))
        self.plot_umaps(adata, '_after_cc1')
        adata_DE = adata.copy()
        sc.tl.rank_genes_groups(adata_DE, 'seurat_clusters', method='wilcoxon', corr_method='bonferroni', key='wilcoxon', pts=True)
        sc.tl.filter_rank_genes_groups(adata_DE, groupby='seurat_clusters', min_in_group_fraction=0.1, min_fold_change=1)
        if self.mono_sample == 'BL_N':
            genes = ['DCX', 'MAP2', 'RBFOX3']
        else:
            genes = ['VIM', 'S100B', 'SOX9']
        clusters = range(len(adata_DE.obs['seurat_clusters'].unique()))   # find all unique clusters
        clusters = [str(x) for x in clusters]
        combined_dp = sc.pl.rank_genes_groups_dotplot(adata_DE, clusters, color_map='Blues',
                                                    var_names=genes, return_fig=True, show=False)
        combined_dp.savefig(os.path.join(self.sample_output, 'DEA', f'dotplot_{self.full_name}.png'))
        combined_dp_df = combined_dp.dot_size_df
        pos = combined_dp_df[(combined_dp_df[genes[0]] > .2) & (combined_dp_df[genes[1]] > .2) & (combined_dp_df[genes[2]] > .2)]
        good_clusters = list(pos.index)
        combined_subset = adata[adata.obs['seurat_clusters'].isin(good_clusters)]
        path_to_h5ad_after_cc = os.path.join(self.sample_output, 'h5ad_RDS_storage', 'combined_after_cc.h5ad')
        # Write adata to 'h5ad_RDS_storage'.
        combined_subset.write(path_to_h5ad_after_cc)
        # Create R variable 'path_to_rds_after_cc', containing the full path to store the new Seurat object. 
        rpy2.robjects.globalenv['path_to_h5ad_after_cc'] = path_to_h5ad_after_cc
        self.path_to_rds_after_cc = os.path.join(self.sample_output, 'h5ad_RDS_storage', 'combined_after_cc.rds')
        rpy2.robjects.globalenv['path_to_rds_after_cc'] = self.path_to_rds_after_cc
        r(
            '''
            # Convert AnnData object to Seurat object, and write results to 'h5ad_RDS_storage'.
            sceasy::convertFormat(
            path_to_h5ad_after_cc,
            from = "anndata",
            to = "seurat",
            outFile = path_to_rds_after_cc
            )
            '''
        )
    
    def re_run_cca(self):
        """
        Runs a second round of Seurat's CCA on co-culture and current mono-culture sample.
        R variables 'path_to_rds_after_cc', 'final_integrated_h5ad' are created, containing the path to the RDS object after cell selection,
          and the final path to store the AnnData object after re-integration with CCA.
        R code is ran with Rpy2. See further documentation below.
        Final integrated AnnData object is then loaded in and self.plot_umaps() is called.
        """        
        print(f'running {self.mono_sample}')
        self.final_integrated_h5ad = os.path.join(self.sample_output, 'AnnData_storage', f'{self.mono_sample}_{self.co_sample}.h5ad')
        rpy2.robjects.globalenv['path_to_rds_after_cc'] = self.path_to_rds_after_cc
        rpy2.robjects.globalenv['final_integrated_h5ad'] = self.final_integrated_h5ad
        r(
            '''
            # Load the Seurat object created after cell selection.
            combined_subset = readRDS(file=path_to_rds_after_cc)

            # First split the combined object based on 'Sample' and store Seurat objects in list.
            combined.list <- SplitObject(combined_subset, split.by = "Sample")

            # Re-perform preprocessing on Seurat objects in list.
            combined.list <- lapply(X = combined.list, FUN = function(x) {
            x <- NormalizeData(x)
            x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
            })

            # Find common genes
            features <- SelectIntegrationFeatures(object.list = combined.list)

            # Re-integrate both samples to create 're_combined'
            anchors <- FindIntegrationAnchors(object.list = combined.list, anchor.features = features)
            re_combined <- IntegrateData(anchorset = anchors)
            DefaultAssay(re_combined) <- "integrated"

            # Run scale, PCA, UMAP, KNN, and Louvain clusters on 're_combined'
            re_combined <- ScaleData(re_combined, verbose = FALSE)
            re_combined <- RunPCA(re_combined, npcs = 30, verbose = FALSE)
            re_combined <- RunUMAP(re_combined, reduction = "pca", dims = 1:30)
            re_combined <- FindNeighbors(re_combined, reduction = "pca", dims = 1:30)
            re_combined <- FindClusters(re_combined, resolution = 0.5, algorithm = 1)

            # Convert Seurat object to AnnData object and write to 'h5ad_RDS_storage'
            sceasy::convertFormat(
            re_combined,
            from = "seurat",
            to = "anndata",
            outFile = final_integrated_h5ad
            )
            '''
        )
        adata_subset = sc.read_h5ad(self.final_integrated_h5ad)
        self.plot_umaps(adata_subset, title='_after_cca2')


    def plot_umaps(self, adata, title):
        """
        Creates four UMAPs of the integrated embedding and plots them in one figure.
        The first one plots the Seurat clusters, the second the different samples,
        the third the mono-culture sample, and the fourth the co-culture sample.

        Args:
            adata (AnnData): The AnnData object that is the integrated sample
            title (String): Title of the .png file
        """        
        title = self.co_sample+'_'+self.mono_sample+title
        fig, axs = plt.subplots(2, 2, figsize=(10,8),constrained_layout=True)
        cocult = adata[adata.obs['Sample'] == self.co_sample]
        sample = adata[adata.obs['Sample'] == self.mono_sample]
        sc.pl.umap(adata, color='seurat_clusters', ax=axs[0,0], title='Both samples seurat_clusters', show=False)
        sc.pl.umap(adata, color='Sample', ax=axs[0,1], title='Both samples colored', show=False)
        sc.pl.umap(sample, color='seurat_clusters', ax=axs[1,0], title=f'{self.mono_sample}', show=False)
        sc.pl.umap(cocult, color='seurat_clusters', ax=axs[1,1], title=f'{self.co_sample}', show=False)
        plt.savefig(os.path.join(self.sample_output, 'UMAPS', title+'.png'))


    def run(self):
        """
        Runs all class functions in chronological order. 
        Created for aesthetic purposes, to avoid running every function in the __init__.
        """        
        self.makedirs()
        self.anndata_to_rds()
        self.run_cca()
        self.perform_cell_selection()
        self.re_run_cca()

