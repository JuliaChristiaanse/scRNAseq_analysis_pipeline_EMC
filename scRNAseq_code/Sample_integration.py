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
    def __init__(self, co_sample, mono_sample, full_name, output_dir, sample_output, markerpath):
        self.co_sample = co_sample
        self.mono_sample = mono_sample
        self.output_dir = output_dir
        self.full_name = full_name
        self.sample_output = sample_output
        self.markerpath = markerpath
        self.run()


    # Directories are created
    def makedirs(self):
        os.makedirs(os.path.join(self.sample_output))
        os.chdir(os.path.join(self.sample_output))
        os.makedirs('DEA')
        os.makedirs('UMAPS')
        os.makedirs('AnnData_storage')
        os.makedirs('h5ad_RDS_storage')

    
    def anndata_to_rds(self):
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
            # Attempt to run the function once
            tryCatch({
                antoseur_Converter(inpath, outpath)
            }, error = function(e) {
                cat("First attempt failed. Retrying...\n")
            })

            # Retry the function with the same parameters
            antoseur_Converter(inpath, outpath)
            }

            initialize_and_run(mono_inpath, mono_outpath)
            initialize_and_run(co_inpath, co_outpath)
            '''
        )


    def run_cca(self):
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
            sample_files <- list(mono_outpath, co_outpath)
            samples.list <- lapply(X=sample_files, FUN=function(x) {
            readRDS(file=x)
            })
            samples.list <- lapply(X=samples.list, FUN=function(x) {
            x <- NormalizeData(x)
            x <- FindVariableFeatures(x, selection.method="vst", nfeatures=3000)
            })
            features <- SelectIntegrationFeatures(object.list = samples.list)
            anchors <- FindIntegrationAnchors(object.list = samples.list, anchor.features = features)
            combined <- IntegrateData(anchorset = anchors)
            DefaultAssay(combined) <- "integrated"
            combined <- ScaleData(combined, verbose = FALSE)
            combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
            combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
            combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
            combined <- FindClusters(combined, resolution = 0.5, algorithm=1)
            saveRDS(combined, file = file.path(combined_rds_path))
            sceasy::convertFormat(
            combined,
            from = "seurat",
            to = "anndata",
            outFile = combined_h5ad_path
            )
            '''
        )

    
    def perform_cell_selection(self):
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
        combined_subset.write(path_to_h5ad_after_cc)
        rpy2.robjects.globalenv['path_to_h5ad_after_cc'] = path_to_h5ad_after_cc
        self.path_to_rds_after_cc = os.path.join(self.sample_output, 'h5ad_RDS_storage', 'combined_after_cc.rds')
        rpy2.robjects.globalenv['path_to_rds_after_cc'] = self.path_to_rds_after_cc
        r(
            '''
            sceasy::convertFormat(
            path_to_h5ad_after_cc,
            from = "anndata",
            to = "seurat",
            outFile = path_to_rds_after_cc
            )
            '''
        )
    
    def re_run_cca(self):
        print(f'running {self.mono_sample}')
        self.final_integrated_h5ad = os.path.join(self.sample_output, 'AnnData_storage', f'{self.mono_sample}_{self.co_sample}.h5ad')
        rpy2.robjects.globalenv['path_to_rds_after_cc'] = self.path_to_rds_after_cc
        rpy2.robjects.globalenv['final_integrated_h5ad'] = self.final_integrated_h5ad
        r(
            '''
            combined_subset = readRDS(file=path_to_rds_after_cc)
            combined.list <- SplitObject(combined_subset, split.by = "Sample")
            combined.list <- lapply(X = combined.list, FUN = function(x) {
            x <- NormalizeData(x)
            x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
            })
            features <- SelectIntegrationFeatures(object.list = combined.list)
            anchors <- FindIntegrationAnchors(object.list = combined.list, anchor.features = features)
            re_combined <- IntegrateData(anchorset = anchors)
            DefaultAssay(re_combined) <- "integrated"
            re_combined <- ScaleData(re_combined, verbose = FALSE)
            re_combined <- RunPCA(re_combined, npcs = 30, verbose = FALSE)
            re_combined <- RunUMAP(re_combined, reduction = "pca", dims = 1:30)
            re_combined <- FindNeighbors(re_combined, reduction = "pca", dims = 1:30)
            re_combined <- FindClusters(re_combined, resolution = 0.5, algorithm = 1)
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
        self.makedirs()
        self.anndata_to_rds()
        self.run_cca()
        self.perform_cell_selection()
        self.re_run_cca()

