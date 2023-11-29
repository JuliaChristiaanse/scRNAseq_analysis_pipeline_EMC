"""
Author:         	Julia Christiaanse
Date:               29-11-2023
Script:             Class Cell_Type_Annotation
Python version:     3.10.9
Imports & settings:
"""
import scanpy as sc
import os
import rpy2
import rpy2.robjects
import rpy2.ipython.html
from rpy2.robjects.packages import importr
from rpy2.robjects import r
import anndata
from matplotlib import pyplot as plt
r('library(Seurat)')
r('library(SingleR)')
r('library(scuttle)')
r('library(scran)')
r('library(reticulate)')
r('use_virtualenv("C:/Users/julia/.virtualenvs/project-VjJne3mB")')


class Cell_Type_Annotation:
    """
    Class Cell_Type_Annotation performs step 6 of the pipeline, which is cell type annotation.
    """

    def __init__(self, sample_name, adata_loc, output_dir, sample_output, reference_data_dir, path_to_ref_data):
        """__init__ of class Cell_Type_Annotation.

        Args:
            sample_name (String): name of the sample
            adata_loc (String): Path to integrated AnnData object
            output_dir (String): Path to master output directory
            sample_output (String): Path to the sample-specific output sub-folder within output_dir
            reference_data_dir (String): Path to all reference data
            path_to_ref_data (String): Path to folder in reference_data_dir where Kriegstein chunks are stored
        """
        self.sample_name = sample_name
        self.adata = sc.read_h5ad(os.path.join(adata_loc, f'{self.sample_name}.h5ad'))
        self.output_dir = output_dir
        self.sample_output = sample_output
        self.reference_data_dir = reference_data_dir
        self.path_to_ref_data = path_to_ref_data
        self.run()


    def makedirs(self):
        """
        Generate necessary sub-folders within the master folder, sample_output.
        """
        os.makedirs(self.sample_output)
        os.chdir(self.sample_output)
        os.makedirs('AnnData_storage')
        os.makedirs('AnnData_raw_storage')
        os.makedirs('Sample_RDS_storage')
        os.makedirs(f'{self.sample_name}_figures_output')
    

    def anndata_to_rds_transformer(self):
        """
        Paths to store the AnnData object and Seurat object are defined here.
        Transforms the AnnData object to a Seurat object, and stores it in 'Sample_RDS_storage'
        """
        self.path_to_raw_sample = os.path.join(self.sample_output, 'AnnData_raw_storage', f'{self.sample_name}_anndata_raw.h5ad')
        self.path_to_sample_rds = os.path.join(self.sample_output, 'Sample_RDS_storage', f'{self.sample_name}.rds')
        self.adata.write(self.path_to_raw_sample)
        # Inherrit Python variables to R.
        rpy2.robjects.globalenv['path_to_raw_sample'] = self.path_to_raw_sample
        rpy2.robjects.globalenv['path_to_sample_rds'] = self.path_to_sample_rds
        r(
            '''# Initialize a flag to track whether the code succeeded
        success <- FALSE
        
        # Define a function to run the code
        run_conversion <- function() {
        sceasy::convertFormat(
            path_to_raw_sample,
            from = "anndata",
            to = "seurat",
            outFile = path_to_sample_rds
            )
        }
        while (!success) {
        tryCatch({
            run_conversion()
            # If no error occurred, set the success flag to TRUE and exit the loop
            success <- TRUE
        }, error = function(err) {
            cat("An error occurred:", conditionMessage(err), "\n")
            cat("Retrying...\n")
            })
        }
        cat("Code executed successfully!\n")
        '''
        )
        


    def annotate_with_singler(self):
        """ 
        **R code obtained from previous R pipeline**
        The R code loops over every Kriegstein chunk and annotates with SingleR.
        The result is 25 annotated RDS objects of the sample.
        """
        kriegstein_annotated_output_dir = os.path.join(self.sample_output, "Annotated_RDS_storage")
        # Inherrit Python variables to R.
        rpy2.robjects.globalenv['kriegstein_annotated_output_dir'] = kriegstein_annotated_output_dir
        rpy2.robjects.globalenv['kriegstein_data_dir'] = self.path_to_ref_data
        rpy2.robjects.globalenv['sample_files'] = self.path_to_sample_rds
        rpy2.robjects.globalenv['sample_names'] = self.sample_name
        rpy2.robjects.globalenv['kriegstein_chunks_input_dir'] = self.reference_data_dir
        r(
            '''
            getGenes <- function(kriegstein_data_dir) {
            genesFile <- file.path(kriegstein_data_dir, "kriegstein_genes.csv")
            if (!file.exists(genesFile)) {
                stop(genesFile, "Does not exist")
            } else {
                genes <- utils::read.table(genesFile, sep = "\\t", col.names = "gene")
                return(genes$gene[2:length(genes$gene)])
                }
            }
            # Create directory to store the annotated sample for all chunks
            dir.create(kriegstein_annotated_output_dir)
            genes <- getGenes(kriegstein_data_dir)

            sample_names <- c(sample_names)
            sample_files <- c(sample_files)
            
            for (j in 1:length(sample_files)) {
                start_iter_time <- Sys.time()
                sample_name <- sample_names[j]
                sample_data <- readRDS(file = sample_files[j])
                SeuratObject::DefaultAssay(sample_data) <- "RNA"
                sample_genes <- rownames(sample_data)
                overlapping_genes <- intersect(sample_genes, genes)
                
                # Create single cell experiment from seurat object
                message("create Seurat sample -> SCE")
                sample_data <- Seurat::as.SingleCellExperiment(sample_data, assay = 'RNA')

                # subset to genes that overlap with Kriegstein chunks
                sample_data <- sample_data[overlapping_genes,]

                # perform transformation to have genes in the same feature space as the reference data
                sample_data <- sample_data[, colSums(counts(sample_data)) > 0]

                # logNormCounts to preprocess the sample_data
                sample_data <- scuttle::logNormCounts(sample_data)
            
                # Loop through Kriegstein chunks and annotate with singleR in this for loop
                for (i in 1:length(list.files(path = kriegstein_chunks_input_dir))) {
                print("loading kriegstein sample")
                load(file.path(kriegstein_chunks_input_dir, paste0("chunk.", i, ".RData"))) # cell_data (R object name)
                cell_data <- cell_data[overlapping_genes,]
                filename_base <- file.path(kriegstein_annotated_output_dir, paste(sample_name, ".iter.", i))  
                
                # run SingleR for each annotation and save RData
                print("Running single R")
                for (annotation in c("custom.clusterv2")) {
                    result <- SingleR::SingleR(
                    test = sample_data,
                    ref = cell_data,
                    labels = SummarizedExperiment::colData(cell_data)[, annotation],
                    clusters = SummarizedExperiment::colData(sample_data)[, "seurat_clusters"],
                    de.method = 'wilcox',
                    aggr.ref = FALSE)
                    save(result, file = paste0(filename_base, ".", annotation, ".RData"))
                    print("Save Seurat result")

                    # Remove chunk before moving on to the next
                    print("Remove Seurat result")
                    rm(result)
                    }
                }
            }
            '''
        )


    def visualize_annotated_data(self):
        """
        **R code taken from previous R pipeline**
        Highest SingleR annotation score is selected out of the 25 Kriegstein Chunks.
        Annotation is visualized and the annotated Seurat object is saved to disk.
        """
        figures_output_dir = os.path.join(self.sample_output, f"{self.sample_name}_figures_output")
        kriegstein_annotated_input_dir = os.path.join(self.sample_output, "Annotated_RDS_storage")
        # Inherrit Python variables to R.
        rpy2.robjects.globalenv['sample_files'] = os.path.join(self.sample_output, 'Sample_RDS_storage', f'{self.sample_name}.rds')
        rpy2.robjects.globalenv['sample_names'] = self.sample_name
        rpy2.robjects.globalenv['output_dir'] = figures_output_dir
        rpy2.robjects.globalenv['kriegstein_data_dir'] = self.path_to_ref_data
        rpy2.robjects.globalenv['kriegstein_annotated_input_dir'] = kriegstein_annotated_input_dir
        rpy2.robjects.globalenv['annotations'] = "custom.clusterv2"
        rpy2.robjects.globalenv['annotations_to_plot'] = "custom.clusterv2"
        rpy2.robjects.globalenv['ref_aggr_strategy'] = "max"
        r(
            '''
            sample_files <- c(sample_files)
            sample_names <- c(sample_names)
            annotations <- c(annotations)
            annotations_to_plot <- c(annotations_to_plot)
            '''
        )
        # Access all Kriegstein metadata
        r(
            '''getMeta <- function(kriegstein_data_dir) {
            metaFile <- file.path(kriegstein_data_dir, "custom.meta.tsv")
            if (!file.exists(metaFile)) {
                stop(metaFile, " does not exist, generate this file with EMC.SKlab.scRNAseq::chunk_kriegstein_data().")
            } else {
                return(utils::read.table(metaFile, header = TRUE, sep = "\t", as.is = TRUE, row.names = 1))
                }
            }
            '''
        )
        # Generate color palette for plotting
        r('''generate_color_palette <- function(type = 'mixed', n = NULL) {
        if (type == 'mixed') {
            palette <- c("#000000", "#DF536B", "#61D04F", "#2297E6", "#28E2E5",
                                "#CD0BBC", "#F5C710", "#9E9E9E", "#F564E3", "#B79F00",
                                "#E69F00", "#009E73", "#0072B2", "#D55E00", "#666666",
                                "#BEAED4", "#FFFF99", "#F0027F", "#A6761D", "#E31A1C",
                                "#16FF32", "#782AB6", "#E4E1E3", "#E5C494", "#CCCCCC",
                                "#B3CDE3", "#CCEBC5", "#DECBE4", "#FDDAEC", "#FFFFCC")
                                names(palette) <- c("black_R4.1", "salmon_R4.2", "lightgreen_R4.3",
                                                    "skyblue_R4.4", "teal_R4.5", "darkpink_R4.6",
                                                    "yellow_R4.7", "grey_R4.8", "pink_ggplot2.6",
                                                    "gold_ggplot2.7", "orange_OI.2", "darkgreen_OI.4",
                                                    "darkblue_OI.6", "darkorange_OI.7", "darkgray_Accent.8",
                                                    "violet_Accent.2", "lightyellow_Accent.4",
                                                    "brightpink_Accent.6", "brown_Dark2.7", "red_Paired.6",
                                                    "brightgeen_Alphabet.7", "purple_Alphabet.4",
                                                    "offwhite_Polychrome36.2", "sand_Set2.7",
                                                    "lightgrey_Pastel2.8", "pastelblue_Pastel1.2",
                                                    "pastelgreen_Pastel1.3", "pastelpurple_Pastel1.4",
                                                    "pastelpink_Pastel1.8", "pastelyellow_1.6")
        }
        if (type == 'colorblind') {
            palette <- grDevices::palette.colors(palette = "Okabe-Ito")
        }
        if (!is.null(n)) {
            if (n > length(palette)) {
            message(paste0('WARNING: asked for ', n, ' colors but only ', length(palette), ' are available and selected from the palette'))
            n <- length(palette)
            }
            palette <- palette[1:n]
        }
        return(palette)
        }
        palette
        '''
        )
        # Create the scores heatmap
        r(
            '''
            visualize_kriegstein_annotated_data <- function(
        sample_names, sample_files, output_dir,
        kriegstein_data_dir, kriegstein_annotated_input_dir,
        annotations = c("custom.clusterv2"),
        annotations_to_plot = c("custom.clusterv2"),
        ref_aggr_strategy = "max") {
    
    library(Seurat)
    dir.create(output_dir, recursive = TRUE)
    names(sample_files) <- sample_names
    meta <- getMeta(kriegstein_data_dir)
    # read and get .rds data
    data.list <- lapply(X = sample_files, FUN = function(x) {
    return(readRDS(file = x))
    })
    names(data.list) <- sample_names
    # initialize list to store SingleR results
    results.list <- list()
    # iterate sample_names and annotations: grab all corresponding files, then load and get corresponding data
    for (sample in sample_names) {
        for (anno in annotations) {
                files <- list.files(path=kriegstein_annotated_input_dir, full.names=T)
                results <- sapply(files, function(file) {
                load(file)
                return(result)
                })
                results.list[[paste(sample, anno)]] <- results
                }
        }
    
        combined.results <- lapply(X=results.list, FUN = function(x){
        combined <- SingleR::combineCommonResults(results = x)
        df <- as.data.frame(combined$scores)

        combined$max.scores <- sapply(unique(names(df)), function(names) {
            sapply(1:nrow(df), function(row){
                col <- max.col(df[names(df) == names]) [row]
                df[names(df) == names][row, col]
            })
        })
        combined$max.labels <- colnames(combined$max.scores)[max.col(combined$max.scores)]
    
        # Write scores to reference figures
        utils::write.csv2(combined$max.scores, file = file.path(output_dir, paste0("Kriegstein_Pearson.correlation.max_", sample, "_", anno ,".csv")))
        
        # get mean score and label of all references
        # for each unique column name select all columns
        combined$mean.scores <- sapply(unique(names(df)), function(names) {
        rowMeans(df[names(df) == names])
        })
        combined$mean.labels <- colnames(combined$mean.scores)[max.col(combined$mean.scores)]
        write.csv2(combined$mean.scores, file = file.path(output_dir, paste0("Kriegstein_Pearson.correlation.max_", sample, "_", anno ,".csv")))
        return(combined)
        })

        # save Kriegstein cluster labels into Seurat object (rds)
        for (sample in sample_names){
        for (anno in annotations) {
            misc_data <- combined.results[[paste(sample, anno)]]
            SeuratObject::Misc(object = data.list[[sample]], slot = paste0("Kriegstein.SingleR.", anno)) <- misc_data
        }
        SeuratObject::Misc(object = data.list[[sample]], slot = "Kriegstein.SingleR.ref_aggr_strategy") <- ref_aggr_strategy
        data.list[[sample]]$kriegstein.seurat.custom.clusters.mean <- data.list[[sample]]$seurat_clusters
        levels(data.list[[sample]]$kriegstein.seurat.custom.clusters.mean) <- paste0(combined.results[[paste(sample, "custom.clusterv2")]]$mean.labels, ".", levels(data.list[[sample]]$kriegstein.seurat.custom.clusters.mean))
        saveRDS(data.list[[sample]], file = sample_files[[sample]])
        }
        # custom visualizations per sample_reference comparison for each annotation
        for (sample in sample_names) {
            annotation_col <- data.frame(row.names = levels(data.list[[sample]]$seurat_clusters))
            annotation_colors <- list()
            for (anno in annotations) {
            annotation_col[ , ncol(annotation_col) + 1] <- data.list[[sample]]@misc[[paste0("Kriegstein.SingleR.", anno)]][[paste0(ref_aggr_strategy, ".labels")]]
            colnames(annotation_col)[ncol(annotation_col)] <- paste0("ref.", anno)
            labels <- unique(meta[which(colnames(meta) == anno)][,1])
            annotation_colors[[length(annotation_colors) + 1]] <- labels[order(labels)]
            names(annotation_colors)[length(annotation_colors)] <- paste0("ref.", anno)
    
            names(annotation_colors[[paste0("ref.", anno)]]) <- generate_color_palette(type = 'mixed', n = length(annotation_colors[[paste0("ref.", anno)]]))

            annotation_colors[[paste0("ref.", anno)]] <- stats::setNames(names(annotation_colors[[paste0("ref.", anno)]]), annotation_colors[[paste0("ref.", anno)]])
            annotation_colors[[paste0("ref.", anno)]] <- annotation_colors[[paste0("ref.", anno)]][unique(annotation_col[which(colnames(annotation_col) == paste0("ref.", anno))][,1])]
            }
        names(annotation_col) <- c("fetal.brain.celltype")
        annotation_col <- annotation_col[, c("fetal.brain.celltype")]

        for (anno in annotations) {
            filename <- file.path(output_dir, paste0("Kriegstein_heatmap_", sample, "_", anno ,".png"))
            message("plotting: ", filename)
            # set rownames for identification of rows during plotting behind <- prints list of clusters
            rownames(combined.results[[paste(sample, anno)]][[paste0(ref_aggr_strategy, ".scores")]]) <- levels((data.list[[sample]]$seurat_clusters))
            p <- pheatmap::pheatmap(t(combined.results[[paste(sample, anno)]][[paste0(ref_aggr_strategy, ".scores")]]),
                                    color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "PiYG"))(100),
                                    main = "Scores",
                                    filename = filename)
            print(levels((data.list[[sample]]$seurat_clusters)))
            }
        }
        }
        '''
        )
        # Call function above
        r(
            '''
            visualize_kriegstein_annotated_data(sample_names=sample_names,
            sample_files=sample_files,
            output_dir=output_dir,
            kriegstein_data_dir=kriegstein_data_dir,
            kriegstein_annotated_input_dir=kriegstein_annotated_input_dir,
            annotations=annotations,
            annotations_to_plot=annotations_to_plot,
            ref_aggr_strategy="max")
            '''
        )
        
    
    def annotated_rds_to_anndata_transformer(self):
        """
        Transforms the now annotated Seurat object back to an AnnData object.
        """
        # Inherrit Python variables to R.
        rpy2.robjects.globalenv['sample_file'] = os.path.join(self.sample_output, 'Sample_RDS_storage', f'{self.sample_name}.rds')
        self.annotated_anndata_storage = os.path.join(self.sample_output, 'AnnData_storage', f'{self.sample_name}_annotated.h5ad')
        rpy2.robjects.globalenv['anndata_storage'] = self.annotated_anndata_storage
        r(
            '''
            sample_data <- readRDS(file = sample_file)
            # Initialize a flag to track whether the code succeeded
            success <- FALSE
            
            run_conversion <- function() {
            sceasy::convertFormat(
                sample_data,
                from = "seurat",
                to = "anndata",
                outFile = anndata_storage
                )
            }

            while (!success) {
            tryCatch({
                run_conversion()
                # If no error occurred, set the success flag to TRUE and exit the loop
                success <- TRUE
            }, error = function(err) {
                cat("An error occurred:", conditionMessage(err), "\n")
                cat("Retrying...\n")
                })
            }
            cat("Code executed successfully!\n")
            '''
        )


    def visualize_annotation(self):
        """
        Loads in the annotated AnnData object.
        Renames the large annotation column name to something smaller.
        Plots a UMAP with the annotations.
        Saves the figure to {self.sample_name}_figures_output.
        """
        adata = sc.read_h5ad(self.annotated_anndata_storage)
        adata.obs['seurat_clusters_annotated'] = adata.obs['kriegstein.seurat.custom.clusters.mean']
        sc.pl.umap(adata, color=['seurat_clusters', 'seurat_clusters_annotated'], legend_loc='right margin', 
                   legend_fontsize='medium', legend_fontweight='normal',show=False)
        plt.savefig(os.path.join(self.sample_output, f"{self.sample_name}_figures_output", f"{self.sample_name}_UMAP.png"))


    def run(self):
        """
        Runs all class functions in chronological order. 
        Created for aesthetic purposes, to avoid running every function in the __init__.
        A final message notifies the user that this step of the pipeline has been completed.
        """
        self.makedirs()
        self.anndata_to_rds_transformer()
        self.annotate_with_singler()
        self.visualize_annotated_data()
        self.annotated_rds_to_anndata_transformer()
        self.visualize_annotation()
        print(f'{self.sample_name} annotated succesfully!')