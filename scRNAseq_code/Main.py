"""
Author:         	Julia Christiaanse
Date:               29-11-2023
Script:             Main script
Python version:     3.10.9
Summary:            Within this script, the pipeline structure is ran. Multiprocessing is systematically applied if possible.
Imports:
"""
import os
import multiprocessing as mp
from Sample_analysis import Sample_Analysis
from Sample_integration import Sample_integration
from Final_dea_GSEA import Final_Dea
from Cell_type_annotation import Cell_Type_Annotation
from TI_RNA_velocity import TI_RNA_Velocity
import time

# Create time variable to time code
s = time.time()

def single_sample_helper(queue, sample_dir, sample, output_dir, markerpath):
    """
    Helper function for class Sample_Analysis.

    Args:
        queue (Queue): Queue for process
        sample_dir (String): Location of the 10X Genomics output folder containing all samples
        sample (String): Name of the sample
        output_dir (String): Path to the master output directory
        markerpath (String): Path to markergenes.txt
    """    
    adata = Sample_Analysis(sample_dir, sample, output_dir, markerpath)
    queue.put(adata)


def integration_helper(queue, mono_h5ad, co_h5ad, full_name, output_dir, sample_output, markerpath):
    """
    Helper function for class Sample_integration.

    Args:
        queue (Queue): Queue for process
        mono_h5ad (String): Name of the co-culture sample
        co_h5ad (String): Name of the mono-culture sample
        full_name (String): Full name of integrated sample
        output_dir (String): Path to the master output directory
        sample_output (String): Path to the sample-specific output sub-folder within output_dir
        markerpath (String): Path to markergenes.txt
    """    
    concat = Sample_integration(co_h5ad, mono_h5ad, full_name, output_dir, sample_output, markerpath)
    queue.put(concat)


def annotation_helper(queue, sample_name, adata_loc, output_dir, sample_output, reference_data_dir, path_to_ref_data):
    """
    Helper function for class Cell_Type_Annotation.

    Args:
        queue (Queue): Queue for process
        sample_name (String): Name of the sample
        adata_loc (String): Path to integrated AnnData object
        output_dir (String): Path to master output directory
        sample_output (String): Path to the sample-specific output sub-folder within output_dir
        reference_data_dir (String): Path to all reference data
        path_to_ref_data (String): Path to folder in reference_data_dir where Kriegstein chunks are stored
    """    
    annotated = Cell_Type_Annotation(sample_name, adata_loc, output_dir, sample_output, reference_data_dir, path_to_ref_data)
    queue.put(annotated)


def tirv_helper(queue, sample_name, annotated_adata_loc, tirv_output_loc, path_to_loomp_C, mono_loomp_files):
    """
    Helper function for TI_RNA_Velocity.

    Args:
        queue (Queue): Queue for process
        sample_name (String): Name of the sample
        annotated_adata_loc (String): Path to annotated AnnData object
        tirv_output_loc (String): Path to subfolder to store output files
        path_to_loomp_C (String): Path to co-culture .loom file
        mono_loomp_files (List): List of paths to mono-culture .loom files
    """    
    tirv_analysis = TI_RNA_Velocity(sample_name, annotated_adata_loc, tirv_output_loc, path_to_loomp_C, mono_loomp_files)
    queue.put(tirv_analysis)

    
def dea_cellrank_helper(queue, sample_name, annotated_adata_loc, final_dea_output_loc, output_dir, markerpath):
    """
    Helper function for class Final_Dea.

    Args:
        queue (Queue): Queue for process
        sample_name (String): Name of the Sample
        annotated_adata_loc (String): Path to annotated AnnData object
        final_dea_output_loc (String): Path to subfolder to store output files
        output_dir (String): Path to the master output directory
        markerpath (String): Path to markergenes.txt
    """    
    deacellrank = Final_Dea(sample_name, annotated_adata_loc, final_dea_output_loc, output_dir, markerpath)
    queue.put(deacellrank)



if __name__ == '__main__':
    """
    Main script. Pathways to necessary input files and the output directory are defined here.
    Programm is given pointers which sample is co-culture, and which are mono-cultures.

    Raises:
        ValueError: Check if the output directory provided does not already exist.
    """    
    # Location of the 10X Genomics output for all three samples
    sample_dir = 'C:/Users/julia/scRNAseq_Analysis_project/samples_dir'
    # Location of chunks
    reference_data_dir = 'C:/Users/julia/Project/data/chunks_25'
    # Location of all Kriegstein data: folder with chunks, metadata, custom annotation from lab, etc.
    path_to_ref_data = 'C:/Users/julia/Project/data'
    # List of sample names in sample_dir
    sample_names = [folder.name for folder in os.scandir(sample_dir) if folder.is_dir()]
    # Main output directory
    output_dir = 'C:/Users/julia/Project/post_documentation_run'
    # Path to markergenes.txt
    markerpath = 'C:/Users/julia/Project/markergenes.txt'
    # Paths to all .loom files
    path_to_loomp_A = 'C:/Users/julia/Project/data/BL_A/velocyto/BL_A.loom'
    path_to_loomp_C = 'C:/Users/julia/Project/data/BL_C/velocyto/BL_C.loom'
    path_to_loomp_N = 'C:/Users/julia/Project/data/BL_N/velocyto/BL_N.loom'
    # Small list with only the .loom files of the mono-cultures
    mono_loomp_files = [path_to_loomp_A, path_to_loomp_N]

    # give programm pointers as to which sample is co culture and which are mono culture
    co_anndata = 'BL_C'
    mono_anndatas = ['BL_N', 'BL_A']
    
    # Check if directory exists
    if os.path.exists(output_dir):
        raise ValueError('Your output directory already exists! Please rename or remove it.')
    
    # # -----------------------------------------------------------------------------------------
    # # Below all steps of the pipeline are performed, a Queue is created, and samples are added to it.
    # # Pipeline steps are run in parallel for samples within the Queue. Once the Queue is done,
    # # the pipeline moves on to the next step.
    # # Step 1-3 Individual Sample Analysis

    q = mp.Queue()
    processes = []
    adatas = []
    
    for sample in sample_names:
        p = mp.Process(target=single_sample_helper, args=[q,
                                                           sample_dir, 
                                                           sample, output_dir, 
                                                           markerpath])
        processes.append(p)
        p.start()
    for p in processes:
        adata = q.get()
        adatas.append(adata)
    for p in processes:
        p.join()

    # # -----------------------------------------------------------------------------------------
    # # Steps 4 and 5 data integration & cell selection 
    
    q = mp.Queue()
    integration_process = []
    integrated_samples = []
    names_of_samples = []

    for solo_sample in mono_anndatas:
        full_name = solo_sample+"_"+co_anndata
        sample_output = os.path.join(output_dir, full_name)
        p = mp.Process(target=integration_helper, args=[q,
                                                        solo_sample,
                                                        co_anndata,
                                                        full_name,
                                                        output_dir,
                                                        sample_output,
                                                        markerpath])
        names_of_samples.append(full_name)
        integration_process.append(p)
        p.start()
    for p in integration_process:
        concat = q.get()
        integrated_samples.append(concat)
    for p in integration_process:
        p.join()

    # # -----------------------------------------------------------------------------------------
    # # Step 6 celltype annotation
    # names_of_samples=['BL_A_BL_C', 'BL_N_BL_C'] # uncomment this line if you want to run step 6-7B seperately.
    q = mp.Queue()
    annotation_process = []
    annotated_anndatas = []

    for sample_name in names_of_samples:
        adata_loc = os.path.join(output_dir, sample_name, "AnnData_storage")
        sample_output = os.path.join(output_dir, "Cell_type_annotation", sample_name)
        p = mp.Process(target=annotation_helper, args=[q, sample_name,
                                                        adata_loc, output_dir,
                                                          sample_output,
                                                            reference_data_dir,
                                                              path_to_ref_data])
        annotation_process.append(p)
        p.start()
    for p in annotation_process:
        annotated = q.get()
        annotated_anndatas.append(annotated)
    for p in annotation_process:
        p.join()
    
    # # -----------------------------------------------------------------------------------------
    # # Step 7 A TI and RNA velocity
    
    q = mp.Queue()
    tirv_process = []
    tirv_analyzed_list = []
    
    for sample_name in names_of_samples:
        annotated_adata_loc = os.path.join(output_dir, "Cell_type_annotation", sample_name, "AnnData_storage")
        tirv_output_loc = os.path.join(output_dir, "TI_RNA_velocity", sample_name)
        p = mp.Process(target=tirv_helper, args=[q, sample_name, annotated_adata_loc, tirv_output_loc, path_to_loomp_C, mono_loomp_files])
        tirv_process.append(p)
        p.start()
    for p in tirv_process:
        tirv_analyzed = q.get()
        tirv_analyzed_list.append(tirv_analyzed)
    for p in tirv_process:
        p.join()

    # # -----------------------------------------------------------------------------------------
    # # Step 7B culture DEA and GSEA

    q = mp.Queue()
    cellrank_dea_process = []
    object_storage = []

    for sample_name in names_of_samples:
        annotated_adata_loc = os.path.join(output_dir, "Cell_type_annotation", sample_name, "AnnData_storage")
        final_dea_output_loc = os.path.join(output_dir, "Final_DEA_GSEA", sample_name)
        p = mp.Process(target=dea_cellrank_helper, args=[q, sample_name, annotated_adata_loc, final_dea_output_loc, output_dir, markerpath])
        cellrank_dea_process.append(p)
        p.start()
    for p in cellrank_dea_process:
        deacellrank = q.get()
        object_storage.append(deacellrank)
    for p in cellrank_dea_process:
        p.join()

    # -----------------------------------------------------------------------------------------

    # end timer
    e = time.time()
    finish = round((e-s),2)

    print("\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n",
          f"\t    PIPELINE FINISHED IN {finish} SECONDS!",
           "\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n")
    print(f'♥ ♥ ♥ ♥   YOUR FILES ARE HERE: {output_dir}  ♥ ♥ ♥ ♥\n')