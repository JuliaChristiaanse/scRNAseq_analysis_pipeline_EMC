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

def single_sample_helper(queue, sample_dir, sample, output_dir, markerpath, reference_data_dir, path_to_ref_data):
    adata = Sample_Analysis(sample_dir, sample, output_dir, markerpath, reference_data_dir, path_to_ref_data)
    queue.put(adata)


def integration_helper(queue, mono_h5ad, co_h5ad, full_name, output_dir, sample_output, markerpath):
    concat = Sample_integration(co_h5ad, mono_h5ad, full_name, output_dir, sample_output, markerpath)
    queue.put(concat)


def annotation_helper(queue, sample_name, adata_loc, output_dir, sample_output, reference_data_dir, path_to_ref_data):
    annotated = Cell_Type_Annotation(sample_name, adata_loc, output_dir, sample_output, reference_data_dir, path_to_ref_data)
    queue.put(annotated)


def tirv_helper(queue, sample_name, annotated_adata_loc, tirv_output_loc, path_to_loomp_C, mono_loomp_files):
    tirv_analysis = TI_RNA_Velocity(sample_name, annotated_adata_loc, tirv_output_loc, path_to_loomp_C, mono_loomp_files)
    queue.put(tirv_analysis)

    
def dea_cellrank_helper(queue, sample_name, annotated_adata_loc, final_dea_output_loc, output_dir, markerpath):
    deacellrank = Final_Dea(sample_name, annotated_adata_loc, final_dea_output_loc, output_dir, markerpath)
    queue.put(deacellrank)



if __name__ == '__main__':
    # Create directories and check if they exist.
    sample_dir = 'C:/Users/julia/scRNAseq_Analysis_project/samples_dir'
    reference_data_dir = 'C:/Users/julia/Project/data/chunks_25'
    path_to_ref_data = 'C:/Users/julia/Project/data'
    sample_names = [folder.name for folder in os.scandir(sample_dir) if folder.is_dir()]
    output_dir = 'C:/Users/julia/Project/15_11_2023_re-run'
    markerpath = 'C:/Users/julia/Project/markergenes.txt'
    path_to_loomp_A = 'C:/Users/julia/Project/data/BL_A/velocyto/BL_A.loom'
    path_to_loomp_C = 'C:/Users/julia/Project/data/BL_C/velocyto/BL_C.loom'
    path_to_loomp_N = 'C:/Users/julia/Project/data/BL_N/velocyto/BL_N.loom'
    mono_loomp_files = [path_to_loomp_A, path_to_loomp_N]

    # give programm pointers as to which sample is co culture and which are mono culture
    co_anndata = 'BL_C'
    mono_anndatas = ['BL_N', 'BL_A']
    
    # if os.path.exists(output_dir):
    #     raise ValueError('Your output directory already exists! Please rename or remove it.')
    
    # # -----------------------------------------------------------------------------------------
    # Step 1-3 Individual Sample Analysis

    q = mp.Queue()
    processes = []
    adatas = []
    
    for sample in sample_names:
        p = mp.Process(target=single_sample_helper, args=[q,
                                                           sample_dir, 
                                                           sample, output_dir, 
                                                           markerpath, 
                                                           reference_data_dir, 
                                                           path_to_ref_data])
        processes.append(p)
        p.start()
    for p in processes:
        adata = q.get()
        adatas.append(adata)
    for p in processes:
        p.join()

    # # -----------------------------------------------------------------------------------------
    # Steps 4 and 5 data integration & cell selection 
    
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
    # names_of_samples=['BL_A_BL_C', 'BL_N_BL_C'] # REMOVE this before running whole pipeline
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
    # Step 7 A TI and RNA velocity
    
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