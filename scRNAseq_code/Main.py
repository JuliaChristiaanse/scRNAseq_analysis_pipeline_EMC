import os
import multiprocessing as mp
from Sample_analysis import Sample_Analysis
from Sample_integration import Sample_integration
from Final_dea import Final_Dea
from Cell_type_annotation import Cell_Type_Annotation
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

    
# def dea_cellrank_helper(queue, subset_adata_path, output_dir, markerpath):
#     deacellrank = Final_Dea(subset_adata_path, output_dir, markerpath)
#     queue.put(deacellrank)



if __name__ == '__main__':

    # Create directories and check if they exist.
    sample_dir = 'C:/Users/julia/scRNAseq_Analysis_project/samples_dir'
    reference_data_dir = 'C:/Users/julia/Project/data/chunks_25'
    path_to_ref_data = 'C:/Users/julia/Project/data'
    sample_names = [folder.name for folder in os.scandir(sample_dir) if folder.is_dir()]
    output_dir = 'C:/Users/julia/Project/final_try'
    markerpath = 'C:/Users/julia/Project/markergenes.txt'

    # give programm pointers as to which sample is co culture and which are mono culture
    co_anndata = 'BL_C'
    mono_anndatas = ['BL_N', 'BL_A']
    
    if os.path.exists(output_dir):
        raise ValueError('Your output directory already exists! Please rename or remove it.')
    
    # # -----------------------------------------------------------------------------------------
    # # Step 1-3 Individual Sample Analysis

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

    
    # q = mp.Queue()
    # annotation_process = []
    
    # for sample_name in integrated_samples:
    #     adata_loc = os.path.join(output_dir, "integrated_"+sample_name, "AnnData_storage")
    #     sample_output = os.path.join(output_dir,  "Cell_type_annotation", sample_name)
    #     p = mp.Process(target=annotation_helper, args=[q, sample_name, adata_loc, output_dir, sample_output, reference_data_dir, path_to_ref_data])
    #     annotation_process.append(p)
    #     p.start()
    # for p in annotation_process:
    #     annotated_adata = q.get()
    # for p in annotation_process:
    #     p.join()


    # -----------------------------------------------------------------------------------------
    # # # Step 7 small DEA, pseudotime & RNA velocity

    # # for i in combis:      # use this for later !
    # #     path = i.path
    # #     full_name = i.full_name
    
    # astroco_subset = os.path.join(output_dir,'integrated_BL_C_BL_A',
    #                               'AnnData_storage',
    #                               'BL_C_BL_A_subset_anndata.h5ad')
    # neuroco_subset = os.path.join(output_dir,'integrated_BL_C_BL_N',
    #                               'AnnData_storage',
    #                               'BL_C_BL_N_subset_anndata.h5ad')
    # subset_adatas_paths = [astroco_subset, neuroco_subset]

    # q = mp.Queue()
    # cellrank_dea_process = []
    # object_storage = []

    # for subset_adata in subset_adatas_paths:
    #     p = mp.Process(target=dea_cellrank_helper, args=[q, subset_adata, output_dir, markerpath])
    #     cellrank_dea_process.append(p)
    #     p.start()
    # for p in cellrank_dea_process:
    #     deacellrank = q.get()
    #     object_storage.append(deacellrank)
    # for p in cellrank_dea_process:
    #     p.join()

    # -----------------------------------------------------------------------------------------

    # end timer
    e = time.time()
    finish = round((e-s),2)

    print("\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n",
          f"\t    PIPELINE FINISHED IN {finish} SECONDS!",
           "\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n")
    print(f'♥ ♥ ♥ ♥   YOUR FILES ARE HERE: {output_dir}  ♥ ♥ ♥ ♥\n')