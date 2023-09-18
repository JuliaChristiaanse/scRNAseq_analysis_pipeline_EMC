
import os
import multiprocessing as mp
from Sample_analysis import Sample_Analysis
from Sample_integration import Sample_integration
from dea_cellrank import DEA_Cellrank
import time

# Create time variable to time code
s = time.time()

# Simple helper function allows for the calling of multiple OOP-imports that represent steps in the pipeline.
def single_sample_helper(queue, sample_dir, sample, output_dir, markerpath):
    # import Sample_Analysis performs individual analysis
    adata = Sample_Analysis(sample_dir, sample, output_dir, markerpath)
    queue.put(adata)

def integration_helper(queue, mono_h5ad, co_h5ad, output_dir, markerpath):
    concat = Sample_integration(co_h5ad, mono_h5ad, output_dir, markerpath)
    queue.put(concat)

# def dea_cellrank_helper(queue, subset_adata_path, output_dir, markerpath):
#     deacellrank = DEA_Cellrank(subset_adata_path, output_dir, markerpath)
#     queue.put(deacellrank)


# Main where users can provide their input and output directories
# Maybe add user input option that asks the user to paste the path and provide an
# output directory name?
if __name__ == '__main__':

    # Create directories and check if they exist.
    sample_dir = 'C:/Users/julia/scRNAseq_Analysis_project/samples_dir'
    sample_names = [folder.name for folder in os.scandir(sample_dir) if folder.is_dir()]
    output_dir = 'C:/Users/julia/Project/Output_directory_test1'
    markerpath = 'C:/Users/julia/Project/markergenes.txt'

    # give programm pointers as to which sample is co culture and which are mono culture
    # !! Might have to be changed, as it is hard coded, not user friendly.
    co_anndata = 'BL_C'
    mono_anndatas = ['BL_N', 'BL_A']
    
    if os.path.exists(output_dir):
        raise ValueError('Your output directory already exists! Please rename or remove it :)')
    
    # -----------------------------------------------------------------------------------------
    # Step 1-3 Individual Sample Analysis

    # Create Queue, and lists to store output of add helper function and processes.
    q = mp.Queue()
    processes = []
    adatas = []
    
    # For loops that starts add_helper() function process and get the output of the process.
    for sample in sample_names:
        p = mp.Process(target=single_sample_helper, args=[q, sample_dir, sample, output_dir, markerpath])
        processes.append(p)
        p.start()
    for p in processes:
        adata = q.get()
        adatas.append(adata)
    for p in processes:
        p.join()

    # -----------------------------------------------------------------------------------------
    # Steps 4 and 5 data integration & cell selection 
    
    q = mp.Queue()
    integration_process = []
    combis = []

    for solo_sample in mono_anndatas:
        p = mp.Process(target=integration_helper, args=[q, solo_sample, co_anndata, output_dir, markerpath])
        integration_process.append(p)
        p.start()
    for p in integration_process:
        concat = q.get()
        combis.append(concat)
    for p in integration_process:
        p.join()

    # -----------------------------------------------------------------------------------------
    # Step 6 celltype annotation



    # -----------------------------------------------------------------------------------------
    # # Step 7 small DEA, pseudotime & RNA velocity

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