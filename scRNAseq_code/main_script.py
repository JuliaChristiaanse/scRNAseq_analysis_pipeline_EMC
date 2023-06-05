
import os
import multiprocessing as mp
from individual_sample_analysis import Sample_Analysis
from sample_integration import Data_integration
import time

# Create time variable to time code
s = time.time()

#Simple helper function allows for the calling of multiple OOP-imports that represent steps in the pipeline.
def single_sample_helper(queue, sample_dir, sample, output_dir):
    # import Sample_Analysis performs individual analysis
    adata = Sample_Analysis(sample_dir, sample, output_dir)
    queue.put(adata)

def integration_helper(queue, mono_h5ad, co_h5ad, output_dir):
    concat = Data_integration(co_h5ad, mono_h5ad, output_dir)
    queue.put(concat)


# Main where users can provide their input and output directories
# Maybe add user input option that asks the user to paste the path and provide an
# output directory name?
if __name__ == '__main__':

    # Create directories and check if they exist.
    sample_dir = 'C:/Users/julia/scRNAseq_Analysis_project/samples_dir'
    sample_names = [folder.name for folder in os.scandir(sample_dir) if folder.is_dir()]
    output_dir = 'C:/Users/julia/Project/TEST2'

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
        p = mp.Process(target=single_sample_helper, args=[q, sample_dir, sample, output_dir])
        processes.append(p)
        p.start()
    for p in processes:
        adata = q.get()
        adatas.append(adata)
    for p in processes:
        p.join()
    

    # -----------------------------------------------------------------------------------------
    # Step 4 data integration 
    
    q = mp.Queue()
    integration_process = []
    combis = []

    for solo_sample in mono_anndatas:
        p = mp.Process(target=integration_helper, args=[q, solo_sample, co_anndata, output_dir])
        integration_process.append(p)
        p.start()
    for p in integration_process:
        concat = q.get()
        combis.append(concat)
    for p in integration_process:
        p.join()

    # -----------------------------------------------------------------------------------------

    # # end timer
    e = time.time()
    finish = round((e-s),2)

    print("\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n",
          f"PRE-PROCESSING, NORMALISATION & CLUSTERING FINISHED IN {finish} SECONDS!",
           "\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n")
    print(f'YOUR FILES ARE HERE: {output_dir}\n')