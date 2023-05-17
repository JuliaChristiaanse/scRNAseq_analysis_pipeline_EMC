
import os
import multiprocessing as mp
from individual_sample_analysis import Sample_Analysis
import time

# Create time variable to time code
s = time.time()

# Simple helper function allows for the calling of multiple OOP-imports that represent steps in the pipeline.
def add_helper(queue, sample_dir, sample, output_dir):
    # import Sample_Analysis performs individual analysis
    adata = Sample_Analysis(sample_dir, sample, output_dir)
    queue.put(adata)

# Main where users can provide their input and output directories
# Maybe add user input option that asks the user to paste the path and provide an
# output directory name?
if __name__ == '__main__':

    # Create directories and check if they exist.
    sample_dir = 'C:/Users/julia/scRNAseq_Analysis_project/samples_dir'
    #sample_names = ['BL_N'] 
    sample_names = [folder.name for folder in os.scandir(sample_dir) if folder.is_dir()]
    output_dir = 'C:/Users/julia/Project/output_dir_doublets_removed'
    
    if os.path.exists(output_dir):
        raise ValueError('Your output directory already exists! Please rename or remove it :)')
    
    # Create Queue, and lists to store output of add helper function and processes.
    q = mp.Queue()
    processes = []
    adatas = []
    
    # For loops that starts add_helper() function process and get the output of the process.
    for sample in sample_names:
        p = mp.Process(target=add_helper, args=[q, sample_dir, sample, output_dir])
        processes.append(p)
        p.start()
    for p in processes:
        adata = q.get()
        adatas.append(adata)
    for p in processes:
        p.join()
    
    # checks if the output of the add_helper() has altered anndata, will be removed.
    for adata in adatas:
        print(adata.adata)

    # end timer
    e = time.time()
    finish = round((e-s),2)

    print("\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n",
          f"PRE-PROCESSING, NORMALISATION & CLUSTERING FINISHED IN {finish} SECONDS!",
           "\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n")
    print(f'YOUR FILES ARE HERE: {output_dir}\n')