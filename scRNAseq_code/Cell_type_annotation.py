from os import listdir
from os.path import isfile, join
import scanpy as sc
import anndata as ad
import symphonypy as symp
import os
from matplotlib import pyplot as plt

# Log 15-8-2023
# all chunks can be loaded into memory one by one with astroco or neuroco in memory as well

# collect all Kriegstein chunks
path_to_data_storage = 'C:/Users/julia/Project/data/Julia Kriegstien Rdata2h5ad'
all_reference_chunks = [f for f in listdir(path_to_data_storage) if isfile(join(path_to_data_storage, f))]
path_to_adata_query = 'C:/Users/julia/Project/Output_directory2/integrated_BL_C_BL_A/AnnData_storage'
adata_query_name = 'BL_C_BL_A_subset_anndata.h5ad'


os.chdir('C:/Users/julia/Project/Annotations_astroco2')
os.makedirs('annotation_TSV')
os.makedirs('annotation_UMAPS')


def annotate_per_chunk(path_to_data_storage, all_reference_chunks, path_to_adata_query, adata_query_name):
    adata_query = sc.read_h5ad(os.path.join(path_to_adata_query, adata_query_name))
    print(adata_query)
    for chunk in all_reference_chunks:
        # start of pre-processing adata-ref
        adata_ref = sc.read(os.path.join(path_to_data_storage, chunk), cache=True)
        sc.pp.normalize_total(adata_ref, target_sum=1e4)
        sc.pp.log1p(adata_ref)
        sc.pp.highly_variable_genes(adata_ref) # edit set no roof to hvg
        adata_ref.raw = adata_ref
        adata_ref = adata_ref[:, adata_ref.var.highly_variable]
        sc.pp.scale(adata_ref, max_value=10)
        sc.pp.pca(adata_ref, n_comps=50)
        sc.pp.neighbors(adata_ref, n_neighbors=10, n_pcs=20)
        sc.tl.umap(adata_ref)
        sc.tl.leiden(adata_ref)
        # end of pre-processing adata_ref

        # start pre-processing adata_query
        adata_query.X = adata_query.layers['rawcounts'] # access raw counts
        sc.pp.normalize_total(adata_query, target_sum=1e4)
        sc.pp.log1p(adata_query)
        # Running Symphony 
        symp.tl.map_embedding(
            adata_query=adata_query,
            adata_ref=adata_ref
            )
        symp.tl.ingest(
            adata_query=adata_query,
            adata_ref=adata_ref,
            use_rep="X_pca_harmony", # choose xpca harmony as the adata_query rep
            )
        # Labels prediction
        symp.tl.transfer_labels_kNN(
            adata_query=adata_query,
            adata_ref=adata_ref,
            ref_labels=["custom.clusterv2", "leiden"],
            ref_basis="X_pca",          # use x_pca as basis for reference # WHY: because we don't use harmony integrate
            query_basis="X_pca_harmony",    # use x_pca_harmony as basis for query
        )
        # break of the Symphonypy_without_harmony_tutorial.ipynb


        # Add this from another vignette to add confidence score.
        symp.tl.per_cell_confidence(
            adata_query=adata_query,
            adata_ref=adata_ref,
            ref_basis_adjusted='X_pca' # no x_pca_harmony for this sample
        )

        # continue of the Symphonypy_without_harmony_tutrotial.ipynb to plot
        fig, axes = plt.subplots(figsize=(8, 4), ncols=2)
        sc.pl.umap(
            adata_query,
            color="custom.clusterv2",
            frameon=False,
            title="Query dataset",
            ax=axes[0],
            show=False,
            legend_loc=None,
        )

        sc.pl.umap(
            adata_ref,
            color="custom.clusterv2",
            frameon=False,
            title="Reference dataset",
            ax=axes[1],
            show=False,
        ) 
        # save the umaps & tsv's
        plt.savefig(os.path.join('C:/Users/julia/Project/Annotations_astroco2', 'annotation_UMAPS', f'annotated_{chunk}_UMAP.png'))
        adata_query.obs[['custom.clusterv2', 'symphony_per_cell_dist']].to_csv(os.path.join('C:/Users/julia/Project/Annotations_astroco2', 'annotation_TSV', f'celltype_symphony_{chunk}.tsv'), sep='\t')
        del(adata_ref)  # remove adata ref to save storage

annotate_per_chunk(path_to_data_storage, all_reference_chunks, path_to_adata_query, adata_query_name)