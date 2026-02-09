from numpy.typing import NDArray
import numpy as np
from collections import namedtuple
from ecifp.utils import (tanimoto_uint32,
                         tanimoto_batch,
                         get_component,
                         generate_ecfp)
from scipy.cluster.hierarchy import linkage, fcluster

ECIFP = namedtuple("ECIFP", "ligand_ifp protein_ifp")

def get_batch_ecifp_sim(
    ligand_ecifp: NDArray[np.uint32],
    protein_ecifp: NDArray[np.uint32],
    ligand_coeff: float,
    protein_coeff: float,
) -> NDArray[np.float32]:
    """
    Compute similarity between ECIFPs.

    Parameters:
        ligand_ecifp: Ligand fingerprint matrix of shape (n_samples, n_words)
            with dtype uint32.
        protein_ecifp: Protein fingerprint matrix of shape (n_samples, n_words)
            with dtype uint32.
        ligand_coeff: Exponent weight for ligand similarity contribution.
        protein_coeff: Exponent weight for protein similarity contribution.

    Returns:
        1D array of combined pairwise similarities with shape
        (n_samples*(n_samples-1)//2,) in upper-triangle row-major order.
    """
    ligand_ecifp_sim = tanimoto_batch(ligand_ecifp)
    protein_ecifp_sim = tanimoto_batch(protein_ecifp)
    ecifp_sim = (ligand_ecifp_sim**ligand_coeff) * (protein_ecifp_sim**protein_coeff)
    return ecifp_sim

def get_intx_ecfp_indices(ccd_id:str, intx_atoms:set[str], fp_radius=2, fp_size=102):
    component = get_component(ccd_id)
    ecfp = generate_ecfp(component.mol, fp_radius, fp_size, atom_to_bits=True)
    on_indices = set()
    for atom in component.mol.GetAtoms():
        atom_name = atom.GetProp("name") if atom.HasProp("name") else None
        if atom_name in intx_atoms:
            on_indices = on_indices.union(set(ecfp[atom.GetIdx()]))

    return on_indices

def get_ligand_ifp(ligand_id:str, intx_atoms:set[str], fp_radius, fp_size):
    ifp = np.zeros(fp_size)
    intx_ecfp_indices = get_intx_ecfp_indices(ligand_id, intx_atoms, fp_radius, fp_size)
    ifp[list(intx_ecfp_indices)] = 1
    return ifp

def get_polymer_ifp(intx_polymer_atoms:dict[str, set[str]], fp_radius, fp_size):
    ifp = np.zeros(fp_size)
    on_indices = set()
    for residue, intx_atoms in intx_polymer_atoms.items():
        residue_on_indices = get_intx_ecfp_indices(residue, intx_atoms, fp_radius, fp_size)
        on_indices = on_indices.union(residue_on_indices)

    ifp[list(on_indices)] = 1
    return ifp

def get_ecifp_sim(ecifp_1, ecifp_2,ligand_coeff, protein_coeff):
    ligand_sim = tanimoto_uint32(ecifp_1.ligand_ifp.view(np.uint32),
                                 ecifp_2.ligand_ifp.view(np.uint32))
    protein_sim = tanimoto_uint32(ecifp_1.protein_ifp.view(np.uint32),
                                  ecifp_2.protein_ifp.view(np.uint32))
    ecifp_sim = (ligand_sim**ligand_coeff) * (protein_sim**protein_coeff)
    return ecifp_sim

def get_optimal_cutoff(linkage_matrix, min_cutoff, max_cutoff):
    heights = linkage_matrix[:, 2]
    gaps = np.diff(heights)
    idxs = np.argsort(gaps)[::-1]
    for _id in idxs:
        if np.max(min_cutoff < heights[_id] <= max_cutoff):
            return round(heights[_id], 3)
    return 0.5


def overlap_similarity(X, Y):
    """
    Compute Overlap similarity between two binary matrices using matrix operations.
    Parameters:
        X: np.ndarray, shape (n, d) - First binary matrix
        Y: np.ndarray, shape (m, d) - Second binary matrix
    Returns:
        np.ndarray, shape (n, m) - Overlap similarity matrix
    """
    # Compute intersection: X @ Y^T (dot product)
    intersection = np.dot(X, Y.T)
    # Compute sums for each row (number of 1s in each vector)
    X_sum = X.sum(axis=1, keepdims=True)  # Shape (n, 1)
    Y_sum = Y.sum(axis=1, keepdims=True).T  # Shape (1, m)
    # Compute minimum fingerprint length
    min_fp = np.minimum(X_sum, Y_sum)
    
    # Compute Overlap similarity
    overlap_sim = intersection / min_fp  # Element-wise division
    # Handle division by zero (if both vectors are zero)
    overlap_sim[np.isnan(overlap_sim)] = 0
    return overlap_sim

def get_bs_clusters(data_matrix):
    distance_matrix = 1 - overlap_similarity(data_matrix, data_matrix)
    distance = 0.4
    N = data_matrix.shape[0]

    # Extract condensed form (upper triangle)
    tri_rows, tri_cols = np.triu_indices(N, k=1)
    condensed_distances = distance_matrix[tri_rows, tri_cols]
    
    # Hierarchical clustering
    Z = linkage(condensed_distances, method='single', metric='precomputed')
    
    # Get cluster assignments
    clusters = fcluster(Z, t=distance, criterion='distance')
    return clusters