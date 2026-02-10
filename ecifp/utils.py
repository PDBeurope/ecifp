import numpy as np
from pathlib import Path
import requests
from requests.exceptions import RequestException, Timeout, HTTPError
from numpy.typing import NDArray
import pyarrow as pa
from numba import njit, uint32, int32
from functools import lru_cache
import os
import tempfile
from pdbeccdutils.core import ccd_reader
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import roc_curve, auc, RocCurveDisplay
import configparser


@njit(int32(uint32))
def popcount32(x: int) -> int:
    """
    Count the number of set bits (1s) in a 32-bit unsigned integer.

    Uses Brian Kernighan's algorithm which iterates only as many times
    as there are set bits.

    Parameters:
        x: A 32-bit unsigned integer.

    Returns:
        The number of bits set to 1 in x.
    """
    count = 0
    while x:
        x &= x - 1
        count += 1
    return count


@njit
def tanimoto_uint32(a: NDArray[np.uint32], b: NDArray[np.uint32]) -> float:
    """
    Compute the Tanimoto similarity between two bit-packed fingerprints.

    The Tanimoto coefficient (Jaccard index) is computed as:
        T(A, B) = |A ∩ B| / |A ∪ B|

    Parameters:
        a: First fingerprint as an array of uint32 words.
        b: Second fingerprint as an array of uint32 words (same length as a).

    Returns:
        Tanimoto similarity coefficient in the range [0.0, 1.0].
    """
    ands = 0
    ors = 0
    for i in range(a.shape[0]):
        ands += popcount32(a[i] & b[i])
        ors  += popcount32(a[i] | b[i])

    return ands / ors


@njit
def tanimoto_batch(fps: NDArray[np.uint32]) -> NDArray[np.float32]:
    """
    Compute all pairwise Tanimoto similarities for a set of fingerprints.

    Parameters:
        fps: Fingerprint matrix of shape (n_fps, n_words) with dtype uint32,
            where each row is a bit-packed fingerprint.

    Returns:
        1D array of shape (n_fps*(n_fps-1)//2,) with dtype float32 containing
        upper-triangle similarities (i < j) in row-major order.
    """
    n = fps.shape[0]
    dist = np.empty((n*(n-1))//2, dtype=np.float32)

    k = 0  # Index in 1D dists array
    for i in range(n):
        for j in range(i+1, n):
            dist[k] = tanimoto_uint32(fps[i], fps[j])
            k += 1

    return dist

def fetch_api_data(url, params=None, timeout=20):
    """
    Fetch data from API with proper error handling
    """
    try:
        response = requests.get(url, params=params, timeout=timeout)

        # Raise exception for 4xx/5xx status codes
        response.raise_for_status()

        # Test content type
        if 'application/json' in response.headers.get('Content-Type', ''):
            return response.json()
        else:
            return response.text

    except Timeout:
        print(f"Request timed out after {timeout} seconds")
        return None
    except HTTPError as e:
        print(f"HTTP error occurred: {e}")
        return None
    except RequestException as e:
        print(f"Request failed: {e}")
        return None

def generate_ecfp(mol, radius, fpSize, atom_to_bits=False):
    """Generates ECFP2 fingerprints for a RDKit mol object
    """
    ecfp = rdFingerprintGenerator.GetMorganGenerator(radius=radius,fpSize=fpSize)
    if not atom_to_bits:
        fp = ecfp.GetFingerprint(mol)
        return fp

    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateAtomToBits()
    fp = ecfp.GetFingerprint(mol,additionalOutput=ao)
    atom_to_bits = ao.GetAtomToBits()
    return atom_to_bits


@lru_cache
def get_ligand_cif(ligand_id):
    base_url = "https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2"
    ligand_dir = None
    if ligand_id.startswith("CLC"):
        ligand_cif = os.path.join("clc", ligand_id[-1], ligand_id, f"{ligand_id}.cif")
    elif ligand_id.startswith("PRD"):
        ligand_cif = os.path.join("prd", ligand_id[-1], ligand_id, f"{ligand_id}.cif")
    else:
        ligand_cif = os.path.join("ccd", ligand_id[0], ligand_id, f"{ligand_id}.cif")

    url = f"{base_url}/{ligand_cif}"
    cif = fetch_api_data(url)
    return cif

def get_component(ligand_id):
    ligand_cif = get_ligand_cif(ligand_id)
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp_file:
            tmp_file.write(ligand_cif)
            tmp_filename = tmp_file.name

        component = ccd_reader.read_pdb_cif_file(tmp_filename).component
        return component
    
    except TypeError:
        print(f"Failed to fetch or parse CIF for ligand ID: {ligand_id}")
        return None
    
    finally:
        os.unlink(tmp_filename)


def get_ecifp(ecifp: pa.Table, ligand_ecifp_name, protein_ecifp_name):
    
    N = ecifp.num_rows

    ligand_ecifp = np.empty((N, 32), dtype=np.uint32)
    protein_ecifp = np.empty((N, 32), dtype=np.uint32)

    for i in range(N):
        ligand_ecifp[i] = np.frombuffer(ecifp[ligand_ecifp_name][i], dtype=np.uint8).view(np.uint32)
        protein_ecifp[i] = np.frombuffer(ecifp[protein_ecifp_name][i], dtype=np.uint8).view(np.uint32)
    
    return ligand_ecifp, protein_ecifp


def get_auc(labels, sim):
    fpr, tpr, thresholds = roc_curve(labels, sim)
    roc_auc = auc(fpr, tpr)
    display = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc)
    values = (fpr, tpr, thresholds)
    return (roc_auc, display, values)

def get_preferred_assembly(entry_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/summary/{entry_id}"
    data = fetch_api_data(url)
    if data and entry_id in data:
        for assembly in data[entry_id][0]['assemblies']:
            if assembly['preferred']:
                return assembly['assembly_id']
    else:
        print(f"Could not retrieve preferred assembly for entry {entry_id}")
        return None


def get_data_dir():
    ## Get the data directory
    config = configparser.ConfigParser()
    config_file = Path().cwd().parent / 'conf.ini'
    config.read(config_file)
    return Path(config['DEFAULT']['Data'])