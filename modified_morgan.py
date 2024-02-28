import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools
import pickle

# Define function to check if two molecules are equivalent
def eq_mol(a, b):
    """Check if two molecules are equivalent."""
    return a.HasSubstructMatch(b) and b.HasSubstructMatch(a)

# Define function to get substructures of radius 2 from a list of molecules
def get_submol_morgan2(mols, r=2, chiral=True):
    """
    Get a list of relevant substructures of radius 2 from a list of molecules.
    
    Parameters:
        mols (list): List of RDKit molecule objects.
        r (int): Radius of substructure search (default is 2).
        chiral (bool): Whether to consider chirality (default is True).
    
    Returns:
        list: List of relevant substructures.
    """
    sublist = []
    for mol in mols:
        bit_info = {}
        fp = AllChem.GetHashedMorganFingerprint(mol, 2, nBits=1024, bitInfo=bit_info, useChirality=chiral)
        for key in bit_info:
            atomidx, radius = bit_info[key][0]
            if radius == r:
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atomidx)
                submol = Chem.PathToSubmol(mol, env)
                if all(eq_mol(submol, s) == False for s in sublist):
                    sublist.append(submol)
    return sublist

# Define function to get modified substructures of radius 2 for a single molecule
def get_mod_morgan2(mol, sublist, binary=False, chiral=True):
    """
    Get modified substructures of radius 2 for a single molecule.
    
    Parameters:
        mol (RDKit Mol object): RDKit molecule object.
        sublist (list): List of relevant substructures.
        binary (bool): Whether to use binary representation (default is False).
        chiral (bool): Whether to consider chirality (default is True).
    
    Returns:
        np.array: Array of modified substructures.
    """
    res = np.zeros(len(sublist))
    for i, submol in enumerate(sublist):
        res[i] = len(mol.GetSubstructMatches(submol, useChirality=chiral))
    if binary:
        return np.where(res >= 1, 1, 0)
    else:
        return res.astype(int)

# Main code

if __name__ == "__main__":
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_dataframe_file>")
        sys.exit(1)

    # Get the path of the dataframe file from command-line argument
    df_path = sys.argv[1]

    try:
        # Read the dataframe from CSV file
        df = pd.read_csv(df_path)

        # Ensure the dataframe has a 'Smiles' column
        if 'Smiles' not in df.columns:
            print("Error: DataFrame must have a 'Smiles' column.")
            sys.exit(1)
        else:
            print("DataFrame loaded")

        # Add structure column to dataframe
        PandasTools.AddMoleculeColumnToFrame(df, "Smiles")

        # Retrieve substructures of radius 1 and 2
        sublist2 = get_submol_morgan2(df['ROMol'], r=2, chiral=True)
        sublist1 = get_submol_morgan2(df['ROMol'], r=1, chiral=True)
        sublist12 = sublist1 + sublist2

        # Make count fingerprint
        x_mod = np.array([np.array(get_mod_morgan2(m, sublist12, chiral=True)) for m in df['ROMol']])

        # Make binary fingerprint
        x_mod_bin = (x_mod > 0).astype(np.int_)

        # Get unique fingerprints
        x_mod_uq, count_indices, _, _ = np.unique(x_mod, axis=1, return_index=True, return_inverse=True, return_counts=True)
        x_mod_bin_uq = x_mod_bin[:, count_indices]
        sublist12_uq = np.array(sublist12)[count_indices]

        # Remove infrequent fingerprints
        min_freq = 10
        selected = x_mod_bin_uq.sum(axis=0) >= min_freq
        x_mod_uq_freq = x_mod_uq[:, selected]
        x_mod_bin_uq_freq = x_mod_bin_uq[:, selected]
        sublist12_freq = sublist12_uq[selected]

        # Save fingerprint and substructures
        with open("x_mod_uq_freq.pkl", "wb") as handle:
            pickle.dump(x_mod_uq_freq, handle)
        with open("x_mod_bin_uq_freq.pkl", "wb") as handle:
            pickle.dump(x_mod_bin_uq_freq, handle)
        with open("struct_uq_freq.pkl", "wb") as handle:
            pickle.dump(sublist12, handle)

    except Exception as e:
        print("An error occurred:", e)