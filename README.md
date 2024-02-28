## Description
This tool generates a pattern-matching fingeprint inspired by Morgan algorithms to extract pertinent substructures from a dataset of molecules, particularly tailored for metabolomic data analysis.

## Dependencies
- Python 3
- pandas
- RDKit
- numpy

## Installation
1. Make sure you have Python 3 installed on your system.
2. Install the required dependencies using pip:

pip install pandas rdkit numpy

## Usage
1. Ensure your dataset is in CSV format with a 'Smiles' column containing the SMILES representations of molecules.
2. Run the script `script.py` from the command line or terminal, passing the path to your dataset file as an argument:

python script.py <path_to_dataframe_file>

Replace `<path_to_dataframe_file>` with the actual path to your dataset file.
3. The script will load the dataset, extract relevant substructures, and perform analysis on their occurrences.
4. The results will be saved in pickle files (`x_mod_uq_freq.pkl`, `x_mod_bin_uq_freq.pkl`, `struct_uq_freq.pkl`) for further analysis.

## Example
python script.py data.csv
