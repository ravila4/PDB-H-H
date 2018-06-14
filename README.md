# PDB-HH
Scripts for finding C-H···H-C contacts in NMR structures.

## Contents
 -  **nmr_finder.py** Filters a local installation of the protein data bank to find NMR sructures containing peptide chains with aliphatic residues.
 - **h2hbond.py** Finds all C-H··H-C pairs within 3 angstroms of each other, and extracts data about their distance, angles, and dihedral.
 - **nmr_files_peptide.pkl** List of NMR peptide files in the PDB as of May 2018. Input for `h2hbond.py`.
 - **hhcontacts.pkl** Raw contact data as a pandas DataFrame
 - **full_pdb_analysis.ipynb** Exploratory analysis of the data.

## Requirements
 - pyspark: `conda install pyspark`
