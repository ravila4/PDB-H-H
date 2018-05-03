from Bio.PDB import *
import glob
import gzip

pdb_dir = "./pdb"

pdb_list = glob.glob(pdb_dir + "/*/*.ent.gz")
parser = PDBParser()
nmr_list = []

# Find nmr structures
for f in pdb_list:
    with gzip.open(f, 'rt') as pdb:
        structure = parser.get_structure("", pdb)
    method = structure.header['structure_method']
    if method == "solution nmr":
        nmr_list.append(f)

# Write to file
with open("nmr_structures.txt", 'w') as outfile:
    for f in nmr_list:    
        outfile.write(f + "\n") 
