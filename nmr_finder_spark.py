from Bio.PDB import *
import glob
import gzip
import pickle
from pyspark.sql import SparkSession


def method_is_nmr(f):
    """Takes a pdb file and finds if it is an NMR structure."""
    with gzip.open(f, 'rt') as pdb:
        header = parse_pdb_header(pdb)
    return header['structure_method'] == "solution nmr"


def contains_alliphatic_aa(f):
    """Takes a pdb file and finds if it contains alliphatic amino acids"""
    aa_list = ['LEU', 'ILE', 'VAL', 'MET', 'ALA', 'PRO', 'GLY']
    with gzip.open(f, 'rt') as pdb:
        parser = PDBParser()
        structure = parser.get_structure("", pdb)
        residues = [res.get_resname() for res in structure.get_residues()]
    for aa in aa_list:
        if aa in residues:
            return True
    return False


if __name__ == "__main__":
    spark = SparkSession.builder.master("local[256]").appName(
                "nmr").getOrCreate()
    sc = spark.sparkContext
    # Directory where the PDB database is located
    PDB_DIR = "./pdb"
    # List of pdbs inside PDB_DIR
    pdb_list = glob.glob(PDB_DIR + "/*/*.ent.gz")
    # Create a parallelized collection
    pdb_data = sc.parallelize(pdb_list)
    # Filter by method
    pdb_data = pdb_data.filter(lambda x: method_is_nmr(x))
    # Filter by residue types
    pdb_data = pdb_data.filter(lambda x: contains_alliphatic_aa(x))
    # Collect
    nmr_files = pdb_data.collect()
    pickle.dump(nmr_files, open("nmr_files_peptide.pkl", 'wb'))
    spark.stop()
