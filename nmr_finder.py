from Bio.PDB import *
import glob
import gzip
import pickle
from pyspark import SparkContext
from pyspark.sql import SparkSession

def find_nmr(f):
    """ Takes a pdb file and finds if it is an NMR structure."""
    print(f)
    with gzip.open(f, 'rt') as pdb:
        structure = parser.get_structure("", pdb)
    method = structure.header['structure_method']
    return (method == "solution nmr")

if __name__ == "__main__":
    spark = SparkSession.builder.master("local[256]").appName("5-MapReduce").getOrCreate()
    sc = spark.sparkContext 
    parser = PDBParser()
    # Directory where the PDB database is located
    PDB_DIR = "./pdb"
    # List of pdbs inside PDB_DIR
    pdb_list = glob.glob(PDB_DIR + "/*/*.ent.gz")
    # Create a parallelized collection
    pdb = sc.parallelize(pdb_list)
    nmr_files = pdb.filter(lambda x: find_nmr(x)).collect()
    pickle.dump(nmr_files, open("nmr_files.pkl", 'wb'))
