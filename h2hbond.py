from Bio.PDB import *
import gzip
import numpy as np
import os
import pandas as pd
import pickle
from pyspark.sql import SparkSession


def h2hanalyze(f):
    # Parse the structure
    parser = PDBParser()
    pdb = gzip.open(f, 'rt')
    pdb_id = os.path.split(f)[1][3:7]
    structure = parser.get_structure(pdb_id, pdb)

    # Get the Chain A from the first model
    chain_A = structure[0]['A']

    # Filter by residues
    aliphatic_res = ["LEU", "ILE", "VAL", "MET", "PRO"]
    residues = (res for res in chain_A if res.get_resname() in aliphatic_res)

    # Get all alliphatic hydrogens
    hydrogens = []
    h_dict = {"LEU": ["HB2", "HB3", "HG", "HD11", "HD12", "HD13", "HD21",
                      "HD22", "HD23"],
              "ILE": ["HB", "HG12", "HG13", "HG21", "HG22", "HG23", "HD11",
                      "HD12", "HD13"],
              "VAL": ["HB", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23"],
              "MET": ["HB2", "HB3", "HG2", "HG3", "HE1", "HE2", "HE3"],
              "PRO": ["HB2", "HB3", "HG2", "HG3"]}
    for res in residues:
        name = res.get_resname()
        # Avoid hetero-residues
        hetero = res.id[0]
        if hetero != " ":
            pass
        for h in h_dict[name]:
            try:
                hydrogens.append(res[h])
            except:
                pass

    # Get pairwise distances and angles for Hydrogen contacts
    hh_contacts = []
    for i in range(len(hydrogens)):
        h1 = hydrogens[i]
        res1 = h1.parent
        for j in range(i + 1, len(hydrogens)):
            h2 = hydrogens[j]
            res2 = h2.parent
            # Check that the two Hs are in different residues
            if res1 is not res2:
                # Get distance
                dist = h1 - h2
                if dist <= 3:
                    # Get IDs
                    res1name = res1.resname
                    res1num = res1.id[1]
                    res2name = res2.resname
                    res2num = res2.id[1]
                    h1name = h1.name
                    h2name = h2.name
                    # Get the parent carbons
                    if len(h1name) <= 3:
                        c1name = 'C' + h1name[1]
                    else:
                        c1name = 'C' + h1name[1:-1]
                    if len(h2name) <= 3:
                        c2name = 'C' + h2name[1]
                    else:
                        c2name = 'C' + h2name[1:-1]
                    c1 = res1[c1name]
                    c2 = res2[c2name]
                    # Find the angles
                    c1vec = c1.get_vector()
                    h1vec = h1.get_vector()
                    h2vec = h2.get_vector()
                    c2vec = c2.get_vector()
                    a1 = np.degrees(calc_angle(c1vec, h1vec, h2vec))
                    a2 = np.degrees(calc_angle(c2vec, h2vec, h1vec))
                    dihed = np.degrees(calc_dihedral(c1vec, h1vec,
                                                     h2vec, c2vec))
                    # Save data
                    hh_contacts.append(
                        {"PDBid": pdb_id, "res1name": res1name,
                         "res1num": res1num, "H1": h1name,
                         "res2name": res2name, "res2num": res2num,
                         "H2": h2name, "distance": dist,
                         "angle1": a1, "angle2": a2, "dihedral": dihed})
    return hh_contacts


if __name__ == "__main__":
    spark = SparkSession.builder.master("local[256]").appName(
                "hhbond").getOrCreate()
    sc = spark.sparkContext
    # Read file list
    pdbs = pickle.load(open("nmr_files_peptide.pkl", 'rb'))
    nmr_files = sc.parallelize(pdbs)
    # Find C-H...H-C contacts
    contacts = nmr_files.flatMap(lambda x: h2hanalyze(x)).collect()
    contacts_df = pd.DataFrame(contacts)
    # Save data
    pickle.dump(contacts_df, open("hhcontacts.pkl", 'wb'))
