 	
#! /usr/bin/env python

__author__ = "Suman Sirimulla"
__authoraffliation__ = "Department of Chemistry & Biochemsitry, Northern Arizona University"
__authoremail___= "suman.sirimulla@nau.edu"

__co_author__='Erik M. Chavez'
__co_authoraffliation__ = "Undergraduate Researcher,Chemistry, Northern Arizona University"
__co_authoremail___ = "emc253@nau.edu"


import re
import math
import numpy as np
import itertools
import sys
import os
import gzip
import pickle
import glob
#from Tkinter import Tk
#from Tkinter import *
#import ttk
#from tkFileDialog import askdirectory

#Tk().withdraw()


nmr_list_file = open("nmr_files.pkl", 'rb') #open('/Users/grant/clayhd/pdb/nmr_files.txt', 'r')

userpdb='./pdb'
pdbfiles = glob.glob(userpdb + "/*/*gz")#os.listdir(userpdb)
useroutpath='.'

useroutname='h2hbond.txt'

#the desired name and filepath of output file
rawoutfile=open(useroutpath+'/'+useroutname,'w')
rawoutfile.close()
rawoutfile=open(useroutpath+'/'+useroutname,'a')


#outward loop
#repeats searh for each file in pdb folder
print("Checking pdb files...")
header='{:4s} {:6s} {:4s} {:3s} {:5s} {:3s} {:1s} {:6s} {:4s} {:3s} {:5s} {:3s} {:6s}'.format('PDB ','  H1# ',
                                                                                  ' H1 ','Res',' Res#',
                                                                                  'Ty ','|','  H2# ',
                                                                                  ' H2 ','Res',' Res#',
                                                                                  'Ty ','H-HDis')
rawoutfile.write(header+ '\n')

nmr_files = pickle.load(nmr_list_file) #nmr_list_file.read().split('\n')
#print nmr_files[1]
wanted_residues = ['LEU', 'ILE', 'VAL', 'MET']
noncarbon_h = [' H  ']

for fil in pdbfiles:
    InFileName = fil
    
    InFile = gzip.open(InFileName, 'rt')
    if not InFileName[-11::] in nmr_files:
        continue
    print(InFileName[-8:-4])
    #OutFile=open(OutFileName,'a')
    dict_hydrogen1 = { }
    dict_hydro1={}
    dict_carbon1= { }
    dict_conect = { }
    dict_atom = { }
    modres=[]
    
    i = 0
    for line in InFile:
            record_type=line[0:6]
            atom_type=line[76:78]
            resname=line[17:20]
            fullname=line[12:16]
            newLine = ' '.join(line.split()).split(' ')
            if record_type == 'MODRES':
                modres.append(newLine[2])
            if line[21] == 'A' and record_type == 'ATOM  ' and (atom_type == ' H' ) and resname in wanted_residues and not fullname in noncarbon_h:
                dict_hydro1[i] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':fullname, 'altloc':line[16:17], 'resname':resname, 'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                #print dict_hydrogen1[i]
                #print line
                i = i + 1
    f = i
    #print (dict_hydrogen1)
    InFile.close()


    i = 0
    d = 0
    s = 0
    for j in range(0,f):
            InFile = open(InFileName, 'r')
            for line in InFile:
                    record_type=line[0:6]
                    serial_number=line[6:11]
                    atom_type=line[76:78]
                    resname=line[17:20]
                    fullname=line[12:16]
                    if line[21] == 'A' and record_type == 'ATOM  ' and atom_type == ' H' and (dict_hydro1[j]['resname'] not in modres) and not serial_number == dict_hydro1[j]['serial_number'] and resname in wanted_residues and not fullname in noncarbon_h:
                        #print line
                        try:
                            distance = math.sqrt(((dict_hydro1[j]['x'] - float(line[30:38]))**2) + ((dict_hydro1[j]['y'] - float(line[38:46]))**2) + ((dict_hydro1[j]['z'] - float(line[46:54]))**2)) # calculate length between hydrogen atoms in protein
                        except(ZeroDivisionError,ValueError):
                                continue
#                        print distance
                        if distance <= 4.0 and distance >= 2.0: # if length is between 2.0 and 4.0 Angstroms
                            dict_atom[s] = {'record_type':line[0:6], 'serial_number':line[6:11],'fullname':line[12:16], 'altloc':line[16:17], 'resname':line[17:20], 'chainid':line[21:22], 'resseq':int(line[22:26]), 'icode':line[26:27], 'x':float(line[30:38]), 'y':float(line[38:46]), 'z':float(line[46:54]), 'atom_type':line[76:78]}
                                #print distance
                                #print dict_atom[s]
                                

                            rawout='{:4s} {:6d} {:4s} {:3s} {:5d} {:3s} {:1s} {:6d} {:3s} {:3s} {:5d} {:3s} {:6.2f}'.format(InFileName[-8:-4],
                                                                                                              int(dict_hydro1[j]['serial_number']),
                                                                                                              dict_hydro1[j]['fullname'],dict_hydro1[j]['resname'],int(dict_atom[s]['resseq']),
                                                                                                              dict_hydro1[j]['atom_type'],'|',int(dict_atom[s]['serial_number']),
                                                                                                              dict_atom[s]['fullname'],dict_atom[s]['resname'],
                                                                                                              int(dict_atom[s]['resseq']),dict_atom[s]['atom_type'],
                                                                                                              round(distance,2))
                            rawoutfile.write(rawout+'\n')


    InFile.close()


rawoutfile.close() 
print('H2H interaction Analysis Complete')
