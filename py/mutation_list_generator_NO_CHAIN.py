from schrodinger import structure, job #, utility
from schrodinger.structutils import build
from schrodinger.job import queue
from schrodinger.job import launchapi
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import multiprocessing as mp # 

import threading as th

# Load the protein structure
protein_structure = structure.StructureReader.read("qcrAB.mae")

# Amino acid dictionary
aas = {"ALA":'A', "ARG":'R', "ASN":'N', "ASP":'D', "CYS":'C', "GLU":'E', "GLN":'Q', "GLY":'G', "HIS":'H', "ILE":'I', "LEU":'L', "LYS":'K', "MET":'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

fhandle = open("mutation_list_qcrB.txt","a")

for residue in protein_structure.residue:
    if residue.chain != 'B':
        continue
    original_residue_name = residue.pdbres.strip()
    if original_residue_name == 'HIE':
        resname = 'H'
    elif original_residue_name == 'ASH':
        resname = 'D'
    elif original_residue_name == 'GLH':
        resname = 'E'
    elif original_residue_name == 'HIP':
        resname = 'H'
    else:
        resname = aas[original_residue_name]
    resnum = residue.resnum
    for new_residue_name in aas.values():
        if new_residue_name == resname:
            continue
        fhandle.write(f'{resname}{resnum}{new_residue_name}\n')
fhandle.close()

