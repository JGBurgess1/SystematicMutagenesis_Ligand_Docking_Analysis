from schrodinger import structure, job #, utility
from schrodinger.structutils import build
from schrodinger.job import queue
from schrodinger.job import launchapi
from schrodinger.forcefield.minimizer import minimize_structure, minimize_substructure
from schrodinger.structutils import measure
from schrodinger.structutils import analyze
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import multiprocessing as mp # doesn't work with schrodinger, use JobDJ

import threading as th # doesnt work with schrodinger - use JobDJ

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Perform mutations, grid generation, and docking for protein-ligand interactions.")
parser.add_argument("-c", "--complex", required=True, help="Path to the protein structure file (mae format)")
parser.add_argument("-l", "--ligand", required=True, help="Path to the ligand structure file (mae format)")
#parser.add_argument("-o", "--output", required=True, help="Output directory for results")
args = parser.parse_args()

"""
# Create the output directory if it doesn't exist
#output_directory = args.output
#if not os.path.exists(output_directory):
 #   os.mkdir(output_directory)
"""

# Load the ligand structure
ligand_structure = structure.StructureReader.read(args.ligand)

# Load the protein structure
complex_structure = structure.StructureReader.read(args.complex)

# Define the amino acids to mutate to (three-letter codes)
amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

# 
counter = 0
num_procs = 0

# Lists to store mutation data
mutations = []
binding_scores = []

# Perform mutations, grid generation, and docking for each mutation
# Initialise the queue object for job submission
jobDJ = queue.JobDJ([("localhost", 24)]) # the key issue that makes the running on multiple procs possible...

for residue in complex_structure.residue:
    if residue.chain != 'B':
        continue
    #if int(residue.resnum) == 14:
    #    continue
    #if int(residue.resnum) > 200:
    #    continue
    original_residue = residue
    original_residue_name = original_residue.pdbres.strip()
    original_chain_id = original_residue.chain
    # Not important to have right now - deals with non standard amino acids - or say Heme??
    if original_residue_name not in amino_acids:
        print(original_residue_name, amino_acids)
        continue

    # For standard amino acids, then loop through each of the possible amino acids at that position:
    for new_residue_name in amino_acids:
        if new_residue_name == original_residue_name:
            continue

        # Make a copy of the protein_structure
        mutated_structure = complex_structure.copy()
        # Mutate the residue for the chosen alternative residue
        minimization_zone_atoms = measure.get_atoms_close_to_subset(mutated_structure, original_residue.atom, 5)
        # Select atoms up to 5A away from the residue to be mutated, for minimization
        build.mutate(mutated_structure, original_residue.atom[1], new_residue_name)
        # minimize the structure, using local minimization around the mutation.
        min_res = minimize_substructure(mutated_structure, minimization_zone_atoms)
        # delete the ligand...
        ligands = analyze.find_ligands(min_res)
        ligand = ligands.pop()
        min_res.deleteAtoms(ligand)

        residue_index = original_residue.resnum
        modified_file_name = f'{original_chain_id}_{original_residue_name}_{residue_index}_{new_residue_name}.mae'
        with structure.StructureWriter(modified_file_name) as writer:
            writer.append(min_res)

# Generate a grid for Glide

        grid_gen_spec = f"JOBNAME   gridgen\nOUTPUTDIR   Output/\nGRID_CENTER   -1.95, -6.80, -13.31\nRECEP_FILE   {modified_file_name}\nGRIDFILE   {modified_file_name[:-4]}_min_grid.zip"
        grid_gen_file_name = f"{modified_file_name[:-4]}_grid_gen.inp"
        with open(grid_gen_file_name, "w") as grid_gen_inp_file:
            grid_gen_inp_file.write(grid_gen_spec)

        grid_gen_job = queue.JobControlJob(["glide", grid_gen_file_name])
        #grid_gen_job.addPrereq(minimize_job)      
        #jobDJ.addJob(grid_gen_job) 

    ####### > glide SP < #######

        glide_SP_inp_spec = f"JOBNAME   glide_sp_dock\nOUTPUTDIR   Output/\nGRIDFILE   {modified_file_name[:-4]}_min_grid.zip\nLIGANDFILE   Q203.mae\nPRECISION   SP"
        glide_SP_inp_file_name = f"{modified_file_name[:-4]}_glide_SP.inp"
        with open(glide_SP_inp_file_name,"w") as glide_SP_inp_file:
            glide_SP_inp_file.write(glide_SP_inp_spec)

        glide_sp_job = queue.JobControlJob(["glide", glide_SP_inp_file_name])
        glide_sp_job.addPrereq(grid_gen_job)
        jobDJ.addJob(glide_sp_job)
    # should run the job for each residue, there should be 19 x 3 jobs to run. 1 CPU at a time. [57 jobs in total]. 
    jobDJ.run()

    jobDJ.wait()
