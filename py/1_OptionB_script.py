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

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Perform mutations, grid generation, and docking for protein-ligand interactions.")
parser.add_argument("-p", "--protein", required=True, help="Path to the protein structure file (mae format)")
parser.add_argument("-l", "--ligand", required=True, help="Path to the ligand structure file (mae format)")
parser.add_argument("-o", "--output", required=True, help="Output directory for results")
args = parser.parse_args()

# Create the output directory if it doesn't exist
output_directory = args.output
if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Load the ligand structure
ligand_structure = structure.StructureReader.read(args.ligand)

# Load the protein structure
protein_structure = structure.StructureReader.read(args.protein)

#Job builder section
"""job_builder = launchapi.JobSpecification()
job_builder.setInputFile(__file__)
job_builder.setInputFile(args.protein)
job_builder.setInputFile(args.ligand)
job_builder.setOutputFile(args.output + "/Output.mae")
job.builder.getJobSpec()
"""
# Define the amino acids to mutate to (single-letter codes)
#amino_acids = "ACDEFGHIKLMNPQRSTVWY"

# Define the amino acids to mutate to (three-letter codes)
amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

# Lists to store mutation data
mutations = []
binding_scores = []

# Perform mutations, grid generation, and docking for each mutation
# Initialise the queue object for job submission
jobDJ = queue.JobDJ([("localhost", 24)]) # the key issue that makes the running on multiple procs possible...

for residue in protein_structure.residue:
    if residue.chain != 'B':
        continue
    #if int(residue.resnum) == 14:
        continue
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
        mutated_structure = protein_structure.copy()
        # Mutate the residue for the chosen alternative residue
        build.mutate(mutated_structure, original_residue.atom[1], new_residue_name)

        residue_index = original_residue.resnum
        modified_file_name = f'{original_chain_id}_{original_residue_name}_{residue_index}_{new_residue_name}.mae'
        with structure.StructureWriter(modified_file_name) as writer:
            writer.append(mutated_structure)


        # 2. Minimize the file, using default parameters, and an implicit membrane. Prime minimization used. 
        minimiz_file_text = f"STRUCT_FILE {modified_file_name}\nPRIME_TYPE  REAL_MIN\nSELECT  asl = all\nUSE_CRYSTAL_SYMMETRY  no\nUSE_RANDOM_SEED yes\nSEED  0\nEXT_DIEL  80.00\nUSE_MEMBRANE  yes"

        # INP file gen for minimization (Prime). 

        minimiz_file_name = f"{modified_file_name[:-4]}_minimiz.inp"
        with open(minimiz_file_name, "w") as minimiz_inp_file:
            minimiz_inp_file.write(minimiz_file_text)

        minimize_job = queue.JobControlJob(["prime", minimiz_file_name])

        # MMGBSA job
        mmgbsa_inp_spec = f""
        mmgbsa_inp_file_name = f""
        with open(mmgbsa_inp_file_name, "w") as mmgbsa_inp_file:
            mmgbsa_inp_file.write(mmgbsa_inp_spec)
        mmgbsa_job = queue.JobControl(["prime_mmgbsa", mmgbsa_inp_file_name])
        mmgbsa_job.addPrereq(minimize_job)
        jobDJ.addJob(mmgbsa_job)


        # Generate a grid for Glide
        """
        grid_gen_spec = f"JOBNAME   gridgen\nOUTPUTDIR   Output/\nGRID_CENTER   -1.45, 7.13, -12.89\nRECEP_FILE   {modified_file_name[:-4]}_minimiz-out.maegz\nGRIDFILE   {modified_file_name[:-4]}_min_grid.zip"
        grid_gen_file_name = f"{modified_file_name[:-4]}_grid_gen.inp"
        with open(grid_gen_file_name, "w") as grid_gen_inp_file:
            grid_gen_inp_file.write(grid_gen_spec)
        
        grid_gen_job = queue.JobControlJob(["glide", grid_gen_file_name])
        grid_gen_job.addPrereq(minimize_job)      
        jobDJ.addJob(grid_gen_job) 

    ####### > glide XP < #######
        
        glide_XP_inp_spec = f"JOBNAME   glide_xp_dock\nOUTPUTDIR   Output/\nGRIDFILE   {modified_file_name[:-4]}_min_grid.zip\nLIGANDFILE   Q203.mae\nPRECISION   XP"
        glide_XP_inp_file_name = f"{modified_file_name[:-4]}_glide_XP.inp"
        with open(glide_XP_inp_file_name,"w") as glide_XP_inp_file:
            glide_XP_inp_file.write(glide_XP_inp_spec)
        
        glide_xp_job = queue.JobControlJob(["glide", glide_XP_inp_file_name])
        glide_xp_job.addPrereq(grid_gen_job)
        jobDJ.addJob(glide_xp_job)
        """
    # should run the job for each residue, there should be 19 x 3 jobs to run. 1 CPU at a time. [57 jobs in total]. 
    #jobDJ.run()
    
    # jobDJ.wait()
    # Wait for all jobs to finish
jobDJ.run()


"""
# Analyze docking results and calculate binding energy for each mutation
for mutation in mutations:
    dock_output_file = os.path.join(output_directory, f"{mutation}_docking_output.sdf")
    docking_scores = utility.analyze_docking(dock_output_file)
    binding_energy = docking_scores[0]["Glide GScore"]
    binding_scores.append(binding_energy)

# Create a dataframe to store mutation data
data = pd.DataFrame({"Mutation": mutations, "Binding Energy": binding_scores})

# Plot the binding energies
plt.figure(figsize=(12, 6))
plt.bar(data["Mutation"], data["Binding Energy"])
plt.xlabel("Mutation")
plt.ylabel("Binding Energy (Glide GScore)")
plt.xticks(rotation=90)
plt.title("Binding Energy of Mutated Proteins")
plt.tight_layout()

# Save the plot to a file
plt.savefig(os.path.join(output_directory, "binding_energy_plot.png"))

# Optionally, display the plot
plt.show()
"""
