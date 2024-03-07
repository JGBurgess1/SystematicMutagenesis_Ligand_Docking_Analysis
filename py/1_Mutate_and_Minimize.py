from schrodinger import structure #, utility
from schrodinger.structutils import build
from schrodinger.forcefield.minimizer import minimize_structure, minimize_substructure
from schrodinger.structutils import measure
from schrodinger.structutils import analyze
import argparse
from schrodinger.job import queue

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Perform mutations, grid generation, and docking for protein-ligand interactions.")
parser.add_argument("-c", "--complex", required=True, help="Path to the protein structure file (mae format)")
parser.add_argument("-l", "--ligand", required=True, help="Path to the ligand structure file (mae format)")
parser.add_argument("-p", "--position", required=True, choices=['0','1','2','3','4','5'], help="Position to use, for parallel running of the script")
args = parser.parse_args()

# Load the ligand structure
ligand_structure = structure.StructureReader.read(args.ligand)

# Load the protein structure
complex_structure = structure.StructureReader.read(args.complex)

position = args.position

# Define the amino acids to mutate to (three-letter codes)
amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

non_standard_aa_conversion = {'ARN':'ARG','ASH':'ASP','GLH':'GLU','LYN':'LYS','HID':'HIS','HIE':'HIS','HIP':'HIS'}

jobDJ = queue.JobDJ([("localhost", 19)], max_retries = 3, max_failures = 10)

for residue in complex_structure.residue:
    if residue.chain != 'B':
        continue

    original_residue = residue
    original_residue_name = original_residue.pdbres.strip()
    original_chain_id = original_residue.chain
    residue_index = original_residue.resnum    
   
    if position == '0':
        pass
    elif position == '1':
        if residue_index < 100:
            continue
    elif position == '2':
        if residue_index < 200:
            continue
    elif position == '3':
        if residue_index < 300:
            continue
    elif position == '4':
        if residue_index < 400:
            continue
    elif position == '5':
        if residue_index < 500:
            continue    
 
    # For non standard amino acid names =the schrodinger alternative protonation state names, for instance, convert these to the normal names.
    if original_residue_name not in amino_acids:
        if original_residue_name in non_standard_aa_conversion.keys():
            original_residue_name = non_standard_aa_conversion[original_residue_name]
        else:
            continue

    # For standard amino acids, then loop through each of the possible amino acids at that position:
    for new_residue_name in amino_acids:
        if new_residue_name == original_residue_name:
            continue

        modified_file_name = f'{original_chain_id}_{original_residue_name}_{residue_index}_{new_residue_name}.mae'
        # 2. Minimize the file, using default parameters, and an implicit membrane. Prime minimization used. 
        minimiz_file_text = f"STRUCT_FILE   {modified_file_name}\nPRIME_TYPE   REAL_MIN\nSELECT   ALL\nUSE_RANDOM_SEED   yes\nUSE_MEMBRANE  yes"

        # INP file gen for minimization (Prime). 

        minimiz_file_name = f"{modified_file_name[:-4]}_minimiz.inp"
        with open(minimiz_file_name, "w") as minimiz_inp_file:
            minimiz_inp_file.write(minimiz_file_text)

        minimize_job = queue.JobControlJob(["prime", minimiz_file_name])
        jobDJ.addJob(minimize_job)
#        """
    ### To add later ####
        ### can also do lig prep of the input ligand for the user ###
        ### ' version 2-> add lig prep, and induced fit docking...
        ### not for now though.
        ###
    jobDJ.run()
    # Wait for all jobs to finish
#jobDJ.run()
