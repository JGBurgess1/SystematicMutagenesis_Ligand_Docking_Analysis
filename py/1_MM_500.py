from schrodinger import structure #, utility
from schrodinger.structutils import build
from schrodinger.forcefield.minimizer import minimize_structure, minimize_substructure
from schrodinger.structutils import measure
from schrodinger.structutils import analyze
import argparse
import multiprocessing as mp # doesn't work with schrodinger, use JobDJ

import threading as th # doesnt work with schrodinger - use JobDJ

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Perform mutations, grid generation, and docking for protein-ligand interactions.")
parser.add_argument("-c", "--complex", required=True, help="Path to the protein structure file (mae format)")
parser.add_argument("-l", "--ligand", required=True, help="Path to the ligand structure file (mae format)")
#parser.add_argument("-o", "--output", required=True, help="Output directory for results")
args = parser.parse_args()

# Load the ligand structure
ligand_structure = structure.StructureReader.read(args.ligand)

# Load the protein structure
complex_structure = structure.StructureReader.read(args.complex)

# Define the amino acids to mutate to (three-letter codes)
amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

non_standard_aa_conversion = {'ARN':'ARG','ASH':'ASP','GLH':'GLU','LYN':'LYS','HID':'HIS','HIE':'HIS','HIP':'HIS'}

for residue in complex_structure.residue:
    if residue.chain != 'B':
        continue

    original_residue = residue
    original_residue_name = original_residue.pdbres.strip()
    original_chain_id = original_residue.chain
    residue_index = original_residue.resnum

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

        # Make a copy of the protein_structure
        mutated_structure = complex_structure.copy()
        # Mutate the residue for the chosen alternative residue
        minimization_zone_atoms = measure.get_atoms_close_to_subset(mutated_structure, original_residue.atom, 5)
        # Select atoms up to 5A away from the residue to be mutated, for minimization
        build.mutate(mutated_structure, original_residue.atom[1], new_residue_name)
        # minimize the structure, using local minimization around the mutation.
        min_res = minimize_substructure(mutated_structure, minimization_zone_atoms)
        # delete the ligand...
        ligand_atoms = analyze.evaluate_asl(min_res, "ligand")
        min_res.deleteAtoms(ligand_atoms)

        modified_file_name = f'{original_chain_id}_{original_residue_name}_{residue_index}_{new_residue_name}.mae'
        with structure.StructureWriter(modified_file_name) as writer:
            writer.append(min_res)
        continue
        #"""
        ## 2. Minimize the file, using default parameters, and an implicit membrane. Prime minimization used. 
        #minimiz_file_text = f"STRUCT_FILE {modified_file_name}\nPRIME_TYPE  REAL_MIN\nSELECT  asl = all\nUSE_CRYSTAL_SYMMETRY  no\nUSE_RANDOM_SEED yes\nSEED  0\nEXT_DIEL  80.00\nUSE_MEMBRANE  yes"

        # INP file gen for minimization (Prime). 

        #minimiz_file_name = f"{modified_file_name[:-4]}_minimiz.inp"
        #with open(minimiz_file_name, "w") as minimiz_inp_file:
        #    minimiz_inp_file.write(minimiz_file_text)

        #minimize_job = queue.JobControlJob(["prime", minimiz_file_name])
#        """
    ### To add later ####
        ### can also do lig prep of the input ligand for the user ###
        ### ' version 2-> add lig prep, and induced fit docking...
        ### not for now though.
        ###

    # Wait for all jobs to finish
#jobDJ.run()
