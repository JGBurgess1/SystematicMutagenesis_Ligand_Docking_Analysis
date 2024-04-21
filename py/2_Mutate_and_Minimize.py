from schrodinger import structure #, utility # type: ignore
from schrodinger.structutils import build # type: ignore
from schrodinger.forcefield.minimizer import minimize_structure, minimize_substructure, MinimizationOptions # type: ignore
from schrodinger.structutils import measure # type: ignore
from schrodinger.structutils import analyze, transform # type: ignore
from schrodinger.job import queue # type: ignore
import argparse
import time
import numpy as np
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Perform mutations, grid generation, and docking for protein-ligand interactions.")
parser.add_argument("-c", "--complex", required=True, default='prime_mmgbsa_test_Nov_16_1-out.mae', help="Path to the protein structure file (mae format)")
parser.add_argument("-l", "--ligand", required=True, default='Q203.mae', help="Path to the ligand structure file (mae format)")
parser.add_argument("-r", "--rank", required=True, help="The rank of the mpi process")
args = parser.parse_args()

# Load the ligand structure
ligand_structure = structure.StructureReader.read(args.ligand)

# Load the protein structure
complex_structure = structure.StructureReader.read(args.complex)

# Define the amino acids to mutate to (three-letter codes)
amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

non_standard_aa_conversion = {'ARN':'ARG','ASH':'ASP','GLH':'GLU','LYN':'LYS','HID':'HIS','HIE':'HIS','HIP':'HIS'}

minimization_options = MinimizationOptions(opls_version=16, nonbond_cutoff=14.0, max_step=1500, energy_convergence=5e-09, gradient_convergence=0.05, line_search_method=1, min_method=0, no_restrain_zob=False, bend_conj_amines=False, perturb=False, restraints=None, constraints=None, debug_outfile=None)

job_dj = queue.JobDJ()

rank = int(args.rank)

def generate_mutations():
    mutation_file = open("mutation_list.txt", "a")
    for residue in complex_structure.residue:
        original_residue = residue
        original_residue_name = original_residue.pdbres.strip()
        original_chain_id = original_residue.chain
        
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
            mutation_file.write(f'{original_residue.atom[1]},{new_residue_name},{original_residue.chain}_{original_residue.resnum}\n')
            # make a list of all amino acid mutations to be made:
    mutation_file.close()

# want to see that you can mutate using the info provided by the text file.
def mutate_and_minimize():
    mutation_arr = []
    with open("mutation_list.txt", "r") as mutation_file:
        lines = mutation_file.readlines()
        for line in lines:
            mutation_arr.append(line.split(","))
    
    
    size = 240
    num_mutations = len(mutation_arr)
    chunk_size = num_mutations//size
    remainder = num_mutations%size

    start_index = rank*chunk_size + min(rank, remainder)
    end_index = start_index + chunk_size + (1 if rank < remainder else 0)

    for new_mutation in mutation_arr[start_index:end_index]:
        atom_num = new_mutation[0][new_mutation[0].find('(')+1:new_mutation[0].find(')')]
        atom_num = int(atom_num)

        new_file_name = f'{new_mutation[2].strip()}_{new_mutation[1].strip()}.pdb'
        if os.path.exists(new_file_name):
            continue
        
        mutated_structure = complex_structure.copy()
        build.mutate(mutated_structure, atom_num, new_mutation[1].strip())
        minimize_structure(mutated_structure, minimization_options)
        ligand_in_structure = analyze.find_ligands(mutated_structure)
        centroid = transform.get_centroid(ligand_in_structure[0].st)
        centroid = centroid[0:3]

        ligand_atoms = analyze.evaluate_asl(mutated_structure, "ligand")
        mutated_structure.deleteAtoms(ligand_atoms)

        print(f'{rank} RANK, Mutation and Minimization complete. {new_mutation}')
        
        with structure.StructureWriter(new_file_name) as writer:
                writer.append(mutated_structure)
        
"""
        # generate a grid for docking
        grid_gen_spec = f"JOBNAME   gridgen\nOUTPUTDIR   April_Output/\nGRID_CENTER   {centroid[0]}, {centroid[1]}, {centroid[2]}\nRECEP_FILE   {new_file_name}\nGRIDFILE   {new_file_name[:-4]}_grid.zip"
        grid_gen_file_name = f"{new_file_name[:-4]}_grid_gen.inp"
        with open(grid_gen_file_name, "w") as grid_gen_inp_file:
            grid_gen_inp_file.write(grid_gen_spec)
        grid_gen_job = queue.JobControlJob(["glide", grid_gen_file_name])
        job_dj.addJob(grid_gen_job)
        job_dj.run()

        # dock the combined ligands to the mutated structure
        glide_XP_inp_spec = f"JOBNAME   glide_xp_dock\nOUTPUTDIR   April_Output/\nGRIDFILE   {new_file_name[:-4]}_grid.zip\nLIGANDFILE   Combined_Ligands.mae\nPRECISION   XP"
        glide_XP_inp_file_name = f"{new_file_name[:-4]}_glide_XP.inp"
        with open(glide_XP_inp_file_name,"w") as glide_XP_inp_file:
            glide_XP_inp_file.write(glide_XP_inp_spec)
        
        glide_xp_job = queue.JobControlJob(["glide", glide_XP_inp_file_name])
        job_dj.addJob(glide_xp_job)
        job_dj.run()
        # summarize results at end.
        # check if worked?!
"""
if __name__ == "__main__":
    if rank == -1:
        generate_mutations()
    else:
        mutate_and_minimize()
