from schrodinger import structure #, utility
from schrodinger.structutils import build
from schrodinger.forcefield.minimizer import minimize_structure, minimize_substructure, MinimizationOptions
from schrodinger.structutils import measure
from schrodinger.structutils import analyze
import argparse
import time
import multiprocessing as mp # doesn't work with schrodinger, use JobDJ
import numpy as np
from contextlib import closing
import threading as th # doesnt work with schrodinger - use JobDJ

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

rank = int(args.rank)
size = 600
num_mutations = 9956

chunk_size = 9956//600 # = 16
remainder = 9956%600 # = 556

start_index = rank*chunk_size + min(rank, remainder)
end_index = start_index + chunk_size + (1 if rank < remainder else 0)

# Define the amino acids to mutate to (three-letter codes)
amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

non_standard_aa_conversion = {'ARN':'ARG','ASH':'ASP','GLH':'GLU','LYN':'LYS','HID':'HIS','HIE':'HIS','HIP':'HIS'}

minimization_options = MinimizationOptions(opls_version=16, nonbond_cutoff=14.0, max_step=1500, energy_convergence=5e-09, gradient_convergence=0.05, line_search_method=1, min_method=0, no_restrain_zob=False, bend_conj_amines=False, perturb=False, restraints=None, constraints=None, debug_outfile=None)

def generate_mutations():
    mutation_file = open("mutation_list.txt", "a")
    for residue in complex_structure.residue:
        if residue.chain != 'B':
            continue

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
            mutation_file.write(f'{original_residue.atom[1]},{new_residue_name},{original_residue}\n')
            # make a list of all amino acid mutations to be made:
    mutation_file.close()

# want to see that you can mutate using the info provided by the text file.
def mutate_and_minimize():
    mutation_arr = []
    with open("mutation_list.txt", "r") as mutation_file:
        lines = mutation_file.readlines()
        for line in lines:
            mutation_arr.append(line.split(","))
    for new_mutation in mutation_arr[start_index:end_index]:
        atom_num = new_mutation[0][new_mutation[0].find('(')+1:new_mutation[0].find(')')]
        atom_num = int(atom_num)
        mutated_structure = complex_structure.copy()
        build.mutate(mutated_structure, atom_num, new_mutation[1].strip())
        minimize_structure(mutated_structure, minimization_options)
        print(f'{rank} RANK, Mutation and Minimization complete. {new_mutation}')
        with structure.StructureWriter(f'{new_mutation[2]}_{new_mutation[1].strip()}.mae') as writer:
                writer.append(mutated_structure)
        pass

if __name__ == "__main__":
    if rank == -1:
        generate_mutations()
    else:
        mutate_and_minimize()
