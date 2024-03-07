from schrodinger import structure #, utility
from schrodinger.structutils import build
from schrodinger.forcefield.minimizer import minimize_structure, minimize_substructure, MinimizationOptions
from schrodinger.structutils import measure
from schrodinger.structutils import analyze
from schrodinger.job import queue
import argparse


# Parse command-line arguments
parser = argparse.ArgumentParser(description="Perform mutations, grid generation, and docking for protein-ligand interactions.")
parser.add_argument("-c", "--complex", required=True, help="Path to the protein structure file (mae format)")
parser.add_argument("-l", "--ligand", required=True, help="Path to the ligand structure file (mae format)")
args = parser.parse_args()

# Load the ligand structure
ligand_structure = structure.StructureReader.read(args.ligand)

# Load the protein structure
complex_structure = structure.StructureReader.read(args.complex)

# Define the amino acids to mutate to (three-letter codes)
amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

non_standard_aa_conversion = {'ARN':'ARG','ASH':'ASP','GLH':'GLU','LYN':'LYS','HID':'HIS','HIE':'HIS','HIP':'HIS'}

jobDJ = queue.JobDJ([("localhost", 600)], max_retries = 3, max_failures = 10) # the key issue that makes the running on multiple procs possible...

jobCounter = 0

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

        # Make a copy of the protein_structure
        mutated_structure = complex_structure.copy()
        # Mutate the residue for the chosen alternative residue     
        
        # Minimize entire protein
        build.mutate(mutated_structure, original_residue.atom[1], new_residue_name)
        # minimize the structure, using local minimization around the mutation.
        
        ### Minimize structure code, load with jobDJ.
	### test with 19 cores, then increase if possible to 600.
        ### job spec...
	### $SCHRODINGER/macromodel [options] <inputfile>        
        ### mmodel_min_job = queue.JobControlJob(['macromodel', mmod_min_file])
        jobDJ.addJob(mmodel_min_job)  
        jobCounter += 1
        
        # to include: write the files, then minimize, then open each file and delete the ligand
        #retain the ligand for now??

        residue_index = original_residue.resnum
        modified_file_name = f'{original_chain_id}_{original_residue_name}_{residue_index}_{new_residue_name}.mae'
        with structure.StructureWriter(modified_file_name) as writer:
            writer.append(mutated_structure)
        continue
 if jobCounter == 600:
    jobDJ.run()
    jobCounter = 0
