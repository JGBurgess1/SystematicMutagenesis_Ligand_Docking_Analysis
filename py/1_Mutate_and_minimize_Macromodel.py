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
args = parser.parse_args()

# Load the protein structure
complex_structure = structure.StructureReader.read(args.complex)

# Define the amino acids to mutate to (three-letter codes)
amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

non_standard_aa_conversion = {'ARN':'ARG','ASH':'ASP','GLH':'GLU','LYN':'LYS','HID':'HIS','HIE':'HIS','HIP':'HIS'}

jobDJ = queue.JobDJ([("localhost", 24)], max_retries = 3, max_failures = 10) # the key issue that makes the running on multiple procs possible...

jobCounter = 0

for residue in complex_structure.residue:
    if residue.chain != 'B':
        continue

    original_residue = residue
    original_residue_name = original_residue.pdbres.strip()
    original_chain_id = original_residue.chain
    residue_index = original_residue.resnum   
 
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
        
        build.mutate(mutated_structure, original_residue.atom[1], new_residue_name)
         
        modified_file_name = f'{original_chain_id}_{original_residue_name}_{residue_index}_{new_residue_name}_macromodel.mae'
         
        with structure.StructureWriter(modified_file_name) as writer:
            writer.append(mutated_structure)
       
        new_comfile_name = f'COM_FILE_{modified_file_name[:-4]}.com'
        # read in the com file.
        with open('mmod_mini_7.com', 'r') as comfile:
            comfile_contents = comfile.readlines()
            comfile_contents[0] = modified_file_name
            comfile_contents[1] = f'{modified_file_name[:-4]}_min_out.mae'
            comfile_string = "\n".join(comfile_contents)
            with open(new_comfile_name, 'w') as new_comfile:
                new_comfile.write(comfile_string)

        #Done reading in the old comfile and writing the new
        mmodel_min_job = queue.JobControlJob(['macromodel', new_comfile_name])
        jobDJ.addJob(mmodel_min_job)  
        jobCounter += 1
        ### Minimize structure code, load with jobDJ.
	### test with 19 cores, then increase if possible to 600.
        ### job spec...
	### $SCHRODINGER/macromodel [options] <inputfile>        
        ### mmodel_min_job = queue.JobControlJob(['macromodel', mmod_min_file])
        
        ### need to write the com file for each instance. [ write from a basic template, modify the name of the structure at the top ] 
        ### or mutate the mmod_mini generated file? And then run that through the minimization protocol? ]
        ### need to make changes to the file, so that it will match the mmod_mini_5.mae file
        ### 
        
        # to include: write the files, then minimize, then open each file and delete the ligand
        #retain the ligand for now??
    if jobCounter >= 24:
        jobDJ.run()
        jobCounter = 0
