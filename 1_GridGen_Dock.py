from schrodinger import structure#, utility
from schrodinger.job import queue
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Perform mutations, grid generation, and docking for protein-ligand interactions.")
parser.add_argument("-c", "--complex", required=True, help="Path to the protein structure file (mae format)")
parser.add_argument("-l", "--ligand", required=True, help="Path to the ligand structure file (mae format)")
parser.add_argument("-pos","--position", required=True, choices = ['0', '1', '2', '3', '4', '5', '6', '7', '8'], help="position within the queue. Allowed values 0, 1, 2, 3, 4, 5, 6, 7, 8")
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

non_standard_aa_conversion = {'ARN':'ARG','ASH':'ASP','GLH':'GLU','LYN':'LYS','HID':'HIS','HIE':'HIS','HIP':'HIS'}

# Perform mutations, grid generation, and docking for each mutation
# Initialise the queue object for job submission
jobDJ = queue.JobDJ([("localhost", 19)], max_retries = 3, max_failures = 10) # the key issue that makes the running on multiple procs possible...

for residue in complex_structure.residue:
    if residue.chain != 'B':
        continue

    # Want to put in different ligands, right? Then to save the file name as part of the ligand.

    if int(residue.resnum) != 222:
        continue 

    original_residue = residue
    original_residue_name = original_residue.pdbres.strip()
    original_chain_id = original_residue.chain
    
    # Not important to have right now - deals with non standard amino acids - or say Heme??

    if original_residue_name not in amino_acids:
        if original_residue_name in non_standard_aa_conversion.keys():
            original_residue_name = non_standard_aa_conversion[original_residue_name]
        else:
            continue
        # This should akip any not covered by the non_standard_aa_conversion dictionary.

    # For standard amino acids, then loop through each of the possible amino acids at that position:
    for new_residue_name in amino_acids:
        if new_residue_name == original_residue_name:
            continue

        # Make a copy of the protein_structure
        # Mutate the residue for the chosen alternative residue
        # Select atoms up to 5A away from the residue to be mutated, for minimization
        residue_index = original_residue.resnum
        modified_file_name = f'{original_chain_id}_{original_residue_name}_{residue_index}_{new_residue_name}.mae'
        #"""
        

        grid_gen_spec = f"JOBNAME   gridgen\nOUTPUTDIR   Output/\nGRID_CENTER   -1.95, -6.80, -13.31\nRECEP_FILE   {modified_file_name}\nGRIDFILE   {modified_file_name[:-4]}_min_grid.zip"
        grid_gen_file_name = f"{modified_file_name[:-4]}_grid_gen.inp"
        with open(grid_gen_file_name, "w") as grid_gen_inp_file:
            grid_gen_inp_file.write(grid_gen_spec)
        
        grid_gen_job = queue.JobControlJob(["glide", grid_gen_file_name])

    ####### > glide SP < #######
        # Three different ligands tried, 

        glide_XP_inp_spec = f"JOBNAME   glide_sp_dock\nOUTPUTDIR   Output/\nGRIDFILE   {modified_file_name[:-4]}_min_grid.zip\nLIGANDFILE   Combined_Ligands.mae\nPRECISION   XP\nWRITE_XP_DESC   True"
        glide_XP_inp_file_name = f"{modified_file_name[:-4]}_glide_XP.inp"
        with open(glide_XP_inp_file_name,"w") as glide_XP_inp_file:
            glide_XP_inp_file.write(glide_XP_inp_spec)
        
        glide_Xp_job = queue.JobControlJob(["glide", glide_XP_inp_file_name])
        glide_Xp_job.addPrereq(grid_gen_job)
        jobDJ.addJob(glide_Xp_job) 
    jobDJ.run() #Will run as 19... i.e. per residue, all possible amino acids, at once.

# chat gpt gave some plot analysis - useful pointers, but likely full of errors that need troubleshooting.
