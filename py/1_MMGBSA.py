from schrodinger import structure #, utility
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

jobDJ = queue.JobDJ([("localhost", 19)], max_retries=3, max_failures=10) # the key issue that makes the running on multiple procs possible...

for residue in complex_structure.residue:
    if residue.chain != 'B':
        continue

    original_residue = residue
    original_residue_name = original_residue.pdbres.strip()
    original_chain_id = original_residue.chain
    residue_index = original_residue.resnum
    # Not important to have right now - deals with non standard amino acids - or say Heme??
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
        #"""
        #    
        glide_PV_file = f"Output/{modified_file_name[:-4]}_glide_XP_pv.maegz"
        ####### > prime_mmgbsa < #######

# No flexibility of active site allowed - calculate the energies from rigid receptor.

        prime_mmgbsa_job = queue.JobControlJob(["prime_mmgbsa", glide_PV_file, "-membrane", "-job_type", "ENERGY"])
        jobDJ.addJob(prime_mmgbsa_job)         
        
    jobDJ.run()

