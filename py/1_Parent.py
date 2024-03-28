from schrodinger import structure #,  utility
from schrodinger.job import queue
from schrodinger.structutils import build, analyze, transform
from schrodinger.forcefield.minimizer import minimize_structure, minimize_substructure, MinimizationOptions


import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Perform mutations, grid generation, and docking for protein-ligand interactions.")
parser.add_argument("-c", "--complex", required=True, help="Path to the protein structure file (mae format)")
parser.add_argument("-l", "--ligand", required=True, help="Path to the ligand structure file (mae format)")
# parser.add_argument("-pos","--position", required=True, choices = ['0', '1', '2', '3', '4', '5', '6', '7', '8'], help="position within the queue. Allowed values 0, 1, 2, 3, 4, 5, 6, 7, 8")
# parser.add_argument("--chain", required=False, default='ALL', help="Specify chains to mutate, or omit for ALL (default)")
parser.add_argument("--cpus", required=True, type=int, help="Enter the number of CPUs available")

args = parser.parse_args()

# Load the protein structure
complex_structure = structure.StructureReader.read(args.complex)

# Load the ligand structure
ligand_structure = structure.StructureReader.read(args.ligand)

# Define the amino acids to mutate to (three-letter codes)
amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
non_standard_aa_conversion = {'ARN':'ARG','ASH':'ASP','GLH':'GLU','LYN':'LYS','HID':'HIS','HIE':'HIS','HIP':'HIS'}

jobDJ = queue.JobDJ([("localhost", int(args.cpus))], max_retries = 3, max_failures = 10) # the key issue that makes the running on multiple procs possible...

job_DJ_counter = 0

Error_Message = ''''''


def Mutate(pos):
    # chain_choice = args.chain.trim().toUpper()
    # pos = int(args.position)
    # chain_number = 0
    #for chain in complex_structure.chain:
    #    print(chain, file=sys.stderr)
    #    chain_number+=1
    #if chain_number>20:
    #    Error_Message+="Too many chains {>20} please use fewer chains in your input structure"
    #    # some exit code here
    #
    #    exit()
    for residue in complex_structure.residue:
        #handle long scripts here.
        original_residue = residue
        original_residue_name = original_residue.pdbres.strip()
        original_chain_id = original_residue.chain
        residue_index = original_residue.resnum
        
        if original_residue_name not in amino_acids:
            if original_residue_name in non_standard_aa_conversion.keys():
                original_residue_name = non_standard_aa_conversion[original_residue_name]
            else:
                continue        
    
        for new_residue_name in amino_acids:
            if new_residue_name == original_residue_name:
                continue
            # each residue here.
            mutated_structure = complex_structure.copy()
            # at each iteration, 
            mutated_file_name = f"{original_chain_id}_{original_residue_name}_{residue_index}_{new_residue_name}.mae"
            build.mutate(mutated_structure, original_residue.atom[1], new_residue_name)
            with structure.StructureWriter(mutated_file_name) as writer:
                writer.append(mutated_structure)
            Minimize_prime(mutated_file_name, mutated_structure)
            #Minimize_macromodel(mutated_file_name, mutated_structure)

def Minimize_prime(mut_file_name, mut_structure):
    minimiz_file_text = f"STRUCT_FILE {mut_file_name}\nPRIME_TYPE  REAL_MIN\nSELECT  asl = all\nUSE_CRYSTAL_SYMMETRY  no\nUSE_RANDOM_SEED yes\nSEED  0\nEXT_DIEL  80.00\nUSE_MEMBRANE  yes"
    minimiz_file_name = f"{mut_file_name[:-4]}_minimiz.inp"
    with open(minimiz_file_name, "w") as minimiz_inp_file:
        minimiz_inp_file.write(minimiz_file_text)
    minimize_job = queue.JobControlJob(["prime", minimiz_file_name])
    jobDJ.addJob(minimize_job)
    #Remember to clean up the unminimized file!?
    #And the minimize_inp_file
    job_DJ_counter += 1
    if job_DJ_counter == 600:
        job_DJ_counter = 0
        jobDJ.run()
    
def Minimize_macromodel(mut_file_name, mut_structure):
    new_comfile_name = f'COM_FILE_{mut_file_name[:-4]}.com'
    # read in the com file.
    with open('mmod_mini_7.com', 'r') as comfile:
        comfile_contents = comfile.readlines()
        comfile_contents[0] = mut_file_name
        comfile_contents[1] = f'{mut_file_name[:-4]}_min_out.mae'
        comfile_string = "\n".join(comfile_contents)
        with open(new_comfile_name, 'w') as new_comfile:
            new_comfile.write(comfile_string)
        mmodel_min_job = queue.JobControlJob(['macromodel', new_comfile_name])
        jobDJ.addJob(mmodel_min_job)
        job_DJ_counter += 1
        if job_DJ_counter == 600:
            job_DJ_counter = 0
            jobDJ.run()
    GenGrid_Dock(mut_file_name, mut_structure)

# GenGrid - mutation_file_name
    # minimized_file_name

def GenGrid_Dock(min_file_name, mut_structure):
    # remove ligands here from min file.
    # need to read the minimized structure file.
    # then remove the ligand>?!

    # Look at minimized structure.
    # (1) Extract the co-ordinates of the ligand.
    # (2) Remove the ligand.
    # (3) Set the grid centre to the co-ordinates of the ligand.

    ligand_atoms = analyze.evaluate_asl(mut_structure, "ligand")

    ligand = analyze.find_ligands(mut_structure)
    centroid = transform.get_centroid(ligand[0].st)

    mut_structure.deleteAtoms(ligand_atoms)
###
    
    # Then want to confirm the grid center of the files that were minimized?
    # Can you get the positions from the ligand?

    #check the grid center for the files minimized
    grid_gen_spec = f"JOBNAME   gridgen\nGRID_CENTER   {centroid[0]}, {centroid[1]}, {centroid[2]}\nRECEP_FILE   {min_file_name}\nGRIDFILE   {min_file_name[:-4]}_grid.zip"
    grid_gen_file_name = f"{min_file_name[:-4]}_grid_gen.inp"
    with open(grid_gen_file_name, "w") as grid_gen_inp_file:
        grid_gen_inp_file.write(grid_gen_spec)  
    grid_gen_job = queue.JobControlJob(["glide", grid_gen_file_name])
    Dock(grid_gen_job)

def Dock(grid_gen_job):
    glide_XP_inp_spec = f"JOBNAME   glide_sp_dock\nOUTPUTDIR   Output/\nGRIDFILE   {modified_file_name[:-4]}_min_grid.zip\nLIGANDFILE   Combined_Ligands.mae\nPRECISION   XP\nWRITE_XP_DESC   True"
    glide_XP_inp_file_name = f"{modified_file_name[:-4]}_glide_XP.inp"
    with open(glide_XP_inp_file_name,"w") as glide_XP_inp_file:
        glide_XP_inp_file.write(glide_XP_inp_spec)
    
    glide_Xp_job = queue.JobControlJob(["glide", glide_XP_inp_file_name])
    glide_Xp_job.addPrereq(grid_gen_job)
    jobDJ.addJob(glide_Xp_job) 

def Analyse():
    pass

if __name__ == '__main__':
    job_DJ_counter = 0
    Mutate(1)