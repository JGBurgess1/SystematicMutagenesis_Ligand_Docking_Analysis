import argparse
from schrodinger import structure, job
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

# Parse command-line arguments
# ... (same as before)

# Create the output directory if it doesn't exist
output_directory = args.output
if not os.path.exists(output_directory):
    os.mkdir(output_directory)

# Load the ligand structure (assuming the ligand is the same for all mutations)
ligand_structure = structure.StructureReader(args.ligand)[0]

# Load the protein structure (each mutation will be performed on the same protein structure)
protein_structure = structure.StructureReader(args.protein)[0]

# Define the amino acids to mutate to (single-letter codes)
amino_acids = "ACDEFGHIKLMNPQRSTVWY"

# Lists to store mutation data
mutations = []
binding_scores = []

# Perform mutations, grid generation, and docking for each mutation
for residue_index in range(protein_structure.residue_count):
    original_residue = protein_structure.residue(residue_index)
    original_residue_name = original_residue.pdbname
    original_chain_id = original_residue.chain

    if original_residue_name not in amino_acids:
        continue

    for new_residue_name in amino_acids:
        if new_residue_name == original_residue_name:
            continue

        # Create a mutated copy of the protein structure
        mutated_structure = protein_structure.copy()
        mutated_structure.mutate_residue(residue_index, new_residue_name)

        # Energy minimize the mutated structure
        em_job = job.ShellJob()
        em_job.prereq = mutated_structure
        em_job.executable = "prime.energy_minimize"
        em_job.arguments = ["-s", mutated_structure, "-o", "em_output.mae"]
        em_job.submit()
        em_job.wait()

        # Generate a grid for Glide
        grid_job = job.ShellJob()
        grid_job.prereq = mutated_structure
        grid_job.executable = "glide"
        grid_job.arguments = [
            "-s", "0.5",
            "-x", "0.0", "-y", "0.0", "-z", "0.0",
            "-n", "20", "-m", "20", "-p", "20",
            "-inner_box", "10.0", "10.0", "10.0",
            "-outer_box", "20.0", "20.0", "20.0",
            "-NOJOBID",
            "-WAIT",
        ]
        grid_job.submit()
        grid_job.wait()

        # Dock the ligand
        dock_job = job.ShellJob()
        dock_job.prereq = [mutated_structure, ligand_structure]
        dock_job.executable = "glide"
        dock_job.arguments = [
            "-XP",
            "-NJOBS", "1",
            "-LIGAND", args.ligand,
            "-s", f"{original_residue_name}{original_chain_id}_{new_residue_name}_docking_output.sdf",
        ]
        dock_job.submit()
        dock_job.wait()

        # Parse docking scores from the output file
        dock_output_file = f"{original_residue_name}{original_chain_id}_{new_residue_name}_docking_output.sdf"
        binding_energy = None
        with structure.StructureReader(dock_output_file) as reader:
            for record in reader:
                if "Glide GScore" in record.property:
                    binding_energy = record.property["Glide GScore"]
                    break
        if binding_energy is not None:
            binding_scores.append(binding_energy)
            mutations.append(f"{original_residue_name}{original_chain_id}{residue_index + 1} to {new_residue_name}")

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

