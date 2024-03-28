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

ligand_atoms = analyze.evaluate_asl(complex_structure, "ligand")

ligand = analyze.find_ligands(complex_structure)
print(ligand)

centroid = transform.get_centroid(ligand[0].st)

print(centroid)
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

ligand_atoms = analyze.evaluate_asl(complex_structure, "ligand")

ligand = analyze.find_ligands(complex_structure)
print(ligand)

centroid = transform.get_centroid(ligand.st)

print(centroid)