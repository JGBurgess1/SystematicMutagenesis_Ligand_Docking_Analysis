from schrodinger import structure
from schrodinger.structutils import analyze, transform
import argparse

parser = argparse.ArgumentParser(description="Check Centroid")
parser.add_argument("-c", "--complex", required=True, help="Path to the protein structure file (mae format)")

args = parser.parse_args()

complex_structure = structure.StructureReader.read(args.complex)
ligand_structure = structure.StructureReader.read(args.ligand)

ligand_atoms = analyze.evaluate_asl(complex_structure, "ligand")
ligand = analyze.find_ligands(complex_structure)
print(ligand)
centroid = transform.get_centroid(ligand[0].st)
print(centroid)