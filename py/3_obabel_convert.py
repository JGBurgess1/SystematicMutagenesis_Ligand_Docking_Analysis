import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description="Converts pdb files to pdbqt files via openbabel.")
parser.add_argument("-r", "--rank", required=True, help="The rank of the mpi process")
args = parser.parse_args()
rank = int(args.rank)

def convert():
    size = 240
    
    arr = os.listdir("1_mutate_minimized")
    
    num_mutations = len(arr)
    chunk_size = num_mutations//size
    remainder = num_mutations%size

    start_index = rank*chunk_size + min(rank, remainder)
    end_index = start_index + chunk_size + (1 if rank < remainder else 0)

    os.chdir("1_mutate_minimized")

    for file in arr[start_index:end_index]:
        new_name = f"{file}qt"
        subprocess.run(['obabel',file,'-O',new_name])

if __name__ == "__main__":
    convert()