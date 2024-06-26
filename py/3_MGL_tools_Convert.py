import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description="Converts pdb files to pdbqt files via mgltools.")
parser.add_argument("-r", "--rank", required=True, help="The rank of the mpi process")
parser.add_argument("-s", "--size", required=True, help="The size of the mpi process pool")
args = parser.parse_args()
rank = int(args.rank)
size = int(args.size)

def convert():
    arr = os.listdir("1_mutate_minimized")
    
    num_mutations = len(arr)
    print(num_mutations)
    chunk_size = num_mutations//size
    remainder = num_mutations%size

    start_index = rank*chunk_size + min(rank, remainder)
    end_index = start_index + chunk_size + (1 if rank < remainder else 0)

    os.chdir("1_mutate_minimized")

    for file in arr[start_index:end_index]:
        subprocess.run(['python2','/home/jburgess1/.conda/envs/open_docking_env/bin/prepare_receptor4.py','-r',file, '-A', 'none','-U','nphs'])

if __name__ == "__main__":
    convert()