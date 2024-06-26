# aweness.
from mpi4py import MPI # type: ignore
import subprocess
import os

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def run_obabel(rank):
    # Command to execute the script with MPI rank as argument
    command = ['python', '3_obabel_convert.py', '-r', str(rank)]
    print(command)
    # Execute the command
    subprocess.run(command)

if __name__ == "__main__":
    # Check if running with MPI
    if size > 1:
        # Distribute work among MPI processes
        run_obabel(rank)
    else:
        # If not using MPI, run the script normally
        print("Running without MPI")
        exit()
