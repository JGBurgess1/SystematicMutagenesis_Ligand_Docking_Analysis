from mpi4py import MPI
import subprocess
import os

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

schrodinger_dir = os.environ.get("SCHRODINGER")

# Function to execute the Schrödinger Python script
def run_schrodinger_script(rank):
    # Path to your Schrödinger Python script
    schrodinger_script_path = "2_Mutate_and_Minimize.py"

    # Command to execute the script with MPI rank as argument
    command = [schrodinger_dir + "/run", schrodinger_script_path, '-c','prime_mmgbsa_test_Nov_16_1-out.mae','-l','Q203.mae',"-r", str(rank)]
    print(command)

    # Execute the command
    subprocess.run(command)

if __name__ == "__main__":
    # Check if running with MPI
    if size > 1:
        # Distribute work among MPI processes
        if rank == 0:
            print(f"Running with {size} MPI processes")

        # Call the function to execute the Schrödinger script
        run_schrodinger_script(rank)

    else:
        # If not using MPI, run the script normally
        print("Running without MPI")
        run_schrodinger_script(0)
