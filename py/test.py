from schrodinger import structure, job #, utility
from schrodinger.structutils import build
from schrodinger.job import queue
from schrodinger.job import launchapi
from schrodinger import jobdistributor
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import multiprocessing as mp
import threading as th
import time

amino_acids = ["ALA", "CYS", "ASP", "GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

jobDJ = queue.JobDJ()

print(mp.cpu_count())

done = False

def worker():
    counter = 0
    while not done:
        time.sleep(1)
        counter += 1
        print(counter)

for amino_acid in amino_acids:
    jobDJ.addJob(["prime", f"B_LEU_14_{amino_acid}_minimiz.inp", "-NJOBS", "24", "-D"])

with open("new.txt", "w") as newFile:
    for job in jobDJ.all_jobs:
        for string in job.getStatusStrings():
            newFile.write(string)

jobDJ.run()
