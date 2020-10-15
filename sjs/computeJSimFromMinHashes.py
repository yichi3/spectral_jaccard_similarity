import argparse, os
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm,trange
import multiprocessing as mp
import shutil, itertools
from tqdm import tqdm

ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--dataset",  help="Folder to dataset eg. /data/MAB_alignment/ecoli_simulated2")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=32 )

args = ap.parse_args()
dataset    = args.dataset        
num_jobs   = args.num_jobs

if dataset[-1] != '/':
	dataset += '/'


reads_lst = []
with open(dataset+"reads.fasta") as handle:
	for values in SimpleFastaParser(handle):
		reads_lst.append(values[1])
n = len(reads_lst)


if not os.path.isdir(dataset+"JSimFromMinHashes"):
	print("Creating JSimFromMinHashes folder")
	os.mkdir(dataset+"JSimFromMinHashes")

## reconstruct JS array from individual files
JSims = np.eye(n)
print("Computing JSim matrix")
"""
To compute the JSim between two read, we need to find the total number of hash collisions out of number of hash
functions. All minhash values can be retrived from dataset/minHashes/minHashes_{}.txt
e.g. minHashArr[i] = np.load(dataset+"minHashes/minHashes_{}.txt".format(i), allow_pickle=True)
"""
for i in trange(n-1):
    minHashArr_i = np.load(dataset+"minHashes/minHashes_{}.txt".format(i), allow_pickle=True)
    length = len(minHashArr_i)
    for j in range(i + 1, n):
        minHashArr_j = np.load(dataset+"minHashes/minHashes_{}.txt".format(j), allow_pickle=True)
        countZero = np.count_nonzero((minHashArr_i - minHashArr_j) == 0)
        val = countZero / length
        JSims[i,j] = val
        JSims[j,i] = val

np.savetxt(dataset+"JSimFromMinHashes/JSims.txt",JSims)
