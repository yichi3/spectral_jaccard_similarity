import argparse, os
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm,trange
import multiprocessing as mp
import shutil, itertools

ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--dataset",  help="Folder to dataset eg. /data/MAB_alignment/ecoli_simulated2", default="NCTC74_filtered/")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=64 )

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
# print("Computing JSim matrix")
"""
To compute the JSim between two read, we need to find the total number of hash collisions out of number of hash
functions. All minhash values can be retrived from dataset/minHashes/minHashes_{}.txt
e.g. minHashArr[i] = np.load(dataset+"minHashes/minHashes_{}.txt".format(i), allow_pickle=True)
"""
minHashMatrix = np.loadtxt(dataset+"minHashes/minHashArr.txt")
_, length = minHashMatrix.shape
# now compute the jsim for each pair
for i in range(n):
    minHashArr_i = minHashMatrix[i]
    diffMatrix = minHashMatrix - minHashArr_i
    countZero = np.count_nonzero(diffMatrix == 0, axis=1)
    JSims[i] = countZero / length

np.savetxt(dataset+"JSimFromMinHashes/JSims_{}.txt".format(length),JSims)
