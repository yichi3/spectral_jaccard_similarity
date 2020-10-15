import numpy as np
import argparse
from tqdm import tqdm, trange
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--dataset",  help="Folder to dataset eg. /data/MAB_alignment/ecoli_simulated2", default="NCTC74_filtered/")
ap.add_argument("--lower",  help="Number of lowest hashes to compute", type=int,default=100)
ap.add_argument("--upper",  help="Number of highest hashes to compute", type=int,default=1000)
ap.add_argument("--step",  help="Steps between two hash number", type=int,default=50)


args = ap.parse_args()
dataset    = args.dataset        
lower      = args.lower
upper      = args.upper
step       = args.step

numHashes = np.arange(lower, upper+1, step)
mse = np.zeros_like(numHashes)
JSims_gt = np.loadtxt("JSims_gt.txt")
for i in trange(len(numHashes)):
    curJSimsMatrix = np.loadtxt(dataset+"JSimFromMinHashes/JSims_{}.txt".format(numHashes[i]))
    mse[i] = np.sum((curJSimsMatrix - JSims_gt)**2)
print(mse)
# plot the mse vs numHash
plt.figure()
plt.plot(numHashes, mse)
plt.title("Plot of MSE vs number of hash functions")
plt.xlabel("Number of hash functions")
plt.ylabel("MSE")
plt.show()
plt.savefig('mse_vs_numHashes.png')