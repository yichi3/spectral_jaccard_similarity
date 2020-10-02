import os
from tqdm import tqdm,trange
import argparse

ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--datasets",  help="Text file with folder to dataset on each line", type = str, default = "Test_ds.txt")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=32 )
ap.add_argument("--lower",  help="Number of lowest hashes to compute", type=int,default=100)
ap.add_argument("--upper",  help="Number of highest hashes to compute", type=int,default=1000)
ap.add_argument("--step",  help="Steps between two hash number", type=int,default=50)

args = ap.parse_args()

num_jobs   = args.num_jobs
datasets   = args.datasets
lower      = args.lower
upper      = args.upper
step       = args.step


# cmd_0 = "python makeFolders.py --datasets {}".format(datasets)
# print(cmd_0)
# os.system(cmd_0)

datasetLst = [line.rstrip('\n') for line in open(datasets)]
# now we only focus on the first dataset
bact = datasetLst[0]
# now conduct the experiment to see the square error between ground truth JS and minHash JS
for i in tqdm(range(lower, upper+1, step)):
    cmd_1 = "python generateMinHashes_pipelined.py --dataset {}_filtered --num_jobs {} --num_hashes {}".format(
        bact, num_jobs, i)
    print(cmd_1)
    os.system(cmd_1)

    cmd_2 = "python computeJSimFromMinHashes.py --dataset {}_filtered --num_jobs {}".format(
        bact, num_jobs)
    print(cmd_2)
    os.system(cmd_2)

# compute the square error vs number of hash functions
