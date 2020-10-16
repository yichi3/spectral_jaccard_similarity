from Bio.SeqIO.FastaIO import SimpleFastaParser
import itertools
import multiprocessing as mp
import numpy as np
import os
import argparse
import shutil
from tqdm import tqdm

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

def runSim(args):

	## avoid one processes starting multiple threads
	os.environ["MKL_NUM_THREADS"] = "1" 
	os.environ["NUMEXPR_NUM_THREADS"] = "1" 
	os.environ["OMP_NUM_THREADS"] = "1" 

	ref_read = args[0]
	dataset,numRandReads,numHashes = args[1]
	# Test use, before set numReads = 1000 to load the first 1000 read from dataset
	# now for test, just load first 3 reads
	numReads = 1000
	# numReads = 3

	## load precomputed minHashes
	minHashArr = np.zeros((numReads,numHashes))
	for i in range(numReads):
		minHashArr[i] = np.load(dataset+"minHashes/minHashes_{}.txt".format(i), allow_pickle=True)
		
	randMinHashArr = np.zeros((numRandReads,numHashes))
	for i in range(numRandReads):
		randMinHashArr[i] = np.load(dataset+"randReads/randMinHashes_{}.txt".format(i), allow_pickle=True)
	# minHashArrExtended = np.vstack((minHashArr[:,:1000],randMinHashArr))
	minHashArrExtended = np.vstack((minHashArr,randMinHashArr))

	## minHash collision matrix
	empiricalMatrix = (minHashArrExtended == minHashArrExtended[ref_read])
	# Test use
	# print(minHashArrExtended)
	# print()
	# print(empiricalMatrix)
	# print()
	empiricalMatrix = np.delete(empiricalMatrix,ref_read,axis=0)
	# Test use
	# print(empiricalMatrix)
	# print()


	## compute SVd
	U,s,VT = np.linalg.svd(empiricalMatrix - np.ones(empiricalMatrix.shape))
	# print(empiricalMatrix - np.ones(empiricalMatrix.shape))
	# print()

	u = U[:,0]
	pHatSVD = 1-np.abs(u[:n-1])/np.abs(np.median(u[n-1:]))
	
	qTemp = VT[0,:]
	qs = 1-np.abs(qTemp)/np.max(np.abs(qTemp))


	## save results
	np.savetxt(dataset+"/SVD/pi_refread_{}.txt".format(ref_read),pHatSVD)
	np.savetxt(dataset+"/SVD/qj_refread_{}.txt".format(ref_read),qs)

	np.savetxt(dataset+"/SVD/raw_pi_refread_{}.txt".format(ref_read),u)
	np.savetxt(dataset+"/SVD/raw_qj_refread_{}.txt".format(ref_read),qTemp)


ap = argparse.ArgumentParser(description="Reproduce the experiments in the manuscript")
ap.add_argument("--dataset",  help="Folder to dataset eg. NCTC5047_filtered")
ap.add_argument("--num_jobs", help="Num of parallel experiments", type=int, default=32 )
ap.add_argument("--num_rand", help="Number of random calibration reads", type = int, default=5)
ap.add_argument("--num_hashes",  help="Number of hashes to compute", type=int,default=1000)

args = ap.parse_args()

num_jobs   = args.num_jobs
dataset    = args.dataset        
numRandReads = args.num_rand
numHashes  = args.num_hashes

if dataset[-1] != '/':
	dataset += '/'

## delete all old SVDs
if os.path.isdir(dataset+"/SVD"):
	print("Deleting all old SVDs")
	shutil.rmtree(dataset+"/SVD")
print("Creating SVD folder")
os.mkdir(dataset+"/SVD")


reads_lst = []
with open(dataset+"reads.fasta") as handle:
	for values in SimpleFastaParser(handle):
		reads_lst.append(values[1])
n = len(reads_lst)

print("Parallelizing SVDs")
pool      = mp.Pool(processes=num_jobs)
arg_tuple =  itertools.product(range(n), [[dataset,numRandReads,numHashes]])
for _ in tqdm(pool.imap_unordered(runSim, arg_tuple), total=n):
	pass
