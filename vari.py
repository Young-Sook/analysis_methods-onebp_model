#!/data/reddylab/software/miniconda2/envs/YoungSook/bin/python2.7

import read
import numpy as np
import math

def initializeGlobalVaris(fragLenValue):
	global N_MGW
	global N_PROT
	global FRAGLEN
	global ENTROPY
	global MIN_TM
	global MAX_TM
	global SCALER
	global PARA1
	global PARA2
	global N_GIBBS
	global CTRLCOEF
	global EXPCOEF
	global EXPSCALERCOEF
	global CTRL0SCALER
	global EXP0SCALER
	global KMER

	KMER = 50

	temp = np.array(read.MGW)
	N_MGW = temp[:,1].astype(float).min()
	temp = np.array(read.PROT)
	N_PROT = temp[:,1].astype(float).min()

	## fragLen
	FRAGLEN = int(fragLenValue)

	##### AMPLIFICATION: GIBBS FREE ENERGY
	ENTROPY = -0.02485
	temp = np.array(read.GIBBS)
	MIN_TM = temp[:,1].astype(float).min()
	MIN_TM = -0.12 / ENTROPY
	MAX_TM = temp[:,1].astype(float).max()
	MAX_TM = -2.7 / ENTROPY
	SCALER = 1
	PARA1 = (math.pow(10, 6) - math.exp(SCALER)) / (math.pow(10, 6) - 1)
	PARA2 =  math.pow(10, -6) / (1-PARA1)
	N_GIBBS = np.median(temp[:,1].astype(float))


	
