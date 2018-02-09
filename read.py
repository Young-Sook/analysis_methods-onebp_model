#!/data/reddylab/software/miniconda2/envs/YoungSook/bin/python2.7
import py2bit
import pyBigWig
import os
import numpy as np

def getDistributions():
	os.chdir("/data/reddylab/YoungSook/kmer/kmer_input/etoh_deep/EtOH_rep/onebp_model/distribution")

	###### SHEARING
        ###  mgw
        input_filename = "mgm_logp"
        input_stream = open(input_filename)
        input_file = input_stream.readlines()
        global MGW

	MGW = []
        for i in range(len(input_file)):
                temp = input_file[i].split()
                temp[1] = float(temp[1])
                MGW.append(temp)

        ### prot
        input_filename = "prot_logp"
        input_stream = open(input_filename)
        input_file = input_stream.readlines()
        global PROT

	PROT = []
        for i in range(len(input_file)):
                temp = input_file[i].split()
                temp[1] = float(temp[1])
                PROT.append(temp)

	#############  AMPLIFICATION_gibbs
        input_filename = "amplification_gibbs"
        input_stream = open(input_filename)
        input_file = input_stream.readlines()
        global GIBBS

	GIBBS = []
        for i in range(len(input_file)):
                temp = input_file[i].split()
                temp[1] = float(temp[1])
                GIBBS.append(temp)


def openInputFiles(ctrlbwFile, expbwFile, tbFile, mappableFile):
	global CTRLBW
	global EXPBW
	global TB

	global CTRLBW_SUM
	global EXPBW_SUM

	global CTRLBW_NUM
	global EXPBW_NUM

	global fileDir
	global MAPPABLE

	MAPPABLE = mappableFile

	dirname = ctrlbwFile[0] + "_" + expbwFile[0]
	fileDir = "/data/reddylab/YoungSook/kmer/kmer_script/modeling/modeling_onebp_DNA_RNA_ver4/" + dirname
	
	CTRLBW_NUM = len(ctrlbwFile)
	EXPBW_NUM = len(expbwFile) 

	CTRLBW = [0] * CTRLBW_NUM
	os.chdir("/data/reddylab/YoungSook/kmer/kmer_input/etoh_deep/EtOH_rep") 
        for i in range(CTRLBW_NUM):
		CTRLBW[i] = pyBigWig.open(ctrlbwFile[i])

	EXPBW = [0] * EXPBW_NUM
	os.chdir("/data/reddylab/YoungSook/kmer/kmer_input/etoh_deep/EtOH_rep/onebp_model/rna")
	for i in range(EXPBW_NUM):
		EXPBW[i] = pyBigWig.open(expbwFile[i])
	
        ## CALCULATE CTRLBW_SUM
	CTRLBW_SUM = float(CTRLBW[0].header().get('sumData'))
        
	#EXPBW_SUM = float(EXPBW.header().get('sumData')) 
	os.chdir("/data/reddylab/YoungSook/ref_genome/hg38")
        TB = py2bit.open(tbFile) 

def closeInputFiles1(ctrlbwFile, expbwFile):
	for i in range(1, CTRLBW_NUM):
		CTRLBW[i].close()
	del CTRLBW[1:CTRLBW_NUM]
	
	for i in range(1, EXPBW_NUM):
                EXPBW[i].close()
	del EXPBW[1:EXPBW_NUM]
	

def closeInputFiles2(ctrlbwFile, expbwFile, tbFile):
	CTRLBW[0].close()
	EXPBW[0].close()
	TB.close()


def getCandiTestRegions(refGenome):
	global SAFE_REG

	os.chdir("/data/reddylab/YoungSook/ref_genome")

	if( refGenome == "hg38"):
		input_filename = "filtered_hg38_ref_genome"
	        input_stream = open(input_filename)
        	input_file = input_stream.readlines()
       	 
        	SAFE_REG = []
        	for i in range(len(input_file)):
                	temp = input_file[i].split()
                	temp[1] = int(temp[1])
			temp[2] = int(temp[2])
			SAFE_REG.append(temp)

		SAFE_REG = np.array(SAFE_REG) 

	#if( refGenome == "hg19"):
		### SHOULD FIND A WAY 


