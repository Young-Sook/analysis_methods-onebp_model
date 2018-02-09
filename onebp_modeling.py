#!/data/reddylab/software/miniconda2/envs/YoungSook/bin/python2.7

import time
import os
import numpy as np
import math
import multiprocessing
import argparse
import py2bit
import pyBigWig
import random
import statsmodels.api as sm
import gc

import read     # distributions
import calculateOnebp 
import vari     # set global variables

def getArgs():
	parser = argparse.ArgumentParser()
	# required
	requiredArgs = parser.add_argument_group('Required Arguments')

	requiredArgs.add_argument('-ctrlbw', help="sorted control bigWig file. [Format: chr start end readcount]", nargs='+', required=True)
	
	requiredArgs.add_argument('-expbw', help="sorted exp bigWig file. [Format: chr start end readcount]", nargs='+', required=True)

	requiredArgs.add_argument('-fa', help="ht38.2bit file", required=True)
	
	requiredArgs.add_argument('-l', help="fragment Length", required=True)

	requiredArgs.add_argument('-rg', help="reference genome: hg38 or hg19", required=True)

	requiredArgs.add_argument('-mappable', help="The bdg file of mappability. [Format: chr start end mappabiality_score]", required=False)
	#optionalArgs = parser.add_argumnet_group('Optional Arguments')

	return parser

def divideGenome(chromo, chromoSize):
	subtask  = []
	## DETERMINE THE NUMBER OF BINS
	binSize = 5000000
	numBin = float(chromoSize) / binSize	
	if(numBin > int(numBin)):
		numBin = int(numBin) + 1

	## DIVIDE THE WHOLE GENOME INTO BINS
	for i in range(numBin):
		start = i*binSize
		end = (i+1)*binSize
		if i == 0:
			start = 1
		if i == (numBin-1): 
			end = chromoSize    # non-included
		subtask.append([chromo, start, end])   # end is not included in the analysis

	return subtask


def selectTrainSet(taskRegion):
	trainSet = []
	taskLen = len(taskRegion)
	trainSetNum = 0 	
	trainSize = 500000

	while trainSetNum < 5:
		idx = random.choice(range(taskLen)) 
		chromo = taskRegion[idx][0]
		start = taskRegion[idx][1]
		end = taskRegion[idx][2]	

		if( (end-trainSize) < start ):
			continue
		testRegion_start = random.choice(range(start, end-trainSize))
		testRegion_end = testRegion_start + trainSize	

		# check whether that region in SAFE or not -> numpy
		check = np.where( (read.SAFE_REG[:,0]==chromo) & (read.SAFE_REG[:,1].astype(int) <= testRegion_start) & (read.SAFE_REG[:,2].astype(int) > testRegion_end))
		if(len(check[0]) == 1 ):
			trainSetNum = trainSetNum + 1
			trainSet.append([chromo, testRegion_start, testRegion_end])

	return trainSet


def performRegression(trainSetResult):
	# the order of the column
	#colName = ['rc', 'group', 'intercect', 'shear1', 'shear2', 'anneal', 'denature', 'map', 'type', 'typeShear1', 'typeShear2', 'typeAnneal', 'typeDenature', 'typeMap', 'rep0', 'rep1' ,,, ]
		
	trainSetResult = np.asarray(trainSetResult)
	Y = trainSetResult[:,0]
	#Y = np.log(Y)
	#group = trainSetResult[:,1]
        X = trainSetResult[:,1:len(trainSetResult[0])]
			
	model = sm.GLM(Y, X, family=sm.families.Gaussian(link=sm.genmod.families.links.log)).fit()
	#model = sm.MixedLM(Y, X, group).fit()
	coef = model.params

	## CORRELATION                          
        corr = np.corrcoef(model.fittedvalues, Y)[0, 1]
	
	## THE MEAN OF REPLICATES FOR CTRL AND EXP
	avg_ctrl = np.mean(coef[ 12 : (12+read.CTRLBW_NUM) ])
	avg_exp = np.mean(coef[ (12+read.CTRLBW_NUM):(12 + read.CTRLBW_NUM + read.EXPBW_NUM) ])

	## COEFFICIENTS
	ctrlCoef = []
        #expCoef = []
        expScalerCoef = []
        ctrlCoef.extend(coef[0:6])
        ctrlCoef[0] = ctrlCoef[0] + avg_ctrl
        #expCoef.extend([ (coef[0]+coef[6] + avg_exp), (coef[1]+coef[7]), (coef[2]+coef[8]), (coef[3]+coef[9]), (coef[4]+coef[10]), (coef[5]+coef[11]) ])  
	expScalerCoef.extend(coef[6:12])
        expScalerCoef[0] = expScalerCoef[0] + avg_exp - avg_ctrl
	
	ctrl0_scaler = np.exp((avg_ctrl - coef[12]))
	exp0_scaler = np.exp((avg_exp - coef[(12+read.CTRLBW_NUM)]))


	print(corr)
	print(coef)	

	## WRITE A FILE TO PLOT OBSERVED VS PREDCITED
	if corr > 0.8:
		#print(corr)
		#print(coef)
	        print(model.summary())

		############### ANOVA TEST
		'''
		del trainSet_Y
		del trainSet_X

		trainSetResult[:,12] = np.log(trainSetResult[:,12]) 

		import pandas as pd
		import statsmodels.formula.api as smf
		from statsmodels.stats.anova import anova_lm

		df = pd.DataFrame(trainSetResult, columns=['intercept', 'shear1', 'shear2', 'anneal', 'denature', 'map', 'type', 'typeShear1', 'typeShear2', 'typeAnneal', 'typeDenature', 'typeMap', 'rc'])
		model = smf.ols(formula='rc ~ 1+shear1+shear2+anneal+denature+map+type+typeShear1+typeShear2+typeAnneal+typeDenature+typeMap', data=df).fit()
		print(model.summary())
		print(model._results.params)
		print(anova_lm(model))	
		'''

		########## WRITE DOWN PREDICTED AND OBSERVED 
		'''
		numRow = len(model.fittedvalues)
		numCol = 3
		result = np.ndarray(shape=(numRow, numCol), dtype=float)
		result[:,0] = model.fittedvalues
		result[:,1] = trainSet_Y
		result[:,2] = trainSetResult[:,6]

		roundcoef = round(coef, 5)
		os.chdir("/data/reddylab/YoungSook/kmer/kmer_script/modeling/modeling_onebp_DNA_RNA_ver2")
		output_filename = "observed_predicted_" + str(roundcoef)
        	output_stream = open(output_filename, "w")

        	for line in result:
                	output_stream.write('\t'.join([str(x) for x in line]) + "\n")

        	output_stream.close()
		'''

	return ctrl0_scaler, exp0_scaler, ctrlCoef, expScalerCoef, corr


def calculateScalerForNorm(resultMeta):
	norm_sum = 0

	for i in range(len(resultMeta)):
		tempFile = pyBigWig.open(resultMeta[i][0])
		norm_sum = norm_sum + float(tempFile.header().get('sumData'))
		tempFile.close()

	scalerNorm = read.CTRLBW_SUM / float(norm_sum)  

	print("Scaler:")

	print(scalerNorm)
	return scalerNorm


def mergeCorrectedBWfiles(resultMeta, scalerNorm, resultBWHeader, ctrlbwFileName, expbwFileName):
	os.chdir(read.fileDir)
	resultFileName1 = "CorrectedCtrlRead_" + ctrlbwFileName[0] + ".bw"  # should be modified (not fixed one)
	resultFileName2 = "CorrectedExpRead_" + expbwFileName[0] + ".bw"
	resultFileName3 = "Acitivity_" + expbwFileName[0] + ".bw"
	resultFile1 = pyBigWig.open(resultFileName1, "w")
	resultFile2 = pyBigWig.open(resultFileName2, "w")
	resultFile3 = pyBigWig.open(resultFileName3, "w")

	resultFile1.addHeader(resultBWHeader)    # put Header	
	resultFile2.addHeader(resultBWHeader)
	resultFile3.addHeader(resultBWHeader)

	for i in range(len(resultMeta)):
		tempFileName1 = resultMeta[i][0]
		tempFileName2 = resultMeta[i][1]
		tempChrom = resultMeta[i][2] 
		tempStart = int(resultMeta[i][3])
		tempEnd = int(resultMeta[i][4])

		tempbw1 = pyBigWig.open(tempFileName1)
		tempbw2 = pyBigWig.open(tempFileName2)
		chroms = np.array([tempChrom] * (tempEnd-tempStart))
		starts = np.array(range(tempStart, tempEnd), dtype=np.int64)
		ends = np.array(range((tempStart+1), (tempEnd+1)), dtype=np.int64)
		ctrlvals = np.array(tempbw1.values(tempChrom, tempStart, tempEnd), dtype=np.float64)

		if len(ctrlvals) == 0 :  
			tempbw1.close()
			tempbw2.close()
			os.remove(tempFileName1)
			os.remove(tempFileName2)
			continue

		ctrlvals = ctrlvals * scalerNorm	
		activals = np.array(tempbw2.values(tempChrom, tempStart, tempEnd), dtype=np.float64)	
		expvals = activals + ctrlvals

		resultFile1.addEntries(chroms, starts, ends, ctrlvals) 
		resultFile2.addEntries(chroms, starts, ends, expvals)
		resultFile3.addEntries(chroms, starts, ends, activals)

		tempbw1.close()
                tempbw2.close()

		# REMOVE TEMP FILES		
		os.remove(tempFileName1)
		os.remove(tempFileName2)

	resultFile1.close()
	resultFile2.close()
	resultFile3.close()

	return resultFileName1, resultFileName2, resultFileName3 


def main():	
	start_time = time.time()

	args = getArgs().parse_args()

	read.getDistributions()                   # read.mgw, read.prot, read.tm
	read.openInputFiles(args.ctrlbw, args.expbw, args.fa, args.mappable)  # read.mappable, read.bed, read.fa
	read.getCandiTestRegions(args.rg)
	vari.initializeGlobalVaris(args.l)    ## all global variables
	
	os.makedirs(read.fileDir)

	####  Divide the genomic regions and append it to task 
	chromosome = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                      "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
		      "chr21", "chr22", "chrX", "chrY"]
	
	task = []   # chromo, analysis coordinates

	resultBWHeader = []

	for i in range(len(chromosome)):
		chromoSize = read.CTRLBW[0].chroms(chromosome[i])
		if(chromoSize != None):
			resultBWHeader.append( (chromosome[i], chromoSize))
			chromoSize = int(chromoSize)
			task.extend(divideGenome(chromosome[i], chromoSize))


	### REPORT MEMORY USAGE
	import psutil
	process = psutil.Process(os.getpid())
	print(process.memory_info().rss)

	####  SELECT A TRAINING SET FROM CTRLBW	
	numProcess = 10
	corr = 0
	while corr < 0.8 :
		trainSet = selectTrainSet(task) #  EXCEP BL REGIONS, REMOVE ZERO COUNTS 
		#trainSet = [ ["chr10", 68885351, 69385351], ["chr11", 92913994, 93412926], ["chr5", 95895849, 96522862] ]
		print(trainSet)
	
		process = psutil.Process(os.getpid())
		print("After selecting trainset %d" % process.memory_info().rss )
		gc.collect()
        	pool = multiprocessing.Pool(numProcess)
		trainSetResult = pool.map_async(calculateOnebp.calculateTrainValues, trainSet).get()
		print("After training %d" % psutil.Process(os.getpid()).memory_info().rss )
		pool.close()
		pool.join()	
		print("After closing and joining pool %d" % psutil.Process(os.getpid()).memory_info().rss )
		del pool
		print("After deleting pool  %d" % psutil.Process(os.getpid()).memory_info().rss )
		gc.collect()

		for i in range(len(trainSetResult)):
			if i == 0:
				temp = list(trainSetResult[i])
			else:
				temp.extend(trainSetResult[i])

		trainSetResult = list(temp)
		del temp
		print(trainSetResult[0:5])

		print("After modifying trainSetResult form %d" % psutil.Process(os.getpid()).memory_info().rss )

		####   DO A REGRESSION USING CTRLBW
		vari.CTRL0SCALER, vari.EXP0SCALER, vari.CTRLCOEF, vari.EXPSCALERCOEF, corr  = performRegression(trainSetResult)
	
		del trainSetResult
		print("After erasing trainSetResult %d" % psutil.Process(os.getpid()).memory_info().rss )

		if( corr < 0.8 ):
			print("TRY A DIFFERENT TRAIN DATA")

		process = psutil.Process(os.getpid())
        	print(process.memory_info().rss)

	print("FINISH TRAINING")
	print("SCALER FOR CTRL AND EXP OBSERVED DATA: ")
	print(vari.CTRL0SCALER)
	print(vari.EXP0SCALER)
	read.CTRLBW_SUM = read.CTRLBW_SUM * vari.CTRL0SCALER

	#### CLOSE OTHER CTRLBW AND EXPBW EXCEPT THE FIRST ONES	
	read.closeInputFiles1(args.ctrlbw, args.expbw)
	print("Closing many ctrlBW and expBW")
	print(process.memory_info().rss)

	####   APPLY THE COEFFICIENTS TO OTHER EGIONS
	#numProcess = multiprocessing.cpu_count() -1 
        pool = multiprocessing.Pool(numProcess)
	resultMeta = pool.map_async(calculateOnebp.calculateTaskValues, task).get()   
	print(resultMeta[0:3]) ## In the result file : bigWig:  [ chr start end fittedReadCount ]
	pool.close()
	pool.join()
	gc.collect()
	print("FINISH WRITING TEMP FILES")

	'''
	## FOR CHECKING : -> WRTIE RESULTMETA
	os.chdir(read.fileDir)
	output_filename = "optimization_ver1_resultMeta" 
	output_stream = open(output_filename, "w")

	for line in resultMeta:
        	output_stream.write('\t'.join([str(x) for x in line]) + "\n")

	output_stream.close()

		
	os.chdir(read.filedir)
	input_filename = "optimization_ver1_resultMeta"
	input_stream = open(input_filename)
	input_file = input_stream.readlines()

	resultMeta = []
	for i in range(len(input_file)):
		temp = input_file[i].split()
		temp[3] = int(temp[3])
		temp[4] = int(temp[4])
		resultMeta.append(temp)
	'''

	## CALCULATE TOTAL READ RATIO OF (OBSERVED) TO  (OBSERVED/PREDICTED) 
	scalerNorm = calculateScalerForNorm(resultMeta)

	## MERGE ALL THE FILES
	print("MERGING ALL TEMP FILES FOR PREDICTED READ COUNTS")
	correctedCtrlBWName, correctedExpBWName, acitivityBWName = mergeCorrectedBWfiles(resultMeta, scalerNorm, resultBWHeader, args.ctrlbw, args.expbw)
	print(correctedCtrlBWName)
	print(correctedExpBWName)
	print(acitivityBWName)

	## should close all the input files
	read.closeInputFiles2(args.ctrlbw, args.expbw, args.fa)

	print("-- RUNNING TIME: %s hour(s)" % ((time.time()-start_time)/3600) )

	
if __name__ == "__main__":
	main()


