#!/data/reddylab/software/miniconda2/envs/YoungSook/bin/python2.7

import sys
import os
import numpy as np
import math
import time
import tempfile

import read
import py2bit
import pyBigWig
import vari


############################################################
#############      USER DEFINED-FUNCTION     ###############
############################################################

def makeMatrix(chromo, start, end):
	result = []
	for i in range(start, end):
		result.append([chromo, i, 0, 0, 0, 0, 0])	
	return result
	

################################## SHEARING

def find5merProb(mer5):
        base_info = 0
	subtract = -1

        for i in range(5):
                if(mer5[i]=='A' or mer5[i]=='a'):
                        base_info = base_info + np.power(4, 4-i) * 0
                elif(mer5[i]=='C' or mer5[i]=='c'):
                        base_info = base_info + np.power(4, 4-i) * 1
                elif(mer5[i]=='G' or mer5[i]=='g'):
                        base_info = base_info + np.power(4, 4-i) * 2
                elif(mer5[i]=='T' or mer5[i]=='t'):
                        base_info = base_info + np.power(4, 4-i) * 3
		
		if(i==0):
			subtract = base_info

        prob_mgw = read.MGW[base_info][1]
	prob_prot = read.PROT[base_info][1]

	base_info = base_info - subtract	

        return base_info, prob_mgw, prob_prot


def edit5merProb(past_mer, oldBase, newBase):
	base_info = past_mer
	base_info = base_info * 4
	subtract = -1

	## newBase
        if(newBase=='A' or newBase=='a'):
                base_info = base_info + 0
        elif(newBase=='C' or newBase=='c'):
                base_info = base_info + 1
        elif(newBase=='G' or newBase=='g'):
                base_info = base_info + 2
        elif(newBase=='T' or newBase=='t'):
                base_info = base_info + 3

        prob_mgw = read.MGW[base_info][1]
        prob_prot = read.PROT[base_info][1]

        ## subtract oldBase
        if(oldBase=='A' or oldBase=='a'):
                subtract = np.power(4, 4) * 0
        elif(oldBase=='C' or oldBase=='c'):
                subtract = np.power(4, 4) * 1
        elif(oldBase=='G' or oldBase=='g'):
                subtract = np.power(4, 4) * 2
        elif(oldBase=='T' or oldBase=='t'):
                subtract = np.power(4, 4) * 3

	base_info = base_info - subtract

	return base_info, prob_mgw, prob_prot


def findComple5merProb(mer5):
        base_info = 0
        subtract = -1

        for i in range(5):
                if(mer5[i]=='A' or mer5[i]=='a'):
                        base_info = base_info + np.power(4, i) * 3
                elif(mer5[i]=='C' or mer5[i]=='c'):
                        base_info = base_info + np.power(4, i) * 2
                elif(mer5[i]=='G' or mer5[i]=='g'):
                        base_info = base_info + np.power(4, i) * 1
                elif(mer5[i]=='T' or mer5[i]=='t'):
                        base_info = base_info + np.power(4, i) * 0

                if(i==0):
                        subtract = base_info

        prob_mgw = read.MGW[base_info][1]
        prob_prot = read.PROT[base_info][1]
        base_info = base_info - subtract

        return base_info, prob_mgw, prob_prot


def editComple5merProb(past_mer, oldBase, newBase):
	base_info = past_mer
	base_info = base_info / 4
	subtract = -1

        # newBase
	if(newBase=='A' or newBase=='a'):
		base_info = base_info + np.power(4, 4) * 3
	elif(newBase=='C' or newBase=='c'):
		base_info = base_info + np.power(4, 4) * 2
	elif(newBase=='G' or newBase=='g'):
		base_info = base_info + np.power(4, 4) * 1
	elif(newBase=='T' or newBase=='t'):
		base_info = base_info + np.power(4, 4) * 0

        prob_mgw = read.MGW[base_info][1]
	prob_prot = read.PROT[base_info][1]

	## subtract oldBase
	if(oldBase=='A' or oldBase=='a'):
                subtract = 3
        elif(oldBase=='C' or oldBase=='c'):
                subtract = 2
        elif(oldBase=='G' or oldBase=='g'):
                subtract = 1
        elif(oldBase=='T' or oldBase=='t'):
                subtract = 0
	
	base_info = base_info - subtract

        return base_info, prob_mgw, prob_prot


################################## GIBBS FREE ENERGY

def findStartGibbs(seq, seqLen):
        gibbs = 0
	subtract = -1

        for i in range(seqLen-1):
                dimer = seq[i:(i+2)].upper()
                if( 'N' in dimer):
                        gibbs = gibbs + vari.N_GIBBS
                else:
                        dimer_idx = 0
                        for j in range(2):
                                if(dimer[j]=='A'):
                                        dimer_idx = dimer_idx + np.power(4, 1-j) * 0
                                elif(dimer[j]=='C'):
                                        dimer_idx = dimer_idx + np.power(4, 1-j) * 1
                                elif(dimer[j]=='G'):
                                        dimer_idx = dimer_idx + np.power(4, 1-j) * 2
                                elif(dimer[j]=='T'):
                                        dimer_idx = dimer_idx + np.power(4, 1-j) * 3
                        gibbs = gibbs + read.GIBBS[dimer_idx][1]

		if(i==0):
			subtract = gibbs

	start_gibbs = gibbs - subtract 

        return start_gibbs, gibbs


def editStartGibbs(oldDimer, newDimer, past_start_gibbs):
        gibbs = past_start_gibbs
	subtract = -1

	# newDimer
        if( 'N' in newDimer):
                gibbs = gibbs + vari.N_GIBBS
        else:
                dimer_idx = 0
                for j in range(2):
                        if(newDimer[j]=='A'):
                                dimer_idx = dimer_idx + np.power(4, 1-j) * 0
                        elif(newDimer[j]=='C'):
                                dimer_idx = dimer_idx + np.power(4, 1-j) * 1
                        elif(newDimer[j]=='G'):
                                dimer_idx = dimer_idx + np.power(4, 1-j) * 2
                        elif(newDimer[j]=='T'):
                                dimer_idx = dimer_idx + np.power(4, 1-j) * 3
                gibbs = gibbs + read.GIBBS[dimer_idx][1]
	
	## remove the old dimer for the next iteration
	if( 'N' in oldDimer):
                subtract = vari.N_GIBBS
        else:
                dimer_idx = 0
                for j in range(2):
                        if(oldDimer[j]=='A'):
                                dimer_idx = dimer_idx + np.power(4, 1-j) * 0
                        elif(oldDimer[j]=='C'):
                                dimer_idx = dimer_idx + np.power(4, 1-j) * 1
                        elif(oldDimer[j]=='G'):
                                dimer_idx = dimer_idx + np.power(4, 1-j) * 2
                        elif(oldDimer[j]=='T'):
                                dimer_idx = dimer_idx + np.power(4, 1-j) * 3
                subtract = read.GIBBS[dimer_idx][1]
 
	start_gibbs = gibbs - subtract

        return start_gibbs, gibbs


def convertGibbs(gibbs):
	tm = gibbs / (vari.ENTROPY*(vari.FRAGLEN-1))
        ## rescaling of tm
        tm = (tm - vari.MIN_TM) / (vari.MAX_TM - vari.MIN_TM)

        # anneal
        anneal_prob = ( math.exp(tm*vari.SCALER) - vari.PARA1 ) * vari.PARA2   # actually this should be anneal. Correct naming
        anneal_prob = math.log(anneal_prob)

        ### denature
        tm = tm -1
        denature_prob = ( math.exp( tm*(-1)*vari.SCALER ) - vari.PARA1 ) * vari.PARA2 # this should be denature. Correct naming
        denature_prob = math.log(denature_prob)

	return anneal_prob, denature_prob


def calculateStartP(start, fasta, mer2, prob_mgw1, prob_prot1, past_mer2, past_start_gibbs):
        mgw_prob = prob_mgw1
        prot_prob = prob_prot1  
     
        if('N' in mer2):
                baseInfo = -1
                mgw_prob = mgw_prob + vari.N_MGW
                prot_prob = prot_prob + vari.N_PROT
        else:
                if(past_mer2 == -1):
                        baseInfo, add1, add2 = findComple5merProb(mer2)
                else:
                        baseInfo, add1, add2 = editComple5merProb(past_mer2, mer2[0], mer2[4])
                mgw_prob = mgw_prob + add1
                prot_prob = prot_prob + add2


        ########  AMPLIFICATION: TM
        if past_start_gibbs == -1:
                start_gibbs, gibbs = findStartGibbs(fasta, vari.FRAGLEN)
        else:
		oldDimer = fasta[0:2].upper()
                newDimer = fasta[(vari.FRAGLEN-2):vari.FRAGLEN].upper()
                start_gibbs, gibbs = editStartGibbs(oldDimer, newDimer, past_start_gibbs)

	anneal_prob, denature_prob = convertGibbs(gibbs)

	prob = [mgw_prob, prot_prob, anneal_prob, denature_prob]
	
        return baseInfo, start_gibbs, gibbs, prob


def addObservedRead(result, chromo, analysis_start, analysis_end):
	ob_ctrlrc = [0] * read.CTRLBW_NUM
	ob_exprc = [0] * read.EXPBW_NUM

	for rep in range(read.CTRLBW_NUM):
		ob_ctrlrc[rep] = read.CTRLBW[rep].values(chromo, analysis_start, analysis_end)

	for rep in range(read.EXPBW_NUM):
		ob_exprc[rep] = read.EXPBW[rep].values(chromo, analysis_start, analysis_end)	

	final = []

	for i in range(len(result)):
		score_ctrl = 0
		score_exp = 0
		pass_ctrl = read.CTRLBW_NUM / float(2)
		pass_exp = read.EXPBW_NUM / float(2)

		ctrl_rep = [0] * read.CTRLBW_NUM
		exp_rep = [0] * read.EXPBW_NUM

		## Quality Check_CTRL
		for rep in range(read.CTRLBW_NUM):
			if (ob_ctrlrc[rep][i] != 0) and (math.isnan(ob_ctrlrc[rep][i]) == False):
				score_ctrl = score_ctrl + 1
				ctrl_rep[rep] = ob_ctrlrc[rep][i]
				
		## Quality Check_EXP
		for rep in range(read.EXPBW_NUM):
                        if (ob_exprc[rep][i] != 0) and (math.isnan(ob_exprc[rep][i]) == False):
                                score_exp = score_exp + 1
				exp_rep[rep] = ob_exprc[rep][i]


		if (score_ctrl >= pass_ctrl) and (score_exp >= pass_exp):
			## the order of the array: ['rc', 'group', 'intercect', 'shear1', 'shear2', 'anneal', 'denature', 'map', 'type', 'typeShear1', 'typeShear2', 'typeAnneal', 'typeDenature', 'typeMap', 'rep0', 'rep1' ,,, ]
			## FOR CONTROL
			for rep in range(read.CTRLBW_NUM):
				temp = [ ctrl_rep[rep], 1  ]
				temp.extend(list(result[i][2:7]))
				temp.extend([ 0, 0, 0, 0, 0, 0])  # 0: binary varialbe 
				# for each replicte
				repArray = [0] * (read.CTRLBW_NUM + read.EXPBW_NUM)
				repArray[rep] = 1				
				temp.extend(repArray)
				final.append(temp)			
				
			## FOR EXP
			for rep in range(read.EXPBW_NUM):
				temp = [ exp_rep[rep], 1  ]
                        	temp.extend(list(result[i][2:7]))
                        	temp.extend([1]) # type variable
				temp.extend(list(result[i][2:7]))
				repArray = [0] * (read.CTRLBW_NUM + read.EXPBW_NUM)
				repArray[(read.CTRLBW_NUM + rep)] = 1
				temp.extend(repArray)
				final.append(temp)

	del result

	return final


def correctRead(result, chromo, analysis_start, analysis_end, chromoEnd):
	print("start CORRECTREAD")
	# subfinal1 : FOR CORRECTED CTRLBW 	
	subfinal1 = tempfile.NamedTemporaryFile(suffix=".bw", dir=read.fileDir, delete=False)
	subfinal1.close()

	## subfinal2 : FOR CORRECTED EXPBW
	subfinal2 = tempfile.NamedTemporaryFile(suffix=".bw", dir=read.fileDir, delete=False)
        subfinal2.close()

	bw1 = pyBigWig.open(subfinal1.name, "w")	
	bw1.addHeader([(chromo, chromoEnd)]) ## gotta work on this part
	bw2 = pyBigWig.open(subfinal2.name, "w")
        bw2.addHeader([(chromo, chromoEnd)]) ## gotta work on this part
	
	ob_ctrlrc = read.CTRLBW[0].values(chromo, analysis_start, analysis_end)
	ob_exprc = read.EXPBW[0].values(chromo, analysis_start, analysis_end)

	for i in range(len(result)):
		temp = list(result[i])
		start = temp[1]
		end = temp[1] + 1
		
		if (ob_ctrlrc[i] == 0) or (math.isnan(ob_ctrlrc[i]) == True) or (ob_exprc[i] == 0) or (math.isnan(ob_exprc[i]) == True):
			continue 
		else:
			### PREDICTED READ COUNT FOR CONTRL DATA
			prdCtrl = vari.CTRLCOEF[0] + vari.CTRLCOEF[1] * temp[2] + vari.CTRLCOEF[2] * temp[3] + vari.CTRLCOEF[3] * temp[4] + vari.CTRLCOEF[4] * temp[5] + vari.CTRLCOEF[5] * temp[6]
			prdCtrl = math.exp(prdCtrl)
			if prdCtrl < 1:
                                prdCtrl= 1.0
			
			### PREDICTED READ COUNT FOR EXP DATA
			'''
			prdExp = vari.EXPCOEF[0] + vari.EXPCOEF[1] * temp[2] + vari.EXPCOEF[2] * temp[3] + vari.EXPCOEF[3] * temp[4] + vari.EXPCOEF[4] * temp[5] + vari.EXPCOEF[5] * temp[6]
			prdExp = math.exp(prdExp)
			if prdExp < 1:
				prdExp = 1.0
			'''			

			scaler = vari.EXPSCALERCOEF[0] + vari.EXPSCALERCOEF[1] * temp[2] + vari.EXPSCALERCOEF[2] * temp[3] + vari.EXPSCALERCOEF[3] * temp[4] + vari.EXPSCALERCOEF[4] * temp[5] + vari.EXPSCALERCOEF[5] * temp[6]	
			scaler = math.exp(scaler)
			if scaler < 1:
                                scaler = 1.0

			correctedCtrl = ( ob_ctrlrc[i] * vari.CTRL0SCALER ) / float(prdCtrl)
			bw1.addEntries([chromo], [start], ends=[end], values=[correctedCtrl])  # rescaled in the next step

			correctedExp =   ((ob_exprc[i] * vari.EXP0SCALER) / float(scaler)) - prdCtrl  # True signal 
			bw2.addEntries([chromo], [start], ends=[end], values=[correctedExp])  # rescaled in the next step

			
	bw1.close()
	bw2.close()

	del result

	return_array = [subfinal1.name, subfinal2.name, chromo, analysis_start, analysis_end]

	return return_array 


def convertMappability(map1, map2):
        ## map1
        if (map1 == 0) or (math.isnan(map1) == True):
                map1 = -6
        else:
                map1 = math.log(map1)

        ## map2
        if (map2 == 0) or (math.isnan(map2) == True):
                map2 = -6
        else:
                map2 = math.log(map2)

        map_prob = map1 + map2

        return map_prob


################################################################
########     CALCULATING PROBABILITY FOR EACH BASE    ##########
################################################################


def calculateTrainValues(args):
	chromo= args[0]
	analysis_start = int(args[1])  # Genomic coordinates(starts from 1)
	analysis_end = int(args[2])
	chromoEnd = int(read.CTRLBW[0].chroms(chromo))

	#### make final result matrix
	frag_start = analysis_start - vari.FRAGLEN + 1  # This is the start position of the first fragment 
	frag_end = analysis_end + vari.FRAGLEN - 1      # NOT INCLUDED: This is the (start position of the last fragment-1)
	shear_start = frag_start - 2 
	shear_end = frag_end + 2		   # NOT INCLUDED

	###### of output_start nad end are out of chromosome size
	if(shear_start < 1): # Position(not index!)
		shear_start = 1
		frag_start = 3    # shear_start + 2
		analysis_start = max(analysis_start, frag_start)
	if(shear_end > chromoEnd): # chromoEnd is included in the real choromosome. chr:1-chromoSize
		shear_end = chromoEnd
		frag_end = shear_end - 2
		analysis_end = min(analysis_end, frag_end)	

	#### mappable file
        os.chdir("/data/reddylab/YoungSook/ref_genome/mappability")
        totalMap = pyBigWig.open(read.MAPPABLE)
	mappable = totalMap.values(chromo, shear_start, shear_end)


	##### CALL FASTA FILE AND MAPPABLE FILE IN THE RANGE
	fa = (read.TB).sequence(chromo, (shear_start-1), (shear_end-1))    # Due to py2bit coordiante system
	#print("problem here? : 1")
	### mappability
        mappable = totalMap.values(chromo, shear_start, shear_end)

	###### MAKE A RESULT MATRIX
	result = makeMatrix(chromo, frag_start, frag_end)

	##### INITILIZATION OF MER1, GIBBS PARAMETERS
	past_mer1 = -1
	past_start_gibbs = -1

	##### IDX CRETERIA: FASTA FILE, MAPPABLE FILE
	start_idx = 2                                    # index in the fasta file (Included in the range)
	end_idx = (frag_end - vari.FRAGLEN) - shear_start + 1   # index in the fasta file (Not included in the range)   

	for i in range(start_idx, end_idx):  # index in the fasta file
		#########  SHEARING_MER1
		mer1 = fa[(i-2):(i+3)]

		if('N' in mer1):
			past_mer1 = -1
			prob_mgw1 = vari.N_MGW
			prob_prot1 = vari.N_PROT
		else:
			if(past_mer1 == -1): # there is no information on past_mer1
				past_mer1, prob_mgw1, prob_prot1 = find5merProb(mer1)
			else:
				past_mer1, prob_mgw1, prob_prot1 = edit5merProb(past_mer1, mer1[0], mer1[4])

		#########  START THE ITERATION WITH VARIOUS WINDOWS
		past_mer2 = -1

		########   FRAGMENT 
		fragEnd_idx = i + vari.FRAGLEN
		fasta = fa[i:fragEnd_idx]
		mer2 = fa[(fragEnd_idx-3):(fragEnd_idx+2)]
		
		########   MAPPABILITY1
                map1 = mappable[i]
                map2 = mappable[(fragEnd_idx - vari.KMER)]
                map_prob = convertMappability(map1, map2)

		past_mer2, past_start_gibbs, past_gibbs, prob  = calculateStartP(i, fasta, mer2, prob_mgw1, prob_prot1, past_mer2, past_start_gibbs)

		####### UPDATE THE RESULT FILE
        	for j in range(0, vari.FRAGLEN):
                	result[i-2+j][2] = result[i-2+j][2] + prob[0]        # shearing_mgm
                	result[i-2+j][3] = result[i-2+j][3] + prob[1]        # shearing_prot
                	result[i-2+j][4] = result[i-2+j][4] + prob[2]        # anneal
                	result[i-2+j][5] = result[i-2+j][5] + prob[3]        # denature
			result[i-2+j][6] = result[i-2+j][6] + map_prob       # mappability
		
	####### select regions
	analysis_start_idx = analysis_start - frag_start
	analysis_end_idx = analysis_end- frag_start
	result = result[analysis_start_idx:analysis_end_idx]

	print("NO")
	final = addObservedRead(result, chromo, analysis_start, analysis_end)

	return final


def calculateTaskValues(args):
	chromo= args[0]
        analysis_start = int(args[1])  # Genomic coordinates(starts from 1)
        analysis_end = int(args[2])
        chromoEnd = int(read.CTRLBW[0].chroms(chromo))

        #### make final result matrix
        frag_start = analysis_start - vari.FRAGLEN + 1  # This is the start position of the first fragment
        frag_end = analysis_end + vari.FRAGLEN - 1      # NOT INCLUDED: This is the (start position of the last fragment-1)
        shear_start = frag_start - 2
        shear_end = frag_end + 2                   # NOT INCLUDED

        ###### of output_start nad end are out of chromosome size
        if(shear_start < 1): # Position(not index!)
                shear_start = 1
                frag_start = 3    # shear_start + 2
                analysis_start = max(analysis_start, frag_start)
        if(shear_end > chromoEnd): # chromoEnd is included in the real choromosome. chr:1-chromoSize
                shear_end = chromoEnd
                frag_end = shear_end - 2
                analysis_end = min(analysis_end, frag_end)

	#### mappable file
        os.chdir("/data/reddylab/YoungSook/ref_genome/mappability")
        totalMap = pyBigWig.open(read.MAPPABLE)
        mappable = totalMap.values(chromo, shear_start, shear_end)

	#print("problem here?????")
        ##### CALL FASTA FILE AND MAPPABLE FILE IN THE RANGE
        fa = (read.TB).sequence(chromo, (shear_start-1), (shear_end-1))    # Due to py2bit coordiante system
	#print("problem here? 3")
	### mappability
        mappable = totalMap.values(chromo, shear_start, shear_end)     


        ###### MAKE A RESULT MATRIX
        result = makeMatrix(chromo, frag_start, frag_end)

        ##### INITILIZATION OF MER1, GIBBS PARAMETERS
        past_mer1 = -1
        past_start_gibbs = -1

        ##### IDX CRETERIA: FASTA FILE, MAPPABLE FILE
        start_idx = 2                                    # index in the fasta file (Included in the range)
        end_idx = (frag_end - vari.FRAGLEN) - shear_start + 1   # index in the fasta file (Not included in the range)

	for i in range(start_idx, end_idx):  # index in the fasta file
                #########  SHEARING_MER1
                mer1 = fa[(i-2):(i+3)]

                if('N' in mer1):
                        past_mer1 = -1
                        prob_mgw1 = vari.N_MGW
                        prob_prot1 = vari.N_PROT
                else:
                        if(past_mer1 == -1): # there is no information on past_mer1
                                past_mer1, prob_mgw1, prob_prot1 = find5merProb(mer1)
                        else:
                                past_mer1, prob_mgw1, prob_prot1 = edit5merProb(past_mer1, mer1[0], mer1[4])

                #########  START THE ITERATION WITH VARIOUS WINDOWS
                past_mer2 = -1

                ########   FRAGMENT
                fragEnd_idx = i + vari.FRAGLEN
                fasta = fa[i:fragEnd_idx]
                mer2 = fa[(fragEnd_idx-3):(fragEnd_idx+2)]

		########   MAPPABILITY1
                map1 = mappable[i]
                map2 = mappable[(fragEnd_idx - vari.KMER)]
                map_prob = convertMappability(map1, map2)

                past_mer2, past_start_gibbs, past_gibbs, prob  = calculateStartP(i, fasta, mer2, prob_mgw1, prob_prot1, past_mer2, past_start_gibbs)

                ####### UPDATE THE RESULT FILE
                for j in range(0, vari.FRAGLEN):
                        result[i-2+j][2] = result[i-2+j][2] + prob[0]        # shearing_mgm
                        result[i-2+j][3] = result[i-2+j][3] + prob[1]        # shearing_prot
                        result[i-2+j][4] = result[i-2+j][4] + prob[2]        # anneal
                        result[i-2+j][5] = result[i-2+j][5] + prob[3]        # denature
			result[i-2+j][6] = result[i-2+j][6] + map_prob

	####### select regions
        analysis_start_idx = analysis_start - frag_start
        analysis_end_idx = analysis_end- frag_start
        result = result[analysis_start_idx:analysis_end_idx]

	final = correctRead(result, chromo, analysis_start, analysis_end, chromoEnd)
	print(final)
	print("Finished Calculating")

	return final


#############################################################
##################       WRITE OUTPUT      ##################
#############################################################

#os.chdir("/data/reddylab/YoungSook/kmer/kmer_input/etoh_deep/EtOH_rep/onebp_model/optimized")
#output_filename = sys.argv[1] + "_optimized_ver1_shear_gc"

#output_stream = open(output_filename, "w")

# for line in final:
#	output_stream.write('\t'.join([str(x) for x in line]) + "\n")

#output_stream.close()

#print("FIXED WINDOW:")
#print(sys.argv[1])
#print("-- RUNNING TIME: %s min" % ((time.time()-start_time)/60) )


###################################
###### FASTA FILE :
# START : ANALYSIS_START - FRAGLen -1
# END: ANALYSIS_START + FRAGLEN + 1


