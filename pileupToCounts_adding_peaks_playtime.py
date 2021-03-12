# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 12:17:59 2020

@author: dal1858
"""

#%%

import pandas as pd

myDir = r'X:\Users\dloftus\allelic_counts_test_expt_comp2\attempt_2_b6_cast_N-masked_genome'
pileupFile = f'{myDir}\SRR7685829Aligned.sortedByCoord.out.pileup'
snpFile = f'{myDir}\\all_SNPs_CAST_EiJ_GRCm38.bed'

pileup = pd.read_csv(pileupFile, 
                     sep = '\t', 
                     names = ['chr', 'end_coord', 'pileup_ref', 'reads', 'bases', 'qual'], 
                     dtype = {'chr': object, 'end_coord': int, 'pileup_ref': object, 'reads': object, 'bases': str, 'qual': str})
print('Pileup file loaded')

snps = pd.read_csv(snpFile, 
                   sep = '\t', 
                   names = ['chr', 'start_coord', 'end_coord', 'ref/alt'], 
                   dtype = {'chr': object, 'start_coord': int, 'end_coord': int, 'ref/alt': str})
print('SNPs file loaded')

#%%

def mergeDFs(myPileup, mySNPs):
    
    ''' merges the pileup and snps dataframes and returns only the needed columns '''
    
    snpsCols = ['chr', 'end_coord', 'ref/alt']
    pileupCols = ['chr', 'end_coord', 'pileup_ref', 'bases']
    
    mergedDF = pd.merge(myPileup[pileupCols], mySNPs[snpsCols], on = ['chr', 'end_coord'])
    
    return mergedDF

df = mergeDFs(pileup, snps)
print('Files merged')

#%%

from pandarallel import pandarallel

pandarallel.initialize()

#%%

#filter out SNPs where the pileup file and the snp file have different reference bases     
def findRefAllele(alleles):
    
    ref = alleles[0].upper()
    return ref

def pileupRefToUpper(base):
    
    baseUpper = base.upper()
    return baseUpper

snpsRef = df['ref/alt'].apply(findRefAllele)
pileupRefUpper = df['pileup_ref'].apply(pileupRefToUpper)
df['snps_ref'] = snpsRef
df['pileup_ref_upper'] = pileupRefUpper

refBool = (df['snps_ref'] == df['pileup_ref_upper'])
df = df[refBool]
df.reindex()



#filter out starts and ends of reads, which tend to show a reference bias 
def filterBases(bases):
    
    ''' 
    removes all bases before the last '$' character and after 
    the first '^' character, which designate the end and start of reads 
    and have a systemic reference bias
    '''
    
    if '$' in bases:
        bases = bases.split('$')[-1]
    if '^' in bases:
        bases = bases.split('^')[0]
        
    return bases

df['bases_filtered'] = df['bases'].apply(filterBases)
print('SNPs filtered')



#count ref and alt bases 
def countRef(bases):
    
    ''' counts and returns the number of ref alleles overlapping a SNP '''
    
    count = bases.count('.') + bases.count(',')
    return count

def countAlt(row):
    
    ''' counts and returns the number of alt alleles overlapping a SNP '''
    
    altAllele = row['ref/alt'][2].upper()
    count = row['bases_filtered'].upper().count(altAllele)
    return count

df['ref_count'] = df['bases_filtered'].apply(countRef)
print('Ref alleles counted')
df['alt_count'] = df.apply(countAlt, axis = 1)
print('Alt alleles counted')






#%%
import pandas as pd
myDir = r'X:\Users\dloftus\allelic_counts_test_expt_comp2\attempt_2_b6_cast_N-masked_genome'

bpFile = f'{myDir}\macs2\SRR7685829_stringent_peaks.broadPeak'

bp = pd.read_csv(bpFile, 
                 sep = '\t', 
                 names = ['chr', 'start_coord', 'end_coord', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue'])




#%%

'''
generate random broadpeaks and counts dataframes to do analysis on
'''

bpSmall = bp.sample(int(1e3))
dfSmall = df.sample(int(1e6))


#%%

import numpy as np
refAllele = 'pat'
def quantPeaks(bpRow):
    
    pileupBoolChr = (df['chr'] == bpRow['chr'])
    pileupBoolCoords1 = (df['end_coord'] <= bpRow['end_coord'])
    pileupBoolCoords2 = (df['end_coord'] >= bpRow['start_coord'])
    peakSNPs = df[pileupBoolChr & pileupBoolCoords1 & pileupBoolCoords2]

    ref_count = sum(np.array(peakSNPs['ref_count']))
    alt_count = sum(np.array(peakSNPs['alt_count']))  
    total_count = ref_count + alt_count
    print(total_count)
    if refAllele == 'mat':
        matFrac = np.divide(ref_count, total_count, where = total_count != 0)
    elif refAllele == 'pat':
        matFrac = np.divide(alt_count, total_count, where = total_count != 0)
        
    return matFrac


#test = bpTest.apply(quantPeaks, axis = 1)
#%%

class peak:
    
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        
    def length(self):
        return(self.end - self.start)

bpClass = []

for i, row in bp.iterrows():
    
    obj = peak(row['chr'], row['start_coord'], row['end_coord'])
    bpClass.append(obj)


#%%
def sizePD(df):
    
    def diff(row):
        myDiff = row['end_coord'] - row['start_coord']
        return myDiff
    
    size = df.apply(diff, axis = 1)
    return size
        


def sizeClass(myList):
    
    myLen = []
    for obj in myList:
        myLen.append(obj.length())
    
    return myLen


#%%

class pileup:
    
    def __init__(self, chrom, end, ref, bases):
        self.chrom = chrom
        self.end = end
        self.ref = ref
        self.bases = bases
        
        
#%%    
import csv

def dummy(file):
    with open(file) as pileupFile:
        
        pileupList = []
        myPileup = csv.reader(pileupFile, delimiter = '\t')
        for row in myPileup:
            pileupObject = pileup(row[0], row[1], row[3], row[4])
            pileupList.append(pileupObject)
    return pileupList






#%%

'''
Try adding multithreading
'''
import numpy as np
import pandas as pd
from multiprocessing import Pool

def parallelize_dataframe(df, func, n_cores = 4):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    myDF = pd.concat(pool.map(func, df_smplit))
    pool.close()
    pool.join()
    return myDF





def meanScore(df):
    
    mean = sum(df.loc[:, 'score'])/len(df)
    return mean






#%%
    
'''
try doing by chromosome so that you only have to find each chromosome one time 
'''

def quantPeaksByChrom(countsDF, bpDF, minSNPs = 1, bedOut = False, bedPrefix = 'out'):
    
    import numpy as np
    
    def quantPeak(bpRow):
        
        pileupBoolCoords1 = (tempDF['end_coord'] <= bpRow['end_coord'])
        pileupBoolCoords2 = (tempDF['end_coord'] >= bpRow['start_coord'])
        peakSNPs = tempDF[pileupBoolCoords1 & pileupBoolCoords2]     
        
        if (len(peakSNPs) >= minSNPs):
            ref_count = sum(np.array(peakSNPs['ref_count']))
            alt_count = sum(np.array(peakSNPs['alt_count']))  
            total_count = ref_count + alt_count
            if refAllele == 'mat':
                matFrac = np.divide(ref_count, total_count, where = total_count != 0)
            elif refAllele == 'pat':
                matFrac = np.divide(alt_count, total_count, where = total_count != 0)
        else:
            matFrac = 'nan'
            
        return matFrac    
    
    peaksResultsPerChrom = []
    
    for chrom in bpDF['chr'].unique():
        tempDFBool = (countsDF['chr'] == chrom)
        tempDF = countsDF[tempDFBool]
        
        tempBPBool = (bpDF['chr'] == chrom)
        tempBP = bpDF[tempBPBool]
        
        peakMatFrac = tempBP.apply(quantPeak, axis = 1)
        tempBP['peakMatFrac'] = peakMatFrac
        peaksResultsPerChrom.append(tempBP)
        
        print(f'Completed chromosome {chrom}')
    
    peaksResultsDF = pd.concat(peaksResultsPerChrom)
    
    if bedOut == True:
        
        filteredResultsBool = (peaksResultsDF['peakMatFrac'] != 'nan')
        filteredResults = peaksResultsDF[filteredResultsBool]
        bedDF = pd.DataFrame()
        bedDF['chr'] = filteredResults['chr']
        bedDF['start'] = filteredResults['start_coord']
        bedDF['end'] = filteredResults['end_coord']
        bedDF['peakMatFrac'] = filteredResults['peakMatFrac']
        
        bedDF.to_csv(f'{bedPrefix}.bedgraph', header = False, index = False, sep = '\t')        
        
    return peaksResultsDF

        

def matFracPeakDist(df):
    '''
    returns the mean maternal fraction of all called peaks
    '''
    peakMatFracs = [x for x in df.peakMatFrac if not isinstance(x, str)]
    mean = sum(peakMatFracs)/len(peakMatFracs)
    return mean


def matFracPeakHist(df, bins = 50):
    
    import matplotlib.pyplot as plt
    peakMatFracs = [x for x in df.peakMatFrac if not isinstance(x, str)]
    plt.hist(peakMatFracs, bins)
    plt.title('Distribution of peak mat prop')
    
#%%

'''test the above histogram for normality'''


