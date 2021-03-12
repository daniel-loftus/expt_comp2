# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 17:39:06 2020

@author: dal1858
"""

#%%



def pileupToCounts(pileupFile, snpFile):
    
    ''' 
    This function takes a pileup file as input (can be generated from a 
    bam file using samtools mpileup) and outputs a table containing allelic 
    counts genome-wide
    '''
#%%

import pandas as pd

pileupFile = 'SRR7685829Aligned.sortedByCoordout.pileup'
snpFile = 'intersect_test.bed'

pileup = pd.read_csv(pileupFile, sep = '\t', names = ['chr', 'end_coord', 'ref', 'reads', 'seq', 'qual'])
snps = pd.read_csv(snpFile, sep = '\t', names = ['chr', 'start_coord', 'end_coord', 'm/p'])


#%%

for chrom in pileup.loc[:, 'chr'].unique():

    '''
    I will perform the main loop chromosome by chromosome in order to improve 
    performance
    '''
    
    chromSNPsBool = snps['chr'].isin([chrom])
    chromPileupBool = pileup['chr'].isin([chrom])  
    
    chromSNPs = snps.loc[chromSNPsBool, ['end', 'm/p']]
    chromPileup = pileup.loc[chromPileupBool, ['coord', 'seq']]
    
    
    for i, row in chromPileup.iterrows():
        
        coord = row['coord']
        print(coord)
        bases = row['seq']



        ''' this much gives me the bases from the pileup file for each snp. now i just need to count them 
        
        try using different strategies to improve performance 
        
        https://engineering.upside.com/a-beginners-guide-to-optimizing-pandas-code-for-speed-c09ef2c6a4d6 ''' 
        
        
#%%
        
'''
optimize subsetting by chromosome
'''
        
def test(myPileup, mySNPs):
    
    chromSNPsBool = snps['chr'].isin([chrom])
    chromPileupBool = pileup['chr'].isin([chrom])  
    
    chromSNPs = snps.loc[chromSNPsBool, ['end', 'm/p']]
    chromPileup = pileup.loc[chromPileupBool, ['coord', 'seq']]

    return chromPileup.head()
    


#%%
    
'''
optimize looping through the dataframes 
'''

def loopTest(myPileup, mySNPs):

    for i, row in myPileup.iterrows():
    
        coord = row['coord']
        print(coord)
        #print(mySNPs.query('end == @coord')['m/p']  )

                
loopTest(pileup, snps)      
        
        
#%%
def test(myPileup, mySNPs):
    
	for chrom in pileup.loc[:, 'chr'].unique():

		chromSNPsBool = mySNPs['chr'].isin([chrom])
		chromPileupBool = myPileup['chr'].isin([chrom])  
		
		chromSNPs = mySNPs.loc[chromSNPsBool, ['end', 'm/p']]
		chromPileup = myPileup.loc[chromPileupBool, ['coord', 'seq']]

		print(chromPileup.head() )







#%%
        
'''try merging the snp and pileup dataframes'''

def mergeDFs(myPileup, mySNPs):
    
    '''
    merges the pileup and snps dataframes and returns only the needed columns
    '''
    
    snpsCols = ['chr', 'end_coord', 'm/p']
    pileupCols = ['chr', 'end_coord', 'reads', 'seq']
    
    mergedDF = pd.merge(myPileup[pileupCols], mySNPs[snpsCols], on = ['chr', 'end_coord'])
    
    return mergedDF
    
test = mergeDFs(pileup, snps)

#%%
refCount = []
altCount = []

for i, row in test.iterrows(): 
    
    basesTemp = row['seq']
    altBase = row['m/p'][2]
    
    ref = basesTemp.count('.') + basesTemp.count(',')
    alt = basesTemp.count(altBase)
    
    if (i % 10000) == 0:
        print(i)
        
#%%

testSub = test.iloc[0:100000, :]
        
#%%

        

def countRef(seq):
    
    ''' counts and returns the number of ref alleles overlapping a SNP '''
    
    count = seq.count('.') + seq.count(',')
    return count

def countAlt(row):
    
    ''' counts and returns the number of alt alleles overlapping a SNP '''
    
    altAllele = row['m/p'][2].upper()
    count = row['seq'].upper().count(altAllele)
    return count

testSub['ref_count'] = testSub['seq'].apply(countRef)
testSub['alt_count'] = testSub.apply(countAlt, axis = 1)
    

#%%




def addAlt(alleles):
    
    ''' adds a column containing the alt allele to the dataframe '''
    
    altAllele = alleles[2].upper()
    return altAllele
    
#test['alt_allele'] = df['m/p'].apply(addAlt) #SUPER FAST 
    


def countAlt2(row):
    
    ''' counts and returns the number of alt alleles overlapping a SNP '''
    
    altAllele = addAlt(row['alt_allele'])
    count = row['seq'].upper().count(altAllele)
    return count







#%%
    
def filterBases(bases):
    
    ''' 
    removes all bases before the last '$' character and after 
    the last '^' character, which designate the end and start of reads 
    and have a systemic reference bias
    '''
    
    if '$' in bases:
        bases = bases.split('$')[-1]
    if '^' in bases:
        bases = bases.split('^')[0]
        
    return bases
        



#%%
        
'''
A different appraoch: SQL instead of pandas 
'''

import sqlite3
from sqlite3 import Error

def create_connection(db_file):
    
    ''' create a database connection to a SQLite database ''' 
    
    conn = None
    
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e: 
        print(e)
    
    return conn
            
def create_table(conn, create_table_sql):
    
    ''' create a table from the create_table_sql statement ''' 
    
    




















#%%

def countsToBed(countsDF, minCount = 10):
    
    '''
    Takes a dataframe as generated by pileupToCounts and outputs a 
    bed file containing the allelic bias at each snp
    '''
    
    
    import pandas as pd 
    
    
    
    countFilteredDFBool = (countsDF['ref_count'] + countsDF['alt_count']) >= minCount
    countFilteredDF = countsDF[countFilteredDFBool]
    
    
    def findRefFrac(row):
        
        ''' calculates the maternal fraction from 0 to 1 at each snp '''
        
        refFrac = row['ref_count']/(row['ref_count'] + row['alt_count'])
        return refFrac
    
    bedDF = pd.DataFrame()  
    bedDF['chr'] = countFilteredDF['chr']
    bedDF['start'] = countFilteredDF['end_coord'].apply(lambda end: end - 1)
    bedDF['end'] = countFilteredDF['end_coord']
    bedDF['ref_frac'] = countFilteredDF.apply(findRefFrac, axis = 1)
    
    return bedDF


#%%
def countsToBed(countsDF, minCount = 10):
    
    '''
    Takes a dataframe as generated by pileupToCounts and outputs a 
    bed file containing the allelic bias at each snp
    '''
    
    
    import pandas as pd 
    import numpy as np
    
    
    countFilteredDFBool = (countsDF['ref_count'] + countsDF['alt_count']) >= minCount
    countFilteredDF = countsDF[countFilteredDFBool]
    
    ref_count = np.array(countFilteredDF['ref_count'])
    alt_count = np.array(countFilteredDF['alt_count'])  
    total_count = ref_count + alt_count
    refFrac = np.divide(ref_count, total_count, where = total_count != 0)
    
    bedDF = pd.DataFrame()  
    bedDF['chr'] = countFilteredDF['chr']
    bedDF['start'] = countFilteredDF['end_coord'].apply(lambda end: end - 1)
    bedDF['end'] = countFilteredDF['end_coord']
    bedDF['ref_frac'] = refFrac
    
    return bedDF

















