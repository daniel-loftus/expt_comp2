# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 15:04:30 2020

@author: dal1858
"""
#%%
def pileupToCounts(pileupFile, snpFile):
    
    
    '''
    NEED TO ADD: 
        Remove bases that come after '^~'. They are (almost?) always reference 
    '''
    ''' 
    This function takes a pileup file (can be generated from a 
    bam file using samtools mpileup) and a bed file containing snps
    as input and outputs a table containing 
    allelic counts genome-wide
    '''
    
    import pandas as pd
    
    
    pileup = pd.read_csv(pileupFile, 
                         sep = '\t', 
                         names = ['chr', 'end_coord', 'ref', 'reads', 'bases', 'qual'], 
                         dtype = {'chr': object, 'end_coord': int, 'ref': str, 'reads': int, 'bases': str, 'qual': str})
    print('Pileup file loaded')
    
    snps = pd.read_csv(snpFile, 
                       sep = '\t', 
                       names = ['chr', 'start_coord', 'end_coord', 'm/p'], 
                       dtype = {'chr': object, 'start_coord': int, 'end_coord': int, 'm/p': str})
    print('SNPs file loaded')
    
    
    def mergeDFs(myPileup, mySNPs):
        
        ''' merges the pileup and snps dataframes and returns only the needed columns '''
        
        snpsCols = ['chr', 'end_coord', 'm/p']
        pileupCols = ['chr', 'end_coord', 'reads', 'bases']
        
        mergedDF = pd.merge(myPileup[pileupCols], mySNPs[snpsCols], on = ['chr', 'end_coord'])
        
        return mergedDF
    
    
    df = mergeDFs(pileup, snps)
    print('Dataframes merged')
    
    
    def countRef(bases):
        
        ''' counts and returns the number of ref alleles overlapping a SNP '''
        
        count = bases.count('.') + bases.count(',')
        return count
    
    
    def countAlt(row):
        
        ''' counts and returns the number of alt alleles overlapping a SNP '''
        
        altAllele = row['m/p'][2].upper()
        count = row['bases'].upper().count(altAllele)
        return count
    
    
    df['ref_count'] = df['bases'].apply(countRef)
    print('Ref alleles counted')
    df['alt_count'] = df.apply(countAlt, axis = 1)
    print('Alt alleles counted')


    return df 