#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse 

parser = argparse.ArgumentParser()

required = parser.add_argument_group('Required arguments')
required.add_argument('--pileup', 
                    type = str, 
                    help = 'Pileup file')
required.add_argument('--snps', 
                    type = str, 
                    help = 'BED file containing SNPs')

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('--outputBed', 
                      help = 'Whether to output a bed file (default = True)',
                      default = True, 
                      choices = [True, False], 
                      required = False)
optional.add_argument('--refAllele', 
                      help = 'Whether the reference genome represents the maternal or paternal genome', 
                      default = 'mat', 
                      choices = ['mat', 'pat'],
                      required = False)
optional.add_argument('--bedMinCount', 
                      type = int, 
                      help = 'Minimum coverage for a SNP to be included in the bed file', 
                      default = 10, 
                      required = False)
optional.add_argument('--bedOutPrefix', 
                      type = str, 
                      help = 'Prefix of output bed file', 
                      default = 'out', 
                      required = False)

args = parser.parse_args()




def main(pileupFile = args.pileup, 
         snpFile = args.snps, 
         outputBed = args.outputBed, 
         refAllele = args.refAllele, 
         bedMinCount = args.bedMinCount, 
         bedOutPrefix = args.bedOutPrefix):
    
    ''' 
    Takes a pileup file (can be generated from a 
    bam file using samtools mpileup) and a bed file containing snps
    as input and outputs a table containing 
    allelic counts genome-wide
    '''
    
    
    import pandas as pd
    
    
    
    #read in pileup and snp-containing bed file 
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
    
    
    
    #merge into a single dataframe 
    def mergeDFs(myPileup, mySNPs):
        
        ''' merges the pileup and snps dataframes and returns only the needed columns '''
        
        snpsCols = ['chr', 'end_coord', 'ref/alt']
        pileupCols = ['chr', 'end_coord', 'pileup_ref', 'bases']
        
        mergedDF = pd.merge(myPileup[pileupCols], mySNPs[snpsCols], on = ['chr', 'end_coord'])
        
        return mergedDF
    
    df = mergeDFs(pileup, snps)
    print('Files merged')
    
 
    
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
    
    
    
    #generate bed-like dataframe 
    def countsToBed(countsDF, minCount = bedMinCount):
        
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
        if refAllele == 'mat':
            matFrac = np.divide(ref_count, total_count, where = total_count != 0)
        elif refAllele == 'pat':
            matFrac = np.divide(alt_count, total_count, where = total_count != 0)
        
        bedDF = pd.DataFrame()  
        bedDF['chr'] = countFilteredDF['chr']
        bedDF['start'] = countFilteredDF['end_coord'].apply(lambda end: end - 1)
        bedDF['end'] = countFilteredDF['end_coord']
        bedDF['mat_frac'] = matFrac
        
        return bedDF

    if outputBed == True:
        print('Writing bed file')
        bed = countsToBed(df)
        bed.to_csv(f'{bedOutPrefix}.bedgraph', header = False, index = False, sep = '\t')
    
    
    
    print('Done')


if __name__ == '__main__':
    main() 