#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

import argparse 

parser = argparse.ArgumentParser()

required = parser.add_argument_group('Required arguments')
required.add_argument('--pileup', 
                    type = str, 
                    help = 'Pileup file', 
                    required = True)
required.add_argument('--snps', 
                    type = str, 
                    help = 'BED file containing SNPs', 
                    required = True)


SNPs_quantification = parser.add_argument_group('SNPs quantification arguments')
SNPs_quantification.add_argument('--quantSNPs', 
                      help = 'Whether to output a bedgraph file containing the maternal proportion at each SNP (default = True)',
                      default = True, 
                      choices = [True, False], 
                      required = False)
SNPs_quantification.add_argument('--refAllele', 
                      help = 'Whether the reference genome represents the maternal or paternal genome', 
                      default = 'mat', 
                      choices = ['mat', 'pat'],
                      required = False)
SNPs_quantification.add_argument('--SNPsMinCount', 
                      type = int, 
                      help = 'Minimum coverage for a SNP to be included in the bed file', 
                      default = 10, 
                      required = False)


features_quantification = parser.add_argument_group('Features quantification arguments')
features_quantification.add_argument('--quantFeatures', 
                                     help = 'Whether to output a bed file containing the maternal proportion at each feature specified in the featuresBed file (default = False)', 
                                     default = False, 
                                     choices = ['True', 'False'], 
                                     required = False)
features_quantification.add_argument('--features', 
                                     type = str,
                                     help = 'Bed file containing genomic features to quantify', 
                                     required = False)
features_quantification.add_argument('--featuresMinSNPs', 
                                     type = int,
                                     help = 'Minimum number of SNPs within a genomic feature for it to be included in quantFeatures', 
                                     default = 1, 
                                     required = False)


other = parser.add_argument_group('Other arguments')
other.add_argument('--outPrefix', 
                      type = str, 
                      help = 'Prefix of output files', 
                      default = 'out', 
                      required = False)

args = parser.parse_args()




def main(pileupFile = args.pileup, 
         snpFile = args.snps, 
         quantSNPs = args.quantSNPs,
         quantFeatures = args.quantFeatures,
         refAllele = args.refAllele, 
         SNPsMinCount = args.SNPsMinCount, 
         featuresMinSNPs = args.featuresMinSNPs,
         outPrefix = args.outPrefix):
    
    ''' 
    Takes a pileup file (can be generated from a 
    bam file using samtools mpileup) and a bed file containing snps
    as input and outputs a table containing 
    allelic counts genome-wide
    '''
    
    
    
    #read in pileup and snp-containing bed file 
    pileup = pd.read_csv(pileupFile, 
                         sep = '\t', 
                         names = ['chr', 'end_coord', 'pileup_ref', 'reads', 'bases', 'qual'], 
                         dtype = {'chr': object, 'end_coord': int, 'pileup_ref': object, 'reads': object, 'bases': str, 'qual': str})
    print('Pileup file loaded', flush = True)
    
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
    def countsToBed(countsDF, minCount = SNPsMinCount):
        
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

    if quantSNPs == True:
        print('Writing bed file')
        bed = countsToBed(df)
        bed.to_csv(f'{outPrefix}.bedgraph', header = False, index = False, sep = '\t')
    
    
    
    def quantFeaturesByChrom(countsDF, bpDF, minSNPs = featuresMinSNPs, bedPrefix = 'out'):
        
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
            
            tempBP['peakMatFrac'] = tempBP.apply(quantPeak, axis = 1)
            peaksResultsPerChrom.append(tempBP)
                    
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

    if quantFeatures == True:
        bed = quantFeaturesByChrom(df)
        bed.to_csv(f'{outPrefix}_features.bedgraph')

    
    print('Done')


if __name__ == '__main__':
    main() 