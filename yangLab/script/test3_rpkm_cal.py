# -*- coding: utf-8 -*-
# The script is used to calculate the RPKM of transcript positions based on 
# the specific bam file. The reads in bam file are only mapped on one position.

"""
Created on Tue Jul 31 15:44:55 2018

@author: liuzhen
"""

import pandas as pd
import numpy as np
import pysam

def parseExon(refgene_line):
    """
    Used to parse the exonStarts and exonEnds in each line.
    Return a list with tuple of each exon region.
    """
    ex_region_list = []
    ex_sslr = refgene_line['exonStarts'].split(',')
    ex_sslr.remove('')
    ex_eslr = refgene_line['exonEnds'].split(',')
    ex_eslr.remove('')
    ex_start_site_list = [int(x) for x in ex_sslr]
    ex_end_site_list = [int(x) for x in ex_eslr]
    for x, y in zip(ex_start_site_list, ex_end_site_list):
        ex_region_list.append((x, y))
    return ex_region_list

def overlap_check(region1, region2):  # region1 or 2 is the list of start site and end site
    if (region1[0] > region2[1]) | (region1[1] < region2[0]):
        return 0
    else:
        return 1

def contain_check(region_list1, region_list2):
    """
    Check whether region_list1 contain region_list2
    The region_list is list of tuples
    """
    for region2 in region_list2:
        flag = 0
        for region1 in region_list1:
            if overlap_check(region2, region1) == 0:
                continue
            else:
                if (region2[0] >= region1[0]) & (region2[1] <= region2[1]):
                    flag = 1
                    break
                else:
                    return 0
        if flag == 0:
            return 0
    return 1
        

def readCount(ex_region_list, samfetch):
    """
    Count the number of reads mapped to the list of exon region.
    The exon region and reads region is 0-based.
    The read is counted when the overlap happens in the exon region.
    The counted reads required to be contained in the exon regions.
    The reads beyond the exon regions will not be counted.
    """
    counts = 0
    for read in samfetch:
        if contain_check(ex_region_list, read.get_blocks()) == 1:
            counts += 1
    return counts

samfile_path = '/Users/liuzhen/intern/data/test/tophat_map_g1/accepted_hits.bam'
refgene_path = '/Users/liuzhen/intern/data/test/refGene.txt'
result_path = '/Users/liuzhen/intern/data/test/RPKM_result.txt'
refgene_column = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 
                  'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 
                  'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
samfile = pysam.AlignmentFile(samfile_path, 'rb')
refgene_df = pd.read_table(refgene_path, sep='\t', header=None, 
                           names=refgene_column)

result_df = refgene_df.loc[:, ['name2', 'name', 'chrom', 'txStart']]
mapped_reads_amount = samfile.count()
RPKM_list = []
for i in range(len(refgene_df)):
    samfetch = samfile.fetch(refgene_df.loc[i, 'chrom'], 
                             refgene_df.loc[i, 'txStart'], 
                             refgene_df.loc[i, 'txEnd'])
    ex_region_list = parseExon(refgene_df.loc[i])
    tmp_matrix = np.array(ex_region_list)
    exon_length = sum(tmp_matrix[:, 1] - tmp_matrix[:, 0])
    read_counts = readCount(ex_region_list, samfetch)
    RPKM = read_counts * 10**9 / (mapped_reads_amount * exon_length)
    RPKM_list.append(RPKM)
result_df.loc[:, 'RPKM'] = RPKM_list
result_df.to_csv(result_path, sep='\t', index=None)





