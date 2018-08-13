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

def merge(region1, region2):  # Require the overlap_check!
    return (min(region1[0], region2[0]), max(region1[1], region2[1]))

def multiRegionMerge(region_list):
    # The rand_df refer to the dataframe contain multi-region; start_site and end_site refer to
    # the column names of start and end sites.
    rand_df = pd.DataFrame(region_list, columns=['start', 'end'])
    if rand_df.empty:
        return []
    sort_df = rand_df.sort_values(by=['start', 'end'])
    sort_df.index = range(len(sort_df))
    merged_region_list = [(sort_df.loc[0, 'start'], sort_df.loc[0, 'end'])]
    id_merged_region = 0
    for i in range(len(sort_df) - 1):
        next_region = (sort_df.loc[i + 1, 'start'], sort_df.loc[i + 1, 'end'])
        if overlap_check(merged_region_list[id_merged_region], next_region) == 1:
            merged_region_list[id_merged_region] = merge(merged_region_list[id_merged_region], next_region)
        else:
            merged_region_list.append(next_region)
            id_merged_region += 1
    return merged_region_list

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

def mergeRefGeneDf(refgene_df):
    """
    Used to merge the refGene dataframe by gene symbel.
    A large dictionary will be returned with the structure like:
        {name2:[['chr1', [transcript_ID], [transcript_len], [exon_region], (TSS, TTS)]
                ['chr2', [transcript_ID], [transcript_len], [exon_region], (TSS, TTS)]
                ...
                ['chrX', [transcript_ID], [transcript_len], [exon_region], (TSS, TTS)]]
         ...
         name2:...}
    """
    sort_df = refgene_df.sort_values(by=['name2', 'chrom'])
    sort_df.index = range(len(sort_df))
    merge_dict = {}
    i = 0
    while i < len(sort_df):
        sym_span = 1
        while (i+sym_span < len(sort_df)) and (sort_df.loc[i+sym_span-1, 'name2'] == sort_df.loc[i+sym_span, 'name2']):
            sym_span += 1
        chrom_ind = i
        gene_merge_list = []
        while chrom_ind < i+sym_span:
            chrom_span = 1
            while (chrom_ind+chrom_span < i+sym_span) and (sort_df.loc[chrom_ind+chrom_span-1, 'chrom'] == sort_df.loc[chrom_ind+chrom_span, 'chrom']):
                chrom_span += 1
            trans_id_list = []
            trans_len_list = []
            ex_region_list = []
            for trans_ind in range(chrom_ind, chrom_ind+chrom_span):
                trans_id_list.append(sort_df.loc[trans_ind, 'name'])
                trans_ex_list = parseExon(sort_df.loc[trans_ind])
                tmp_matrix = np.array(trans_ex_list)
                ex_region_list.extend(trans_ex_list)
                trans_len_list.append(sum(tmp_matrix[:, 1] - tmp_matrix[:, 0]))
            merge_ex_list = multiRegionMerge(ex_region_list)
            if '_' not in sort_df.loc[chrom_ind, 'chrom']:
                gene_merge_list.append([sort_df.loc[chrom_ind, 'chrom'],
                                        trans_id_list, trans_len_list, merge_ex_list, 
                                        (merge_ex_list[0][0], merge_ex_list[-1][-1])])
            chrom_ind += chrom_span
        if len(gene_merge_list) != 0:
            merge_dict[sort_df.loc[i, 'name2']] = gene_merge_list
        i += sym_span
    return merge_dict

    
samfile_path = '/Users/liuzhen/intern/data/test/tophat_map_g1/accepted_hits.bam'
refgene_path = '/Users/liuzhen/intern/data/test/refGene.txt'
result_path = '/Users/liuzhen/intern/data/test/RPKM_gene_dpsfd_avglen_result.txt'
refgene_column = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 
                  'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 
                  'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
samfile = pysam.AlignmentFile(samfile_path, 'rb')
refgene_df = pd.read_table(refgene_path, sep='\t', header=None, 
                           names=refgene_column)

mapped_reads_amount = samfile.count()
merge_dict = mergeRefGeneDf(refgene_df)
result_list = []
for gene in merge_dict:
    read_counts = 0
    trans_len = 0
    trans_num = 0
    trans_id_list = []
    tss_list = []
    chrom_list = []
    for each_chrom in merge_dict[gene]:
        samfetch = samfile.fetch(each_chrom[0], each_chrom[4][0], each_chrom[4][1])
        read_counts += readCount(each_chrom[3], samfetch)
        trans_num += len(each_chrom[2])
        trans_len += sum(each_chrom[2])
        trans_id_list.extend(each_chrom[1])
        for each_trans in each_chrom[1]:
            chrom_list.append(each_chrom[0])
    avg_trans_len = trans_len / trans_num
    RPKM = read_counts * 10**9 / (mapped_reads_amount * avg_trans_len)
    result_list.append([gene, RPKM, ','.join(trans_id_list), ','.join(chrom_list)])
result_df = pd.DataFrame(result_list, columns=['gene', 'RPKM', 'transcript_ID', 'chrom'])
result_df.to_csv(result_path, sep='\t', index=None, )


    
    





