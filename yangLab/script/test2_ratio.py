# The script is used to calculate the ratio of each elements on hg19 genome

import pandas as pd
import numpy as np


def overlap_check(region1, region2):  # region1 or 2 is the list of start site and end site
    if (region1[0] > region2[1]) | (region1[1] < region2[0]):
        return 0
    else:
        return 1


def merge(region1, region2):  # Require the overlap_check!
    return [min(region1[0], region2[0]), max(region1[1], region2[1])]


def inter(region1, region2):  # Require the overlap_check!
    int_list = region1 + region2
    int_list.sort()
    return [int_list[1], int_list[2]]


def MultiRegionMerge(rand_df, start_site, end_site):
    # The rand_df refer to the dataframe contain multi-region; start_site and end_site refer to
    # the column names of start and end sites.
    if rand_df.empty:
        return rand_df
    sort_df = rand_df.sort_values(by=[start_site, end_site])
    sort_df.index = range(len(sort_df))
    merged_region_list = [[sort_df.loc[0, start_site], sort_df.loc[0, end_site]]]
    id_merged_region = 0
    for i in range(len(sort_df) - 1):
        next_region = [sort_df.loc[i + 1, start_site], sort_df.loc[i + 1, end_site]]
        if overlap_check(merged_region_list[id_merged_region], next_region) == 1:
            merged_region_list[id_merged_region] = merge(merged_region_list[id_merged_region], next_region)
        else:
            merged_region_list.append(next_region)
            id_merged_region += 1
    new_chrom_df = pd.DataFrame(np.array(merged_region_list), columns=[start_site, end_site])
    return new_chrom_df


def lenth_ratio(chrom_df, region_df, start_site, end_site):
    whole_lenth = sum(chrom_df['lenth'])
    region_lenth = sum(region_df[end_site] - region_df[start_site])
    return [region_lenth, region_lenth / whole_lenth]


def parseExon(refflat_line):  # Used to parse the exonStarts and exonEnds in each line. Return a dataframe
    ex_region_list = []
    ex_sslr = refflat_line['exonStarts'].split(',')
    ex_sslr.remove('')
    ex_eslr = refflat_line['exonEnds'].split(',')
    ex_eslr.remove('')
    ex_start_site_list = [int(x) for x in ex_sslr]
    ex_end_site_list = [int(x) for x in ex_eslr]
    for x, y in zip(ex_start_site_list, ex_end_site_list):
        ex_region_list.append([x, y])
    ex_region_df = pd.DataFrame(np.array(ex_region_list), columns=['ExonStart', 'ExonEnd'])
    # mer_ex_region_df = MultiRegionMerge(ex_region_df, 'ExonStart', 'ExonEnd')
    return ex_region_df


def parseIntron(refflat_line):  # Used to parse the intronStarts and intronEnds in each line. Return a dataframe
    in_region_list = []
    empty_df = pd.DataFrame()
    ex_region_df = parseExon(refflat_line)
    tss = refflat_line['txStart']
    tts = refflat_line['txEnd']
    in_start_site_list = [tss]
    in_start_site_list.extend(list(ex_region_df['ExonEnd']))
    in_end_site_list = list(ex_region_df['ExonStart'])
    in_end_site_list.append(tts)
    for x, y in zip(in_start_site_list, in_end_site_list):
        if x != y:
            in_region_list.append([x, y])
    if len(in_region_list) == 0:
        return empty_df
    in_region_df = pd.DataFrame(np.array(in_region_list), columns=['IntronStart', 'IntronEnd'])
    return in_region_df


def parseUTR5(refflat_line):  # Used to parse the 5UTR for each line in the refflat file. Return a dataframe.
    utr5_region_list = []
    ex_region_df = parseExon(refflat_line)
    cds_start = refflat_line['cdsStart']
    cds_end = refflat_line['cdsEnd']
    empty_df = pd.DataFrame()
    if cds_start == cds_end:
        return empty_df
    for i in range(len(ex_region_df)):
        exon_start = ex_region_df.loc[i, 'ExonStart']
        exon_end = ex_region_df.loc[i, 'ExonEnd']
        if exon_start < cds_start:
            if exon_end < cds_start:
                utr5_region_list.append([exon_start, exon_end])
            else:
                utr5_region_list.append([exon_start, cds_start])
                break
        else:
            return empty_df
    utr5_region_df = pd.DataFrame(np.array(utr5_region_list), columns=['UTR5Start', 'UTR5End'])
    return utr5_region_df


def parseUTR3(refflat_line):  # Used to parse the 5UTR for each line in the refflat file. Return a dataframe.
    utr3_region_list = []
    ex_region_df = parseExon(refflat_line)
    revsort_ex_region_df = ex_region_df.sort_values(by=['ExonEnd', 'ExonStart'], ascending=False)
    revsort_ex_region_df.index = range(len(revsort_ex_region_df))
    cds_start = refflat_line['cdsStart']
    cds_end = refflat_line['cdsEnd']
    empty_df = pd.DataFrame()
    if cds_start == cds_end:
        return empty_df
    for i in range(len(revsort_ex_region_df)):
        exon_start = revsort_ex_region_df.loc[i, 'ExonStart']
        exon_end = revsort_ex_region_df.loc[i, 'ExonEnd']
        if exon_end > cds_end:
            if exon_start > cds_end:
                utr3_region_list.append([exon_start, exon_end])
            else:
                utr3_region_list.append([cds_end, exon_end])
                break
        else:
            return empty_df
    utr3_region_df = pd.DataFrame(np.array(utr3_region_list), columns=['UTR3Start', 'UTR3End'])
    return utr3_region_df


def parseCDS(refflat_line):
    cds_region_list = []
    empty_df = pd.DataFrame()
    cds_start = refflat_line['cdsStart']
    cds_end = refflat_line['cdsEnd']
    if cds_start == cds_end:
        return empty_df
    ex_region_df = parseExon(refflat_line)
    for i in range(len(ex_region_df)):
        cds_region = [cds_start, cds_end]
        exon_region = [ex_region_df.loc[i, 'ExonStart'], ex_region_df.loc[i, 'ExonEnd']]
        if overlap_check(cds_region, exon_region) == 1:
            cds_region_list.append(inter(cds_region, exon_region))
    cds_region_df = pd.DataFrame(np.array(cds_region_list), columns=['CDSStart', 'CDSEnd'])
    return cds_region_df

def work_flow(refflat_df, chrom_df, element, parseFun):
    # The work flow contain the intervals integration and merging, the dataframe generation and ratio cal
    integrated_region_df = pd.DataFrame()
    for chrom in chrom_df['chrom']:
        one_chrom_df = pd.DataFrame(refflat_df[refflat_df['chrom'].isin([chrom])])
        if one_chrom_df.empty:
            continue
        one_chrom_df.index = range(len(one_chrom_df))
        ele_region_df = pd.DataFrame()
        for i in range(len(one_chrom_df)):
            tmp_ele_region_df = parseFun(one_chrom_df.loc[i])
            if not tmp_ele_region_df.empty:
                ele_region_df = pd.concat([ele_region_df, tmp_ele_region_df])
        if not ele_region_df.empty:
            mer_region_df = MultiRegionMerge(ele_region_df, '{0}Start'.format(element), '{0}End'.format(element))
            mer_region_df.loc[:, 'chrom'] = chrom
            integrated_region_df = pd.concat([integrated_region_df, mer_region_df])
    integrated_region_df.to_csv(r'/home/liuzhen2018/data/test/result/test2/{0}s.txt'.format(element),
                                sep='\t', index=None)
    lr_list = lenth_ratio(chrom_df, integrated_region_df, '{0}Start'.format(element), '{0}End'.format(element))
    print(['{0}s'.format(element), lr_list[0], lr_list[1]])
    return ['{0}s'.format(element), lr_list[0], lr_list[1]]


def firstExonsCal(refflat_df, chrom_df):  # Return a list [Region_of_hg19, Length, Ratio]
    integrated_region_df = pd.DataFrame()
    for chrom in chrom_df['chrom']:
        one_chrom_df = pd.DataFrame(refflat_df[refflat_df['chrom'].isin([chrom])])
        if one_chrom_df.empty:
            continue
        one_chrom_df.index = range(len(one_chrom_df))
        firExonStart_list = []
        firExonEnd_list = []
        for i in range(len(one_chrom_df)):
            firExonStart_list.append(int(one_chrom_df.loc[i, 'exonStarts'].split(',')[0]))
            firExonEnd_list.append(int(one_chrom_df.loc[i, 'exonEnds'].split(',')[0]))
        one_chrom_df.loc[:, 'firExonStart'] = firExonStart_list
        one_chrom_df.loc[:, 'firExonEnd'] = firExonEnd_list
        new_chrom_df = MultiRegionMerge(one_chrom_df, 'firExonStart', 'firExonEnd')
        new_chrom_df.loc[:, 'chrom'] = chrom
        integrated_region_df = pd.concat([integrated_region_df, new_chrom_df])
    integrated_region_df.to_csv(r'/home/liuzhen2018/data/test/result/test2/firstExons.txt', sep='\t',
                                index=None)
    lr_list = lenth_ratio(chrom_df, integrated_region_df, 'firExonStart', 'firExonEnd')
    print(['FirstExons', lr_list[0], lr_list[1]])
    return ['FirstExons', lr_list[0], lr_list[1]]


def exonsCal(refflat_df, chrom_df):
    return work_flow(refflat_df, chrom_df, 'Exon', parseExon)


def intronsCal(refflat_df, chrom_df):
    return work_flow(refflat_df, chrom_df, 'Intron', parseIntron)


def utr5Cal(refflat_df, chrom_df):
    return work_flow(refflat_df, chrom_df, 'UTR5', parseUTR5)


def utr3Cal(refflat_df, chrom_df):
    return work_flow(refflat_df, chrom_df, 'UTR3', parseUTR3)


def cdsCal(refflat_df, chrom_df):
    return work_flow(refflat_df, chrom_df, 'CDS', parseCDS)


refflat_columns = ['geneName', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
                   'exonStarts', 'exonEnds']
refflat_df = pd.read_table(r'/home/liuzhen2018/data/test/refFlat.txt', sep='\t', header=None,
                           names=refflat_columns)
chrom_df = pd.read_table(r'/home/liuzhen2018/data/test/chrom.size', sep='\t', header=None,
                         names=['chrom', 'lenth'])
result_list = []
result_list.append(firstExonsCal(refflat_df, chrom_df))
result_list.append(exonsCal(refflat_df, chrom_df))
result_list.append(intronsCal(refflat_df, chrom_df))
result_list.append(utr5Cal(refflat_df, chrom_df))
result_list.append(utr3Cal(refflat_df, chrom_df))
result_list.append(cdsCal(refflat_df, chrom_df))
result_df = pd.DataFrame(np.array(result_list), columns=['Region_of_hg19', 'Length', 'Ratio'])
result_df.to_csv(r'/home/liuzhen2018/data/test/result/test2/ratio_result.txt', sep='\t', index=None)
