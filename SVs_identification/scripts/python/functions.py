import pandas as pd
import io
import csv
import matplotlib.pyplot as plt
import math


def read_delly_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    df[['PRECISE', 'EXTRA']] = df['INFO'].str.split(';', 1, expand=True)
    df = df.drop(['EXTRA'], axis=1)
    return df


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def read_gtf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('#!')]
    genes = pd.read_csv(io.StringIO(''.join(lines)), header=None,
                        names=['seqname', 'source', 'feature', 'start', 'end', 'score',
                               'strand', 'frame', 'attributes'],
                        dtype={'seqname': str, 'source': str, 'feature': str, 'start': int, 'end': int, 'score': str,
                               'strand': str,
                               'frame': str, 'attributes': str}, sep='\t')

    # genes_split = genes['attributes'].str.split(";", expand=True)
    # df = pd.DataFrame(genes_split)
    # df2 = pd.DataFrame()
    # df2[0,['col_name', 'value']] = df.iloc[0, 0].split(' ')
    # print(df2.iloc[0,:])

    #for row in range(len(df)):
    #    for column in range(len(df.columns)):
    #        df2[row][['col_name', 'value']] = df.iloc[row,column].split(' ')

    #print(df2)
    return genes


def filter_precise(df):
    return df.loc[df['PRECISE'] == 'PRECISE']


def filter_qual(df):
    return df.loc[df['FILTER'] == 'PASS']


def read_survivor_stats(path):
    file = open(path)
    tsv = csv.reader(file, delimiter="\t")
    df = pd.DataFrame(tsv)
    df.columns = df.loc[0, :]
    new_df = df.loc[1:5, :]
    final_df = new_df.astype({'Len': str, 'Del': int, 'Dup': int, 'Inv': int, 'INS': int, 'TRA': int, 'UNK': int})
    return final_df


def get_genome_name(path):
    split = path.split('/')
    split2 = split[-1].split('.')
    sample_name = split2[0]
    return sample_name


def summarize_stats(df):
    return [df['Del'].sum(), df['Dup'].sum(), df['Inv'].sum(), df['INS'].sum(), df['TRA'].sum(), df['UNK'].sum()]


def plot_stats(sum_stats, sample_name):
    x_axis = ['DEL', 'DUP', 'INV', 'INS', 'TRA', 'UNK']
    plt.plot(x_axis, sum_stats, '.', label=sample_name)


def count_sv_type(alt_column):
    count_del = 0
    count_dup = 0
    count_ins = 0
    count_inv = 0
    count_unk = 0
    count_tra = 0
    for j in alt_column:
        if j == '<DEL>':
            count_del = count_del + 1
        elif j == '<DUP>':
            count_dup = count_dup + 1
        elif j == '<INS>':
            count_ins = count_ins + 1
        elif j == '<INV>':
            count_inv = count_inv + 1
        elif j == '<TRA>':
            count_tra = count_tra + 1
        else:
            count_unk = count_unk + 1

    return [count_del, count_dup, count_inv, count_ins, count_tra, count_unk]


# def find_overlap(vcf_list, vcf_files):
#     joined = pd.DataFrame(columns=['CHR_POS'])
#     for i in range(len(vcf_list)):
#         chrom_col = vcf_list[i][['CHROM']]['CHROM'].astype(str)
#         pos_col = vcf_list[i][['POS']]['POS'].astype(str)
#         combined = pd.DataFrame()
#         combined['CHR_POS'] = chrom_col.str.cat(pos_col, sep=":")
#         if i == 0:
#             colname = 'CHR_POS' + get_genome_name(vcf_files[i])
#             joined = pd.DataFrame({colname: combined['CHR_POS']})
#         else:
#             vcf = pd.DataFrame({'CHR_POS': combined['CHR_POS']})
#             joined = joined.join(vcf, how='outer', lsuffix=get_genome_name(vcf_files[i - 1]),
#                                  rsuffix=get_genome_name(vcf_files[i]))
#
#     print(joined.to_string())
#     index = joined.columns
#     genome = []
#     for row in index:
#         split = row.split('CHR_POS')
#         genome.append(split[-1])
#
#     # rows = len(df)
#     # cols = len(df.columns)
#     # results = pd.DataFrame(columns=['CHR', 'POS', 'COUNT', 'SAMPLES'])
#     # samples = []
#     # for row in range(rows):
#     #     count = 0
#     #     chr = df[row]['CHROM']
#     #     pos = df[row]['POS']
#     #     print(df.iloc[1,:])
#     #     for col in range(cols):
#     #         if df.iloc[row, col] != 'NaN':
#     #             count = count + 1
#     #             #samples.append(sample)
#     #
#     #     results[row]['CHR'] = chr
#     #     results[row]['POS'] = pos
#     #     results[row]['COUNT'] = count
#     #     results[row]['SAMPLES'] = samples


def count_svs(vcf_list, vcf_files):
    counts = {}
    genomes_pos = {}
    all_info = {}
    for i in range(len(vcf_list)):
        genome = get_genome_name(vcf_files[i])
        for row in range(len(vcf_list[i])):
            chrom = vcf_list[i].iloc[row]['CHROM']
            pos = vcf_list[i].iloc[row]['POS'].astype(str)
            chrom_pos = chrom + ':' + pos
            if chrom_pos in counts.keys():
                counts[chrom_pos] += 1
                genomes_pos[chrom_pos] += ';' + genome
                all_info[chrom_pos] += ';' + vcf_list[i].iloc[row].astype(str)
            else:
                counts[chrom_pos] = 1
                genomes_pos[chrom_pos] = genome
                all_info[chrom_pos] = vcf_list[i].iloc[row].astype(str)
    return counts, genomes_pos, all_info

