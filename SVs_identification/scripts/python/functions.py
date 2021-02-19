import pandas as pd
import io
import csv
import matplotlib.pyplot as plt


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


def get_sample_name(path):
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