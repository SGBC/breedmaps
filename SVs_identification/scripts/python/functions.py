import pandas as pd
import io
import csv

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


def read_survivor_stats(path):
    file = open(path)
    tsv = csv.reader(file, delimiter="\t")
    df = pd.DataFrame(tsv)
    df.columns = df.loc[0, :]
    new_df = df.loc[1:5, :]
    final_df = new_df.astype({'Len': str, 'Del': int, 'Dup': int, 'Inv': int, 'INS': int, 'TRA': int, 'UNK': int})
    return final_df


def summarize_stats(df):
    return [df['Del'].sum(), df['Dup'].sum(), df['Inv'].sum(), df['INS'].sum(), df['TRA'].sum(), df['UNK'].sum()]

