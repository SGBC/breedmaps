import pandas as pd
import io
import matplotlib.pyplot as plt


def read_delly_vcf(path):
    """Reads VCF files but specifically modified for VCF files created by DELLY. Returns a dataframe"""
    # Extracts the PRECISE/IMPRECISE flag and the END column from the INFO field and saves it as seperate columns
    with open(path, 'r') as f:
        # Read lines that is not a header line
        lines = [l for l in f if not l.startswith('##')]
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    # Save the PRECISE/IMPRECISE flag as an additional column
    df[['PRECISE', 'EXTRA']] = df['INFO'].str.split(';', 1, expand=True)
    df[['1', '2', 'END_full', '3']] = df['EXTRA'].str.split(';', 3, expand=True)
    df[['4', 'END']] = df['END_full'].str.split('=', 1, expand=True)
    df = df.drop(['EXTRA', '1', '2', '3', '4', 'END_full'], axis=1)
    return df


def read_gtf(path):
    """Reads a gene annotation file, in the format GTF, and returns it as a data frame"""
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('#!')]
    genes = pd.read_csv(io.StringIO(''.join(lines)), header=None,
                        names=['seqname', 'source', 'feature', 'start', 'end', 'score',
                               'strand', 'frame', 'attributes'],
                        dtype={'seqname': str, 'source': str, 'feature': str, 'start': int, 'end': int, 'score': str,
                               'strand': str,
                               'frame': str, 'attributes': str}, sep='\t')
    return genes


def filter_precise(df):
    """Filter VCF data frame on the flag PRECISE(will remove IMPRECISE). Specific for DELLY"""
    # Saves only the variants with the PRECISE flag. IMPRECISE are filtered out
    return df.loc[df['PRECISE'] == 'PRECISE']


def filter_qual(df):
    """Filter on the FILTER coulmn in VCF files. Common for VCF files"""
    # Saves only the variants that passed the vcf filters. LowQual are filtered out
    # Poor quality and insufficient number of PEs(paired ends) and SRs(split reads).
    return df.loc[df['FILTER'] == 'PASS']


def read_survivor_stats(path):
    """A custom function for importing SURVIVOR stats files. Returns a data frame"""
    file = open(path)
    df = pd.read_csv(file, sep="\t")
    df.columns = df.loc[0, :]
    new_df = df.loc[1:5, :]
    final_df = new_df.astype({'Len': str, 'Del': int, 'Dup': int, 'Inv': int, 'INS': int, 'TRA': int, 'UNK': int})
    return final_df


def get_file_name(path):
    """Takes a path and extracts the file name."""
    # Split the path and extract the final part of the path(which will be the name)
    # Adjusted to shorten the names of the files starting with RDC
    split = path.split('/')
    split2 = split[-1].split('.')
    split3 = split2[0].split('_aligned')
    file_name = split3[0]
    return file_name


def summarize_stats(df):
    """Summerize the counts create in the SURVIVOR stats. Returns a list of counts"""
    return [df['Del'].sum(), df['Dup'].sum(), df['Inv'].sum(), df['INS'].sum(), df['TRA'].sum(), df['UNK'].sum()]


def plot_stats(sum_stats, sample_name):
    """Plots the summerized stats. Add plt.show to be able to show the graph"""
    x_axis = ['DEL', 'DUP', 'INV', 'INS', 'TRA']
    plt.plot(x_axis, sum_stats, '.', label=sample_name)


def count_sv_type(sample):
    """Counts how many variants of each SV type(DEL, DUP, INV, INS and TRA). Returns a list of counts"""
    alt_column = sample['ALT']
    count_del = 0
    count_dup = 0
    count_ins = 0
    count_inv = 0
    count_tra = 0
    for j in alt_column:
        if j == '<DEL>':
            count_del += 1
        elif j == '<DUP>':
            count_dup += 1
        elif j == '<INS>':
            count_ins += 1
        elif j == '<INV>':
            count_inv += 1
        else:
            count_tra += 1

    return [count_del, count_dup, count_inv, count_ins, count_tra]


def count_svs(vcf_list, vcf_files):
    """Counts how many SV is in each location. Returns one dict object with positions as keys and counts as values and
    one dict with positions and what genomes have SVs there"""
    # The scripts loops each file and looks for all positions and summerizes the counts of each location
    counts = {}
    genomes_pos = {}
    for i in range(len(vcf_list)):
        genome = get_file_name(vcf_files[i])
        for row in range(len(vcf_list[i])):
            # Combine the chromosome and the start position to be able compare positions
            # The breakends will only handle the initial position, not the position on another chromosome
            chrom = vcf_list[i].iloc[row]['CHROM']
            pos = vcf_list[i].iloc[row]['POS']
            end = vcf_list[i].iloc[row]['END']
            sv_id = vcf_list[i].iloc[row]['ID']
            chrom_pos_end = chrom + ':' + pos + ':' + end

            if chrom_pos_end in counts.keys():
                counts[chrom_pos_end] += 1
                genomes_pos[chrom_pos_end] += ';' + genome + "(" + sv_id + ")"
            else:
                counts[chrom_pos_end] = 1
                genomes_pos[chrom_pos_end] = genome + "(" + sv_id + ")"
    return counts, genomes_pos
