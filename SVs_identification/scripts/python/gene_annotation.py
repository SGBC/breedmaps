import pandas as pd
import glob
import functions
import numpy as np

# gene_ann_path = '/Users/jj/breedmaps/SVs_identification/data/annotation/Bos_taurus.ARS-UCD1.2.102.gtf'
# genes = functions.read_gtf(gene_ann_path)
top_svs = pd.read_csv('top_svs.csv', sep=',')
top_svs.columns = ['CHR:POS', 'GENOMES']
top_svs[['CHR', 'POS']] = top_svs['CHR:POS'].str.split(':', expand=True)

# Read in vcf files to get all info
vcf_files = glob.glob("/Users/jj/breedmaps/SVs_identification/data/vcf/*.vcf")
vcf = {}
genome_names = []
for i in range(len(vcf_files)):
    file = functions.read_delly_vcf(vcf_files[i])
    # Filtered by the FILTER column(FITLER == PASS) and PRECISE flag
    filtered = functions.filter_qual(file)
    precise = functions.filter_precise(filtered)
    genome_name = functions.get_genome_name(vcf_files[i])
    genome_names.append(genome_name)
    vcf[genome_name] = precise

# print(functions.get_genome_name(vcf_files[0]))
# print(vcf[functions.get_genome_name(vcf_files[0])].columns)

# Seperate the genomes
genomes = top_svs['GENOMES'].str.split(';', expand=True)
genome_df = pd.DataFrame(np.nan, index=range(len(top_svs)), columns=genome_names)
genome_dict_list = []
for i in range(len(genome_names)):
    genome_dict = {}
    for row in range(len(top_svs)):
        if genome_names[i] in genomes[row].values:
            # genome_df.loc[genome_df.iloc[row][genome_names[i]]] = True
            genome_dict[genome_names[i]] = True
        else:
            # genome_df.loc[genome_df.iloc[row][genome_names[i]]] = False
            genome_dict[genome_names[i]] = False
    genome_dict_list.append(genome_dict)
# Choose one location and one genome to go forward with
df1 = vcf['BTA127_L4_SV']

df1_split = df1['INFO'].str.split(';', expand=True)
# Remove precise column since it is already in the previous data frame
df1_split = df1_split.iloc[:, 1:17]
df1_split.columns = ['SVTYPE', 'SVMETHOD', 'CHR2', 'END', 'PE', 'MAPQ', 'CT', 'CIPOS', 'CIEND', 'SRMAPQ', 'INSLEN', 'HOMLEN', 'SR', 'SRQ', 'CONSENSUS', 'CE']

for row in range(len(df1_split)):

    for col in range(len(df1_split.columns)):
        column = df1_split.iloc[row, col]
        df1_split_col = column.split('=')
        df1_split.iloc[row, col] = df1_split_col[1]

df1 = df1.join(df1_split)

print(df1.iloc[0,:])