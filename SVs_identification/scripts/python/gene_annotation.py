import functions
import pandas as pd

#gene_ann_path = '/Users/jj/breedmaps/SVs_identification/data/annotation/Bos_taurus.ARS-UCD1.2.102.gtf'
#genes = functions.read_gtf(gene_ann_path)
top_svs = pd.read_csv('top_svs.csv', sep=',')
top_svs.columns = ['CHR:POS', 'GENOMES', 'INFO']
split = top_svs['INFO'].str.split(";", expand=True)
print(len(split))
print(split.iloc[0,:])


