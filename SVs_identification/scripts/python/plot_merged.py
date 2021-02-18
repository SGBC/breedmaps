import pandas as pd
import matplotlib.pyplot as plt
from functions import read_survivor_stats, summarize_stats
stats = read_survivor_stats('/Users/jj/breedmaps/SVs_identification/SURVIVOR/stats.txt')
x_axis = ['DEL', 'DUP', 'INV', 'INS', 'TRA', 'UNK']
y_axis = summarize_stats(stats)
plt.plot(x_axis,y_axis, '*', label = 'All files merged')
plt.legend()
plt.title('Summary of number of structural variants found')
plt.show()
#summary = pd.DataFrame({'DEL': final_df['Del'].sum(), 'DUP': final_df['Dup'].sum(), 'INV': final_df['Inv'].sum(), 'INS': final_df['INS'].sum(), 'TRA': final_df['TRA'].sum(), 'UNK': final_df['UNK'].sum()})

