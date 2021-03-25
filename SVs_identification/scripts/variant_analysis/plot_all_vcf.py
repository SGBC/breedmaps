import glob
import matplotlib.pyplot as plt
from functions import read_survivor_stats, summarize_stats, get_genome_name, plot_stats

stats_files = glob.glob("/Users/jj/breedmaps/SVs_identification/results/SURVIVOR/*.stats")

for i in range(len(stats_files)):
    file_path = stats_files[i]
    stats = read_survivor_stats(file_path)
    print(stats)
    sample_name = get_genome_name(file_path)
    sum_stats = summarize_stats(stats)
    plot_stats(sum_stats, sample_name)

plt.legend(loc='upper center', ncol=5)
plt.title('Number of SVs found in each genome')
plt.show()
