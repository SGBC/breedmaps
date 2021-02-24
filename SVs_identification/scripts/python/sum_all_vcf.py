import functions
import glob
import matplotlib.pyplot as plt
import pandas as pd

# Read in all vcf files in your directory
vcf_files = glob.glob("/Users/jj/breedmaps/SVs_identification/data/vcf/*.vcf")
vcf_filtered = [None] * len(vcf_files)
vcf_raw = [None] * len(vcf_files)

# # # -------- Load the files -------- # # #
# Read in all files and save both raw and filtered.
for i in range(len(vcf_files)):
    file = functions.read_delly_vcf(vcf_files[i])
    vcf_raw[i] = file
    # Filtered by the FILTER column(FITLER == PASS) and PRECISE flag
    filtered = functions.filter_qual(file)
    precise = functions.filter_precise(filtered)
    vcf_filtered[i] = precise


# # # -------- Plot -------- # # #
# Plot the number of each type of structural variation
plt.figure(1)
# Figure 1 is the raw identified SVs.
for i in range(len(vcf_raw)):
    sample = vcf_raw[i]
    counts = functions.count_sv_type(sample['ALT'])
    functions.plot_stats(counts, functions.get_genome_name(vcf_files[i]))

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of all SVs found in each genome')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("nonfiltered_svs.png")
#plt.show()

# Figure 2 plots the filtered identified SVs
plt.figure(2)
for i in range(len(vcf_filtered)):
    sample = vcf_filtered[i]
    counts = functions.count_sv_type(sample['ALT'])
    functions.plot_stats(counts, functions.get_genome_name(vcf_files[i]))

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of filtered SVs found in each genome')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("filtered_svs.png")
#plt.show()

# # -------- Read in the files -------- # # #
## Counting how many svs are repeated in multiple genomes
# First we want to join all vcfs to get all positions
plt.figure(3)
for i in range(len(vcf_filtered)):
    count = functions.count_svs(vcf_filtered[i])
    sorted_count = sorted(count.items())
    chr_pos, num = zip(*sorted_count)
    plt.plot(num, label=functions.get_genome_name(vcf_files[i]))
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.show()
