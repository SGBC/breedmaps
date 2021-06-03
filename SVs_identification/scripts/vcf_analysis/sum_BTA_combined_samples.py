import functions
import glob
import matplotlib.pyplot as plt
import pandas as pd
from os.path import expanduser

# Read in all vcf files in your directory
home_path = expanduser("~")
vcf_files_raw = glob.glob(home_path + "/breedmaps/SVs_identification/data/vcf/combined*")
vcf_files_filt = glob.glob(home_path + "/breedmaps/SVs_identification/results/filtered_variants/precise_combined*")
vcf_filtered = [None] * len(vcf_files_filt)
vcf_raw = [None] * len(vcf_files_raw)

# # # -------- Load the files -------- # # #
# Read in all files and save both raw and filtered.
for i in range(len(vcf_files_raw)):
    file = functions.read_delly_vcf(vcf_files_raw[i])
    vcf_raw[i] = file

for i in range(len(vcf_files_filt)):
    precise = pd.read_csv(
        vcf_files_filt[i],
        sep='\t',
        dtype=str
    )
    vcf_filtered[i] = precise

# # # -------- Plot -------- # # #
# Plot the number of all variants in each file
genome_names = []
num_var_raw = []
plt.figure(1, figsize=(19.20,10.80))
print("Number of variants in each genome [#SVs]")

for i in range(len(vcf_raw)):
    sample = vcf_raw[i]
    genome_name = functions.get_file_name(vcf_files_raw[i])
    genome_names.append(genome_name)
    num_var_raw.append(len(sample))
    print(genome_name, len(sample))

plt.plot(genome_names, num_var_raw, '.')
plt.title('Number of SVs found in each genome [#SVs]')
plt.xticks(fontsize=4, rotation=45, ha='right')
plt.ylabel('# SVs')
plt.savefig('combined_BTA_count_per_genome.png')

# Plot the number of filtered variants in each file
num_var_filt = []
plt.figure(2, figsize=(19.20,10.80))
print("Number of filtered variants in each genome")
for i in range(len(vcf_filtered)):
    sample = vcf_filtered[i]
    num_var_filt.append(len(sample))
    print(genome_names[i], len(sample))

plt.plot(genome_names, num_var_filt, '.')
plt.title('Number of filtered SVs found in each genome')
plt.xticks(fontsize=4, rotation=45, ha='right')
plt.ylabel('# SVs')
plt.savefig('combined_BTA_filt_count_per_genome.png')

# Plot the number of each type of structural variation
plt.figure(3, figsize=(19.20,10.80))
# Figure 3 is the raw identified SVs.
counts = [None] * len(vcf_raw)
print("The counts of SV types in raw files['DEL', 'DUP', 'INV', 'INS', 'TRA']")
for i in range(len(vcf_raw)):
    sample = vcf_raw[i]
    counts[i] = functions.count_sv_type(sample)
    functions.plot_stats(counts[i], genome_names[i])
    print(genome_names[i], counts[i])

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of each SV type in all files')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("combined_BTA_nonfiltered_svs.png")



# Figure 4 plots the filtered identified SVs
plt.figure(4, figsize=(19.20,10.80))
counts = [None] * len(vcf_filtered)
print("The counts of SV types in filtered files ['DEL', 'DUP', 'INV', 'INS', 'TRA']")
for i in range(len(vcf_filtered)):
    sample = vcf_filtered[i]
    counts[i] = functions.count_sv_type(sample)
    functions.plot_stats(counts[i], genome_names[i])
    print(genome_names[i], counts[i])

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of each SV type in the filtered files')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("combined_BTA_filtered_svs.png")

plt.figure(5, figsize=(19.20,10.80))
count_raw, genome_pos_raw = functions.count_svs(vcf_raw, vcf_files_raw)
sorted_count_raw = sorted(count_raw.items())
chr_pos_raw, num_raw = zip(*sorted_count_raw)
plt.plot(num_raw)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV has been found in the raw VCF files')
plt.savefig("combined_BTA_counted_nonfiltered_svs.png")

plt.figure(6, figsize=(19.20,10.80))
count_filt, genome_pos_filt = functions.count_svs(vcf_filtered, vcf_files_filt)
sorted_count_filt = sorted(count_filt.items())
chr_pos_filt, num_filt = zip(*sorted_count_filt)
plt.plot(num_filt)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV has been found in the filtered VCF files')
plt.savefig("combined_BTA_counted_filtered_svs.png")

top_svs = {k: v for (k, v) in count_filt.items() if v > 6}
df_count = pd.Series(top_svs)
top_svs_genomes = {k: v for (k, v) in genome_pos_filt.items() if k in top_svs.keys()}
df_genomes = pd.Series(top_svs_genomes)
merged_top_svs = pd.concat([df_count, df_genomes], axis=1)
merged_top_svs.to_csv('combined_BTA_top_svs.csv', sep=',')

all_svs = {k: v for (k, v) in count_filt.items()}
df_all = pd.Series(all_svs)
all_svs_genomes = {k: v for (k, v) in genome_pos_filt.items() if k in all_svs.keys()}
df_all_genomes = pd.Series(all_svs_genomes)
merged_all_svs = pd.concat([df_all, df_all_genomes], axis=1)
merged_all_svs.to_csv('combined_BTA_all_svs.csv', sep=',')

# Uncomment this line to output the plots instead of just saving them
#plt.show()