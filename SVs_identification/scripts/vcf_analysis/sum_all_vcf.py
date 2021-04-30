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
# Plot the number of all variants in each file
genome_names = []
num_var_raw = []
plt.figure(1)
for i in range(len(vcf_raw)):
    sample = vcf_raw[i]
    genome_name = functions.get_file_name(vcf_files[i])
    genome_names.append(genome_name)
    num_var_raw.append(len(sample))

plt.plot(genome_names, num_var_raw, '.')
plt.title('Number of SVs found in each genome')
plt.xticks(fontsize=4, rotation=45, ha='right')
plt.ylabel('# SVs')
plt.savefig('count_per_genome.png')
plt.show()
# Plot the number of filtered variants in each file
num_var_filt = []
plt.figure(2)
for i in range(len(vcf_filtered)):
    sample = vcf_filtered[i]
    num_var_filt.append(len(sample))

plt.plot(genome_names, num_var_filt, '.')
plt.title('Number of filtered SVs found in each genome')
plt.xticks(fontsize=4, rotation=45, ha='right')
plt.ylabel('# SVs')
plt.savefig('filt_count_per_genome.png')
plt.show()

# Plot the number of each type of structural variation
plt.figure(3)
# Figure 3 is the raw identified SVs.
counts = [None] * len(vcf_raw)
for i in range(len(vcf_raw)):
    sample = vcf_raw[i]
    counts[i] = functions.count_sv_type(sample)
    functions.plot_stats(counts[i], genome_names[i])

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of all times each SVs occurs')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("nonfiltered_svs.png")



# Figure 4 plots the filtered identified SVs
plt.figure(4)
counts = [None] * len(vcf_filtered)
for i in range(len(vcf_filtered)):
    sample = vcf_filtered[i]
    counts[i] = functions.count_sv_type(sample)
    functions.plot_stats(counts[i], genome_names[i])

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of all times each filtered SVs occurs')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("filtered_svs.png")

#summary = pd.DataFrame(
#    {"File": genome_names, "Number of non-filtered variants": num_var_raw, "Number of filtered variants": num_var_filt,
#     "filt_DEL": counts[:][1], "filt_DUP": counts[:][2], "filt_INV": counts[:][3], "filt_INS": counts[:][4], "filt_TRA": counts[:][5]})


plt.figure(5)

count_raw, genome_pos_raw = functions.count_svs(vcf_raw, vcf_files)
sorted_count_raw = sorted(count_raw.items())
chr_pos_raw, num_raw = zip(*sorted_count_raw)
plt.plot(num_raw)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV has been found in the raw VCF files')
plt.savefig("counted_nonfiltered_svs.png")

plt.figure(6)
count_filt, genome_pos_filt = functions.count_svs(vcf_filtered, vcf_files)
sorted_count_filt = sorted(count_filt.items())
chr_pos_filt, num_filt = zip(*sorted_count_filt)
plt.plot(num_filt)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV has been found in the filtered VCF files')
plt.savefig("counted_filtered_svs.png")

top_svs = {k: v for (k, v) in count_filt.items() if v > 10}
df_count = pd.Series(top_svs)
top_svs_genomes = {k: v for (k, v) in genome_pos_filt.items() if k in top_svs.keys()}
df_genomes = pd.Series(top_svs_genomes)
merged_top_svs = pd.concat([df_count, df_genomes], axis=1)
merged_top_svs.to_csv('top_svs.csv', sep=',')

all_svs = {k: v for (k, v) in count_filt.items()}
df_all = pd.Series(all_svs)
all_svs_genomes = {k: v for (k, v) in genome_pos_filt.items() if k in all_svs.keys()}
df_all_genomes = pd.Series(all_svs_genomes)
merged_all_svs = pd.concat([df_all, df_all_genomes], axis=1)
merged_all_svs.to_csv('all_svs.csv', sep=',')

plt.show()
