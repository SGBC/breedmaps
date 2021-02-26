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
genome_names = []
for i in range(len(vcf_raw)):
    sample = vcf_raw[i]
    counts = functions.count_sv_type(sample['ALT'])
    genome_name = functions.get_genome_name(vcf_files[i])
    genome_names.append(genome_name)
    functions.plot_stats(counts, genome_name)

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of all SVs found in each genome')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("nonfiltered_svs.png")

# Figure 2 plots the filtered identified SVs
plt.figure(2)
for i in range(len(vcf_filtered)):
    sample = vcf_filtered[i]
    counts = functions.count_sv_type(sample['ALT'])
    functions.plot_stats(counts, genome_names[i])

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of filtered SVs found in each genome')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("filtered_svs.png")

plt.figure(3)

count_raw, genome_pos_raw, all_info_raw = functions.count_svs(vcf_raw, vcf_files)
sorted_count_raw = sorted(count_raw.items())
chr_pos_raw, num_raw = zip(*sorted_count_raw)
plt.plot(num_raw)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV(chr and start pos) has been found in the raw VCF files')
plt.savefig("counted_nonfiltered_svs.png")

plt.figure(4)
count_filt, genome_pos_filt, all_info_filt = functions.count_svs(vcf_filtered, vcf_files)
sorted_count_filt = sorted(count_filt.items())
chr_pos_filt, num_filt = zip(*sorted_count_filt)
plt.plot(num_filt)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV(chr and start pos) has been found in the filtered VCF files')
plt.savefig("counted_filtered_svs.png")
plt.show()

top_svs = {k:v for (k,v) in count_filt.items() if v > 10}
df_count = pd.Series(top_svs)
print(f"There were {len(df_count)} SVs found in 10 genomes or more.")
top_svs_genomes = {k:v for (k,v) in genome_pos_filt.items() if k in top_svs.keys()}
df_genomes = pd.Series(top_svs_genomes)
print("CHR:POS              COUNT   GENOMES")
merged_top_svs = pd.concat([df_count, df_genomes], axis=1)
print(merged_top_svs.to_string())

top_svs_all_info = {k:v for (k,v) in all_info_filt.items() if k in top_svs.keys()}
print(top_svs_all_info)
df_all_info = pd.Series(top_svs_all_info)
merged_top_svs_all_info = pd.concat([df_genomes, df_all_info], axis=1)
merged_top_svs_all_info.to_csv('top_svs.csv', sep=',')
