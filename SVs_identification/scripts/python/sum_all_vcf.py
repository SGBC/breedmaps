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

# # -------- Read in the files -------- # # #
## Counting how many svs are repeated in multiple genomes
# First we want to join all vcfs to get all positions
plt.figure(3)

count_raw = functions.count_svs(vcf_raw)
sorted_count_raw = sorted(count_raw.items())
chr_pos_raw, num_raw = zip(*sorted_count_raw)
plt.plot(num_raw)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV(chr and start pos) has been found in the raw VCF files')
plt.savefig("counted_nonfiltered_svs.png")

plt.figure(4)
count_filt = functions.count_svs(vcf_filtered)
sorted_count_filt = sorted(count_filt.items())
chr_pos_filt, num_filt = zip(*sorted_count_filt)
plt.plot(num_filt)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV(chr and start pos) has been found in the filtered VCF files')
plt.savefig("counted_filtered_svs.png")
plt.show()

x = {k:v for (k,v) in count_filt.items() if v > 10}
df = pd.Series(x)
print(f"There were {len(df)} SVs found in 10 genomes or more.")
print("CHR:POS              COUNT")
print(df.to_string())