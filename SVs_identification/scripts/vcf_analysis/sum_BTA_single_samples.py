import functions
import glob
import matplotlib.pyplot as plt
import pandas as pd
from os.path import expanduser

# Read in all vcf files in your directory
home_path = expanduser("~")
vcf_files_raw = glob.glob(home_path + "/breedmaps/SVs_identification/data/vcf/combined_BTA*")
vcf_files_filt = glob.glob(home_path + "/breedmaps/SVs_identification/results/filtered_variants/precise_combined_BTA*")
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
file_names = []
num_var_raw = []
plt.figure(1, figsize=(19.20,10.80))
print("Number of variants in each file [#SVs]")

for i in range(len(vcf_raw)):
    sample = vcf_raw[i]
    file_name = functions.get_file_name(vcf_files_raw[i])
    file_names.append(file_name)
    num_var_raw.append(len(sample))
    print(file_name, len(sample))

plt.plot(file_names, num_var_raw, '.')
plt.title('Number of SVs found in each file [#SVs]')
plt.xticks(fontsize=4, rotation=45, ha='right')
plt.ylabel('# SVs')
plt.savefig('single_BTA_count_per_file.png')

# Plot the number of filtered variants in each file
num_var_filt = []
plt.figure(2, figsize=(19.20,10.80))
for i in range(len(vcf_filtered)):
    sample = vcf_filtered[i]
    num_var_filt.append(len(sample))

plt.plot(file_names, num_var_filt, '.')
plt.title('Number of filtered SVs found in each file')
plt.xticks(fontsize=4, rotation=45, ha='right')
plt.ylabel('# SVs')
plt.savefig('single_BTA_filt_count_per_file.png')

# Plot the number of each type of structural variation
plt.figure(3, figsize=(19.20,10.80))
# Figure 3 is the raw identified SVs.
counts = [None] * len(vcf_raw)
for i in range(len(vcf_raw)):
    sample = vcf_raw[i]
    counts[i] = functions.count_sv_type(sample)
    functions.plot_stats(counts[i], file_names[i])

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of each SV type in all files')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("single_BTA_nonfiltered_svs.png")


# Figure 4 plots the filtered identified SVs
plt.figure(4, figsize=(19.20,10.80))
counts = [None] * len(vcf_filtered)
for i in range(len(vcf_filtered)):
    sample = vcf_filtered[i]
    counts[i] = functions.count_sv_type(sample)
    functions.plot_stats(counts[i], file_names[i])

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of each SV type in the filtered files')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("single_BTA_filtered_svs.png")

##################
# Counting how many times the SVs occur in the dataset(all vcf files loadedd)
##################
plt.figure(5, figsize=(19.20,10.80))
# Count_svs counts how many times each location(chr:start:end) occurs in all files
# It returns both counts and the position as two dict objects where the keys are the locations
count_raw, file_pos_raw = functions.count_svs(vcf_raw, vcf_files_raw)
sorted_count_raw = sorted(count_raw.items())
chr_pos_raw, num_raw = zip(*sorted_count_raw)
plt.plot(num_raw)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV has been found in the raw VCF files')
plt.savefig("single_BTA_counted_nonfiltered_svs.png")

plt.figure(6, figsize=(19.20,10.80))
# Count_svs counts how many times each location(chr:start:end) occurs in all files
# It returns both counts and the position as two dict objects where the keys are the locations
count_filt, file_pos_filt = functions.count_svs(vcf_filtered, vcf_files_filt)
sorted_count_filt = sorted(count_filt.items())
chr_pos_filt, num_filt = zip(*sorted_count_filt)
plt.plot(num_filt)
plt.xlabel('SVs')
plt.ylabel('Number of times each SV is found')
plt.title('The number of times each SV has been found in the filtered VCF files')
plt.savefig("single_BTA_counted_filtered_svs.png")

# Extract the counts for each position to be outputted into a csv file
# Only the filtered SVs are outputted into a csv file
all_svs = {pos: count for (pos, count) in count_filt.items()}
df_all = pd.Series(all_svs)
all_svs_files = {pos: file for (pos, file) in file_pos_filt.items() if pos in all_svs.keys()}
df_all_files = pd.Series(all_svs_files)
merged_all_svs = pd.concat([df_all, df_all_files], axis=1)
merged_all_svs.to_csv('single_BTA_counted_svs.csv', sep=',')

# Do the same as above but filter out the low count SVs
# Extract the SVs which occurs in 70 % of the files/individuals ( in this case 11 files/individuals)
top_svs = {pos: count for (pos, count) in count_filt.items() if count > 11}
df_count = pd.Series(top_svs)
# Fetch the files for each high count SV
top_svs_files = {pos: file for (pos, file) in file_pos_filt.items() if pos in top_svs.keys()}
df_files = pd.Series(top_svs_files)
# Paste together the counts and the files/individuals
merged_top_svs = pd.concat([df_count, df_files], axis=1)
merged_top_svs.to_csv('single_BTA_top_svs.csv', sep=',')



