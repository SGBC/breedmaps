from functions import read_delly_vcf, filter_precise, filter_qual, plot_stats, get_sample_name, count_sv_type
import glob
import matplotlib.pyplot as plt
vcf_files = glob.glob("/Users/jj/breedmaps/SVs_identification/data/vcf/*.vcf")
vcf_filtered = [None] * len(vcf_files)
vcf_raw = [None] * len(vcf_files)
for i in range(len(vcf_files)):
    file = read_delly_vcf(vcf_files[i])
    vcf_raw[i] = file

    filtered = filter_qual(file)
    precise = filter_precise(filtered)
    vcf_filtered[i] = precise


plt.figure(1)
for i in range(len(vcf_raw)):
    sample = vcf_raw[i]
    counts = count_sv_type(sample['ALT'])
    plot_stats(counts, get_sample_name(vcf_files[i]))

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of all SVs found in each genome')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("nonfiltered_svs.png")
plt.show()

plt.figure(2)
for i in range(len(vcf_filtered)):
    sample = vcf_filtered[i]
    counts = count_sv_type(sample['ALT'])
    plot_stats(counts, get_sample_name(vcf_files[i]))

plt.legend(loc='upper center', fontsize='xx-small', ncol=3)
plt.title('Number of filtered SVs found in each genome')
plt.xlabel('SV type')
plt.ylabel('# SVs')
plt.savefig("filtered_svs.png")
plt.show()
