from functions import read_delly_vcf, filter_precise, filter_qual
import glob
vcf_files = glob.glob("/Users/jj/breedmaps/SVs_identification/data/vcf/*.vcf")
vcf = [None] * len(vcf_files)
for i in range(len(vcf_files)):
    file = read_delly_vcf(vcf_files[i])
    filtered = filter_qual(file)
    precise = filter_precise(filtered)
    vcf[i] = precise
    print(f"The file has gone from {len(file)} to {len(vcf[i])} by filtering by FILTER and PRECISE flag\n")



