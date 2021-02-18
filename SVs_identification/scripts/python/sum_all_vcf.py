from functions import read_delly_vcf, filter_precise
import glob
vcf_files = glob.glob("/Users/jj/breedmaps/SVs_identification/data/vcf/*.vcf")
vcf = [None] * len(vcf_files)
for i in range(len(vcf_files)):
    new_vcf = read_delly_vcf(vcf_files[i])
    vcf[i] = filter_precise(new_vcf)

