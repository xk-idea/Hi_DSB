import os

# bin number
NUM = {
    5000 : [49850, 48639, 39604, 38230, 36183, 34223, 31827, 29272, 28242, 27106, 27001, 26770, 23033, 21469, 20506, 18070, 16239, 15615, 11825, 12605, 9625, 10260],
    10000 : [24926, 24320, 19803, 19116, 18092, 17112, 15914, 14637, 14122, 13554, 13501, 13386, 11517, 10735, 10254, 9036, 8120, 7808, 5913, 6303, 4813, 5131, 15528, 5938],
    50000 : [4986, 4864, 3961, 3824, 3619, 3423, 3183, 2928, 2825, 2711, 2701, 2678, 2304, 2147, 2051, 1808, 1624, 1562, 1183, 1261, 963, 1027]
}

# get each bin peak_density
def calculate_peak_density(chromosome, start, end, bin_size=10000):
    bin_start = start // bin_size
    bin_end = end // bin_size
    bins = range(bin_start, bin_end + 1)
    return [(bin, 1) for bin in bins]


def calculate_value_density(chromosome, start, end, value, bin_size=10000):
    bin_start = start // bin_size
    bin_end = end // bin_size
    bins = range(bin_start, bin_end + 1)
    overlap_length = end - start
    density = value * overlap_length / bin_size
    return [(bin, density) for bin in bins]

for filename in ['CTCF_Broad_NHEK_hg19', 'DNase_Uw_NHEK_hg19', 'H3K4me3_Broad_NHEK_hg19', 'H3K27ac_Broad_NHEK_hg19']:
    
    peak_densities = {}
    with open("../../data/NHEK/Epigenome/" + filename+ ".bed", "r") as file:
        for line in file:
            columns = line.strip().split("\t")
            chr, start, end = columns[:3]

            
            if chr in ["chrX", "chrY"]:
                continue

            start = int(start)
            end = int(end)
            chr_num = int(chr[3:])  
            
            
            bin_counts = calculate_peak_density(chr_num, start, end)
            
            
            if chr not in peak_densities:
                peak_densities[chr] = {}
            for bin, count in bin_counts:
                if bin not in peak_densities[chr]:
                    peak_densities[chr][bin] = 0
                peak_densities[chr][bin] += count

   
    for chr_num in range(1, 23):
        chr = f"chr{chr_num}"
        num_bins = NUM[chr_num - 1]
        densities = [peak_densities[chr].get(bin, 0) for bin in range(num_bins)]
        
        
        folder_path = f"Epigenome/{filename}"
        os.makedirs(folder_path, exist_ok=True)
        
        with open(f"Epigenome/" + filename + "/" + chr + "_dense.txt", "w") as file:
            file.write("\n".join(map(str, densities)))