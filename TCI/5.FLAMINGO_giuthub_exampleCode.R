library(FLAMINGOr)
library(GenomicFeatures)
library(Matrix)
all_size <- getChromInfoFromUCSC("hg19")
path <- 'W:/DSB_TAD_density/NHEK_run_FLAMINGO/1.FLAMINGO'
aim_flod <- 'K:/2024_NC/HiC_normalize/KR_MCFS/NHEK/7.FLAMINGO'
dir.create('K:/2024_NC/HiC_normalize/KR_MCFS/NHEK/7.FLAMINGO')
dirs <- list.files(path)
for (chr in dirs){
  save_flod <- paste0(aim_flod,'/',chr)
  dir.create(paste0(aim_flod,'/',chr))
  setwd(save_flod)
  print(chr)
  chr_name = chr
  chr_size = all_size[all_size$chrom==chr_name,2]
  res = flamingo.main_func_large(hic_data_low=paste0('K:/2024_NC/HiC_normalize/KR_MCFS/NHEK/4.sparse_2_hic/',chr,'/',chr,'.hic'),
                                 file_format='hic',
                                 domain_res=1e6,frag_res=10e3,
                                 chr_size=chr_size,
                                 chr_name=chr_name,
                                 normalization='KR',
                                 downsampling_rates=0.75,
                                 lambda=10,max_dist=0.01,nThread=30,n_row=500000)
  write.table(res, paste0(save_flod,"/",chr),quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  write.vtk(points=res[,-1],lookup_table=rep(1,dim(res)[1]),name=paste0(chr,' 10kb 3D structure'),opt_path=paste0(save_flod,"/",chr,'_10kb.vtk'))
}


