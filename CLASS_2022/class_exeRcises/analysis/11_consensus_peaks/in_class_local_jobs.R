library(tidyverse)
library(GenomicRanges)

source("/scratch/Shares/rinnclass/CLASS_2022/alhi9177/bchm5631_alex/CLASS_2022/util/class_2022_functions.R")

consensus_peaks <- create_consensus_peaks(broadpeakfilepath = "/scratch/Shares/rinnclass/CLASS_2022/data/peaks/") 

for(i in 1:length(consensus_peaks)) {
  rtracklayer::export(consensus_peaks[[i]], 
                      paste0("/scratch/Shares/rinnclass/CLASS_2022/alhi9177/bchm5631_alex/CLASS_2022/class_exeRcises/analysis/11_consensus_peaks/consensus_peaks/", 
                             names(consensus_peaks)[i], 
                             "_consensus_peaks.bed"))
}

