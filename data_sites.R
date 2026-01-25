T1 <- list.files("/home/kunqi.chen/human/")[-57]
for(i in 1:length(T1)){
  if(i==1){
    T2 <- readRDS(paste0("/home/kunqi.chen/human/",T1[i]))
  }else{
    T2 <- c(T2,readRDS(paste0("/home/kunqi.chen/human/",T1[i])))
  }
}

T1 <- readRDS("/data/kunqidir/AC16_EV/exomePeak2_T/Mod.rds")
C1 <- readRDS("/data/kunqidir/AC16_EV/exomePeak2_con/Mod.rds")
T1gr <- unlist(attr(T1,"rowRanges"))
C1gr <- unlist(attr(C1,"rowRanges"))
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
hg19gtf <- makeTxDbFromGFF("/home/kunqi.chen/reference/gtf/hg19.refGene.gtf")

Over1 <- findOverlaps(T1gr,C1gr)
T1gr_p <- T1gr[-unique(queryHits(Over1))]
T1gr_p_s <- subsetByOverlaps(unique(T2),T1gr_p)


C1gr_p <- T1gr[-unique(subjectHits(Over1))]

T1gr_exons <- unique(subsetByOverlaps(transcripts(hg19gtf), T1gr_p))

###T1
library(m6ALogisticModel)
T1gr_exons_m6A <- m6ALogisticModel::sample_sequence("RRACH",T1gr_exons,BSgenome.Hsapiens.UCSC.hg19)
T1gr_p_m6A <- unique(subsetByOverlaps(T1gr_exons_m6A, T1gr_p))
T1gr_n_m6A <-T1gr_exons_m6A[-queryHits(findOverlaps(T1gr_exons_m6A,c(T1gr,C1gr)))]

saveRDS(T1gr_p_s,"/home/kunqi.chen/m6AHD/data/AC16EV_UP_p_m6A_s.rds")
saveRDS(T1gr_n_m6A,"/home/kunqi.chen/m6AHD/data/AC16EV_UP_n_m6A.rds")
###C1
C1gr_p_m6A <- m6ALogisticModel::sample_sequence("RRACH",C1gr_p,BSgenome.Hsapiens.UCSC.hg19)
C1gr_n_m6A <- subsetByOverlaps(sample_sequence("RRACH",T1gr,BSgenome.Hsapiens.UCSC.hg19),
                               sample_sequence("RRACH",C1gr,BSgenome.Hsapiens.UCSC.hg19))

saveRDS(C1gr_p_m6A,"/home/kunqi.chen/m6AHD/data/AC16EV_do_p_m6A.rds")
saveRDS(C1gr_n_m6A,"/home/kunqi.chen/m6AHD/data/AC16EV_do_n_m6A.rds")
