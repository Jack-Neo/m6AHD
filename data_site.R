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
T1gr_exons <- unique(subsetByOverlaps(transcripts(hg19gtf), T1gr_p))

###T1
library(m6ALogisticModel)
T1gr_exons_m6A <- m6ALogisticModel::sample_sequence("RRACH",T1gr_exons,BSgenome.Hsapiens.UCSC.hg19)
T1gr_n_m6A <-T1gr_exons_m6A[-queryHits(findOverlaps(T1gr_exons_m6A,c(T1gr,C1gr)))]

#saveRDS(T1gr_p_s,"/home/lujiajie/ml/new/hypertrophy/up/hypertrophy_UP_p_m6A_s.rds")
#saveRDS(T1gr_n_m6A,"/home/lujiajie/ml/new/hypertrophy/up/hypertrophy_UP_n_m6A.rds")


##C1
# 筛选对照组特有peak（不与实验组peak重叠）
Overlap_control <- findOverlaps(C1gr, T1gr)
C1gr_p <- C1gr[-unique(queryHits(Overlap_control))]  # 对照组特有peak

# 获取正样本：对照组特有peak中的单碱基位点（且不在实验组数据中）
C1gr_p_s <- subsetByOverlaps(unique(T2), C1gr_p)  # 对照组特有修饰位点


# 获取负样本：两组共有peak中的修饰位点
Common_peaks <- findOverlaps(T1gr, C1gr)  # 重叠区域
Common_gr <- pintersect(T1gr[queryHits(Common_peaks)], 
                        C1gr[subjectHits(Common_peaks)])  # 精确交集

common_motif=m6ALogisticModel::sample_sequence("RRACH",Common_gr,BSgenome.Hsapiens.UCSC.hg19)

# 取实验组和对照组在共有区域中的交集位点
Negative_gr <- subsetByOverlaps(unique(T2), common_motif)  # 实验组数据


# 保存结果
saveRDS(C1gr_p_s, "/home/lujiajie/ml/new/AC16EV/down/AC16EV_do_p_m6A_s.rds")
saveRDS(Negative_gr, "/home/lujiajie/ml/new/AC16EV/down/AC16EV_do_n_m6A.rds")