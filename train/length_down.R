set.seed(123)
#C=n,T=p
# 加载需要的库
library(pROC)
library(openxlsx)  # 用于写入xlsx文件
# 初始化一个空的数据框来存储每次训练的AUC
auc_results <- data.frame(length = character(0), AUC = numeric(0))
for (max_length in c(21,31,41,51,61,71,81)){
C <- readRDS("D:\\Users\\Administrator\\Desktop\\machine learning\\new\\TKIs\\down\\TKIs_do_n_m6A.rds")
T<- readRDS("D:\\Users\\Administrator\\Desktop\\machine learning\\new\\TKIs\\down\\TKIs_do_p_m6A_s.rds")

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

# 提取起始和结束位置
start_positions <- start(C)
end_positions <- end(C)

# 创建新的GRanges对象，调整范围
adjusted_granges <- GRanges(
  seqnames = seqnames(C),
  ranges = IRanges(
    start = start_positions - (max_length-1)/2,  # 将起始位置调整为原位置前个碱基
    end = end_positions + (max_length-1)/2       # 将结束位置调整为原位置后个碱基
  ),
  strand = strand(C)
)


# 获取调整后的序列
adjusted_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, adjusted_granges)

# 转换为字符类型并创建数据框
sequence_strings <- as.character(adjusted_sequences)
output_df_n <- data.frame(
  seqnames = seqnames(C),
  start = start_positions,
  end = end_positions,
  strand = strand(C),
  adjusted_start = start(adjusted_granges),
  adjusted_end = end(adjusted_granges),
  sequence = sequence_strings
)



start_positions <- start(T)
end_positions <- end(T)

# 创建新的GRanges对象，调整范围
adjusted_granges <- GRanges(
  seqnames = seqnames(T),
  ranges = IRanges(
    start = start_positions - (max_length-1)/2,  # 将起始位置调整为原位置前2个碱基
    end = end_positions + (max_length-1)/2       # 将结束位置调整为原位置后23个碱基
  ),
  strand = strand(T)
)

# 获取调整后的序列
adjusted_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, adjusted_granges)

# 转换为字符类型并创建数据框
sequence_strings <- as.character(adjusted_sequences)
output_df_p <- data.frame(
  seqnames = seqnames(T),
  start = start_positions,
  end = end_positions,
  strand = strand(T),
  adjusted_start = start(adjusted_granges),
  adjusted_end = end(adjusted_granges),
  sequence = sequence_strings
)


library(dplyr)
# 假设 df1 和 df2 已经存在
output_df_p$type <- "T"  # 在 df1 中添加 type 列，并将其值设置为 1
output_df_n$type <- "C"  # 在 df2 中添加 type 列，并将其值设置为 0

output_df_n <- output_df_n[!grepl("N", output_df_n$sequence), ]
output_df_p <- output_df_p[!grepl("N", output_df_p$sequence), ]

set.seed(123)
# 随机抽取 5000 行
output_df_p_sample <- output_df_p %>%
  sample_n(1000)
output_df_n_sample <- output_df_n %>%
  sample_n(1000)
result <- rbind(output_df_n_sample, output_df_p_sample)

# 载入需要的库
library(e1071)  # 支持向量机库
library(dplyr)  # 数据操作库
library(tidyr)  # 数据整理库

data=result

pad_sequence <- function(seq) {
  chars <- strsplit(seq, "")[[1]]
  if (length(chars) > max_length) {
    chars <- chars[1:max_length] # 截断长序列
  } else {
    chars <- c(chars, rep("N", max_length - length(chars))) # 填充"N"
  }
  return(paste(chars, collapse = ""))
}

data$sequence <- sapply(data$sequence, pad_sequence)

#OH编码
one_hot_encode <- function(base) {
  bases <- c("A", "T", "C", "G")
  if (base == "N") {
    return(rep(0, 4))
  } else {
    return(as.numeric(bases == base))
  }
}

# 将序列转换为独热编码矩阵,G0001,T0100,A1000,C0010
encoded_seqs_oh <- lapply(data$sequence, function(seq) {
  chars <- strsplit(seq, "")[[1]]
  matrix(sapply(chars, one_hot_encode), ncol = 4, byrow = TRUE)
})

# 展平为二维特征矩阵（样本数 × (max_length × 4)）
feature_matrix_oh <- t(sapply(encoded_seqs_oh, function(mat) as.vector(t(mat))))
#创建列名，50可随序列长度改变
# 创建一个空的字符向量，用于存储列名
col_names_oh <- c()

# 循环创建列名
for (i in 1:max_length) {
  col_names_oh <- c(col_names_oh, paste("A_oh", i, sep = "_"), 
                    paste("T_oh", i, sep = "_"), 
                    paste("C_oh", i, sep = "_"), 
                    paste("G_oh", i, sep = "_"))
}
# 假设 feature_matrix_oh 是一个 50 行的数据框
colnames(feature_matrix_oh) <- col_names_oh
data=cbind(data,feature_matrix_oh)


# EIIP编码,EIIP值定义
eiip_values <- list(
  A = c(0.1260),
  T = c(0.1335),
  C = c(0.0806),
  G = c(0.1340)
)

# 将序列转换为EIIP值向量
get_eiip_vector <- function(seq) {
  # 将每个碱基转换为对应的EIIP值
  chars <- strsplit(seq, "")[[1]]
  eiip_vector <- sapply(chars, function(base) {
    if (base %in% names(eiip_values)) {
      return(eiip_values[[base]])
    } else {
      return(0) # 如果是'N'，用0填充
    }
  })
  return(unlist(eiip_vector))
}
# 将序列转换为EIIP值向量
eiip_seqs <- lapply(data$sequence, get_eiip_vector)

# 展平为二维特征矩阵（样本数 × max_length）
feature_matrix_eiip <- t(sapply(eiip_seqs, function(seq) {
  # 通过EIIP值向量展开成一维向量
  return(seq)
}))
#创建列名，50可随序列长度改变
# 创建一个空的字符向量，用于存储列名
col_names_eiip <- c()
# 循环创建列名
for (i in 1:max_length) {
  col_names_eiip <- c(col_names_eiip, paste("eiip", i, sep = "_"))
}
# 假设 feature_matrix_oh 是一个 50 行的数据框
colnames(feature_matrix_eiip) <- col_names_eiip
data=cbind(data,feature_matrix_eiip)

# 使用NAC编码
# 定义函数：计算DNA序列的四维向量（A, T, C, G的相对频率）
extract_four_dimensional_features <- function(seq) {
  total_length <- nchar(seq)  # 序列的总长度
  
  # 计算每个碱基的出现次数
  A_count <- sum(strsplit(seq, NULL)[[1]] == "A")
  T_count <- sum(strsplit(seq, NULL)[[1]] == "T")
  C_count <- sum(strsplit(seq, NULL)[[1]] == "C")
  G_count <- sum(strsplit(seq, NULL)[[1]] == "G")
  
  # 计算每种碱基的相对频率
  A_freq <- A_count / total_length
  T_freq <- T_count / total_length
  C_freq <- C_count / total_length
  G_freq <- G_count / total_length
  
  # 返回一个四维向量
  return(c(A_freq, T_freq, C_freq, G_freq))
}

# 将序列转换为NAC编码（相对频率特征）
encoded_seqs_nac <- lapply(data$sequence, extract_four_dimensional_features)

# 将NAC编码结果合并为特征矩阵
feature_matrix_nac <- do.call(rbind, encoded_seqs_nac)
#创建列名
col_names_nac <- c("A_freq","T_freq","C_freq","G_freq")


colnames(feature_matrix_nac) <- col_names_nac
data=cbind(data,feature_matrix_nac)


#DNC编码：提取二核苷酸频率的函数
extract_dinucleotide_frequencies <- function(sequence) {
  # 定义所有的二联体
  dinucleotides <- c("AA", "AC", "AT", "AG", "TA", "TT", "TC", "TG", 
                     "CA", "CT", "CC", "CG", "GA", "GT", "GC","GG")
  
  # 初始化一个向量，用来存储二联体的频率
  freq_vector <- rep(0, length(dinucleotides))
  
  # 序列长度
  seq_length <- nchar(sequence)
  
  # 如果序列长度小于 2，返回零频率
  if (seq_length < 2) {
    return(freq_vector)
  }
  
  # 获取序列的所有二联体
  for (i in 1:(seq_length - 1)) {
    dinucleotide <- substr(sequence, i, i + 1)
    if (dinucleotide %in% dinucleotides) {
      freq_vector[which(dinucleotides == dinucleotide)] <- freq_vector[which(dinucleotides == dinucleotide)] + 1
    }
  }
  
  # 计算相对频率
  relative_freq <- freq_vector / (seq_length - 1)
  
  return(relative_freq)
}

# 对所有DNA序列应用DAC编码
encoded_seqs_dnc <- lapply(data$sequence, extract_dinucleotide_frequencies)

# 将所有序列的DAC特征合并成一个矩阵
feature_matrix_dnc <- do.call(rbind, encoded_seqs_dnc)
#创建列名
col_names_dnc <- c("freq_AA", "freq_AC", "freq_AT", "freq_AG", "freq_TA", "freq_TT", "freq_TC", "freq_TG", 
                   "freq_CA", "freq_CT", "freq_CC", "freq_CG", "freq_GA", "freq_GT", "freq_GC","freq_GG")


colnames(feature_matrix_dnc) <- col_names_dnc
data=cbind(data,feature_matrix_dnc)


# 基于化学物理性质的编码,CP编码
CP_encoding <- list(
  A = c(1, 1, 1),
  T = c(0, 0, 1),
  C = c(0, 1, 0),
  G = c(1, 0, 0)
)

# 修改编码函数，使用CP_encoding
cp_encode <- function(base) {
  if (base %in% names(CP_encoding)) {
    return(CP_encoding[[base]])
  } else {
    return(c(0, 0, 0))  # 对于'N'或者无法识别的碱基返回[0, 0, 0]
  }
}

# 将序列转换为基于CP_encoding的矩阵
encoded_seqs_cp <- lapply(data$sequence, function(seq) {
  chars <- strsplit(seq, "")[[1]]
  matrix(sapply(chars, cp_encode), ncol = 3, byrow = TRUE)  # 每个碱基有3个特征
})

# 展平为二维特征矩阵（样本数 × (max_length × 3)）
feature_matrix_cp <- t(sapply(encoded_seqs_cp, function(mat) as.vector(t(mat))))

#创建列名
col_names_dnc <- c("ring structure", "functional groups", "hydrogen bonds")
#创建列名，50可随序列长度改变
# 创建一个空的字符向量，用于存储列名
col_names_cp <- c()

# 循环创建列名
for (i in 1:max_length) {
  col_names_cp <- c(col_names_cp, paste("ring structure", i, sep = "_"), 
                    paste("functional groups", i, sep = "_"), 
                    paste("hydrogen bonds", i, sep = "_"))
}
# 假设 feature_matrix_oh 是一个 50 行的数据框
colnames(feature_matrix_cp) <- col_names_cp
data=cbind(data,feature_matrix_cp)


# 计算累计频率,anf编码
calculate_cumulative_frequencies <- function(sequence) {
  n <- nchar(sequence)
  cumulative_frequencies <- matrix(0, nrow = n, ncol = 4)
  colnames(cumulative_frequencies) <- c("A", "G", "C", "T")
  
  seq_vector <- strsplit(sequence, NULL)[[1]]
  for (i in 1:n) {
    sub_seq <- seq_vector[1:i]
    cumulative_frequencies[i, "A"] <- sum(sub_seq == "A") / i
    cumulative_frequencies[i, "G"] <- sum(sub_seq == "G") / i
    cumulative_frequencies[i, "C"] <- sum(sub_seq == "C") / i
    cumulative_frequencies[i, "T"] <- sum(sub_seq == "T") / i
  }
  
  return(cumulative_frequencies)
}

# 将序列转换为累计频率特征矩阵
encoded_seqs_anf <- lapply(data$sequence, function(seq) {
  cumulative_frequencies <- calculate_cumulative_frequencies(seq)
  # 展平矩阵，将每行变成一个特征向量
  return(as.vector(t(cumulative_frequencies)))
})

# 生成特征矩阵（样本数 × (max_length × 4)）
feature_matrix_anf <- do.call(rbind, encoded_seqs_anf)
#创建列名，50可随序列长度改变
# 创建一个空的字符向量，用于存储列名
col_names_anf <- c()

# 循环创建列名
for (i in 1:max_length) {
  col_names_anf <- c(col_names_anf, paste("A_anf", i, sep = "_"), 
                     paste("G_anf", i, sep = "_"), 
                     paste("C_anf", i, sep = "_"), 
                     paste("T_anf", i, sep = "_"))
}
# 假设 feature_matrix_oh 是一个 50 行的数据框
colnames(feature_matrix_anf) <- col_names_anf
data=cbind(data,feature_matrix_anf)

#PseDNC
dnacheck = function (x) {
  
  DNADict = c("A","G","C","T")
  return(all(strsplit(x, split = "")[[1]]%in%DNADict)) 
  
}
make_kmer_index = function (k, alphabet = 'ACGT') {
  
  dict = unlist(strsplit(alphabet, ''))
  make_index = list()
  temp = rep(dict, k)
  make_index = sort(unique(apply(expand.grid(split(temp,
                                                   rep(1:k,each=4))), 1, paste, collapse = '')))
  return(make_index)
  
}

extrPseDNC = function (x, lambda = 2, w = 0.05, 
                       normalize = FALSE, customprops = NULL) {
  
  if (dnacheck(x) == FALSE) 
    stop("x has unrecognized  type !")
  
  Didx = data.frame(AA = c(0.06, 0.5, 0.27, 1.59, 0.11, -0.11),
                    AC = c(1.50, 0.50, 0.80, 0.13, 1.29, 1.04),
                    AG = c(0.78, 0.36, 0.09, 0.68, -0.24, -0.62),
                    AT = c(1.07, 0.22, 0.62, -1.02, 2.51, 1.17),
                    CA = c(-1.38, -1.36, -0.27, -0.86, -0.62, -1.25),
                    CC = c(0.06, 1.08, 0.09, 0.56, -0.82, 0.24),
                    CG = c(-1.66, -1.22, -0.44, -0.82, -0.29, -1.39),
                    CT = c(0.78, 0.36, 0.09, 0.68, -0.24, -0.62),
                    GA = c(-0.08, 0.5, 0.27, 0.13, -0.39, 0.71),
                    GC = c(-0.08, 0.22, 1.33, -0.35, 0.65, 1.59),
                    GG = c(0.06, 1.08, 0.09, 0.56, -0.82, 0.24),
                    GT = c(1.50, 0.50, 0.80, 0.13, 1.29, 1.04),
                    TA = c(-1.23, -2.37, -0.44, -2.24, -1.51, -1.39),
                    TC = c(-0.08, 0.5, 0.27, 0.13, -0.39, 0.71),
                    TG = c(-1.38, -1.36, -0.27, -0.86, -0.62, -1.25),
                    TT = c(0.06, 0.5, 0.27, 1.59, 0.11, -0.11))
  
  if (!is.null(customprops))
    Didx = rbind(Didx, customprops)
  
  n = dim(Didx)[1]
  Ddict = make_kmer_index(k = 2) 
  if (normalize) { 
    H = matrix(ncol = 16, nrow = n)
    Didx = as.matrix(Didx)
    for (i in 1:n) H[i, ] = (Didx[i, ] - mean(Didx[i, ]))/(sqrt(sum((Didx[i,] - 
                                                                       mean(Didx[i, ])) ^ 2) / 16))
    H = round(H, 3)
    colnames(H) = Ddict
  }
  else H = Didx
  
  Theta = vector("list", lambda)
  for (i in 1:lambda) Theta[[i]] = vector("list", n)
  xSplit = strsplit(x, split = "")[[1]]
  xPaste = paste0(xSplit[1:length(xSplit) - 1], xSplit[2:length(xSplit)])
  N = length(xPaste)
  for (i in 1:lambda) {
    temp=c()
    for (j in 1:(N - i)) {
      temp = append(temp, mean((H[, xPaste[j]] - H[, xPaste[j + i]]) ^ 2))
    }
    Theta[[i]] = temp
  }
  theta = sapply(Theta, mean)
  
  fc = summary(factor(xPaste, levels = Ddict), maxsum = 16)
  fc = fc/sum(fc)
  
  Xc1 = fc/(1 + (w * sum(theta)))
  names(Xc1) = paste("Xc1.", names(Xc1), sep = "")
  Xc2 = (w * theta)/(1 + (w * sum(theta)))
  names(Xc2) = paste("Xc2.lambda.", 1:lambda, sep = "")
  Xc = c(Xc1, Xc2)
  
  Xc = round(Xc, 3)
  return (Xc)
  
}
# 将序列转换为PseDNC特征矩阵
encoded_seqs_psednc <- lapply(data$sequence, function(seq) {
  psednc <- extrPseDNC(seq)
  return(as.vector(t(psednc)))
})

# 合并所有特征向量为一个矩阵
feature_matrix_psednc <- do.call(rbind, encoded_seqs_psednc)

# 设置列名为每个向量的名字
# 假设每个 psednc 向量的名字相同，可以取第一个向量的名字作为列名
colnames(feature_matrix_psednc) <- names(extrPseDNC(data$sequence[1]))
data=cbind(data,feature_matrix_psednc)

# 假设 data 是你原始的数据框
data_n <- data[1:1000, ]  # 取前 1000 行
data_p <- data[1001:2000, ]  # 取后 1000 行

# 加载必要的包
library(e1071)
library(caret)  # 用于数据拆分和评估模型

# 假设 data_sample 已经加载
# 提取特征和目标变量
#features <- data_sample[, 9:208]
#target <- as.factor(data_sample$type)  # 确保目标变量是因子类型（分类任务）
features_n <- data_n[, -c(1:8)]
features_p <- data_p[, -c(1:8)]

oh_n=features_n[, c(1:(max_length*4))]
eiip_n=features_n[, c((max_length*4+1):(max_length*5))]
nac_n=features_n[, c((max_length*5+5):(max_length*5+20))]
cp_n=features_n[, c((max_length*5+21):(max_length*8+20))]
anf_n=features_n[, c((max_length*8+21):(max_length*12+20))]
pseknc_n=features_n[, -c(1:max_length*12+20)]

oh_p=features_p[, c(1:(max_length*4))]
eiip_p=features_p[, c((max_length*4+1):(max_length*5))]
nac_p=features_p[, c((max_length*5+5):(max_length*5+20))]
cp_p=features_p[, c((max_length*5+21):(max_length*8+20))]
anf_p=features_p[, c((max_length*8+21):(max_length*12+20))]
pseknc_p=features_p[, -c(1:max_length*12+20)]

f <- list(oh_n, eiip_n, nac_n, cp_n, anf_n, pseknc_n)
feature_combinations_n <- list()
for (i in 1:length(f)) {
  # 获取长度为 i 的所有特征组合
  combs <- combn(f, i, simplify = FALSE)
  
  # 将每个组合应用 cbind，存储在 feature_combinations 中
  for (comb in combs) {
    feature_combinations_n <- append(feature_combinations_n, list(do.call(cbind, comb)))
  }
}
f_p <- list(oh_p, eiip_p, nac_p, cp_p, anf_p, pseknc_p)
feature_combinations_p <- list()
for (i in 1:length(f_p)) {
  # 获取长度为 i 的所有特征组合
  combs <- combn(f_p, i, simplify = FALSE)
  
  # 将每个组合应用 cbind，存储在 feature_combinations 中
  for (comb in combs) {
    feature_combinations_p <- append(feature_combinations_p, list(do.call(cbind, comb)))
  }
}

m <- list("oh", "eiip", "nac", "cp", "anf", "pseknc")
m_n <- list()
for (i in 1:length(m)) {
  # 获取长度为 i 的所有特征组合
  combs <- combn(m, i, simplify = FALSE)
  
  # 将每个组合应用 cbind，存储在 feature_combinations 中
  for (comb in combs) {
    m_n <- append(m_n, list(do.call(cbind, comb)))
  }
}


# 获取当前特征组合
cb_n <- feature_combinations_n[1]
cb_n <- as.data.frame(cb_n)
cb_p <- feature_combinations_p[1]
cb_p <- as.data.frame(cb_p)

# 随机分割训练集和测试集
set.seed(123)
train_index <- sample(1:nrow(cb_n), size = 0.8 * nrow(cb_n))
train_n <- cb_n[train_index, ]
test_n <- cb_n[-train_index, ]
train_p <- cb_p[train_index, ]
test_p <- cb_p[-train_index, ]

# 添加 'type' 列标记
train_n$type <- "n"
train_p$type <- "p"
test_n$type <- "n"
test_p$type <- "p"

# 合并训练集和测试集
train_data = rbind(train_n, train_p)
test_data = rbind(test_n, test_p)

# 准备训练数据
train_features <- train_data[, !(colnames(train_data) == "type")]
train_labels <- train_data$type
test_features <- test_data[, !(colnames(test_data) == "type")]
test_labels <- test_data$type

set.seed(123)
# 加载所需的包
library(caret)
library(randomForest)
library(pROC)

set.seed(123)
# 训练控制参数
train_control <- trainControl(
  method = "cv",             # 交叉验证
  number = 5,                # 5 折交叉验证
  classProbs = TRUE,         # 启用类别概率估计
  summaryFunction = twoClassSummary,  # 使用二分类问题的评估指标（如AUC, Sensitivity等）
  verboseIter = TRUE         # 输出训练过程
)

# 使用随机森林训练模型
rf_model <- train(
  x = train_features,        # 输入特征
  y = train_labels,          # 输入标签
  method = "rf",             # 使用随机森林
  trControl = train_control, # 训练控制
  probability = TRUE         # 启用概率估计
)

# 使用训练好的模型进行预测
predictions <- predict(rf_model, test_features, type = "raw")
predictions_prob <- predict(rf_model, test_features, type = "prob")

# 计算ROC曲线
roc_curve <- roc(test_labels, predictions_prob[, 2])
auc_results <- rbind(auc_results, data.frame(Combination = paste("length", max_length), AUC = auc(roc_curve)))
}
write.xlsx(auc_results, "D:\\Users\\Administrator\\Desktop\\machine learning\\new\\TKIs\\down\\auc_length.xlsx")