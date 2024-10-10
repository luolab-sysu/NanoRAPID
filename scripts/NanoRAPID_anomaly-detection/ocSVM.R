#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(parallel))


# 命令行参数解析
option_list = list(
  make_option(c("-u", "--unmod"), type="character", default=NULL, help="unmodified RData", metavar="character"),
  make_option(c("-m", "--mod"), type="character", default=NULL, help="Modified RData", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$unmod) | is.null(opt$output)|is.null(opt$mod) ){
  print_help(opt_parser)
  stop("Input files and output name must be supplied.", call.=FALSE)
}


# 函数：处理单个位点
process_position <- function(b, modified_data, unmodified_data, output_file) {
  modified_data_at_position_b <- modified_data[position == b, .(read_name, event_level_mean, event_stdv, count, reference_kmer)]
  modified_data_at_position_b[, `:=`(event_level_mean = (event_level_mean), event_stdv = log(event_stdv), count = log(count))]

  kmer <- modified_data_at_position_b$reference_kmer[1]

  unmodified_data_at_position_b <- unmodified_data[position == b, .(event_level_mean, event_stdv, count)]
  unmodified_data_at_position_b[, `:=`(event_level_mean = (event_level_mean), event_stdv = log(event_stdv), count = log(count))]

  if (nrow(unmodified_data_at_position_b) < 5) {
    return(NULL)
  }

  svm_model <- svm(unmodified_data_at_position_b[, .(event_level_mean, event_stdv)], y = NULL,
                   type = 'one-classification', nu = 0.01, gamma = 0.0009, kernel = "radial")

  svm_predtest <- predict(svm_model, modified_data_at_position_b[, .(event_level_mean, event_stdv)])

  outdata <- modified_data_at_position_b[, predicted_label := ifelse(svm_predtest, 0, 1)]
  outdata[, position := b]

  fwrite(outdata, file = output_file, sep = "\t", append = TRUE, col.names = !file.exists(output_file)) # 使用 fwrite 提高效率

  return(list(predictions = svm_predtest, kmer = kmer, n_unmod = nrow(unmodified_data_at_position_b)))
}



# 加载数据
modified_data <- as.data.table(get(load(opt$mod)))
setnames(modified_data, "site", "read_name") # 直接修改列名
modified_data[event_stdv == 0, event_stdv := 0.01]

unmodified_data <- as.data.table(get(load(opt$unmod)))
setnames(unmodified_data, "site", "read_name") # 直接修改列名
unmodified_data[event_stdv == 0, event_stdv := 0.01]


# 预处理和初始化
pos <- max(modified_data$position)
pos_m <- sort(unique(modified_data$position))
gene_name <- unique(modified_data$contig)

output_file <- paste(opt$output,"_test.txt",sep="")


# 并行处理
numCores <- detectCores() - 1  # 使用大部分核心，留一个核心给系统
cl <- makeCluster(numCores)
clusterExport(cl, c("process_position", "modified_data", "unmodified_data", "output_file", "opt", "svm")) # 导出需要的变量和函数

results <- parLapply(cl, pos_m, process_position, modified_data = modified_data, unmodified_data = unmodified_data, output_file = output_file)

stopCluster(cl)


# 后续处理
# 初始化结果矩阵
n_unique_reads <- uniqueN(modified_data$read_name)
mod_mat <- matrix(NA, nrow = n_unique_reads, ncol = pos + 1)
rownames(mod_mat) <- unique(modified_data$read_name)
colnames(mod_mat) <- 0:pos
unmod_mat <- matrix(0, ncol = 1, nrow = pos + 1)
nt_mat <- matrix("N", ncol = 1, nrow = pos + 1)

# 汇总结果
for (i in seq_along(results)) {
  if (!is.null(results[[i]])) {
    b <- pos_m[i]
    predictions <- results[[i]]$predictions
    read_names_at_b <- modified_data[position == b, read_name]
    mod_mat[read_names_at_b, (b + 1)] <- 0
    mod_mat[read_names_at_b[!predictions], (b + 1)] <- 1
    unmod_mat[(b + 1), ] <- results[[i]]$n_unmod
    nt_mat[(b + 1), ] <- substr(results[[i]]$kmer, 3, 3)
  }
}

# 转换为数值矩阵
mode(mod_mat) <- "numeric"

# 计算修饰位点数量和比例
mod_strands <- apply(mod_mat, 2, function(x) sum(!is.na(x)))
dat.results <- colSums(mod_mat, na.rm = TRUE)
dat.percentage <- dat.results / mod_strands

# 创建输出数据框
gene_mat <- matrix(gene_name, ncol = 1, nrow = (pos + 1))
dat.output <- data.table(
  Gene = gene_mat,
  Position = 3:(pos + 3),
  NT = nt_mat,
  Mod_percentage = dat.percentage,
  Mod_number = dat.results,
  Mod_strands = mod_strands,
  Training_strands = unmod_mat
)

# 输出结果
output_csv_file <- paste0("Output_", gene_name, "_", opt$output, ".csv")
fwrite(dat.output, output_csv_file, sep = "\t")


print("done")