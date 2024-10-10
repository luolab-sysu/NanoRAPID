
# 获取所有命令行参数
args <- commandArgs(trailingOnly = TRUE)
#print(args)
filename <- args[1]
print(filename)

# 读取数据
data <- read.table(file = filename,header = F)
head(data)

#长宽表转换
library(reshape2)
trans_df <- data
head(trans_df)
wide_df <- dcast(trans_df,V1 ~ V2)

df <- wide_df[c(-1)]
rownames(df) <- wide_df$V1
head(df[c(1:5)])

## 计算基因整体修饰比例
reads_counts <- dim(df)[1]
non_na_elements <- sum(!is.na(df))
print(non_na_elements)
mod_elements <- sum(df, na.rm = TRUE)
print(mod_elements)
mod_ratio <- mod_elements/non_na_elements
print(mod_ratio)


## 筛选去除背景位点
# 创建一个空列表存储符合条件的列名
selected_columns <- list()
# 循环遍历数据框的列
for (col_name in colnames(df)) {
  # 计算每列非NA元素的和
  sum_non_na <- sum(df[[col_name]], na.rm = TRUE)
  # 计算每列非NA元素的数量
  count_non_na <- sum(!is.na(df[[col_name]]))
  # 计算比例
  ratio <- sum_non_na / count_non_na
  # 如果比例大于mod_ratio，则将该列名存储到列表中
  if (ratio > mod_ratio && count_non_na > reads_counts*0.75) {
    selected_columns[[col_name]] <- df[[col_name]]
  }
}


# 将符合条件的列提取出来
result_df <- as.data.frame(selected_columns)
row.names(result_df) <- row.names(df)
colnames(result_df) <- gsub("X","site_",colnames(result_df))


# 填充缺省值
result_df <- replace(result_df, is.na(result_df), 0)
print(head(result_df[c(1:10)]))

write.table(
  result_df,file = gsub("out_","filter_wide_",filename),
  row.names = T,col.names = T,quote = F,sep = "\t"
)


