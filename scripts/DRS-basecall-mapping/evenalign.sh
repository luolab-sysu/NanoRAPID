#!/bin/bash

check_program() {
    if !(command -v "$1" &> /dev/null); then
        echo "Error: $1 not found. Make sure it is installed and in your PATH."
        exit 1
    fi
}
# 检查程序是否可运行
check_program "nanopolish"

wd="$1"
data="$wd/data"
fa="$2"
# 检查参数是否为空
if [ -z $wd ] && [ -z $fa ]; then
    echo "missing argvs" # 打印错误信息
    exit                # 退出脚本
fi
# 如果 $wd 不为空，则执行后续代码
echo "Working directory: $wd"

# 生成输出路径
mkdir -p $wd/result/evenalign

echo "----- nanopolish index ------"
date
echo "-----------"
cd $data
# 建立index
nanopolish index --directory=$data/fast5 $data/merge.fastq

echo "----- nanopolish eventalign ------"
date
echo "-----------"
for file in $wd/result/bam/tmp/*.bam
do
{
echo $file
info=(${file//// })
echo ${info[-1]}
nanopolish eventalign --reads $data/*.fastq --bam $file --genome $fa -t 15 --scale-events --samples --signal-index --summary $wd/result/evenalign/${info[-1]%%.bam}_summary.txt --print-read-names > $wd/result/evenalign/${info[-1]%%.bam}_evenalign.txt
} &
done

echo "----- END ------"
date
echo "-----------"
