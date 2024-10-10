#!/bin/bash

check_program() {
    if !(command -v "$1" &> /dev/null); then
        echo "Error: $1 not found. Make sure it is installed and in your PATH."
        exit 1
    fi
}
# 检查程序是否可运行
#check_program "NanoPlot"
check_program "minimap2"
#check_program "java -jar picard.jar"
check_program "samtools"

# 设置输入路径

wd="$1"
data="$wd/data"
fa="/disk2/work/zehui/NBT_ONT_result/Tetra_test/ref/ribosomal_RNA_self-splicing_intron.fa"
# 检查参数是否为空
if [ -z $wd ] && [ -z $fa ]; then
    echo "missing argvs" # 打印错误信息
    exit                # 退出脚本
fi

# 生成输出路径
#mkdir -p $wd/result/qc
mkdir -p $wd/result/bam

echo "----- nanoplot ------"
date
echo "-----------"
# 合并fastq
cat $data/pass/*fastq $data/fail/*fastq > $data/merge.fastq 
# 数据质控
#NanoPlot --summary $data/sequencing_summary.txt -t 16 -o $wd/result/qc -p sample_name 

echo "----- minimap2 ------"
date
echo "-----------"
minimap2 -ax map-ont -k 14 $fa -t 25 --secondary=no $data/merge.fastq -o $wd/result/bam/sample_name.sam
samtools view -@ 30 -F 2048 -F 4 -b $wd/result/bam/sample_name.sam | samtools sort -O BAM -@ 20  -o $wd/result/bam/sample_name.bam
samtools index -@ 16 $wd/result/bam/sample_name.bam
samtools flagstat -@ 16 $wd/result/bam/sample_name.sam > $wd/result/bam/sample_name.maplog

echo "----- PicardSplit ------"
date
echo "-----------"
#把bam文件拆成多个（25），让后续步骤可以并行操作
cd $wd/result/bam
mkdir -p tmp
java -jar /public/work/zehui/software/picard.jar SplitSamByNumberOfReads --INPUT sample_name.bam --SPLIT_TO_N_FILES 25 --OUTPUT tmp 
for bamfile in tmp/*bam
do
{
samtools index $bamfile
} &
done

echo "----- END ------"
date
echo "-----------"
