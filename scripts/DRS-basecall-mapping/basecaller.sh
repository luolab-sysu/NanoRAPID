#!/bin/bash

check_program() {
    if !(command -v "$1" &> /dev/null); then
        echo "Error: $1 not found. Make sure it is installed and in your PATH."
        exit 1
    fi
}
# 检查程序是否可运行
check_program "guppy_basecaller"


# 设置输入路径

wd="$1"
data="$wd/data"
gpu="${2:-3}"
# 检查参数是否为空
if [ -z $wd ] && [ -z $gpu ]; then
    echo "missing argvs" # 打印错误信息
    exit                # 退出脚本
fi


# 运行 Guppy Basecaller
echo "---- basecaller ----"
date
echo "--------"
guppy_basecaller -i $data/fast5 -s $data -c rna_r9.4.1_70bps_hac.cfg -x "cuda:${gpu}"

echo "---- END ----"
date
echo "--------"
