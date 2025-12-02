#!/usr/bash
set -e

# 设置输入路径
wd=$1
py=rawsignal_extra_singlemotif.py
data=$wd/data
code=./DRS-Features
fa=reference

# 检查 $wd 变量是否已设置
if [ -z $wd ] && [ -z $py ]; then
  echo "missing wd"
  exit 1  # 使用非零退出码表示错误
fi

# 如果 $wd 不为空，则执行后续代码
echo "Working directory: $wd"


echo 'load features' `date`
mkdir -p $wd/result/temp_single_features
cd $wd/result/temp_features/
batchs=(shard_0001 shard_0002 shard_0003 shard_0004 shard_0005 shard_0006 shard_0007 shard_0008 shard_0009 shard_0010 shard_0011 shard_0012 shard_0013 shard_0014 shard_0015 shard_0016 shard_0017 shard_0018 shard_0019 shard_0020 shard_0021 shard_0022 shard_0023 shard_0024 shard_0025)
for batch in ${batchs[@]}
do
{
python -u $code/$py $wd/result/evenalign/${batch}_evenalign.txt $wd/result/temp_single_features $batch
} &
done
wait
ls *_extra.txt|sed "s/_extra.txt//g" > sample_list.txt
