wd=$1
list=$2
gpu1=${3:-2}
gpu2=$(($gpu1+1))

code=/NanoRAPID/scripts/NanoRAPID_CNN_denoise/split_data_train

# 检查 $wd 变量是否已设置
if [ -z $wd ] || [ -z $list ]; then
  echo "missing wd"
  exit 1  # 使用非零退出码表示错误
fi
# 如果 $wd 不为空，则执行后续代码
echo "Working directory: $wd; GPU: ${gpu}"
#exit

data=$wd/data
for name in $(cat $list)
do
  out=$wd/train/$name
  mkdir -p $out
  echo `date` "|" "$name train start"
  {
  python -u $code/split_part_train.py $data part1_rep1_${name} $out $gpu1 > $out/part1_rep1.log &
  python -u $code/split_part_train.py $data part2_rep1_${name} $out $gpu2 > $out/part2_rep1.log &
  python -u $code/split_part_train.py $data part1_rep2_${name} $out $gpu1 > $out/part1_rep2.log &
  python -u $code/split_part_train.py $data part2_rep2_${name} $out $gpu2 > $out/part2_rep2.log &
  python -u $code/split_part_train.py $data part1_rep3_${name} $out $gpu1 > $out/part1_rep3.log &
  python -u $code/split_part_train.py $data part2_rep3_${name} $out $gpu2 > $out/part2_rep3.log &
  }
  wait
  echo `date` "|" "$name train finish"
done

