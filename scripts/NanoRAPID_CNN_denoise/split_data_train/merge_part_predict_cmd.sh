code=/disk1/work/zehui/gpu_public/public_data_result/NBT_PORE-cupine-ONT_result/model/rawsignal_train/result/denoise/code/split_data_train


wd=$1
sample=$2
outdir=$wd/filter
mkdir -p $outdir


for name in $(cat $sample);do
  out=$outdir/${name}.filter.txt
  echo $wd/train/$name $out
  python -u $code/merge_part_predict.py $wd/test/$name > $out

done

