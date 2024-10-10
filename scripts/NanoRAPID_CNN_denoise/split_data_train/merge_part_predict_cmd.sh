code=/NanoRAPID/scripts/NanoRAPID_CNN_denoise/split_data_train


wd=$1
sample=$2
outdir=$wd/filter
mkdir -p $outdir


for name in $(cat $sample);do
  out=$outdir/${name}.filter.txt
  echo $wd/train/$name $out
  python -u $code/merge_part_predict.py $wd/test/$name > $out

done

