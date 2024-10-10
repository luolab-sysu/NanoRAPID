code=/NanoRAPID/scripts/NanoRAPID_CNN_denoise/split_data_train

wd=$1
sample=$2
gpu=${3:-1}
outdir=$wd/test

for name in $(cat $sample ); do
  echo $name
  out=$outdir/$name
  mkdir -p $out
  python -u $code/predict_test.py -d $wd/data/part1_rep1_$name  -m $wd/train/$name/part2_rep1_${name}_model_10.pth.tar -g $gpu -o $out/part1_rep1_test.txt
  python -u $code/predict_test.py -d $wd/data/part2_rep1_$name  -m $wd/train/$name/part1_rep1_${name}_model_10.pth.tar -g $gpu -o $out/part2_rep1_test.txt
  python -u $code/predict_test.py -d $wd/data/part1_rep2_$name  -m $wd/train/$name/part2_rep2_${name}_model_10.pth.tar -g $gpu -o $out/part1_rep2_test.txt
  python -u $code/predict_test.py -d $wd/data/part2_rep2_$name  -m $wd/train/$name/part1_rep2_${name}_model_10.pth.tar -g $gpu -o $out/part2_rep2_test.txt
  python -u $code/predict_test.py -d $wd/data/part1_rep3_$name  -m $wd/train/$name/part2_rep3_${name}_model_10.pth.tar -g $gpu -o $out/part1_rep3_test.txt
  python -u $code/predict_test.py -d $wd/data/part2_rep3_$name  -m $wd/train/$name/part1_rep3_${name}_model_10.pth.tar -g $gpu -o $out/part2_rep3_test.txt
done

