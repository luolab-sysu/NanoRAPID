import numpy as np
import os, sys, gzip
import time

input_dir = sys.argv[1]
output_dir = sys.argv[2]
sample_list_file = sys.argv[3]
outname = sys.argv[4]

def nowtime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def split_motif_write2txt(batch,idx_list):
    with open(f"{input_dir}/{batch}_sequence.txt", "r") as f_seq, \
         open(f"{input_dir}/{batch}_signal.txt", "r") as f_sig, \
         open(f"{input_dir}/{batch}_extra.txt", "r") as f_extra:
        for n, (line_seq, line_sig, line_extra) in enumerate(zip(f_seq, f_sig, f_extra)):
            motif = idx_list[n]
            with open(f"{output_dir}/{outname}_{motif}_sequence.txt", "a") as fh_seq, \
                    open(f"{output_dir}/{outname}_{motif}_signal.txt", "a") as fh_sig, \
                    open(f"{output_dir}/{outname}_{motif}_extra.txt", "a") as fh_extra:
                fh_seq.write(line_seq)
                fh_sig.write(line_sig)
                fh_extra.write(line_extra)
            if n % 500000 == 0:
                print(f'{nowtime()} | {batch} load {n} lines ')
                

def read_batch_data():
    chunks = {}
    idxs = {}
    motifs = set()
    with open(f'{sample_list_file}', "r") as f:
        for line in f:
            print(line.strip())
            info = line.strip().split(" ")
            #batch_chunk = int(info[1])
            batch_chunk = 1
            batch = info[0]
            chunks[batch] = batch_chunk

            idx = []
            # 提取序列号和对应的motif
            with open(f"{input_dir}/{batch}_extra.txt", "r") as f_extra:
                for line_extra in f_extra:
                    infoL = line_extra.strip().split("|")
                    motif = infoL[-2]
                    idx.append(motif)
                    motifs.add(motif)
            idxs[batch] = idx

    batchs = chunks.keys()
    motif_list = list(motifs) 

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(batchs)
    print(len(motif_list),motif_list)
    for batch in batchs:
        idx_list = idxs[batch]
        split_motif_write2txt(batch,idx_list)
        print(f'{nowtime()} | {batch} finish')
    
    print(f"{nowtime()}|All finish")

if __name__ == '__main__':
    read_batch_data()
