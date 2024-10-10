import sys,os,time,random
import numpy as np

def nowtime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def get_embed_index_dict(file):
    embed_index_dict = {}
    for line in open(file):
        lineL = line.strip().split(' ')
        motif = lineL[-2]
        index = lineL[1]
        embed_index_dict[index] = motif
    return embed_index_dict
#----

def write_table(input_name, output_name, embed_dict, counts):
    sig_memmap_filename=f"{output_name}_signal.npy"
    seq_memmap_filename=f"{output_name}_sequence.npy"
    ext_memmap_filename=f"{output_name}_extra.npy"
    sig_memmap = np.memmap(sig_memmap_filename, dtype="float32", mode='w+', shape=(counts, 400))
    seq_memmap = np.memmap(seq_memmap_filename, dtype="int8", mode='w+', shape=(counts, 400, 5))
    ext_memmap = np.memmap(ext_memmap_filename, dtype="<U80", mode='w+', shape=(counts))

    sig_data_list,seq_data_list,ext_data_list = [],[],[]
    line_number = 0
    memmap_start = 0
    with open(f"{input_name}_sequence.txt", "r") as f_seq, \
         open(f"{input_name}_signal.txt", "r") as f_sig, \
         open(f"{input_name}_extra.txt", "r") as f_extra:
        for n, (line_seq, line_sig, line_extra) in enumerate(zip(f_seq, f_sig, f_extra)):
            line_number += 1
            site_name = line_extra.strip()
            ext_data_list.append(site_name)

            row = []
            sig_feature_list = line_sig.strip().split("\t")[1:]
            row = [float(x) for x in sig_feature_list]
            sig_data_list.append(row)

            row = []
            seq_feature_list = line_seq.strip().split("\t")[1:]
            for index in seq_feature_list:
                k = embed_dict[index]
                kmer = [int(x) for x in k]
                row.append(kmer)
            seq_data_list.append(row)

            if line_number % 500000 == 0:
                memmap_end = memmap_start + len(ext_data_list)
                seq_data_array = np.array(seq_data_list, dtype="int8")
                seq_memmap[memmap_start:memmap_end] = seq_data_array[:]
                sig_data_array = np.array(sig_data_list, dtype="float32")
                sig_memmap[memmap_start:memmap_end] = sig_data_array[:]
                ext_data_array = np.array(ext_data_list, dtype="<U80")
                ext_memmap[memmap_start:memmap_end] = ext_data_array[:]
                print(f"{nowtime()} | loaded {memmap_end}/{counts} lines")
                #refresh
                memmap_start = memmap_end
                sig_data_list,seq_data_list,ext_data_list = [],[],[]

    if ext_data_list:
        memmap_end = memmap_start + len(ext_data_list)
        seq_data_array = np.array(seq_data_list, dtype="int8")
        seq_memmap[memmap_start:memmap_end] = seq_data_array[:]
        sig_data_array = np.array(sig_data_list, dtype="float32")
        sig_memmap[memmap_start:memmap_end] = sig_data_array[:]
        ext_data_array = np.array(ext_data_list, dtype="<U80")
        ext_memmap[memmap_start:memmap_end] = ext_data_array[:]
        print(f"{nowtime()} |loaded {memmap_end}/{counts} lines")

    seq_memmap.flush()
    sig_memmap.flush()
    ext_memmap.flush()


#----------
inputdir = sys.argv[1]
outputdir = sys.argv[2]
motif = sys.argv[3]
counts = int(sys.argv[4])
#---------
embed_file="/disk1/work/zehui/index/data/embed.txt"
embed_motif_dict = get_embed_index_dict(embed_file)
input_name=f'{inputdir}/{motif}'
output_name=f'{outputdir}/{motif}'
write_table(input_name, output_name, embed_motif_dict, counts)
print(f"{nowtime()} | finish txt2npy")
