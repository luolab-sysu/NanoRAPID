import sys,os,time,random
import numpy as np

def get_embed_index_dict(file):
    embed_index_dict = {}
    for line in open(file):
        lineL = line.strip().split(' ')
        motif = lineL[-2]
        index = lineL[1]
        embed_index_dict[index] = motif
    return embed_index_dict
#----
def pos_get_filter(filter_file):
    filter_dic = {}
    #--
    with open(filter_file) as f:
        for line in f:
            lineL = line.strip().split("\t")
            reads_site_index = lineL[0]
            infoL = lineL[-1].split("_")
            seq_motif = infoL[-1]
            motif = infoL[-2]
            filter_dic.setdefault(motif,{})
            filter_dic[motif].setdefault(seq_motif,[])
            filter_dic[motif][seq_motif].append(reads_site_index)
    print("load postive site")
    #---
    select_dic = {}
    for motif in filter_dic:
        select_list = []
        remian_list = []
        for seq_motif in filter_dic[motif]:
            my_list = filter_dic[motif][seq_motif]
            # 打乱列表顺序
            random.shuffle(my_list)
            select_list += my_list[:1000]
            remian_list += my_list[1000:]
        if len(select_list) < 16*1000:
            #print("append remian site")
            cutline = 16*1000 - len(select_list)
            random.shuffle(remian_list)
            select_list += remian_list[:cutline]
        select_dic[motif] = select_list
    print("select postive site")
    #print(select_dic["GGACT"][:5])
    return select_dic
#----
def get_feature_dic(filename,flag=0,filter=False):
    feature_d = {}
    with open(filename) as f:
        for line in f:
            lineL = line.strip().split("\t")
            site_name = lineL[0]
            #print(site_name)
            if filter:
                infoL = site_name.split("|")
                motif = infoL[-2]
                if site_name not in select_dic[motif]:
                    continue
            index = f'{site_name}|{flag}'
            feature_d[index] = flag
    return feature_d
#-----
def merge_ctrl_test(test_extra_file, ctrl_extra_file):
    test_flag = 1
    test_feature_d = get_feature_dic(test_extra_file,test_flag,True)
    test_counts = len(test_feature_d.keys())
    print(f'{test_extra_file}\t{test_counts}')

    ctrl_flag = 0
    ctrl_feature_d = get_feature_dic(ctrl_extra_file,ctrl_flag)
    ctrl_counts = len(ctrl_feature_d.keys())
    print(f'{ctrl_extra_file}\t{ctrl_counts}')

    random_len = max([test_counts*5,50000])
    random_len = min([random_len,ctrl_counts])
    random_ctrl_name = random.sample(ctrl_feature_d.keys(), random_len)
    name_list = list(test_feature_d.keys()) + random_ctrl_name
    counts = len(name_list)
    print(f"out_counts\t{counts}")
    print(name_list[:5])
    return name_list
#----
def write_table(input_name_L, out_list, embed_dict):
    counts = len(out_list)
    sig_memmap = np.memmap(sig_memmap_filename, dtype="float32", mode='w+', shape=(counts, 400))
    seq_memmap = np.memmap(seq_memmap_filename, dtype="int8", mode='w+', shape=(counts, 400, 5))
    ext_memmap = np.memmap(ext_memmap_filename, dtype="<U80", mode='w+', shape=(counts))

    sig_data_list,seq_data_list,ext_data_list = [],[],[]
    line_number = 0
    memmap_start = 0
    for input_name in input_name_L:
        with open(f"{input_name}_sequence.txt", "r") as f_seq, \
             open(f"{input_name}_signal.txt", "r") as f_sig, \
             open(f"{input_name}_extra.txt", "r") as f_extra:
            for n, (line_seq, line_sig, line_extra) in enumerate(zip(f_seq, f_sig, f_extra)):
                site_name = line_extra.strip()
                if f'{site_name}' not in out_list:
                    continue
                site_name_label = out_list[site_name] 
                line_number += 1
                ext_data_list.append(site_name_label)

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
                    print(f"loaded {memmap_end}/{counts} lines")
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
            print(f"loaded {memmap_end}/{counts} lines")

    seq_memmap.flush()
    sig_memmap.flush()
    ext_memmap.flush()


def write_part_npy(name_list,out_flag):
    half_len = int(len(name_list)/2)
    print(f"half_len\t{half_len}")
    random.shuffle(name_list)

    out_list_part1 = {i[:-2]:i for i in name_list[:half_len]}
    out_list_part2 = {i[:-2]:i for i in name_list[half_len:]}

    global sig_memmap_filename,seq_memmap_filename,ext_memmap_filename
    sig_memmap_filename = f"{outdir}/part1_{out_flag}_{motif}_signal.npy"
    seq_memmap_filename = f"{outdir}/part1_{out_flag}_{motif}_sequence.npy"
    ext_memmap_filename = f"{outdir}/part1_{out_flag}_{motif}_extra.npy"
    write_table(input_name_L,out_list_part1, embed_motif_dic)

    sig_memmap_filename = f"{outdir}/part2_{out_flag}_{motif}_signal.npy"
    seq_memmap_filename = f"{outdir}/part2_{out_flag}_{motif}_sequence.npy"
    ext_memmap_filename = f"{outdir}/part2_{out_flag}_{motif}_extra.npy"
    write_table(input_name_L,out_list_part2, embed_motif_dic)

def write_npy():
    #----------
    input_name_L = [f'{ctrl_path}_{motif}',f'{test_path}_{motif}']
    sig_memmap_filename = f"{outdir}/{motif}_signal.npy"
    seq_memmap_filename = f"{outdir}/{motif}_sequence.npy"
    ext_memmap_filename = f"{outdir}/{motif}_extra.npy"
    write_table(input_name_L,out_name_list, embed_motif_dic)
    print(f"{motif} finish")


#---------
outdir = sys.argv[1]
motif =sys.argv[2]
#---------
embed_file="/disk1/work/zehui/index/data/embed.txt"
embed_motif_dic = get_embed_index_dict(embed_file)
#----------
filter_file = "/disk1/work/zehui/gpu_public/public_data_result/NBT_PORE-cupine-ONT_result/model/PORE-cupine-SVM/split_data/pos_site.txt"
select_dic = pos_get_filter(filter_file)
#----------
wd = "/disk1/work/zehui/gpu_public/public_data_result/NBT_PORE-cupine-ONT_result/model/rawsignal_train/data"
ctrl_path = f"{wd}/unmod/motif_single_rawsignal/unmod_train"
test_path = f"{wd}/mod/motif_single_rawsignal/mod"
test_extra_file = f'{test_path}_{motif}_extra.txt'
ctrl_extra_file = f'{ctrl_path}_{motif}_extra.txt'
out_name_list = merge_ctrl_test(test_extra_file, ctrl_extra_file)
#----------
input_name_L = [f'{ctrl_path}_{motif}',f'{test_path}_{motif}']
write_part_npy(out_name_list,"rep1")
write_part_npy(out_name_list,"rep2")
write_part_npy(out_name_list,"rep3")

