
import sys,os,time
import numpy as np

def nowtime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def get_embed_index_dict(file):
    embed_index_dict = {}
    for line in open(file):
        lineL = line.strip().split(' ')
        motif = lineL[-1]
        index = lineL[1]
        embed_index_dict[motif] = index
    return embed_index_dict
file="/public/work/zehui/index/data/embed.txt"
embed_motif_dic = get_embed_index_dict(file)

eventalignfile = sys.argv[1]
output_dir = sys.argv[2]
sample_name = sys.argv[3]

feature_dic = {}
motif_dic = {}
with open(eventalignfile) as f:
    next(f)
    for line in f:
        lineL = line.strip().split("\t")
        site = int(lineL[1])
        extra_info = f'{lineL[0]}|{lineL[3]}'
        feature_dic.setdefault(extra_info,{})
        feature_dic[extra_info].setdefault(site,[])
        samples_list = lineL[15].split(",")
        feature_dic[extra_info][site] += samples_list
        
        chr = lineL[0]
        motif_dic.setdefault(chr,{})
        if site not in motif_dic[chr]:
            motif_dic[chr][site] = lineL[2]

sig_file = open(output_dir + "/" + sample_name +"_signal.txt","w")
seq_file = open(output_dir + "/" + sample_name +"_sequence.txt","w")
ext_file = open(output_dir + "/" + sample_name +"_extra.txt","w")
 
for extra_info in feature_dic:
    chr = extra_info.split("|")[0]
    reads_site_counts = len(feature_dic[extra_info].keys())
    reads_site_max = max(feature_dic[extra_info].keys())
    if reads_site_counts < 200: continue

    for site in feature_dic[extra_info]:
        samples_list = feature_dic[extra_info][site]
        
        motif = motif_dic[chr][site]
        seq_list = []
        ran = len(samples_list)
        seq_list +=  [motif for n in range(ran)]
        
        temp_site = site
        if len(samples_list) < 400:
            samples_list.extend(["0"]*400)
            seq_list.extend(["NNNNN"]*400)
        if len(samples_list) < 400:
            print(f"{extra_info}|{site}")
            continue
        out_samples_list = samples_list[:400]
        out_seq_list = seq_list[:400]
        #print(f"{extra_info}|{site}\t{len(out_samples_list)}")
        kmer_slice_array = [embed_motif_dic[i] for i in out_seq_list]
        seq_out = "\t".join(kmer_slice_array)
        
        signal_out = "\t".join(out_samples_list)
        
        can_motif = out_seq_list[0]
        extra_out = f"{extra_info}|{can_motif}|{site}"
        
        print(extra_out,file=ext_file)
        print(f"{extra_out}\t{seq_out}",file=seq_file)
        print(f"{extra_out}\t{signal_out}",file=sig_file)
        
sig_file.close()
seq_file.close()
ext_file.close()

print(f"{nowtime()} | {sample_name} all finish")     
        
