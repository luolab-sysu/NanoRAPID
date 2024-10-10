import sys
import re
import numpy as np
import time

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

def read_process(database,read_name,strand,idx):
    sites, kmers, signals = database
    cans = set([sites[i] for i, kmer in enumerate(kmers) if motif.match(kmer)])
    idx = idx
    for can in cans:
        first_index = sites.index(can)
        last_index = len(sites) - 1 - sites[::-1].index(can)
        middle_index = (first_index + last_index) // 2

        # +200 -200
        start_index = max(0, middle_index - 199)
        end_index = min(len(kmers), middle_index + 201)

        can_motif = kmers[middle_index]
        kmer_slice = kmers[start_index:end_index]
        signal_slice = signals[start_index:end_index]

        #save to txt
        if len(signal_slice) == 400:
            
            kmer_slice_array = [motif_dic[i] for i in kmer_slice]
            seq_out = " ".join(kmer_slice_array)
            signal_out = "\t".join(signal_slice)
            extra_out = can+"|"+read_name+"|"+can_motif
            
            print(extra_out,file=ext_file)
            print(seq_out,file=seq_file)
            print(signal_out,file=sig_file)
            #print(can,read_name,can_motif,file=extra_file,sep="\t")
            idx = idx + 1
    return idx


bedfile = sys.argv[1]
eventalignfile = sys.argv[2]
output_dir = sys.argv[3]
sample_name = sys.argv[4]

file="/public/work/zehui/work/public_data_result/NBT_ONT_result/h9_mod_rep1/code/rawsignal/embed.txt"
motif_dic = get_embed_dict()
#base to onehot
encoding_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

strands = {}
with open(bedfile,"r") as f:
    for line in f:
        info = line.strip().split("\t")
        strands[info[3]] = info[5]

#candidate motif
#motif = re.compile(r'.[ATG]A[TC].')
motif = re.compile(r'[ATCG][ATCG][ATCG][ATCG][ATCG]')


sig_file = open(output_dir + "/" + sample_name +"_signal.txt","w")
seq_file = open(output_dir + "/" + sample_name +"_sequence.txt","w")
ext_file = open(output_dir + "/" + sample_name +"_extra.txt","w")
#extra_file = open(output_dir + "/" + sample_name +"_extra_info.txt","w")

with open(eventalignfile,"r") as f:
    
    # header
    next(f)
    #first line
    line = f.readline().strip()
    info = line.strip().split("\t")
    database = [[],[],[]]
    read_name = info[3]
    strand = strands.get(read_name)
    if strand == "+" and info[9] != "NNNNN":
        site = info[0]+"|"+info[1]+"|+"
        ran = int(info[14]) - int(info[13])
        database[0][:0] = [site] * ran
        database[1][:0] = [info[9]] * ran
        database[2][:0] = list(map(float, info[15].split(",")))
    elif strand == "-" and info[9] != "NNNNN":
        site = info[0]+"|"+info[1]+"|-"
        ran = int(info[14]) - int(info[13])
        database[0].extend([site] * ran)
        database[1].extend([info[9]] * ran)
        database[2].extend(map(float, info[15].split(",")))
    else:
        pass

    n = 2
    idx = 0
    for line in f:
        info = line.strip().split("\t")
        if info[3] == read_name:
            strand = strands.get(info[3])
            if strand == "+" and info[9] != "NNNNN":
                site = info[0]+"|"+info[1]+"|+"
                ran = int(info[14]) - int(info[13])
                database[0][:0] = [site] * ran
                database[1][:0] = [info[9]] * ran
                database[2][:0] = list(map(float, info[15].split(",")))
            elif strand == "-" and info[9] != "NNNNN":
                site = info[0]+"|"+info[1]+"|-"
                ran = int(info[14]) - int(info[13])
                database[0].extend([site] * ran)
                database[1].extend([info[9]] * ran)
                database[2].extend(map(float, info[15].split(",")))
            else:
                pass
            
            n = n+1
            if n % 500000 == 0:
                print(f'{nowtime()} | {sample_name} {n} line finish')
        else:
            #previous read process
            strand = strands.get(read_name)
            idx = read_process(database,read_name,strand,idx)

            database = [[],[],[]]
            read_name = info[3]
            strand = strands.get(read_name)
            if strand == "+" and info[9] != "NNNNN":
                site = info[0]+"|"+info[1]+"|+"
                ran = int(info[14]) - int(info[13])
                database[0][:0] = [site] * ran
                database[1][:0] = [info[9]] * ran
                database[2][:0] = list(map(float, info[15].split(",")))
            elif strand == "-" and info[9] != "NNNNN":
                site = info[0]+"|"+info[1]+"|-"
                ran = int(info[14]) - int(info[13])
                database[0].extend([site] * ran)
                database[1].extend([info[9]] * ran)
                database[2].extend(map(float, info[15].split(",")))
            else:
                pass
            
            n = n+1
            if n % 500000 == 0:
                print(f'{nowtime()} | {sample_name} {n} line finish')
    
    #last read process
    strand = strands.get(read_name)
    idx = read_process(database,read_name,strand,idx)

sig_file.close()
seq_file.close()
ext_file.close()
#extra_file.close()

print(f"{nowtime()} | {sample_name} all finish")

