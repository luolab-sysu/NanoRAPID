import sys
import numpy as np

def get_embed_dict(file):
    """从文件中读取嵌入字典"""
    embed_dict = {}
    with open(file, 'r') as f:
        for line in f:
            lineL = line.strip().split()
            embed_dict[lineL[1]] = lineL[2] 
    return embed_dict

def process_txt_files(sample):
    embed_dict = get_embed_dict("/public/work/zehui/work/public_data_result/NBT_ONT_result/h9_mod_rep1/code/rawsignal/embed.txt")
    data_list = []

    # 读取三个文本文件
    for txt_filename in [f"part1_{sample}_sequence.txt", f"part2_{sample}_sequence.txt", f"part3_{sample}_sequence.txt"]:
        print(f"load {txt_filename} start")
        line_number = 0
        with open(txt_filename, 'r') as f:
            for line in f:
                row = []
                for index in line.strip().split():
                    k = embed_dict[index]
                    kmer = [int(x) for x in k]
                    row.append(kmer)
                data_list.append(row)
                line_number += 1
                if line_number % 100000 == 0:
                    print(f"load {txt_filename} {line_number} lines")

    counts = len(data_list)
    shape = (counts, 400, 5) # 根据实际数据量调整
    dtype="int8"  # 根据实际数据类型调整
    print(f"{sample} : {counts} site")
    data_array = np.array(data_list, dtype=dtype)
    print(f"{sample} : write npy")
    memmap_filename = f"{sample}_sequence.npy"  # 定义内存映射文件名
    memmap = np.memmap(memmap_filename, dtype=dtype, mode='w+', shape=shape)

    memmap[:] = data_array[:]
    memmap.flush()

if __name__ == '__main__':
    sample = sys.argv[1]
    process_txt_files(sample)

