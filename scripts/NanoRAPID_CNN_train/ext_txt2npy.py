import sys
import numpy as np

def process_txt_files(sample):
    data_list = []

    # 读取三个文本文件
    for txt_filename in [f"part1_{sample}_extra.txt", f"part2_{sample}_extra.txt", f"part3_{sample}_extra.txt"]:
        print(f"load {txt_filename}")
        with open(txt_filename, 'r') as f:
            for line in f:
                data_list.append(line.strip())

    counts = len(data_list)
    shape = (counts,)  # 根据实际数据量调整
    dtype = "<U70"  # 根据实际数据类型调整
    print(f"{sample} : {counts} site")
    data_array = np.array(data_list, dtype=dtype)
    print(f"{sample} : write npy")
    memmap_filename = f"{sample}_extra.npy"  # 定义内存映射文件名
    memmap = np.memmap(memmap_filename, dtype=dtype, mode='w+', shape=shape)

    memmap[:] = data_array[:]
    memmap.flush()

if __name__ == '__main__':
    sample = sys.argv[1]
    process_txt_files(sample)

