import sys
import numpy as np

def process_txt_files(sample):
    data_list = []

    # 读取三个文本文件
    for txt_filename in [f"part1_{sample}_signal.txt", f"part2_{sample}_signal.txt", f"part3_{sample}_signal.txt"]:
        print(f"load {txt_filename} start")
        line_number = 0
        with open(txt_filename, 'r') as f:
            for line in f:
                row = [float(x) for x in line.strip().split()]
                data_list.append(row)
                line_number += 1
                if line_number % 100000 == 0:
                    print(f"load {txt_filename} {line_number} lines")

    counts = len(data_list)
    shape = (counts, 400)  # 根据实际数据量调整
    dtype = "float32"  # 根据实际数据类型调整
    print(f"{sample} : {counts} site")
    data_array = np.array(data_list, dtype=dtype)
    print(f"{sample} : write npy")
    memmap_filename = f"{sample}_signal.npy"  # 定义内存映射文件名
    memmap = np.memmap(memmap_filename, dtype=dtype, mode='w+', shape=shape)

    memmap[:] = data_array[:]
    memmap.flush()

if __name__ == '__main__':
    sample = sys.argv[1]
    process_txt_files(sample)
