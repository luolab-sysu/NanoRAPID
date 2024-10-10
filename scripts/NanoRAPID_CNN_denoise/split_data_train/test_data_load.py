import sys,os,random
import collections
import numpy as np
import scipy.stats as ss
from datetime import datetime

import torch
import torch.nn as nn
from torch.nn import functional as F
import torch.multiprocessing as mp
from torch.nn.parallel import DistributedDataParallel as DDP
import torch.distributed as dist
from torch.autograd import Function
from torch.utils.data import Dataset, Subset, DataLoader, random_split
from torch.cuda.amp import autocast, GradScaler


def get_time():
    now = datetime.now()
    return now.strftime("%Y-%m-%d %H:%M:%S")


class MilDataset(Dataset):
    def __init__(self, seq, sig, extra_info):
        self.seq = seq
        self.sig = sig
        self.extra_info = extra_info
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, index):
        extra_info_list = self.extra_info[index].split("|")
        label = torch.tensor(int(extra_info_list[-1]), dtype=torch.long)
        extra = "|".join(extra_info_list[:-1])
        return torch.from_numpy(self.seq[index].copy()), torch.from_numpy(self.sig[index].copy()), label, extra
    def n_features(self):
        return 2


def main():
    #command
    print("-------------------------------------------")
    print(f"{get_time()} | python"," ".join(sys.argv),sep="\t")

    npy_dir = sys.argv[1]

    motif_extra = np.memmap(npy_dir + "_extra.npy", mode='r', dtype="<U80")
    #load data
    #total chunks number
    lengths = motif_extra.shape[0]
    motif_seq = np.memmap( npy_dir + "_sequence.npy", mode='r', shape=(lengths,400,5), dtype="int8")
    motif_sig = np.memmap( npy_dir + "_signal.npy", mode='r', shape=(lengths,400), dtype="float32")
    print(f'{get_time()} | training within {npy_dir} start {lengths}')


    dataset = MilDataset(motif_seq, motif_sig, motif_extra)
    test_dl = DataLoader(dataset,shuffle = False,batch_size = 64)
    for i,datas in enumerate(test_dl):
        seq,signal,label,extra = datas
        print(seq)
        print(signal)
        print(label)
        print(extra)
        break


if __name__ == '__main__':
    main()

