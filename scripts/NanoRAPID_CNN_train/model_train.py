import sys,os,random
import collections
import numpy as np
import scipy.stats as ss
from datetime import datetime

import torch
import torch.nn as nn
import torch.multiprocessing as mp
import torch.distributed as dist
from torch.nn import functional as F
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.autograd import Function
from torch.utils.data import Dataset, Subset, DataLoader, random_split
from torch.cuda.amp import autocast, GradScaler


def get_time():
    now = datetime.now()
    return now.strftime("%Y-%m-%d %H:%M:%S")

class MilDataset(Dataset):
    def __init__(self, seq, sig, extra):
        self.seq = seq
        self.sig = sig
        self.extra = extra
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, index):
        extra_info_list = self.extra[index].split("|")
        label = torch.tensor(int(extra_info_list[-1]), dtype=torch.long)
        return torch.from_numpy(self.seq[index].copy()), torch.from_numpy(self.sig[index].copy()), label, self.extra[index]
    def n_features(self):
        return 2

class Model(nn.Module):
    def __init__(self):
        super(Model,self).__init__()
        #add layers
        #signal
        self.conv1 = nn.Sequential(
            nn.Conv1d(1, out_channels=16, kernel_size=3,padding = 1),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(in_channels=16, out_channels=32, kernel_size=3,padding = 1),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(in_channels=32, out_channels=64, kernel_size=3,padding = 1),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(in_channels=64, out_channels=128, kernel_size=3,padding = 1),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Conv1d(in_channels=128, out_channels=256, kernel_size=3,padding = 1),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.MaxPool1d(2)
        )

        #sequence
        self.conv2 = nn.Sequential(
            nn.Conv1d(25, out_channels=16, kernel_size=3,padding = 1),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(in_channels=16, out_channels=32, kernel_size=3,padding = 1),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(in_channels=32, out_channels=64, kernel_size=3,padding = 1),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(in_channels=64, out_channels=128, kernel_size=3,padding = 1),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Conv1d(in_channels=128, out_channels=256, kernel_size=3,padding = 1),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.MaxPool1d(2)
        )

        #dropout
        self.dropout = nn.Dropout(p=0.5)

        #merge
        self.conv3 = nn.Sequential(
            nn.Conv1d(512, out_channels=64, kernel_size=3,padding = 1),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(in_channels=64, out_channels=64, kernel_size=3),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Conv1d(in_channels=64, out_channels=64, kernel_size=5),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Conv1d(in_channels=64, out_channels=64, kernel_size=7),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Conv1d(in_channels=64, out_channels=64, kernel_size=11),
            nn.BatchNorm1d(64),
            nn.ReLU()
        )

        self.fc = nn.Linear(64 * 3, 2)

        #activate
        self.softmax = nn.Softmax(dim=1)

    def forward(self,sigs,seqs):
        #forward
        #print(sigs.size(),seqs.size())

        #module 1
        sigs_x = self.conv1(sigs)
        #module 2
        seqs_x = self.conv2(seqs)

        #print(sigs_x.size(),seqs_x.size())
        #merge
        z = torch.cat((self.dropout(sigs_x), self.dropout(seqs_x)), 1)
        z = self.conv3(z)
        z = torch.flatten(z, start_dim=1)
        z = self.dropout(z)
        z = self.fc(z)

        #combine
        probs = self.softmax(z)[:,1]

        return z,probs


def train(train_dl,test_dl,epochs,out_dir,device):
        #define model, loss and optimizer
        model = Model().to(device)
        criterion = nn.CrossEntropyLoss()

        lr = 1e-3
        weight_decay = 1e-5
        optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
        torch.set_num_threads(2)

        # 定义 StepLR 学习率调度器
        step_size = 50  # 每 50 个 epoch 调整一次学习率
        gamma = 0.65  # 衰减因子
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=step_size, gamma=gamma)

        time = {}

        #onehot coding
        base_to_onehot = torch.eye(5)
        base_to_onehot = base_to_onehot.to(device)

        for epoch in range(epochs):
            model.train()

            train_loss, test_loss = 0, 0
            for i,datas in enumerate(train_dl):
                seq,signal,label,extra = datas
                seq = seq.long().to(device)
                seq = torch.clamp(seq, max=4)
                seq = base_to_onehot[seq.flatten()].view(-1, 400, 25).permute(0, 2, 1) #seq2onehot
                sig = signal.reshape(-1,1,400).to(device)
                label = label.to(device)
                out,prob  = model(sig,seq)
                #print(prob,label)
                loss = criterion(out,label)

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                # 保存每个batch的loss
                train_loss += loss
            mean_train_loss = train_loss/i

            scheduler.step()

            model.eval()
            with torch.no_grad():
                for i,datas in enumerate(test_dl):
                    seq,signal,label,extra = datas
                    seq = seq.long().to(device)
                    seq = torch.clamp(seq, max=4)
                    seq = base_to_onehot[seq.flatten()].view(-1, 400, 25).permute(0, 2, 1) #seq2onehot
                    sig = signal.reshape(-1,1,400).to(device)
                    label = label.to(device)
                    out,prob  = model(sig,seq)
                    loss = criterion(out,label)
                    # 保存每个batch的loss
                    test_loss += loss
            mean_test_loss = test_loss/i

            epoch_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            time[epoch+1] = f'{epoch+1}\t{epoch_time}\t{mean_train_loss}\t{mean_test_loss}'
            torch.save({'state_dict': model.state_dict()}, out_dir + f'/{sys.argv[2]}_model_{str(epoch+1)}.pth.tar', _use_new_zipfile_serialization=False)
            #if epoch == 0 or (epoch+1) % 5 == 0:
            if (epoch+1) > 0:
                for index in time:
                    print(time[index])
                time = {}

def main():
    #command
    print("-------------------------------------------")
    print(f"{get_time()} | python"," ".join(sys.argv),sep=" ")
    wd = sys.argv[1]
    motif = sys.argv[2]
    npy_dir = f'{wd}/{motif}'
    extra_file = npy_dir + "_extra.npy"
    out_dir = sys.argv[3]
    if not os.path.exists(extra_file):
        print(f"{get_time()} | npy_dir:{npy_dir} doesnt exit")
        return 
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    #load data
    motif_extra = np.memmap(extra_file, mode='r', dtype="<U80")
    print(motif_extra[:5])
    #total chunks number
    lengths = motif_extra.shape[0]
    motif_seq = np.memmap( npy_dir + "_sequence.npy", mode='r', shape=(lengths,400,5), dtype="int8")
    motif_sig = np.memmap( npy_dir + "_signal.npy", mode='r', shape=(lengths,400), dtype="float32")
    print(f'{get_time()} | training within {motif} motif start {lengths}')
    #print(motif_seq[:5])

    if lengths > 1000000 :
        trian_size = 256
    elif lengths > 500000:
        trian_size = 128
    elif lengths > 50000:
        trian_size = 64
    elif lengths > 5000:
        trian_size = 32
    else:
        trian_size = 8

    dataset = MilDataset(motif_seq, motif_sig, motif_extra)

    # 创建随机数生成器对象并设置种子
    rng = torch.Generator().manual_seed(72)
    # 划分数据集
    train_size, test_size = 0.9, 0.1
    train_dataset, test_dataset = random_split(dataset, [train_size, test_size], generator=rng)

    train_dl = DataLoader(dataset,shuffle = True,batch_size = trian_size)
    test_dl = DataLoader(dataset,shuffle = False,batch_size = 64)

    if len(sys.argv) > 4:
        gpu = sys.argv[4]
    else:
        gpu = "3"
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    if torch.cuda.is_available() :
        print(f"{get_time()} | cuda {gpu} is available {torch.cuda.is_available()}; batchs: {trian_size}")
    else:
        print(f"{get_time()} | GPU is not available")
        return
    epochs = 50
    train(train_dl,test_dl,epochs,out_dir,device)
    print(f'{get_time()} | trianing within {motif} motif finish')

if __name__ == '__main__':
    main()
