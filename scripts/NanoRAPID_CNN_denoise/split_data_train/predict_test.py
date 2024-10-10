import sys
import os
import torch
import torch.nn as nn
from torch.nn import functional as F
import random
import torch.multiprocessing as mp
from torch.nn.parallel import DistributedDataParallel as DDP
import torch.distributed as dist
from torch.autograd import Function
import numpy as np
from torch.utils.data import Dataset, Subset,DataLoader
import scipy.stats as ss
import collections
from torch.cuda.amp import autocast, GradScaler
from datetime import datetime
import argparse

class MilDataset(Dataset):
        def __init__(self, seq,sig,extra):
                self.seq = seq
                self.sig = sig
                self.extra = extra
        def __len__(self):
                return len(self.seq)
        def __getitem__(self, index):
                return torch.from_numpy(self.seq[index].copy()), torch.from_numpy(self.sig[index].copy()), self.extra[index]
                #return torch.from_numpy(self.seq[index]),torch.from_numpy(self.sig[index]),self.extra[index]
        def n_features(self):
                return 2

class Model(nn.Module):
        def __init__(self):
                super(Model,self).__init__()
                # add layers
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
                #module 1
                sigs_x = self.conv1(sigs)

                #module 2
                seqs_x = self.conv2(seqs)

                #merge
                z = torch.cat((self.dropout(sigs_x), self.dropout(seqs_x)), 1)
                z = self.conv3(z)
                z = torch.flatten(z, start_dim=1)
                z = self.dropout(z)
                z = self.fc(z)

                #combine
                probs = self.softmax(z)[:,1]

                return probs

def inference(test_dl,model_pare,out_file,device):
        #load model
        model = Model().to(device)
        checkpoint = torch.load(model_pare,map_location=device)
        model.load_state_dict({k.replace('module.',''):v for k,v in checkpoint['state_dict'].items()})

        #open file to save predict result
        f1 = open(out_file,"w")

        #onehot coding
        base_to_onehot = torch.eye(5)
        base_to_onehot = base_to_onehot.to(device)
        model.eval()
        with torch.no_grad():
                for i,datas in enumerate(test_dl):
                        seq, sig, extra = datas
                        seq = seq.long().to(device)
                        seq = torch.clamp(seq, max=4)
                        seq = base_to_onehot[seq.flatten()].view(-1, 400, 25).permute(0, 2, 1)
                        sig = sig.reshape(-1,1,400).to(device)
                        prob  = model(sig,seq)
                        test_prob = prob.tolist()
                        test_extra = list(extra)

                        for site,score in zip(test_extra,test_prob):
                                #print(site,score)
                                print(site,score,sep="\t",file=f1,flush=True)
        f1.close()

def main():
        #command
        print("-------------------------------------------")
        print("python"," ".join(sys.argv),sep="\t")

        #parameters setting
        parser = argparse.ArgumentParser(description='Prediction using SingleMod model within specific motif')
        parser.add_argument('-d','--npy_dir',required = True, help="the directory to npy files")
        parser.add_argument('-m','--model',required = True, help="the used model")
        parser.add_argument('-g','--gpu',type=int,default=0,required = False, help="cuda index, default is 0")
        parser.add_argument('-o','--out_file',required = True, help="the file to output prediction result")
        args = parser.parse_args(sys.argv[1:])

        motif_extra = np.memmap(args.npy_dir + "_extra.npy", mode='r', dtype="<U80")
        print(f'predicting start')
        #load data
        #total chunks number
        lengths = motif_extra.shape[0]
        motif_seq = np.memmap(args.npy_dir + "_sequence.npy", mode='r', shape=(lengths,400,5), dtype="int8")
        motif_sig = np.memmap(args.npy_dir + "_signal.npy", mode='r', shape=(lengths,400), dtype="float32")

        dataset = MilDataset(motif_seq, motif_sig,motif_extra)
        test_dl = DataLoader(dataset,shuffle = False,batch_size = 30000)

        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(device,torch.cuda.is_available())
        inference(test_dl,args.model,args.out_file,device)
        print(f'predicting finish')

if __name__ == '__main__':
        main()
