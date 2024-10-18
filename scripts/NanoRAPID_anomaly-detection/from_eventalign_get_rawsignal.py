#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import time

def nowtime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

eventalignfile = sys.argv[1]
outfile = sys.argv[2]
reads_num = 0
pre_read_name = ""
read_dict = {}
with open(eventalignfile,"r") as f:
    # header
    next(f)
    for line in f:
        line = line.strip()
        lineL = line.split("\t")
        read_name = lineL[3]
        if read_name != pre_read_name:
            for (gene_name,site,motif),singal_list in read_dict.items():
                #outlist = [f"{gene_name}|{site}|{motif}|{pre_read_name}",np.mean(singal_list),np.std(singal_list),np.median(singal_list),len(singal_list)]
                outlist = [gene_name,site,motif,pre_read_name,np.mean(singal_list),np.std(singal_list),np.median(singal_list),len(singal_list)]
                print("\t".join([str(n) for n in outlist]),file=open(outfile,"a"))
                #print("----")
            #--
            pre_read_name = read_name
            read_dict = {}
            reads_num += 1
            if reads_num % 1000 == 0:
                print("%s\t%d"%(nowtime(),reads_num))
        gene_name = lineL[0]
        site = int(lineL[1])
        motif = lineL[2]
        read_dict.setdefault((gene_name,site,motif),[])
        singal_list = [float(n) for n in lineL[15].split(",")]
        read_dict[(gene_name,site,motif)] += singal_list



