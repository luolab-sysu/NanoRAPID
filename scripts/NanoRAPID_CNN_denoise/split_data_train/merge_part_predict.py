import os,sys
import numpy as np

filenameL = [f'{sys.argv[1]}/{file}' for file in os.listdir(sys.argv[1]) if file.find("rep")>0]

adic = {}
for filename in filenameL:
    for line in open(filename):
        lineL = line.strip().split("\t")
        name = lineL[0]
        score = float(lineL[1])
        adic.setdefault(name,[])
        adic[name].append(score)

for name in adic:
    mean = np.mean(adic[name])
    outlist = [str(n) for n in adic[name]]
    if len(outlist) < 3:
        print("error",name,outlist)
        continue
    print(name,mean,' '.join(outlist))
