import sys,os
ctrl_d = {}
test_d = {}
filenameL = [f"{sys.argv[1]}/{n}" for n in os.listdir(sys.argv[1]) if n.find(".filter.txt")>0]
for file in filenameL:
    for line in open(file):
        lineL = line.split(' ')
        name_info_l = lineL[0].split("|")
        try:
            motif = name_info_l[-3]
        except:
            print(file,line.strip())
            sys.exit()
        index = motif
        test_d.setdefault(index,{})
        ctrl_d.setdefault(index,[])
        score = float(lineL[1])
        site = lineL[0][:-2]
        if name_info_l[-1] == "1":
            test_d[index][site] = score
        else:
            ctrl_d[index].append(score)
#
for index in ctrl_d:
    ctrl_list = sorted(ctrl_d[index])
    cutline = int(len(ctrl_list) * 0.95)
    cutline = min([cutline,len(ctrl_list)-5])
    try:
        cutscore = ctrl_list[cutline]
    except:
        print(ctrl_list,cutline)
        sys.exit()
    for site in test_d[index]:
        score = test_d[index][site]
        if score > cutscore:
            print(site,score,cutscore)


