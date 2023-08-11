#!usr/bin/python3
#coding:utf-8

import sys,re
if len(sys.argv)<2:
    print("usage: python3 combine.reads.count.py out.file input.files.....")
    sys.exit(-1)
else:
    file_dict={}
    out=open(sys.argv[1],'w')
    keys={}
    files=[]
    for i in sys.argv[2:len(sys.argv)]:
        fh=open(i,'r')
        i=i.replace(".all.rRNA.sort.bam.umapped.bam.fq.gz.map.to.lncRNA.bam.unmapped.bam.fq.gz.map.to.pc.sort.bam.dedup.bam.read.count","")
        files.append(i)
        for line in iter(fh):
            if line.startswith("ENSG"):
                line_info=re.split("\t",line.strip())
                keys[line_info[0]]=1
                file_dict[(i,line_info[0])]=line_info[1]
    out.write("genes")
    for f in files:
        out.write("\t"+f)
    out.write("\n")
    for k in keys.keys():
        out.write(k)
        for f in files:
            if (f,k) in file_dict.keys():
                out.write("\t"+str(file_dict[(f,k)]))
            else:
                out.write("\t0")
        out.write("\n")
