#!usr/bin/python3
#coding:utf-8
import sys,re,pysam
if len(sys.argv)<2:
    print("""usage: python3 merge.tailfindr.and.nannopolish.results.py transcript.bam tailfindr.results nanopolish.results combined.results""")
    sys.exit(-1)
else:
    samfile = pysam.AlignmentFile(sys.argv[1],'rb')
    tailfindr=open(sys.argv[2],'r')
    nanopolish=open(sys.argv[3],'r')
    out=open(sys.argv[4],'w')
    tailfindr_result,nanopolish_result={},{}
    for line in iter(tailfindr):
        line_info=re.split(",",line.strip())
        if line_info[4] =="": line_info[4] = "-"
        tailfindr_result[line_info[0]]=line_info[4]
    for line in iter(nanopolish):
        line_info=re.split("\t",line.strip())
        if line_info[9] =="PASS":
            nanopolish_result[line_info[0]]=line_info[8]
        else:
            nanopolish_result[line_info[0]]="-"

    for read in samfile.fetch(until_eof=True):
        output=[read.qname,"-","-","-"]
        if read.is_mapped:
            if read.qname in nanopolish_result.keys():
                output[2]=nanopolish_result[read.qname]
            if read.qname in tailfindr_result.keys():
                output[1]=tailfindr_result[read.qname]
            output[3]=read.reference_name
        out.write("\t".join(output)+"\n")
