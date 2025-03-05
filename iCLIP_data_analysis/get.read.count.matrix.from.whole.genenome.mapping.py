#!usr/bin/python
#coding:utf-8
import sys,re
if len(sys.argv)<2:
    print("""usage: python3 get.read.count.matrix.from.whole.genenome.mapping.py normalized.bw gene.ID.table gtf full.matrix
    this script count the signal for a gene based on the gtf annotation:
    #  -------|----------------|--------------------|---------------------------|----------------
    # up50   TSS     UTR5  50+cds_start-1   CDS     50+cds_end      UTR3       TES  down 50
        """)
    sys.exit(-1)
else:
    import pyBigWig
    import numpy as np
    bw=pyBigWig.open(sys.argv[1])
    id_file=open(sys.argv[2],'r')
    longest_protein={}
    for line in iter(id_file):
        line_info=re.split("\t",line.strip())
        longest_protein[line_info[1]] = line_info[3:]
    gtf=open(sys.argv[3],'r')
    out=open(sys.argv[4],'w')
    transcript_signals={}
    last_transcript=""
    downstream=[]
    transcript_id=""
    for line in iter(gtf):
        if not line.startswith("#"):
            line_info=re.split("\t",line.strip())
            transcript_id=re.split("\"",line_info[8])[3]
            if transcript_id in longest_protein.keys():
                if line_info[2] == "transcript":
                # first line of a transcript
                        if last_transcript !="":
                        # add the downstream signal to the last record
                            downstream_signal=bw.values(downstream[0],downstream[1],downstream[2])
                            if downstream[3] == "-":
                                downstream_signal=downstream_signal[::-1]
                            transcript_signals[last_transcript]=transcript_signals[last_transcript] + downstream_signal
                            # put the upstream signal to the dict & store downstream region
                            if line_info[6] =="+":
                                # update the downstream list
                                downstream=[line_info[0],int(line_info[4])+1,int(line_info[4])+51,line_info[6]]
                                # TES +1 ----> TES + 51
                                transcript_signals[transcript_id]=bw.values(line_info[0],int(line_info[3])-51,int(line_info[3])-1)
                                # TSS - 51 -----> TSS -1  upstream signal
                            else:
                                up_stream_signal=bw.values(line_info[0],int(line_info[4])+1,int(line_info[4])+51)
                                # end + 1 -----> end + 51
                                downstream=[line_info[0],int(line_info[3])-51,int(line_info[3])-1,line_info[6]]
                                # start - 51 ----> start -1 
                                transcript_signals[transcript_id]=up_stream_signal[::-1]
                        else:
                        # the first record, put the upstream signal to the dict & store downstream region
                            if line_info[6] =="+":
                                downstream=[line_info[0],int(line_info[4])+1,int(line_info[4])+51,line_info[6]]
                                # TES +1 ----> TES + 51
                                transcript_signals[transcript_id]=bw.values(line_info[0],int(line_info[3])-51,int(line_info[3])-1)
#                               transcript_signals[transcript_id]= downstream + bw.values(line_info[0],int(line_info[3])-51,int(line_info[3])-1)
                                # TSS - 51 -----> TSS -1
                            else:
                                up_stream_signal=bw.values(line_info[0],int(line_info[4])+1,int(line_info[4])+51)
                                # end + 1 -----> end + 51
                                downstream=[line_info[0],int(line_info[3])-51,int(line_info[3])-1,line_info[6]]
                                # start - 51 ----> start -1 
                                #transcript_signals[transcript_id]=downstream + up_stream_signal[::-1]
                                transcript_signals[transcript_id]=up_stream_signal[::-1]
                        last_transcript = transcript_id
                if line_info[2] == "exon":
                    signal=bw.values(line_info[0],int(line_info[3])-1,int(line_info[4]))
                    if line_info[6] == "-":
                        signal=signal[::-1]
                    transcript_signals[transcript_id]=transcript_signals[transcript_id]+signal
    # add the downstream signal to the last record
    if transcript_id in transcript_signals.keys():
        downstream_signal=bw.values(downstream[0],downstream[1],downstream[2])
        if downstream[3] == "-":
            downstream_signal=downstream_signal[::-1]
        transcript_signals[transcript_id]=transcript_signals[transcript_id]+downstream_signal

    for transcript_id in transcript_signals:
        protein_info=longest_protein[transcript_id]
        # transcript_id => [gene_length, cds_region]
        gene_length=int(protein_info[0])
        [cds_start,cds_end]=re.split("-",protein_info[1])
        cds_start,cds_end=int(cds_start),int(cds_end)
        if cds_start >1 and gene_length > cds_end:
            transcript_signals_values=np.array(transcript_signals[transcript_id])
            transcript_signals_values=list(np.where(np.isnan(transcript_signals_values),0,transcript_signals_values))
            #  -------|----------------|--------------------|---------------------------|---------------------
            #  up    50    UTR5  50+cds_start-1   CDS     50+cds_end      UTR3   gene_length+50(-50)   down
            
            upstream_signal=" ".join(map(str,transcript_signals_values[0:50]))
            utr5_signal=" ".join(map(str,transcript_signals_values[50:50+cds_start-1]))

            cds_signal=" ".join(map(str,transcript_signals_values[50+cds_start-1 : 50+cds_end-3]))
            utr3_signal=" ".join(map(str,transcript_signals_values[50+cds_end-3 : gene_length+50]))

            downstream_signals=" ".join(map(str,transcript_signals_values[-50:]))
            out.write(transcript_id + "  " + str(len(transcript_signals_values))+"  " + upstream_signal + "  " + utr5_signal + "  " + cds_signal + "  " + utr3_signal +"  "+ downstream_signals +"\n")
