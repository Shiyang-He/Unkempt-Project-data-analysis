#!usr/bin/python3
#coding:utf-8

import sys,re
from optparse import OptionParser
def parseOptions():
    desc='''this script is to fish out some line from pool files using the columns specified in bait file.'''
    epilog='''example: python3 fish.column.py bait,tab,1,2 pool,space,2,3 out.xls which use 1st and 2nd column (\\t separated) in the bait file to fish out the common lines of pool file (space separated 2nd and 3rd column)'''
    parser=OptionParser(description=desc,epilog=epilog)
    parser.add_option("-H","--Head",type=str,help="print out head line in pool file",default="no")
    parser.add_option("-r","--reverse",type=str,help="print out lines in pool file but not in bait file",default="no")
    parser.add_option("-o","--out",type=str,help="print the output to file otherwise to STDOUT")
    return parser.parse_args()


if len(sys.argv)<2:
    print("""usage: python3 fish.column.py bait,separator,columns pools,separator,columnss -o out.xls
    separator can be vertical:|, comma:, tab  
    example: python3 fish.column.py bait,tab,1,2 pool,space,2,3 out.xls
    which use 1st and 2nd column (\\t separated) in the bait file to fish
    out the common lines of pool file (space separated 2nd and 3rd column)
    """)
    sys.exit(-1)
else:
    (options,args)=parseOptions()
    header=options.Head
    
    bait_info=re.split(",",sys.argv[1])
    bait_file=bait_info[0];del(bait_info[0])
    bait_file_fh=open(bait_file,'r')
    bait_sep=bait_info[0];del(bait_info[0])
    bait_col=bait_info
    bait_dict={}
    pool_info=re.split(",",sys.argv[2])
    pool_file=pool_info[0];del(pool_info[0])
    pool_file_fh=open(pool_file,'r')
    pool_sep=pool_info[0];del(pool_info[0])
    pool_col=pool_info
    if options.out:
        out=open(options.out,'w')
    if bait_sep=="vertical":bait_sep="\|"
    if bait_sep=="comma":bait_sep=","
    if bait_sep=="tab":bait_sep="\t"
    if bait_sep=="space":bait_sep=" "
    if pool_sep=="vertical":pool_sep="\|"
    if pool_sep=="comma":pool_sep=","
    if pool_sep=="tab":pool_sep="\t"
    if pool_sep=="space":pool_sep=" "
    for line in iter(bait_file_fh):
        store_columns=[]
        line_info=re.split(bait_sep,line.strip())
        for c in bait_col:
            c=int(c)-1
            store_columns.append(line_info[c])
        bait_dict[tuple(store_columns)]=1
    if options.reverse=="no":
        i=1
        for line in iter(pool_file_fh):
            if i==1 and header =="yes":
                if options.out:
                    out.write(line)
                else:
                    print(line.strip())
                i+=1
                continue
            else:   
                i+=1
                line_info=re.split(pool_sep,line.strip())
                search_columns=[]
                for c in pool_col:
                    c=int(c)-1
                    search_columns.append(line_info[c])
                if tuple(search_columns) in bait_dict:
                    if options.out:
                        out.write(line)
                    else:
                        print(line.strip())
    else:
        i=1
        for line in iter(pool_file_fh):
            if i==1 and header =="yes":
                if options.out:
                    out.write(line)
                else:
                    print(line.strip())
                i+=1
                continue
            else:   
                i+=1
                line_info=re.split(pool_sep,line.strip())
                search_columns=[]
                for c in pool_col:
                    c=int(c)-1
                    search_columns.append(line_info[c])
                if tuple(search_columns) not in bait_dict:
                    if options.out:
                        out.write(line)
                    else:
                        print(line.strip())
