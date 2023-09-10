import os
from sys import argv
import re
cellline=argv[1]
#"/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/fragnumNormal_CL81_HT29.fragnum.txt"
translate="/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/00.result/all_barcodeTranslate.tsv"
junctionsite="/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/result/08.junctionsite/all_final.txt"
outfile=argv[2]
out=open(outfile,"wt")
myhash={}
myhash2={}

with open(cellline,"r") as fh1:
	line=fh1.readline().strip()
	while(line):
		myhash[line.split(" ")[1]]=0
		line=fh1.readline().strip()

with open(translate,"r") as fh2:
	line=fh2.readline().strip()
	while(line):
		if line.split("\t")[1] in myhash:
			myhash2[line.split("\t")[0]]=line.split("\t")[1]
		line=fh2.readline().strip()

with open (junctionsite,"r") as fh3:
	line=fh3.readline().strip()
	while(line):
		seqid=line.split("\t")[0].split("|||CB:Z:")[1]
		if seqid in myhash2:
			out.write(line + "\t"+ myhash2[seqid] + "\n")
		line=fh3.readline().strip()
