import os
from sys import argv
mypos=argv[1]
g38 = "/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/gene.pos_4.txt"
h19 = "/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/ref/gene.pos_hg19.txt"
mm10 = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/ref/mm10/gene_pos.txt"
outfile=argv[2]
myhash={}
out=open(outfile,"wt")
with open (mypos,"r") as fh:
	line= fh.readline().strip()
	while (line):
		myhash.setdefault(line.split("\t")[0],[]).append(line.split("\t")[1]+"_"+line.split("\t")[2])
		line= fh.readline().strip()
fh.close()
with open (mm10,"r") as fh:
	line= fh.readline().strip()
	while (line):
		if "chr"+line.split("\t")[1] in myhash:
			for i in myhash["chr"+line.split("\t")[1]]:
				if int(line.split("\t")[2]) > int(i.split("_")[0]) and int(line.split("\t")[3]) < int(i.split("_")[1]):
					out.write("chr"+line.split("\t")[1]+":"+i+"\t"+line.split("\t")[0]+"\n")
					break
		line= fh.readline().strip()

