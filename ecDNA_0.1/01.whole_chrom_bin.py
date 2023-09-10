import os
import gzip
from sys import argv
import math
fragments = argv[1]
fragnum=argv[2]
#information=argv[3]
information="/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/ref/gene_s_e.txt"
outfile=argv[3]
#interval=int(argv[5])
interval=100000
myhash2={}
myhash1={}
myhash={}
out=open(outfile,"wt")
with open (information,"r") as fh:
	line=fh.readline().strip()
	while(line):
		myhash1[line.split("\t")[0]]=line.split("\t")[1]+":"+line.split("\t")[2]
		line=fh.readline().strip()
fh.close()
with open (fragnum,"r") as fh:
	line=fh.readline().strip()
	while(line):
		myhash[line.split(" ")[1]]=line.split(" ")[0]
		line=fh.readline().strip()
fh.close()


with gzip.open (fragments,"r") as fh:
	line=fh.readline().strip().decode()
	while(line):
		chrom= line.split("\t")[0]
		if chrom in myhash1:
			start=int(myhash1[chrom].split(":")[0])-interval
			if start < 0 :
				start = 0
			end=int(myhash1[chrom].split(":")[1])
			mylist = range(start,end,interval)
			if int(line.split("\t")[1]) >=start and int(line.split("\t")[2]) <= end:
				for i in range(len(mylist)-1):
					s=mylist[i]
					e=mylist[i+1]
					if int(line.split("\t")[1]) >=s and int(line.split("\t")[2]) <= e:
						myhash2.setdefault(chrom+":"+str(s)+"_"+str(e),[]).append(line.split("\t")[3])
						break
		line=fh.readline().strip().decode()
fh.close()

for binn in myhash2:
	uniqlist=list(set(myhash2[binn]))
	for CB in uniqlist:
		num=myhash2[binn].count(CB)
		prelog=num/int(myhash[CB])*10000
		#final=math.log(prelog)
		out.write(binn+"\t"+CB+"\t"+str(prelog)+"\n")


