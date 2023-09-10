import os
import gzip
from sys import argv
import math
fragments = argv[1]
fragnum=argv[2]
outfile=argv[3]
interval=100000
myhash2={}
myhash={}
out=open(outfile,"wt")

with open (fragnum,"r") as fh:
	line=fh.readline().strip()
	while(line):
		myhash[line.split("\t")[1]]=line.split("\t")[0]
		line=fh.readline().strip()
fh.close()


with open (fragments,"r") as fh:
	line=fh.readline().strip()
	while(line):
		chrom= line.split("\t")[0]
		s=int (int(line.split("\t")[1])/100000)*100000
		e=(int (int(line.split("\t")[1])/100000)+1)*100000
		myhash2.setdefault(chrom+":"+str(s)+"_"+str(e),[]).append(line.split("\t")[3])
		line=fh.readline().strip()
fh.close()

for binn in myhash2:
	uniqlist=list(set(myhash2[binn]))
	for CB in uniqlist:
		num=myhash2[binn].count(CB)
		prelog=num/int(myhash[CB])*10000
		#final=math.log(prelog)
		out.write(binn+"\t"+CB+"\t"+str(prelog)+"\n")


