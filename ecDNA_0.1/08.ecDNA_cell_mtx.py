import os
import gzip
from sys import argv
import math
fragments = argv[1]
fragnum=argv[2]
#information=argv[3]
information=argv[3]
outfile=argv[4]
#interval=int(argv[5])
interval=100000
myhash2={}
myhash1={}
myhash={}
out=open(outfile,"wt")
with open (information,"r") as fh:
	fh.readline().strip()
	line=fh.readline().strip()
	while(line):
		myhash1.setdefault(line.split("\t")[0].split(":")[0],[]).append(line.split("\t")[0].split(":")[1])
		line=fh.readline().strip()
fh.close()
#print(myhash1["chr1"])
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
		if chrom in myhash1:
			s= int(line.split("\t")[1])
			e= int(line.split("\t")[2])
			for i in myhash1[chrom]:
				if s >=int(i.split("_")[0]) and e <= int(i.split("_")[1]):
					myhash2.setdefault(chrom+":"+i,[]).append(line.split("\t")[3])
					break
		line=fh.readline().strip()
fh.close()

for binn in myhash2:
	uniqlist=list(set(myhash2[binn]))
	for CB in uniqlist:
		num=myhash2[binn].count(CB)
		prelog=num/int(myhash[CB])*10000
		#final=math.log(prelog)
		out.write(binn+"\t"+CB+"\t"+str(prelog)+"\n")


