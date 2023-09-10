import gzip
import os
from sys import argv
infile1="/jdfssz1/ST_SUPERCELLS/P21Z10200N0134/Project/22.ecDNA/02.ArchR/result/cell_matedata.txt"
infile2=argv[1]
outfile=argv[2]
myhash={}
out=open(outfile,"wt")
with open (infile1,"r") as fh:
	fh.readline().strip()
	line = fh.readline().strip()
	while(line):
		if line.split("\t")[16] == "C2" or line.split("\t")[16] == "C3":
			cell=line.split("\t")[0].split(".")[1]
			myhash[cell]=1
		line = fh.readline().strip()
with gzip.open (infile2,"r") as fh:
	line = fh.readline().strip().decode()
	while(line):
		if line.split("\t")[3] in myhash:
			out.write(line+"\n")
		line = fh.readline().strip().decode()
