import os
from sys import argv
inputfile = argv[1]
outfile=argv[2]
myhash={}
out=open(outfile,"wt")
with open (inputfile,"r") as fh:
	line= fh.readline().strip()
	while(line):
		chrr =line.split("\t")[0].split(":")[0]
		start = line.split("\t")[0].split(":")[1].split("_")[0]
		end = line.split("\t")[0].split(":")[1].split("_")[1]
		out.write(chrr+"\t"+start+"\t"+end+"\n")
		line= fh.readline().strip()

