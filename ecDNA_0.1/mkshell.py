import os
from sys import argv
ID=argv[1]
outfile=argv[2]
out=open(outfile,"wt")
out.write(
"#less ../data/"+ID+".fragments.tsv.gz|awk '{print $4}'|sort|uniq -c >../fragnum/"+ID+".fragnum.txt\n"
"#/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python ../bin/01.whole_chrom_bin.v2.py ../data/"+ID+".fragments.tsv.gz ../fragnum/"+ID+".fragnum.txt ../result/01.bin_distribution/"+ID+".preheatmap.txt\n"
"/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/Rscript ../bin/02.bin_density_2.R ../result/01.bin_distribution/"+ID+".preheatmap.txt ../result/02.filter_7_peak/"+ID+"\n"
"/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/Rscript ../bin/03.candidate_cell_density.R ../result/02.filter_7_peak/"+ID+"_candidate.txt ../result/03.filtercell_15cells/"+ID+"\n" 
"/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python ../bin/04.get_pos.py ../result/03.filtercell_15cells/"+ID+"_candidate_bin.txt ../result/04.region/"+ID+"_candidate_pos.txt\n" 
"less ../result/04.region/"+ID+"_candidate_pos.txt |sort -k1,1 -k2,2n > ../result/04.region/"+ID+"_candidate_pos_sort.txt\n"
"/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python ../bin/04.hebing_bin_step2.v3.py ../result/02.filter_7_peak/"+ID+"_candidate.txt ../result/04.region/"+ID+"_candidate_pos_sort.txt ../result/04.region/"+ID+"_candidate_pos.txt\n"
"/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python ../bin/05.find_gene.py  ../result/04.region/"+ID+"_candidate_pos.txt ../result/05.getgene/"+ID+"_gene.txt\n" 
"/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python ../bin/06.tssportion_region_v3.py ../result/05.getgene/"+ID+"_gene.txt ../data/"+ID+".fragments.tsv.gz ../result/06.celltssproportion/"+ID+"_tssportion_region.txt\n" 
"/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/Rscript ../bin/07.tss_filtered.R ../result/06.celltssproportion/"+ID+"_tssportion_region.txt ../result/05.getgene/"+ID+"_gene.txt ../result/07.filterTSS/"+ID+"\n"
"/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python ../bin/08.ecDNA_cell_mtx.py ../data/"+ID+".fragments.tsv.gz ../fragnum/"+ID+".fragnum.txt ../result/07.filterTSS/"+ID+"_filterd_region.txt ../result/08.ecDNA_cell_mtx/"+ID+".ecDNA_cell_mtx.txt"
)
