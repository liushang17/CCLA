use warnings;
use strict;

my $infile = $ARGV[0];
my $outdir = $ARGV[1];
my $pre = $ARGV[2];

print "perl /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/00.fragnum.pl $infile $outdir/$pre.fragnum.txt $outdir/$pre.fragment.new.txt\n";
print "/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/01.whole_chrom_bin.v2.py $outdir/$pre.fragment.new.txt $outdir/$pre.fragnum.txt $outdir/$pre.pheatmap.txt\n";
print "/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/Rscript /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/02.bin_density_2.R $outdir/$pre.pheatmap.txt $outdir/$pre\n";
print "/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/Rscript /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/03.candidate_cell_density.R $outdir/$pre.candidate.txt $outdir/$pre\n";
print "/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/04.get_pos.py $outdir/$pre.candidate_bin.txt $outdir/$pre.candidate_pos.txt\n";
print "less $outdir/$pre.candidate_pos.txt |sort -k1,1 -k2,2n > $outdir/$pre.candidate_pos_sort.txt\n";
print "/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/04.hebing_bin_step2.v3.py $outdir/$pre.candidate.txt $outdir/$pre.candidate_pos_sort.txt $outdir/$pre.candidate_pos.txt\n";
print "/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/05.find_gene.py $outdir/$pre.candidate_pos.txt $outdir/$pre.gene.txt\n";
print "/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/06.tssportion_region_v3.py $outdir/$pre.gene.txt $outdir/$pre.fragment.new.txt $outdir/$pre.tssportion_region.txt\n";
print "/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/Rscript /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/07.tss_filtered.R $outdir/$pre.tssportion_region.txt $outdir/$pre.gene.txt $outdir/$pre\n";
print "/hwfssz1/ST_SUPERCELLS/Reference/software/miniconda3_py3.9/bin/python /hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/version.01/08.ecDNA_cell_mtx.py $outdir/$pre.fragment.new.txt $outdir/$pre.fragnum.txt $outdir/$pre.filterd_region.txt $outdir/$pre.ecDNA_cell_mtx.txt\n";
