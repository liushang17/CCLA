 ls *.fragments.tsv.gz |while read p ;do echo "less $p|awk '{print $4}'|sort|uniq -c >$(basename $p .fragments.tsv.gz).fragnum.txt";done
