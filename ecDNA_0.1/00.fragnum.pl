use warnings;
use strict;

open I,$ARGV[0] or die "$!";
my %cell;
while(<I>){
	chomp;
	my @word = split " ",$_;
	$cell{$word[3]} += $word[4];
}
close I;

open O1,">$ARGV[1]" or die "$!";
foreach(keys %cell){
	next if($cell{$_} < 1000);
	print O1 "$cell{$_}\t$_\n";
}
close O1;

open I,$ARGV[0] or die "$!";
open O2,">$ARGV[2]" or die "$!";
while(<I>){
        chomp;
        my @word = split " ",$_;
        if($cell{$word[3]} >= 1000){
		 print O2 "$_\n";
	}
}
close I;
close O2;
