#!/usr/bin/perl
use strict;
use warnings;
die "usage: perl $0 infile.fas outfile.fas
# 
# 2016-12-1
# RSdetector_v.0.2.pl is used to extract RS domain containing proteins from the wanted proteome
# Note: dipeptide region is defined according to the definition of RS domain in plants 

###########################################################################################
#The required parameters:
infile.fas indicates the input proteome fasta file from such as the wanted species.
outfile.fas indicates the output file including RS domain containing proteins

#The optional parameter:
Didpepdide:RS, default is RS to search RS domain containing proteins in the proteome
Density:20, default is 20 indicating RS domain with at least 20% RS content
Width:50, default is 50 indicating RS domain with at least 50 amino acids
############################################################################################
 " unless (@ARGV>=2);

my $infile=$ARGV[0];
my $outfile=$ARGV[1];

my $dipeptide="(RS)|(SR)"; my $win_width=50; my $content=20;  
if (@ARGV >=3) {
  for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i]=~/Didpepdide/i) { my (undef,$tmp1) = split/:/,$ARGV[$i]; my $tmp2=reverse $tmp1; $dipeptide="($tmp1)|($tmp2)";  }
	if ($ARGV[$i]=~/Density/i) { (undef,$content) = split/:/,$ARGV[$i];  }
	if ($ARGV[$i]=~/Width/i) { (undef,$win_width) = split/:/,$ARGV[$i];  }
  }
}


my %hash;my $seq;my $id;my $key; my $density;my $get=0;my $window=3;
open IN, "$infile" or die "cannot find $infile\n";
open OUT,">$outfile";



while (<IN>) {
	chomp;

if ($_=~/^>/) {
	($id,undef)=split/\s/,$_;    
		 $id=~s/>//;

	}
	else {
		$hash{$id}.="$_";	
	}
}
foreach $key (keys %hash ) {
    $seq=$hash{$key};
	my $fragment;
	for (my $i=0;$i<(length($seq)-$win_width+1);$i++) {
	 $fragment=substr($seq,$i,$win_width);
my @a= split//,$fragment;
my $count=0;  my $flag=-1;
foreach  (my $i=0;$i<@a-1;$i++) {
		my $j=$i+1; 
		my $substr=$a[$i].$a[$j];
		if ( $substr=~/$dipeptide/i) {
			 if ($i==$flag) { $count ++;   
			 } else {$count +=2; }
			 $flag = $j;
		}
}
$density=$count/$win_width*100;
if ($density>=$content) {
		$get++;
		 print OUT ">$key Density:$density dipeptide:$dipeptide\n$seq\n";
		last;
	}
          
}
}
	print  "We found $get\n";

