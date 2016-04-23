#!usr/bin/perl -w
# GJR 1/21/2015
#get only longest isoform from trinity output

use strict;
my $infile='Trinity.fasta';
my $outfile='TrinityLongestIso.fasta';

my %fastaHash;
my @order;
my $term=$/;
$/=">";
open IN, "<$infile";
while (<IN>) {
  if (/(^c\d+_g\d+)/) {
    my $geneId=$1;
    $_ =~ s/>//;
    my $seq=$_;
    $seq =~ s/^c\d+.*\n//;
     if (my $storedSeq=$fastaHash{$geneId}) {
       $storedSeq =~ s/^c\d+.*\n//;
       $fastaHash{$geneId}=$_ if length($seq) > length($storedSeq);
    } else {
      $fastaHash{$geneId}=$_;
      push @order, $geneId;
    }
  }
}

close IN;

open OUT, ">$outfile";
for my $geneId (@order) {
  print OUT ">$fastaHash{$geneId}";
}
close OUT;
