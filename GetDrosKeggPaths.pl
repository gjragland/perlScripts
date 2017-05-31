#!/usr/bin/perl -w
#GJR 10/1/2015
#query KEGG database to map drosophila ids to kegg ids


use strict;
use LWP::Simple;

# File containin all kegg pathways, downloaded from 'http://www.genome.jp/kegg-bin/get_htext'
my $infile="br08901.keg";
my $outfile="DrosKeggPaths.txt";

my %KeggPath;
open IN, "<$infile";
while (<IN>) {
  chomp;
  if (/^C\s+(\d+)\s+(\S.*)/) {
    $KeggPath{$1} = $2;
  }
}
close IN;

open OUT, ">$outfile";
print OUT "KeggPathway\tPathwayDescription\tDrosId\tKeggId\\n";
for my $kegg (keys(%KeggPath)) {
  my $desc = $KeggPath{$kegg};
  my @genes = retrieve_kegg($kegg);
  next if $genes[0] =~ m/NoHit/;
  for my $gene (@genes) {
    print OUT "$kegg\t$desc\t$gene\n"
  }
}
close OUT;

#-----subs-------------

#query kegg database via http
sub retrieve_kegg {
  my $query=shift;
  my $base='http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:dme';
  #form http query
  $query="$base$query";
  my @genes;
  #perform query through LWP::Simple function 'get'
  my $content = get $query;
  if (defined $content and $content !~ m/<b>No link information was found/) {
    while ($content =~ m/<a.*?(CG\d+).*?\;\s+(K\d+)/g) {
      push @genes, "$1\t$2";
      my $dum;
    }
  }
  push @genes, "NoHit" if @genes==0;
  return @genes;
}


#http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:dme00020
