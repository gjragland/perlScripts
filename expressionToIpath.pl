#!/usr/bin/perl -w
#GJR 10/1/2015
# create data lists for input into iPath
#usage CreateIpath.pl -keggFile <kegg map file> -eFile <expression file> -outFile <output file>

use strict;
use Getopt::Long;

#specify default dros cg to kegg id map file
my $keggMapFile="DrosKeggPaths.txt";
#specify default colors for up and down regulation
my $up="RGB(0,255,0)";
my $down="RGB(255,0,0)";  



my $keggMapFileTemp;
my $expressionFile;
my $outfile;
my $invert;
  
GetOptions ("keggFile=s" => \$keggMapFileTemp,    
            "eFile=s"   => \$expressionFile,      
            "outFile=s"  => \$outfile,
           "invert" => \$invert);

if ($keggMapFileTemp) {$keggMapFile=$keggMapFileTemp}

my %DrosToKegg;
open IN, "<$keggMapFile";
while (<IN>) {
  next if $.==1;
  chomp;
  my @vals=split "\t";
  $DrosToKegg{$vals[2]}=$vals[3];
}
close IN;

my $annoIndex;
my $fdrIndex;
my $valIndex;
open IN, "<$expressionFile";
open OUT, ">$outfile";
while (<IN>) {
  chomp;
  my @vals=split "\t";
  if ($.==1) {
    #$annoIndex = 
    ($annoIndex) = grep { $vals[$_] eq "ANNOTATION_SYMBOL" } 0..$#vals;
    ($fdrIndex) = grep { $vals[$_] eq "FDR" } 0..$#vals;
    ($valIndex) = grep { $vals[$_] eq "logFC" } 0..$#vals;
    next;
  }
  if ($vals[$fdrIndex] < 0.1 and my $ko = $DrosToKegg{$vals[$annoIndex]}) {
    
    my $color=$up;
    my $lfc=$vals[$valIndex];
    $lfc = $lfc*(-1) if $invert;
    $color=$down if $lfc < 0;
    my $width = assign_width($lfc);
    print OUT "$ko\t$color\tW$width\n";
  }
}
close OUT;
close IN;

#subs

sub assign_width {
  my $lfc=shift;
  #scale: 0-.6 = 1, .6 - 2 = 2, 2 - 4 = 3, >4 = 4;
  my $width;
  $lfc=abs($lfc);
  if ($lfc <= 0.6) {
    $width=1;
  } elsif ($lfc > 0.6 and $lfc <= 2 ) {
    $width=3;
  } elsif ($lfc > 2 and $lfc <= 4) {
    $width=5;
  } elsif ($lfc > 4) {
    $width=7;
  }
  return $width;
}
