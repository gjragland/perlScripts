#!/usr/bin/perl -w
#GJR 12/07/2012
# mod 5/24/2017 GJR
#filter vcf files based on MAF, depth, missingness, SNP quality score, Genotype Quality score, and binomial test for overassembly 

# specify packages and specify filter variables
use strict;
use Math::CDF;
use List::MoreUtils qw(first_index);

$"="\t";
my $minAltFreq=0.05;
my $minDepth=50;
my $maxMissing=100;
my $minGQ=15;
my $minQual=21;
my $binomThresh=0.05;
die 'must specify infile and outfile' unless @ARGV == 2;
my $infile="<".shift @ARGV;
my $outfile=">".shift @ARGV;

#my $infile="<snps.raw.vcf";
#my $outfile=">snps.05-300-100.vcf";
my $nSnps=0;





open INFILE, "$infile";
open OUTFILE, "$outfile";
while (<INFILE>) {
  if (/^\#/) {
    print OUTFILE "$_";
    next;
  }
  chomp;
  my @vals=split "\t";

  #primary filters for minimum depth, minimum quality, bi-allelic loci, and minimum alternate allele frequency
  next if $vals[4] =~ m/\S\S/;
  my $qual=$vals[5];
  next if $qual < $minQual;
  /DP=(\d+)/;
  next if $1 < $minDepth;
  /AF=([\d\.e\-]+)/;
  next if $1 < $minAltFreq and $1 < (1-$minAltFreq);
  my $nMissing=0;
  my @fields = split ":", $vals[8];
  my @GQind = grep { $fields[$_] =~ /GT/ } 0..$#fields; #get index of GQ field
  my $GQind = first_index { $_ eq 'GQ' } @fields;
  my $valsInd=8;
  for my $info (@vals[9..$#vals]) {
    $valsInd++;
    my @ivals = split ":", $info; #individual genotype-associated values
    if ($info =~ m/\.\/\./) {
      $nMissing++;
      next;
    }
    my $GQ = $ivals[$GQind];
    if ($GQ < $minGQ) {
      $vals[$valsInd] = "./.";
      $nMissing++;
    }
  }
  next unless $nMissing < $maxMissing;

  #binomial filter
  my $contig=$vals[0];
  my $allele1=0;
  my $allele2=0;
  for my $info (@vals[9..$#vals]) {
    next if $info =~ m/\.\/\./;
    my @AllInfo=split ":", $info;
    my @Counts=split ",", $AllInfo[1];
    if ($AllInfo[0] =~ m/0\/1/) {
      $allele1=$allele1+$Counts[0];
      $allele2=$allele2+$Counts[1];
    }
  }
  my @acounts=sort {$a <=> $b} ($allele1,$allele2);
  my $prob = 0.5;
  if ($allele1 > 0 and $allele2 > 0) {
    $prob = &Math::CDF::pbinom($acounts[0], ($acounts[0]+$acounts[1]), 0.5);
  }
  if ( $prob < (1-$binomThresh) and $prob > $binomThresh) {
    print OUTFILE "@vals\n";
    $nSnps++;
  }


}
close INFILE;
close OUTFILE;

print "Number SNPs Retained: $nSnps\n";



# sub getField {
#   my $fname=shift; # target field name
#   my $fline=shift; #vcf column with field names
#   my $line=shift; #vcf column with filed values
#   my @fields=split ":", $fline;
#   my @vals=split ":", $line;
#   #my $ind=grep $fname, @fields;
#   my $ind = grep { $fields[$_] =~ /$fname/ } 0..$#fields;
#   my $dum;
  
# }
