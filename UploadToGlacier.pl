#!/usr/bin/perl 
#GJR 3/9/2016
#use glacier-cmd utility to upload all R. pomonella RNAseq fastq files to Amazon Glacier Vault: 'PomonellaCerasiRNAseqDataFeb2016'

use strict;
use File::Find;

my $vault='PomonellaCerasiRNAseqDataFeb2016';
my $startDir='RawData/HO_PE100_20160112_1mismatch';

opendir my ($dh), $startDir;
my @dirNames=readdir $dh;
closedir $dh;
@dirNames=grep(/Project_Cevos/,@dirNames);


for my $dir (@dirNames) {
  find(\&upToGlacier,"$startDir/$dir");
}




  
sub upToGlacier {
  next unless /.*.fastq.gz/;
  print "$File::Find::fullname\n";
  #print "uploading $_ to amazon glacier vault $vault\n";
  #system("glacier-cmd upload PomonellaCerasiRNAseqDataFeb2016 c_Hi_2M_01_CTGAAGCT-TATAGCCT_L001_R1_001.fastq.gz --description "RCerasiCTGAAGCT-TATAGCCT_forwardPE"");
  my $dum;
}
