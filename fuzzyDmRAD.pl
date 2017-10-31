#!usr/bin/perl -w
#GJR 3/20/2017
#Modified 10/31/2017, added support for gzipped files

# Demultiplex RADseq data with in-line barcodes
# Neither Lauren's 'trimmer.py' nor the stacks demultiplexer appear to check whether allowing barcode mismathches produce ambiguous mappings
# Here, we use fuzzy matching, tossing any putative barcodes in the sequence that fuzzy match to more than one of the user-supplied list of actual barcodes
# Current version assumes 1) sequence starts with barcode, 2) a 1bp protector base between the barcode and cutsite, 3) cutsite
#    (hard coded as 'AATT' for now) following the protector base
# Current version also hardcodes max/min barcode lengths ($m and $n)
# Finally, current version only handles single-end data
# uses the two column, user-supplied barcode file to map ids to reads
# REQUIRES the String::Approx package

# logic flow:
# 1.	Find cutsite with exact match, allowing for protector base
# 2.	If exact barcode match, retain read and map to id
# 3.	If unambiguous fuzzy match, reatin read (try fuzzy match for 1 - max mismatches, substitutions only) and map to id
# 4.	If no exact cutsite match, attempt to fuzzy match cutsite (max one substitution), starting at max barcode length, working backwards to min barcode length
#       Also requires exact match to protector base in proper position
# 5.	If fuzzy match to the cutsite
# 6.	If exact match to barcode, retain read and map to id
# 7.	If unambiguous fuzzy match, reatin read (try fuzzy match for 1 - max mismatches, substitutions only) and map to id
# 8.	END




use strict;
use warnings;
use String::Approx qw(amatch);
use Getopt::Long;



my $outfile="fuzzyDmRAD.out.fq";
#my $infile="/media/raglandlab/ExtraDrive2/IpsPE_RADseq/Project_Martin64/Plate6_16_Plate4_712_NoIndex_L008_R1_001.fastq";
my $infile;
my $logfile="fuzzyDmRAD.log";
#my $bcodeFile="barcodes_64.txt";
my $bcodeFile;
my $help;
my %barcodes;
my %barcodeCounts;


GetOptions ('fastq=s' => \$infile,
	    'bcodefile=s' => \$bcodeFile,
	    'logfile:s' => \$logfile,
	    'outfq:s' => \$outfile,
	    'help' => \$help
	   );
if ($help) {
  print "USAGE: perl fuzzyDmRAD.pl --fastq='yourFastqFile' --bcodefile='yourBcodeFile'
     --logfile='yourLogFile' --outfq='yourOutputFqFile'
--logfile and --outfq are optional, defaults = 'fuzzyDmRAD.log' and 'fuzzyDmRAD.out.fq'
--help: print usage info
";
exit;    
}

die "\nError: User must specify barcode file and input fastq file: try 'perl fuzzyDmRAD.pl --help'\n\n" unless $infile and $bcodeFile;


open my $outfileHandle, ">$outfile";

open IN, "<$bcodeFile";
while (<IN>) {
  chomp;
  my @vals=split "\t";
  $vals[0] = uc ($vals[0]);
  $barcodes{$vals[0]}=$vals[1];
  $barcodeCounts{$vals[0]}=0;
}
close IN;


#min length of barcode
#my $n=5;
my $n=8;
#max length of barcode
#my $m=6;
my $m=10;
#max mismatches
my $maxMismatch=3;
my $cutsite="AATT";
#accounts for protector base b/t barcode and cutsite
my $protectorSeqLen=1;
my $len = length( $cutsite) + $protectorSeqLen;

my $nReads=0;
my $nPerfectMatchCutSite=0;
my $nPerfectMatchNoCutSite=0;
my $nFuzzyMatchCutsite=0;
my $nFuzzyMatchNoCutsite=0;

#added support for zipped files, Halloween, 2017
if ($infile =~ /.gz$/) {
  open(IN, "gunzip -c $infile |") || die "can’t open pipe to $infile\n";
} else {
  open(IN, "<$infile") || die "can’t open $infile";
}


while (<IN>) {
  #read 4 lines at a time into array
  my @lines = ($_,'','','');
  for (my $i=1; $i <= 3; $i++) {
    $lines[$i]= <IN>;
  }
  if ($.==3968) {
    my $dum;
  }
  $nReads++;
  my $sequence = $lines[1];
  #if ($sequence =~ m/\S{$n,$m}/) {
  if ($sequence =~ m/(^\S{$n,$m})\S{$protectorSeqLen}${cutsite}(\S+)/) { # if cutsite is found, perform fuzzy match allowing only substitutions (length of bcode is known)
    #putative barcode
    my $bcode = $1;
    #biological sequence following cutsite
    my $bioSeq=$2;
    if (exists $barcodes{$bcode}) {
      mapIdtoRead($barcodes{$bcode},$bioSeq,$outfileHandle,@lines);
      $barcodeCounts{$bcode}++;
      $nPerfectMatchCutSite++;
    } else {
      #fuzzy match if no exact match, for all possible bcode lengths and allowing only subs
      my $matchedBcode = fuzzy ($maxMismatch, $bcode,\%barcodes);
      if ($matchedBcode) {
	mapIdtoRead($barcodes{$matchedBcode},$bioSeq,$outfileHandle,@lines);
	$nFuzzyMatchCutsite++;
	$barcodeCounts{$matchedBcode}++;
      }
    }
  } else { #if can't find cutsite, try all leading sequences up to the max barcode length)
    for (my $i=$m; $i >= $n ; $i--) {
      my ($bcode, $bioSeq) = $sequence =~ m/(^\S{$i})(\S+)/;
      #increase stringency by requiring the correct protector base and a default fuzzy match to the cutsite
      #next unless $bioSeq =~ m/^C/ and amatch(substr($bioSeq,1,4),[ 'i',2,"S2","I0","D0" ],$cutsite);
      #next unless amatch(substr($bioSeq,1,4),[ 'i',2,"S2","I0","D0" ],$cutsite);
      next unless amatch(substr($bioSeq,0,5),[ 'i',1,"S1","I0","D0" ],"C".$cutsite);
      if (exists $barcodes{$bcode}) {
	$bioSeq =~ s/^\S{$len}//;
	mapIdtoRead($barcodes{$bcode},$bioSeq,$outfileHandle,@lines);
	$barcodeCounts{$bcode}++;
	$nPerfectMatchNoCutSite++;
	#print "$.,";
	last;
      } else {
	my $matchedBcode = fuzzy ($maxMismatch, $bcode,\%barcodes);
	if ($matchedBcode) {
	  $bioSeq =~ s/^\S{$len}//;
	  mapIdtoRead($barcodes{$matchedBcode},$bioSeq,$outfileHandle,@lines);
	  $barcodeCounts{$matchedBcode}++;
	  $nFuzzyMatchNoCutsite++;
	  last;
	}
      }
    }
  }
}

close IN;
close $outfileHandle;


my $retained = $nPerfectMatchCutSite + $nPerfectMatchNoCutSite + $nFuzzyMatchCutsite + $nFuzzyMatchNoCutsite;
my $perRetained = $retained/$nReads;
print "
Total reads: $nReads
Retained reads: $retained
Percent reads retained: $perRetained
Retained read categories ...
Perfect match with cutsite: $nPerfectMatchCutSite
Perfect match without cutsite: $nPerfectMatchNoCutSite
Fuzzy match with cutsite: $nFuzzyMatchCutsite
Fuzzy match without cutsite: $nFuzzyMatchNoCutsite
";

open OUT, ">$logfile";
print OUT "Barcode\tReadCount\n";
while (my ($bcode,$count) = each %barcodeCounts) {
  print OUT "$bcode\t$count\n";
}
close OUT;




#subs

#atcgc
#tcgc

# e.g., mapIdtoRead($barcodes{$bcode},$bioSeq,$outfileHandle,@lines)
sub mapIdtoRead {
  my ($id,$seq,$fileHandle,@lines) = @_;
  my $l1=length($lines[1]);
  my $l2=length($lines[3]);
  $lines[1] = "$seq\n";
  chomp $lines[3];
  #$lines[0] = "${id}_".$lines[0];
  $lines[0] =~ s/^\@/\@${id}_/;
  my $start = length( $lines[3]) - length ($seq);
  my $len = length ($seq); 
  $lines[3] = substr($lines[3],$start,$len)."\n";
  $l1=length($lines[1]);
  $l2=length($lines[3]);
  print $fileHandle @lines;
}

#e.g,  fuzzy ($maxMismatch, $bcode,\%barcodes)
sub fuzzy {
  my $maxMismatch = shift;
  my $bcode = shift;
  my %barcodes = %{$_[0]};
  my $i=1;
  my @matched;
  my @MODS;
  my $matchedBcode;
  while ($i <= $maxMismatch and ! $matchedBcode) {
    #specify $i possible substitutions (no indels)
    #@MODS = ('i',0,"S$i");
    @MODS = ('i',$i,"S$i","I0","D0");
    @matched = amatch($bcode, [ @MODS ], keys(%barcodes));
#    @matched = amatch($bcode, keys(%barcodes));
    if (@matched == 1) {
      $matchedBcode = $matched[0];
    }
    $i++;
  }
  return($matchedBcode);
}

