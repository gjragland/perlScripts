#!usr/bin/perl -w
# GJR 4/21/2012
#add read groups to sam file corresponding to barcode ids
# ALSO: changes overall quality scores of 255 to VAL to make compatible with GATK if 'changeQual' is selected
#mod 4/14/2017 adding command line options and tweaking regex

use strict;
use Getopt::Long;



my $idfile;
my $outfile;
my $idcol=1;
my $samfile;
my $newqual;
my $sepChar="_";
my $help;


GetOptions ('idfile=s' => \$idfile,
	    'idcol:i' => \$idcol,
	    'sepChar:s' => \$sepChar,
	    'samfile=s' => \$samfile,
	    'outfile:s' => \$outfile,
	    'changequal:i' => \$newqual,
	    'help' => \$help
	   );
if ($help) {
  print "USAGE: perl sampToSam.pl --idfile='yourIdsFile' --idcol='yourIdCol' --sepchar='yourSepChar' --samfile='yourSamFile'
     --outfile='yourOutputSamFile' --changequal <yourReplacemntQualityScore>
--idcol, --outfile, and --changeQual are optional, defaults = 1, 'yourSamFile.samples.sam', and NULL
--help: print usage info

Assign read groups to reads based on identifiers in the headers. Searches for
identifier information in fastq headers bookended (e.g., '_yourid_', '_' is
the sepchar) by the 'sepchar', which should appear exactly twice in the header
and NOT be present in the identifier. This conforms to output from the hotRAD
demultiplexer. Will also match the case where the identifier is the first
string of the header and terminates with 'sepchar' to conform to the output of
fuzzyDmRAD. If an identifier is found, adds an 'RG:Z:ID' field and RG headers.
Automatically adds a 'none' readgroup, and 'RG:Z:none' to reads in which
identifiers can not be located.

Required:
--idfile       tab-separated file with one column containing identifiers
--samfile      input sam file
Optional:
--idcol        column of tabbed file containing idendifiers (default=1)
--sepChar      character 'setting off' the identifier in the header. '_' is the
               default, non-alphanumeric characters may cause problems, though most should work
--outfile      name for output file (defaut=yourfile.samples.sam). Will write
               to the same directory as that containing the input file by default.
--changequal   integer value to replace all mapping quality scores of '255'
               (we use '99'). This is occasionally handy, e.g., older versions
               of bwa produced 255 values not compatible with GATK.

";
exit;    
}

die "\nError: User must specify file containing ids and input sam file: try 'perl sampToSam.pl --help'\n\n" unless $idfile and $samfile;



#decrement id column number to work with perl 0-indexing
$idcol--;
#escape special characters
$sepChar=quotemeta $sepChar;
unless ($outfile) {
  $outfile=$samfile;
  $outfile =~ s/\.sam//;
  $outfile = $outfile.'.samples.sam';
}
my %barcodeIds;
my @barcodeArray;
open INFILE, "<$idfile" or die "\nCan't locate idfile: check filename or path\n";
while (<INFILE>) {
  chomp;
  my @vals = split "\t";
  my $id=$vals[$idcol];
  #optional, add regex to remove extraneous info when ids from idfile don't perfectly match ids in samfile
  #$id =~ s/_sequence//;
  $barcodeIds{$id} = 0;
  push(@barcodeArray,$id);
}
close INFILE;

open OUTFILE, ">$outfile";
open INFILE, "<$samfile" or die "\nCan't locate samfile: check filename or path\n";
my $switch=0;
while (<INFILE>) {
  if (/^\@/) {
    print OUTFILE "$_";
    next;
  }
  if ($_ !~ m/^\@/ and $switch == 0) {
    $switch = 1;
    for my $id (@barcodeArray) {
      print OUTFILE "\@RG\tID\:$id\tSM\:$id\n";
    }
    print OUTFILE "\@RG\tID\:none\tSM\:none\n";
  }
  chomp;
  #replace '255' with new quality score
  if ($newqual) {s/^\S+\s255\s/\t$newqual\t/}
  print OUTFILE "$_";
  #modify regex to custom-grab sequence id
  my $id;
  if (/${sepChar}([\w\d\.]+?)${sepChar}/) {$id=$1}
  elsif (/^([\w\d\.]+?)${sepChar}/) {$id=$1}
  if ($id and exists($barcodeIds{$id})) {
    print OUTFILE "\tRG\:Z\:$id\n";
  } else {
    print OUTFILE "\tRG\:Z\:none\n";
  }
}
close INFILE;
close OUTFILE;
