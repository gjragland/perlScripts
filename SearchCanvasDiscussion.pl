#!/usr/bin/perl
#GJR 5/13/2019
#script to search list of student names against pdf versions of discussions from Canvas

use  CAM::PDF;
use strict;
my %students;
my @ids;

#tabbed list of student names and ids, including header
my $students<-open IN, "<studentId.txt";
while (<IN>) {
  next if $.==1;
  chomp;
  my @vals=split "\t";
  $students{$vals[1]}=$vals[0];
  push @ids, $vals[1];
}
close IN;

# array of pdf file names
my @filenames = ("Discussion1.pdf","Discussion2.pdf","Discussion3.pdf","Discussion4.pdf","Discussion5.pdf","Discussion6.pdf","Discussion7.pdf","Discussion8.pdf","Discussion9.pdf","Discussion10.pdf","Discussion11.pdf");



my %counts;
for my $id (@ids) {
  $counts{$id} = [];
}

for my $filename (@filenames) {


  my $text;
  my $doc = CAM::PDF->new($filename) || die "$CAM::PDF::errstr\n";
  for my $pagenum (1 .. $doc->numPages()) {
    $text = $text.$doc->getPageText($pagenum);
  }

  for my $id (@ids) {
    my ($number) = scalar( @{[ $text=~/$id/gi ]} );
    push @{ $counts{$id} }, $number;
  }

}


open OUT, ">StudentDiscussionCounts.txt";
$"="\t";
print OUT "Student\tid\t@filenames\n";
while (my ($id, $student) = each (%students)) {
  print OUT "$student\t$id\t@{$counts{$id}}\n";
}

close OUT;

