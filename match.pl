#equivlent of 'match' function in R
sub match {
  my $arrayRefa=shift;
  my $arrayRefb=shift;
  my @indexs;
  for my $val (@$arrayRefa) {
    push @indexs,  grep { $arrayRefb->[$_] eq $val } 0..$#$arrayRefb;
  }
  return @indexs;
}
