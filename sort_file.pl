#!/usr/bin/perl

$cmd="cat $ARGV[0] | awk \'{if(\$1==\"X\" || \$1==\"Y\" || \$1==\"MT\"){print;}}\' > $ARGV[0].XY";
print STDERR "$cmd\n";
system($cmd);
$cmd="cat $ARGV[0] | awk \'{if(\$1!=\"X\" && \$1!=\"Y\" && \$1!=\"MT\"){print;}}\' > $ARGV[0].num";
print STDERR "$cmd\n";
system($cmd);
$cmd="sort -k 1,1d -k 2,2n $ARGV[0].XY > $ARGV[0].XY.sorted";
print STDERR "$cmd\n";
system($cmd);
$cmd="sort -k 1,1n -k 2,2n $ARGV[0].num > $ARGV[0].num.sorted";
print STDERR "$cmd\n";
system($cmd);
$cmd="cat $ARGV[0].num.sorted $ARGV[0].XY.sorted > $ARGV[0].sorted";
print STDERR "$cmd\n";
system($cmd);
$cmd="rm $ARGV[0].num* $ARGV[0].XY*";
print STDERR "$cmd\n";
system($cmd);
