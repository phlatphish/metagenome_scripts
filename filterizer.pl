#!/usr/bin/env perl

#
# Open a krakenuniq report and eliminate rows falling below a kmer threshold
#

use strict;

my($report,$kmers_per_million_reads) = @ARGV;

open (IN,"<$report");

my $header;

my $h = <IN>; $header = $header . $h;
$h = <IN>; $header = $header . $h;
$h = <IN>; $header = $header . $h;
$h = <IN>; $header = $header . $h;

my $classified_reads = -1;

print $header;

# get each line from report  
LINE: while(my $line = <IN>)
{
   chomp $line;
             
   if($line =~ /root/)
   {
      my ($pct,$reads1,$reads2,$kmers,$dup,$cov,$taxid,$rank,$taxon) = split /\t/,$line;
      $classified_reads = $reads1;
   }
   
   my $threshold = ($classified_reads/1000000)*$kmers_per_million_reads;
   
   my ($pct,$reads1,$reads2,$kmers,$dup,$cov,$taxid,$rank,$taxon) = split /\t/,$line;
   $taxon =~ s/^\s+//;
      
   #print "$rank,$taxon [$kmers] $threshold\n";
             
   if($kmers > $threshold)
   {
      print "$line\n";
   }
      
}
close(IN);



