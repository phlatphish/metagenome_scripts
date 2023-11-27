#!/usr/bin/env perl

#? create a folder tree and run krakenuniq on some metagenome sequences

use strict;

# edit carefully--------------------------------------------------------------------
my $data_dir   = "/fasta/file/location";
my $output_dir = "/output/directory";
my $database   = "/krakenuniq/index/location";
# ----------------------------------------------------------------------------------

# map the fasta file names to their (new) output folder names
my %name_map = qw(Clone_639_2009.fasta	S_GC_2
                  B11.fasta				S_CF_B11
                  B15.fasta				S_CF_B15
                  B16.fasta				S_CF_B16s
                  T4.fasta					S_CF_T4
                  N-1-L.fasta				S_OP_lobe_1
                  N-1-N.fasta				S_OP_nod_1
                  N-2-L.fasta				S_OP_lobe_2
                  N-2-N.fasta				S_OP_nod_2
                  Aberd_2008.fasta		S_OP_3
                  2010.fasta				S_OP_4
                  1989.fasta				S_GC_1                 
                  T15.fasta				S_CF_T15);

while ( my($seq,$name) = each %name_map )
{
   chdir($output_dir);
   mkdir($name);
   chdir($name);
    
   my $fasta_file = $data_dir . "/" . $seq;
  
   print STDERR "name:       $name\n";
   print STDERR "fasta_file: $fasta_file\n";
   print STDERR "data_dir:   $data_dir\n";
   print STDERR "database:   $database\n";
   print STDERR "output_dir: $output_dir\n";
   print "\n";
       
   system("krakenuniq --db $database --threads 16  --report-file report.txt  --unclassified-out unclassified.fa --classified-out classified.fa $fasta_file > output.txt");   
             
}

