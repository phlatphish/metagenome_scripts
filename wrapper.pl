#!/usr/bin/env perl

use strict;

# edit carefully--------------------------------------------------------------------
my $data_dir   = "/home/projects/redgillite/cjb_metagenome/seq";
my $output_dir = "/home/projects/redgillite/cjb_metagenome/krakenuniq/analysis_090722";
my $database   = "/home/projects/redgillite/cjb_metagenome/krakenuniq/db_052422";
# ----------------------------------------------------------------------------------


my %name_map = qw(Clone_639_2009.fasta	S_GC_2
                  B11.fasta				S_CF_B11
                  B15.fasta				S_CF_B15
                  B16.fasta				S_CF_B16
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
#   mkdir($name);
   chdir($name);
    
   my $fasta_file = $data_dir . "/" . $seq;
  
   print STDERR "name:       $name\n";
   print STDERR "fasta_file: $fasta_file\n";
   print STDERR "data_dir:   $data_dir\n";
   print STDERR "database:   $database\n";
   print STDERR "output_dir: $output_dir\n";
   print "\n";
       
   system("krakenuniq --db $database --threads 16  --report-file report  --unclassified-out unclassified --classified-out classified $fasta_file > output");   
             
}

