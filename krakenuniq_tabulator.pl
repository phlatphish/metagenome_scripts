#!/usr/bin/env perl

#? Flexibly collate read counts from a set of krakenuniq reports

use strict;
use FileHandle;
use GD::Graph::hbars;
use List::Util qw(sum min max);

#
# specify these as args
#
my ($data_dir,$tax_rank,$infile,$outid,$threshold,$pattern) = @ARGV;


#
# some input checking
#
if(((scalar @ARGV) < 5) || ((scalar @ARGV) > 6))
{
   print "Usage: $0 data_dir_path {species|genus} infile outid threshold [pattern]\n"; 
   exit -1;
}

if (! -d $data_dir)
{
   print "Missing directory: $data_dir\n";
   exit -1;
}

#
# CONFIGURE
#
# Directories in which to find krakenuniq reports
#
#my @dir_list = qw(S_CF_B11
#                  S_CF_B15
#                  S_CF_B16
#                  S_CF_T15
#                  S_CF_T4
#                  S_GC_1
#                  S_GC_2
#                  S_OP_3
#                  S_OP_4
#                  S_OP_nod_1
#                  S_OP_lobe_1
#                  S_OP_nod_2
#                  S_OP_lobe_2);

#
# Labels for graph samples
#
my @nod_list = qw(CF_B11
                  CF_B15
                  CF_B16
                  CF_T15
                  CF_T4
                  GC_1
                  GC_2
                  OP_3
                  OP_4
                  OP_nod_1
                  OP_lobe_1
                  OP_nod_2
                  OP_lobe_2);

my @dir_list = @nod_list;

#
# get the number of samples for later calculations
#
my $num_samples = scalar @nod_list;
 
#
# Check things are there
#
DIR: for my $dir ( @dir_list )
{
   chdir($data_dir);

   if (! -d $dir)
   {
      print "Missing directory: $dir\n";
      exit -1;
   }
   
   chdir($dir);
   
   if (! -e $infile)
   {
      print "Missing file: $data_dir   $dir   $infile\n";
      exit -1;      
   }
}


#
# main data structure: read_data->dir->taxon_name = reads;
#
my $read_data;


#
# main data structure: kmer_data->dir->taxon_name = kmers
#
my $kmer_data;


#
# totals of kmers across nodules for each taxon: taxon_read_totals{taxon} = kmers
#
my %taxon_kmer_totals;


#
# totals of reads across nodules for each taxon: taxon_read_totals{taxon} = reads
#
my %taxon_read_totals;


#
# totals of kmers across nodules for each taxon: taxon_kmer_totals{taxon} = kmers
#
my %taxon_kmer_totals;


#
#  totals of reads for each nodule: sample_read_totals{sample} = reads;
#
my %sample_read_totals;


#
#  totals of kmers for each nodule: sample_read_totals{sample} = kmers;
#
my %sample_kmer_totals;

#
# classified reads and kmers for each sample
#
my %classified_reads;
my %classified_kmers;

#
# gather all the data in useful form
#
DIR: for my $dir ( @dir_list )
{
   chdir($data_dir);
   chdir($dir);
   open (IN,"<$infile");
      
   #print "$dir\n";

   # get each line from report  
   LINE: while(my $line = <IN>)
   {
      chomp $line;
            
      if($line =~ /root/)
      {
         my ($pct,$reads1,$reads2,$kmers,$dup,$cov,$taxid,$rank,$taxon) = split /\t/,$line;
         $classified_reads{$dir} = $reads1;
         $classified_kmers{$dir} = $kmers;
      }
      
      # consider only lines with chosen taxonomic rank
      if($line =~ /\d\t$tax_rank\t/)
      {
         my ($pct,$reads1,$reads2,$kmers,$dup,$cov,$taxid,$rank,$taxon) = split /\t/,$line;
         $taxon =~ s/^\s+//;                        
         $read_data->{$dir}->{$taxon} = $reads1;
         $kmer_data->{$dir}->{$taxon} = $kmers;
         $taxon_read_totals{$taxon}   += $reads1;
         $taxon_kmer_totals{$taxon}   += $kmers;
         $sample_read_totals{$dir}    += $reads1;
         $sample_kmer_totals{$dir}    += $kmers;

      }
   }
   close(IN);
}


#
# Write read and kmer tables
#
tabulate($data_dir,
         \@dir_list,
         \@nod_list,
         $read_data,
         $kmer_data,
         \%classified_reads,
         $outid,
         $threshold,
         \%taxon_read_totals,
         \%taxon_kmer_totals,
         $pattern);

       
#
# write two tables, one for reads and one for kmers, applying the kmer threshold
#
sub tabulate
{
   my ($s_data_dir,
       $s_dir_list,
       $s_nod_list,
       $s_read_data,
       $s_kmer_data,
       $s_classified_reads,
       $s_outid,
       $s_threshold,
       $s_taxon_read_totals,
       $s_taxon_kmer_totals,
       $s_pattern) = @_;
         
   #
   # reverse order taxa by read totals across samples
   #   
   my @ordered_taxa = (sort { $s_taxon_read_totals->{$b} <=> $s_taxon_read_totals->{$a} } keys %$s_taxon_read_totals);
   
   #
   # filter on the pattern, if there is one
   #
   if(defined $s_pattern)
   {
      @ordered_taxa = grep(/$s_pattern/, @ordered_taxa);
   }   
   
   #
   # filter on kmer threshold and fill in the taxa having no counts
   #
   for my $dir (@$s_dir_list)
   {
      # get millions of classified reads * kmer threshold. That is the bar each taxon needs to clear
      my $test_kmers = int ($s_classified_reads->{$dir}/1000000)*$s_threshold;
   
      for my $taxon(keys %$s_taxon_read_totals)
      {
         if($s_kmer_data->{$dir}->{$taxon} < $test_kmers)
         {
            $s_read_data->{$dir}->{$taxon} = 0;
            $s_kmer_data->{$dir}->{$taxon} = 0;
         }
                  
         unless (exists $s_read_data->{$dir}->{$taxon})
         {
             $s_read_data->{$dir}->{$taxon} = 0;
         }
         unless (exists $s_kmer_data->{$dir}->{$taxon})
         {
             $s_kmer_data->{$dir}->{$taxon} = 0;
         }
      }
   }
   
   #
   # set up two output files, one for reads and one for fractions
   #
   my $outfile_path1 = $s_data_dir . "/" . $s_outid . "_reads.csv";
   my $outfile_path2 = $s_data_dir . "/" . $s_outid . "_kmers.csv";
   my $fh_out1 = FileHandle->new(">$outfile_path1");
   my $fh_out2 = FileHandle->new(">$outfile_path2");


   print $fh_out1 ",";
   $, = ",";
   print $fh_out1 @$s_nod_list;
   print $fh_out1 "\n";

   #format nicely
   LINE: for my $taxon (@ordered_taxa)
   {
      print $fh_out1 "$taxon,";
 
      my @reads;
 
      for my $dir (@$s_dir_list)
      {
         my $reads = $s_read_data->{$dir}->{$taxon};
         push @reads,$reads;
      }
  
      print $fh_out1 @reads;
      print $fh_out1 "\n";
 
   }
   $fh_out1->close();
 
   
   print $fh_out2 ",";
   $, = ",";
   print $fh_out2 @$s_nod_list;
   print $fh_out2 "\n";

   #format nicely
   LINE: for my $taxon (@ordered_taxa)
   {
      print $fh_out2 "$taxon,";
 
      my @kmers;
 
      for my $dir (@$s_dir_list)
      {
         my $kmers = $s_kmer_data->{$dir}->{$taxon};
         push @kmers,$kmers;
      }
  
      print $fh_out2 @kmers;
      print $fh_out2 "\n";
 
   }
   $fh_out2->close();
   
}

  



















