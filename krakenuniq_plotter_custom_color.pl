#!/usr/bin/env perl

#? Flexibly plot read counts from a set of krakenuniq outputs

use strict;
use FileHandle;
use GD::Graph::hbars;
use List::Util qw(sum min max);

#
# specify these as args
#
my ($data_dir,$tax_rank,$infile,$outid,$title,$threshold,$pattern) = @ARGV;

#
# max number of taxa in each chart bar
#
my $num_taxa = 20;

#
# some input checking
#
if(((scalar @ARGV) < 6) || ((scalar @ARGV) > 7))
{
   print "Usage: $0 data_dir_path {S|G} infile outid title [pattern]\n"; 
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
# Directories in which to find reports
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
   open (IN,"<$infile") || die "Unable to open $infile\n";
      
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



plot  ($data_dir,
       \@dir_list,
       \@nod_list,
       $read_data,
       $kmer_data,
       \%classified_reads,
       \%classified_kmers,
       $num_taxa,
       $outid,
       $title,
       \%taxon_read_totals,
       \%taxon_kmer_totals,
       $threshold,
       $num_samples,
       $pattern);

      
        
       
sub plot
{
   my ($s_data_dir,
       $s_dir_list,
       $s_nod_list,
       $s_read_data,
       $s_kmer_data,
       $s_classified_reads,
       $s_classified_kmers,
       $s_num_taxa,
       $s_outid,
       $s_title,
       $s_taxon_read_totals,
       $s_taxon_kmer_totals,
       $s_threshold,
       $s_num_samples,
       $pattern) = @_;

  
   #
   # filter the read and kmer data based on kmers per million reads threshold
   #
   
   # top level, list retained taxa
   my @retained_taxa;
   
   # for each dir/sample
   while( my ($dir,$dir_data) = each %$s_read_data)
   {
      # get millions of classified reads * kmer threshold. That is the bar each taxon needs to clear
      my $test_kmers = int ($s_classified_reads->{$dir}/1000000)*$s_threshold;
   
      # for each taxon
      while( my ($taxon,$reads) = each %$dir_data)
      {
         # get number of kmers
         #print "$dir\t$taxon\t$s_kmer_data->{$dir}->{$taxon}\n";
 
         # pass test, add to retained taxa
         if($s_kmer_data->{$dir}->{$taxon} > $test_kmers)
         {
            push @retained_taxa,$taxon;
            #print "$dir,$taxon,$s_kmer_data->{$dir}->{$taxon},[$test_kmers]\n";
         }
         # fails test, delete from read_data and kmer_data
         else
         {
            delete $s_read_data->{$dir}->{$taxon};
            delete $s_kmer_data->{$dir}->{$taxon};
         }
      }       
   }
   
         
   # go through edited read_data and add zeroes to missing points
   for my $dir (keys %$s_read_data)
   {
      for my $taxon (@retained_taxa)
      {
          if(!exists $s_read_data->{$dir}->{$taxon})
         {
            $s_read_data->{$dir}->{$taxon} = 0;
            $s_kmer_data->{$dir}->{$taxon} = 0;
         }       
      }
   }
    
   #
   # reverse order taxa by read totals across samples
   #   
   my @ordered_taxa = (sort { $s_taxon_read_totals->{$b} <=> $s_taxon_read_totals->{$a} } keys %$s_taxon_read_totals);
   
   my @ordered_retained_taxa;
   
   #
   # get out the filtered/retained taxa in the same order
   #
   for my $ordered_taxon (@ordered_taxa)
   {
      if(grep /$ordered_taxon/,@retained_taxa)
      {
          push @ordered_retained_taxa,$ordered_taxon;
      }
   }
   
   #print "@ordered_retained_taxa\n";
   
   #
   # We plot just a subset of the taxa
   #
   my @plot_taxa; 
   my @label_taxa;
   
   #
   # If a pattern is suggested, select for it or take the top N
   #  
   if(defined $pattern)
   {
      @plot_taxa = grep(/$pattern/, @ordered_retained_taxa);
      @plot_taxa = @plot_taxa[0 .. ($s_num_taxa - 1)];
   }
   else
   {
      @plot_taxa = @ordered_retained_taxa[0 .. ($s_num_taxa - 1)];
   }
  
   #
   # CONFIGURE. May comment out
   #
   # For select plots, redefine plot taxa and legend taxa here.
   #
   
   print "KrakenUniq\n";
   
   if ($tax_rank eq "genus" || $tax_rank eq "G")
   {   
      print "Condition 1\n";
      print "tax_ranks: $tax_rank\n";
      print "pattern:   $pattern\n";
      
      
      @plot_taxa = qw(Frankia
                      Streptomyces
                      Bradyrhizobium
                      Cupriavidus
                      Pseudomonas
                      Kutzneria
                      Nocardia
                      Mycolicibacterium
                      Micromonospora
                      Mycobacterium
                      Nocardioides
                      Amycolatopsis
                      Actinoplanes
                      Rhodococcus
                      Mesorhizobium
                      Pseudonocardia
                      Microbacterium
                      Staphylococcus
                      Actinomadura
                      Variovorax); 
       
        @label_taxa = @plot_taxa;
   }                   
   
   elsif (($tax_rank eq "species" || $tax_rank eq "S") && $pattern eq "Frankia")
   {
      print "Condition 2\n";
      print "tax_rank: $tax_rank\n";
      print "pattern:  $pattern\n";
      
      
       @plot_taxa = ("Frankia sp. QA3",
                     "Frankia sp. ACN1ag",
                     "Frankia canadensis",
                     "Frankia sp. AvcI1",
                     "Frankia alni",
                     "Frankia sp. EAN1pec",
                     "Frankia sp. EI5c",
                     "Frankia sp. Cc1.17",
                     "Frankia inefficax",
                     "Frankia sp. BMG5.36",
                     "Frankia asymbiotica",
                     "Frankia discariae",
                     "Frankia sp. ArI3",
                     "Candidatus Frankia californiensis",
                     "Frankia sp. DC12",
                     "Frankia elaeagni",
                     "Frankia irregularis",
                     "Frankia sp. EUN1f",
                     "Candidatus Frankia datiscae",
                     "Frankia sp. R43");

       @label_taxa = ("Frankia sp. QA3 (1A)",
                     "Frankia sp. ACN1ag (1A)",
                     "Frankia canadensis (1A)",
                     "Frankia sp. AvcI1 (1A)",
                     "Frankia alni (1A)",
                     "Frankia sp. EAN1pec (3)",
                     "Frankia sp. EI5c (3)",
                     "Frankia sp. Cc1.17 (3)",
                     "Frankia inefficax (4)",
                     "Frankia sp. BMG5.36 (4)",
                     "Frankia asymbiotica (4)",
                     "Frankia discariae (3)",
                     "Frankia sp. ArI3 (1A)",
                     "Frankia californiensis (2)",
                     "Frankia sp. DC12 (4)",
                     "Frankia elaeagni (3)",
                     "Frankia irregularis (3)",
                     "Frankia sp. EUN1f (3)",
                     "Frankia datiscae (2)",
                     "Frankia sp. R43 (3)");
    }
                  
    
   elsif ($tax_rank eq "species" || $tax_rank eq "S")
   {
      print "Condition 3\n";
      print "tax_rank: $tax_rank\n";
      print "pattern:  $pattern\n";
      
      
       @plot_taxa = ("Frankia sp. QA3",
                     "Frankia sp. ACN1ag",
                     "Frankia canadensis",
                     "Kutzneria sp. CA-103260",
                     "Bradyrhizobium erythrophlei",
                     "Frankia inefficax",
                     "Bradyrhizobium diazoefficiens",
                     "Catenulispora acidiphila",
                     "Bradyrhizobium lablabi",
                     "Afipia sp. GAS231",
                     "Rhodopseudomonas palustris",
                     "Candidatus Phytoplasma ziziphi",
                     "Bradyrhizobium canariense",
                     "Bradyrhizobium sp. CCBAU 051011",
                     "Bradyrhizobium paxllaeri",
                     "Bradyrhizobium icense",
                     "Bradyrhizobium genosp. B",
                     "Rhodoplanes sp. Z2-YC6860",
                     "Bradyrhizobium sp. CCBAU 51753",
                     "Bradyrhizobium guangzhouense");
                     
       @label_taxa = @plot_taxa;
   }
   # END CONFIGURE

                    
#    #
#    # Colors borrowed from google sheets
#    #
#    my @hues = qw(#4285f4
#                  #ea4335
#                  #fbbc04
#                  #34a853
#                  #ff6d01
#                  #46bdc6
#                  #7baaf7
#                  #f07b72
#                  #71c287
#                  #fcd04f
#                  #ff994d
#                  #7ed1d7
#                  #b3cefb
#                  #f7b4ae
#                  #fde49b
#                  #aedcba
#                  #ffc599
#                  #b5e5e8
#                  #ecf3fe
#                  #fdeceb);

   #
   # Colors borrowed from Paul Tol (https://personal.sron.nl/~pault/) light and muted palettes together
   #
   my @hues = ('#88CCEE', 
               '#44AA99', 
               '#117733', 
               '#332288', 
               '#DDCC77', 
               '#999933',
               '#CC6677', 
               '#882255', 
               '#AA4499', 
               '#DDDDDD',
               '#BBCC33', 
               '#AAAA00', 
               '#77AADD', 
               '#EE8866', 
               '#EEDD88', 
               '#FFAABB', 
               '#99DDFF', 
               '#44BB99', 
               '#DDDDDD');

   
   #
   # sum the reads for each nodule for the selected subset of taxa
   #
   # nodule_read_totals->$dir = total reads for nodule
   #
   my $nodule_read_totals;
   
   for my $dir(@$s_dir_list)
   {
      for my $taxon (@plot_taxa)
      {
         $nodule_read_totals->{$dir} += $s_read_data->{$dir}->{$taxon};;
      }
   }
       
   #
   # structure for dir nodule read percentages instead of reads
   #
   # For stack plotting nodule_read_percents->dir->taxon = percent_of_nodule_reads
   #         
   my $nodule_read_percents;
   
   for my $dir (@$s_dir_list)
   {
      for my $taxon (@plot_taxa)
      {
         my $reads = $s_read_data->{$dir}->{$taxon};
         
         $nodule_read_totals->{$dir} = 1 if ($nodule_read_totals->{$dir} == 0);
         
         my $percent = ($reads/$nodule_read_totals->{$dir}) * 100;
         $nodule_read_percents->{$dir}->{$taxon} = $percent;
      }      
   }
    
    
   #
   # Table of proportions in the nodule list order
   #
   my $all_plot_percents;
   
   # 20 taxa
   for my $taxon (@plot_taxa)
   {     
      my $nodule_plot_percents;
           
      # 13 dirs
      for my $dir (@$s_dir_list)
      {    
         push @$nodule_plot_percents,$nodule_read_percents->{$dir}->{$taxon};
      }    
      
      # should get 20 arrays, each with 13 values
      push @$all_plot_percents,$nodule_plot_percents;
   } 


   my $graph = new GD::Graph::hbars(1000,1000);
   
   
   # orientation is flipped so y label is on the x
   
   $graph->set(y_label           => "Percent of reads",
               y_max_value       => 100,
               y_tick_number     => 4,
               bar_spacing       => 5,
               bargroup_spacing  => 10,
               long_ticks        => 0,
               l_margin          => 20,
               r_margin          => 20,
               t_margin          => 20,
               b_margin          => 20,
               title             => $s_title,
	            cumulate          => 1,
	            box_axis          => 0,
               );
               
   #
   # CONFIGURE, may comment out and use system fonts
   #
   # non standard location, but nice way to specify
   #
   #my $font = '/path/to/OpenSans-Regular.ttf';
   my $font = '/usr/share/fonts/open-sans/OpenSans-Regular.ttf';
   
   if ( ! -f $font )
   {
      die "Missing font file: $font\n";
   }
   # END CONFIGURE

   $graph->set( dclrs => \@hues ); 
   $graph->set( borderclrs => [undef] );
   $graph->set_title_font($font,18);      
   $graph->set_x_label_font($font,14);
   $graph->set_y_label_font($font,14);
   $graph->set_x_axis_font($font,14);
   $graph->set_y_axis_font($font,14);
   $graph->set_legend_font($font,12);
        
   $graph->set(legend_placement => 'BC', legend_marker_width => 15, legend_marker_height => 15);
   $graph->set_legend( @label_taxa );
        
   my @data = (\@nod_list,@$all_plot_percents);
      
   my $splotfile_path = $s_data_dir . "/" . $outid . ".png";
   
   print "\$splotfile_path: $splotfile_path\n";
  
   open  OUT,">$splotfile_path" || die "Unable to open $splotfile_path\n";
   print OUT $graph->plot(\@data)->png;
   close OUT;  
}


