#!/usr/bin/env perl

#? Flexibly plot read counts from a set of krakenuniq outputs

use strict;
use FileHandle;
use GD::Graph::hbars;
use List::Util qw(sum min max);

#
# specify these as args
#
my ($data_dir,$tax_rank,$infile,$outid,$title,$threshold,$pattern,$hack) = @ARGV;

# print "----------------------------\n";
# print "data_dir: $data_dir\n";
# print "hack: $hack\n";
# print "infile: $infile\n";
# print "outid: $outid\n";
# print "pattern: $pattern\n";
# print "tax_rank: $tax_rank\n";
# print "threshold: $threshold\n";
# print "title: $title\n";
# print "----------------------------\n";
# exit 0;


# args
#
# data_dir:  where to find the folders containing krakenuniq output
# tax_rank:  S or G
# infile:    name of the report file, usually "report"
# outid:     string to use to name the output of this script e.g. krakenuniq_genus_0
# title:     Title of the plotted graph, e.g. "Krakenuniq Genus Proportions (No kmer threshold)"
# threshold  kmer threshold per million reads. usually 0 or 2000
# pattern    Pattern to filter, like "Frankia", or 0 for no filter
# hack       0 or 1, to force taxon naming using a list

# examples:
#
# /home/projects/redgillite/cjb_metagenome/krakenuniq/metagenome_scripts/krakenuniq_plotter_hack.pl \
# /home/projects/redgillite/cjb_metagenome/krakenuniq/analysis_080124 \
# genus \
# report \
# krakenuniq_genus_0 \
# "Krakenuniq Genus Proportions (No kmer threshold)" \
# 0 \
# 0 \
# 0

#
# /home/projects/redgillite/cjb_metagenome/krakenuniq/metagenome_scripts/krakenuniq_plotter_hack.pl \
# /home/projects/redgillite/cjb_metagenome/krakenuniq/analysis_080124 \
# genus \
# report \
# krakenuniq_genus_0 \
# "Krakenuniq Genus Proportions (No kmer threshold)" \
# 0 \
# 0 \
# 1

#
# hack and pattern cannot both be true
#
if ($hack != 0  && $pattern != 0)
{
   print "hack and pattern are mutually exclusive\n"; 
   exit -1;   
}


#
# max number of taxa in each chart bar
#
my $num_taxa = 20;

#
# some input checking
#
if( scalar @ARGV != 8 )
{
   print "Usage: $0 data_dir_path {S|G} infile outid title threshold pattern hack\n"; 
   exit -1;
}

if (! -d $data_dir)
{
   print "Missing directory: $data_dir\n";
   exit -1;
}

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
# optionally, force taxon list for plotting
#
my @plot_taxa  = ();
my @label_taxa = ();


#
# CONFIGURE BELOW
#
if ($hack == 1)
{
      #print STDERR "The hack == 1\n";
      
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
                
#         @plot_taxa = ("Frankia sp. QA3",
#                      "Frankia sp. ACN1ag",
#                      "Frankia canadensis",
#                      "Frankia sp. AvcI1",
#                      "Frankia alni",
#                      "Frankia sp. EAN1pec",
#                      "Frankia sp. EI5c",
#                      "Frankia sp. Cc1.17",
#                      "Frankia inefficax",
#                      "Frankia sp. BMG5.36",
#                      "Frankia asymbiotica",
#                      "Frankia discariae",
#                      "Frankia sp. ArI3",
#                      "Candidatus Frankia californiensis",
#                      "Frankia sp. DC12",
#                      "Frankia elaeagni",
#                      "Frankia irregularis",
#                      "Frankia sp. EUN1f",
#                      "Candidatus Frankia datiscae",
#                      "Frankia sp. R43");
# 
#        @label_taxa = ("Frankia sp. QA3 (1A)",
#                      "Frankia sp. ACN1ag (1A)",
#                      "Frankia canadensis (1A)",
#                      "Frankia sp. AvcI1 (1A)",
#                      "Frankia alni (1A)",
#                      "Frankia sp. EAN1pec (3)",
#                      "Frankia sp. EI5c (3)",
#                      "Frankia sp. Cc1.17 (3)",
#                      "Frankia inefficax (4)",
#                      "Frankia sp. BMG5.36 (4)",
#                      "Frankia asymbiotica (4)",
#                      "Frankia discariae (3)",
#                      "Frankia sp. ArI3 (1A)",
#                      "Frankia californiensis (2)",
#                      "Frankia sp. DC12 (4)",
#                      "Frankia elaeagni (3)",
#                      "Frankia irregularis (3)",
#                      "Frankia sp. EUN1f (3)",
#                      "Frankia datiscae (2)",
#                      "Frankia sp. R43 (3)");
#     
     
}


#
# CONFIGURE ABOVE
#


#
# main data structure: read_data->dir->taxon_name = reads;
#
my $read_data;


#
# main data structure: kmer_data->dir->taxon_name = kmers
#
my $kmer_data;


#
# totals of kmers ACROSS nodules for each taxon: taxon_read_totals{taxon} = kmers
#
my %taxon_kmer_totals;


#
# totals of reads ACROSS nodules for each taxon: taxon_read_totals{taxon} = reads
#
my %taxon_read_totals;


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
my %classified_reads;  # classified_reads{sample} = reads
my %classified_kmers;  # classified_kmers{sample} = kmers


#
# gather all the data in useful form
#
DIR: for my $dir ( @dir_list )
{
   chdir($data_dir);
   chdir($dir);
   open (IN,"<$infile") || die "Unable to open $infile\n";
      
   #print "$dir\n";

   LINE: while(my $line = <IN>)
   {
      chomp $line;
      
      # common to the whole sample      
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
# ready to plot
#

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
       $pattern,
       $hack,
       \@plot_taxa,
       \@label_taxa
       );


#
# clean up
#
      
        

#
# subs below 
#

       
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
       $s_pattern,
       $s_hack,
       $s_plot_hack_taxa,
       $s_label_hack_taxa) = @_;
  
   # list taxa retained after filtering out the false positives
   my @retained_taxa;
   
   # taxa for plotting and labeling
   my @s_plot_taxa = ();
   my @s_label_taxa = ();
   
   # Apply kmer threshold per million reads
   while( my ($dir,$dir_data) = each %$s_read_data)
   {
      # get millions of classified reads * kmer threshold. That is the bar each taxon needs to clear. It can be zero.
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
         # fails test, set values to zero
         else
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
   
   #
   # taxa retained after filtering falsies, in reverse order by read totals
   #
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
        
   #
   # If a pattern is suggested, select for it or take the top N
   #  
   #print STDERR "pattern: >$s_pattern<\n";

   if($s_pattern ne 0)
   {
      #print STDERR"$s_pattern must not be 0\n";
      
      @s_plot_taxa = grep(/$pattern/, @ordered_retained_taxa);
      @s_plot_taxa = @s_plot_taxa[0 .. ($s_num_taxa - 1)];
   }
   else
   {
      @s_plot_taxa = @ordered_retained_taxa[0 .. ($s_num_taxa - 1)];
      
      #print STDERR "$s_pattern must be 0\n";
   }
  
   @s_label_taxa = @s_plot_taxa;
  

   #
   # If a hack of the taxon names is authorized, reassign plot taxa and label taxa
   #
   if ($hack != 0)
   {
      @s_plot_taxa = @$s_plot_hack_taxa;
      @s_label_taxa = @$s_label_hack_taxa;
   } 

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
      for my $taxon (@s_plot_taxa)
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
      for my $taxon (@s_plot_taxa)
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
   for my $taxon (@s_plot_taxa)
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
   $graph->set_legend( @s_label_taxa );
        
   my @data = (\@nod_list,@$all_plot_percents);
      
   my $splotfile_path = $s_data_dir . "/" . $outid . ".png";
   
   print "\$splotfile_path: $splotfile_path\n";
  
   open  OUT,">$splotfile_path" || die "Unable to open $splotfile_path\n";
   print OUT $graph->plot(\@data)->png;
   close OUT;  
}


