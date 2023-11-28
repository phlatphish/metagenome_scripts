# Scripts for tabulating and plotting KrakenUniq read and kmer counts 

These scripts are highly specific for KrakenUniq and Bracken reports and output files. All code is written in perl. Modules used are:

- GD::Graph::hbars
- List::Util
- FileHandle

They have been only been tested on perl 5.32.1 running on CentOS 7. The ability to read and write perl is an advantage. Some editing of configurable variables is done in the code where indicated by CONFIGURE.

## The scripts:

### filterizer.pl

Eliminates rows from a KrakenUniq file falling below a specific kmer threshold. Typically this is 2000 kmers per million reads.

Usage: filterizer.pl krakenuniq\_report threshold
Args:
- KrakenUniq report filename
- kmer hreshold (e.g. 2000)
Output: STDOUT

### krakenuniq\_tabulator.pl

Prints read an kmer count tables for all samples and all taxa. Applies a kmer threshold if specified.

Usage: krakenuniq\_tabulator.pl data\_dir {species|genus} infile outid threshold [pattern]
Args
- data\_dir (path to dir containing sequence files)
- taxon level (genus or species)
- input file (KrakenUniq report)
- outid (string used to name output files)
- threshold (minimum kmers per million reads)
- pattern (only treats rows matching pattern, e.g. genus name)
Output: read and kmer count tables in CSV format

### bracken\_tabulator.pl

Usage: bracken\_tabulator.pl data_dir {S|G} infile outid [pattern]\n"; 
Args:
- data\_dir (path to dir containing sequence files)
- taxon level (G or S)
- input file (Bracken report)
- outid (string used to name output files)
- pattern (only treats rows matching pattern, e.g. genus name)
Output: read and kmer count tables in CSV format

### krakenuniq\_plotter\_custom\_color.pl

Create a horizontal bar plot for taxa in a krakenuniq report passing a kmer threshold

Usage: krakenuniq\_plotter\_custom\_color.pl data\_dir {G|S} infile outid title threshold [pattern]
Args:
- data\_dir (path to dir containing sequence files)
- taxon level (G or S)
- input file (KrakenUniq report)
- outid (string used to name output files)
- title (literal title for the output graph, in quotes)
- threshold (minimum kmers per million reads)
- pattern (only treats rows matching pattern, e.g. genus name)
Output: PNG file

### bracken\_plotter\_custom\_color.pl

Create a horizontal bar plot for taxa in a bracken report

Usage: bracken\_plotter\_custom\_color.pl data\_dir {G|S} infile outid title [pattern]
Args:
- data\_dir (path to dir containing sequence files)
- taxon level (G or S)
- input file (KrakenUniq report)
- outid (string used to name output files)
- title (literal title for the output graph, in quotes)
- pattern (only treats rows matching pattern, e.g. genus name)
Output: PNG file



