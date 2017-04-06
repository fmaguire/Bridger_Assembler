#!/usr/bin/env perl
use strict;

use threads;
use FindBin;
use File::Basename;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $PACKAGE = "Bridger";
my $VERSION = "2014-12-01";
my $ROOTDIR = "$FindBin::Bin";
my $SRC_DIR = "$ROOTDIR/build/src";
my $PERL_DIR = "$ROOTDIR/perllib";
my $FASTOOL_DIR = "$ROOTDIR/plugins/fastool";
my $BIN_DIR = "$ROOTDIR/bin";

# defaults:
my $kmer_length = 25;
my $CPU = 2;
my $no_run_pathsearch = 0;
my $output_directory = "bridger_out_dir";

my $Cmd_line = "A typical command might be:\n\nperl $PACKAGE.pl --seqType fq --left reads1.fq --right reads2.fq --SS_lib_type RF --CPU 6\n\n"
               . "Use option --help for more information.\n";


my $usage = <<_EOUSAGE_;


#######################################################################################
#			 
# $PACKAGE : An Efficient De novo Transcriptome Assembler For RNA-Seq Data
#
# version : $VERSION
#
## USAGE ##
#		   
# ** Required **
#
#  --seqType <string>  : type of reads: (fa, fq, cfa, cfq)
#
#  If paired reads:
#     --left  <string>  : left reads
#     --right <string>  : right reads
#
#  (Compressed files ending with .gz are able to be read directly, it 
#    will save your time to uncompress them.)
#
#  Or, if unpaired reads:
#     --single <string>  : single reads
#
# 
# ** Optional **
#
#  if strand-specific data, set:
#  --SS_lib_type <string>        : Strand-specific RNA-Seq reads orientation.
#				   if paired: RF or FR,  if single: F or R.
#                                  (dUTP method = RF)
#  --kmer_length/-k <int>        : length of kmer, default: $kmer_length. 
#  --output/-o <string>          : name of directory for output, default: $output_directory.
#  --CPU <int>                   : number of CPUs for PathSearch, default: $CPU.
#  --pair_gap_length             : gap length of paired reads, default: 200.
#  --min_seed_coverage <int>     : minimum coverage of kmer as a seed, default: 2 .
#  --min_seed_entropy <float>    : minimum entropy of kmer as a seed, default: 1.5.
#  --min_kmer_coverage <int>     : minimum coverage of kmer used to extend, default: 1.
#  --min_kmer_entropy <float>    : minimum entroy of kmer used to extend, default: 0.0.
#  --min_junction_coverage <int> : minimum of the coverage of a junction, default: 2.
#  --min_ratio_non_error <float> : min ratio for low/high alternative extension that is 
#                                  not an error, default: 0.05.
#  --min_reads_span_junction<int>: minimum number of reads supporting a junction.
#  --clean                       : clear all intermediate files
#  --debug                       : display more information for debugging
#  --version                     : report current version and exit.    
#  --help/-h                     : show help information.
#  
# 
# ** Note **
#
#  A typical command might be:
#    perl $PACKAGE.pl --seqType fq --left reads1.fq --right reads2.fq --CPU 6
#  (If your data are strand-strand, it is recommended to set --SS_lib_type option.)
#  
#  For more details, visit : http://to_be_done.net 
#
###########################################################################################

_EOUSAGE_
	;


unless (@ARGV) {
    die "$usage\n";
}


# option list:
my ($seqType, $left_file, $right_file, $single_file, $SS_lib_type, 
    $min_seed_coverage, $min_seed_entropy, $min_kmer_coverage, $min_kmer_entropy, 
    $min_junction_coverage, $min_ratio_non_error, $double_strand, $pair_gap_length,
    $fr_strand, $show_citation, $show_version, $show_help, $debug, $clean);

&GetOptions( 
    ## required
    "seqType=s" => \$seqType,
    "left|l=s" => \$left_file,
    "right|r=s" => \$right_file,
    "single=s" => \$single_file,
    "SS_lib_type=s" => \$SS_lib_type,
    # optional 
    "kmer_length|k=i" => \$kmer_length,
    "output|o=s" => \$output_directory,
    "CPU=i" => \$CPU,
    "min_seed_coverage=i" => \$min_seed_coverage,
    "min_seed_entropy=f" => \$min_seed_entropy,
    "min_kmer_coverage=i" => \$min_kmer_coverage,
    "min_kmer_entropy=f" => \$min_kmer_entropy,
    "min_junction_coverage=i" => \$min_junction_coverage,
    "min_ratio_non_error=f" => \$min_ratio_non_error,
    "pair_gap_length|g=i" => \$pair_gap_length,
    "no_run_pathsearch"=> \$no_run_pathsearch,
    "clean" => \$clean,
    "debug" => \$debug,
    "cite" => \$show_citation,
    "version|v" => \$show_version,
    "help|h" => \$show_help,
    );


if (@ARGV) {  # if still have other options
    die "[Error] unknown option $ARGV[0]!\n$usage\n";
}

if ($show_help) {die "$usage\n";}

if ($show_version) {
    print "current version : $VERSION\n";
    exit;
}

if ($show_citation) {
    &show_lit_citation();
    exit;
}

## Check options set:
if ($kmer_length < 19) {
    die "[Error] Length of kmer is too small, range [19,32] is acceptable!\n";
} elsif ($kmer_length > 32) {
    die "[Error] Length of kmer is too big, range [19,32] is acceptable!\n";
}


if ($seqType) {
    unless ($seqType =~ /^(fa|fq|cfa|cfq)$/) {
        die "[Error] option '--seqType' is not one of (fa, fq, cfa, cfq).\n";
    }
} else {
    die "[Error] '--seqType' is a required option.\n\n$Cmd_line\n";
}

if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(R|F|RF|FR)$/) {
	die "[Error] unrecongnized value of '--SS_lib_type', set it RF or FR if paired, F or R if single\n";
    }
    if ($SS_lib_type =~ /^RF$/) {
        $fr_strand = 1;
    } elsif ($SS_lib_type =~ /^FR$/) {
        $fr_strand = 2;
    }
}

unless (($left_file && $right_file) || $single_file ) {
    die "[Error] need either options '--left' and '--right' or option '--single'.\n\n$Cmd_line\n";
}

if ($single_file && ($left_file || $right_file)) {
    die "[Error] option '--single' can not be setted simultaneously with option '--left' or '--right'.\n\n";
}

unless ($output_directory) {
    $output_directory = "bridger_out_dir";
}



main: {

    my $start_dir = cwd();

    ## create complete paths for input files:
    $left_file = &create_full_path($left_file) if $left_file;
    $right_file = &create_full_path($right_file) if $right_file;
    $single_file = &create_full_path($single_file) if $single_file;
    $output_directory = &create_full_path($output_directory);	

    # create output directory
    if (-d $output_directory) {
	# clear
        print "[Warning] $output_directory exists already. It will be overwritten.\n";
        &process_cmd("rm -r $output_directory/*");
        
    } else {
	mkdir $output_directory or die "Error, cannot mkdir $output_directory!\n";
    } 
	
    chdir ($output_directory) or die "Error, cannot cd to $output_directory!\n";
	
    my $target_fa = ($single_file) ? "single.fa" : "both.fa";
	
    if ($left_file && $right_file) { 	
	#if ($run_ALLPATHSLG_error_correction_flag) {
	#    &process_cmd("$ROOTDIR/util/run_AllpathsLG_error_correction.pl $left_file $right_file");
	#    $left_file = "$left_file.ErrCor.fq";
	#    $right_file = "$right_file.ErrCor.fq";
	#}
        unless (-s "both.fa") {
 	    my ($left_SS_type, $right_SS_type);
	    if ($SS_lib_type) {
	        ($left_SS_type, $right_SS_type) = split(//, $SS_lib_type);
    	    }
		
            print "\nConverting input files... (in parallel) \n";
            my $thr1;
            my $thr2;
            if (!(-s "left.fa")) {
                $thr1 = threads->create('prep_seqs', $left_file, $seqType, "left", $left_SS_type);
            } else {
                $thr1 = threads->create(sub { print ("left file exists, nothing to do\n"); });
            }
            if (!(-s "right.fa")) {
                $thr2 = threads->create('prep_seqs', $right_file, $seqType, "right", $right_SS_type);
            } else {
                $thr2 = threads->create(sub { print ("right file exists, nothing to do\n"); });
            }
            $thr1->join();
            $thr2->join();
            print "\nDone converting input files. \n\n" ;
	    &process_cmd("cat left.fa right.fa > $target_fa") unless (-s $target_fa && (-s $target_fa == ((-s "left.fa") + (-s "right.fa"))) );		
  	    unlink ("left.fa", "right.fa");  # delete these two file, no longer need them
	}
    } elsif ($single_file) {
	#if ($run_ALLPATHSLG_error_correction_flag) {
	#    &process_cmd("$ROOTDIR/util/run_AllpathsLG_error_correction.pl $single_file");
	#    $single_file = "$single_file.ErrCor.fq";
	#}
        
        &prep_seqs($single_file, $seqType, "single", $SS_lib_type) unless (-s "single.fa");
    } else {
 	die "not sure what to do. "; # should never get here.
    }
	
    ##================
    ## Assemble step:
    print "\n### Splicing Graphs Reconstruction ###\n\n";

    my $assemble_cmd = "$SRC_DIR/Assemble --reads $target_fa -k $kmer_length ";
    #$assemble_cmd .= " --CPU $CPU" if ($CPU);
    $assemble_cmd .= " --pair_end" unless ($target_fa eq "single.fa");
    $assemble_cmd .= " --double_stranded_mode" unless ($SS_lib_type);
    $assemble_cmd .= " --fr_strand $fr_strand" if ($fr_strand); 
    $assemble_cmd .= " --min_seed_coverage $min_seed_coverage" if ($min_seed_coverage);
    $assemble_cmd .= " --min_seed_entropy $min_seed_entropy" if ($min_seed_entropy);
    $assemble_cmd .= " --min_kmer_coverage $min_kmer_coverage" if ($min_kmer_coverage);
    $assemble_cmd .= " --min_kmer_entropy $min_kmer_entropy" if ($min_kmer_entropy);
    $assemble_cmd .= " --min_junction_coverage $min_junction_coverage" if ($min_junction_coverage);
    $assemble_cmd .= " --pair_gap_length $pair_gap_length" if ($pair_gap_length && ($target_fa eq "both.fa"));
    $assemble_cmd .= " --debug" if ($debug);
    $assemble_cmd .= " 2>Assemble.log";
    #$assemble_cmd .=  " 2>&1 | tee Assemble.log";
    &process_cmd($assemble_cmd);
    
    ##=================
    ## prepare the list of file
    my @raw_graphs;
    
    open(Graph,"$output_directory/RawGraphs/raw_graph.list") || die "Can not open directory/RawGraphs/raw_graph.list!\n";
    @raw_graphs = <Graph>;
    chomp(@raw_graphs);
    close Graph;

    open(Cmd_list, ">$output_directory/path_search.commands") || die "can not open $output_directory/path_search.commands!\n";
    foreach my $raw_graph (@raw_graphs) {
        my $cmd = "$SRC_DIR/PathSearch -i $raw_graph";
        $cmd .= " -k $kmer_length";
        $cmd .= " --pair_end" unless ($target_fa eq "single.fa");
        $cmd .= " --double_stranded_mode" unless ($SS_lib_type);
        $cmd .= " --pair_gap_length $pair_gap_length" if ($pair_gap_length); 
        $cmd .= " --debug" if ($debug);
	$cmd .= "  >$output_directory/transcripts/$raw_graph.transcripts.fasta";
	print Cmd_list "$cmd\n";
    }
    close Cmd_list;

    ##===================
    ## Path search step:
    exit if ($no_run_pathsearch);
    print "\n### Search paths from Splicing Graphs ###\n\n";
    

    unless (-s "$output_directory/path_search.commands") {
	die "No commands for path search, Could it be that your read data is too sparse to generate minimal length contigs?\n";
    }

    unless (-d "$output_directory/transcripts") {
        mkdir "$output_directory/transcripts" || die "$!\n" ;
    }

    my $path_search_cmds = "$PERL_DIR/cmd_process_forker.pl -c $output_directory/path_search.commands --CPU $CPU --shuffle";
    &process_cmd($path_search_cmds);

    ##======================
    ## collect all transcripts
    # no longer scan the file system... we know which files should exist
    #my $cmd = 'find ./transcripts -regex ".*transcripts.fasta" -exec cat {} \; Bridger.fasta';
    &collect_transcripts(\@raw_graphs, 'Bridger.fasta');

    print "\n\n";
    print "############################################################################\n";
    print "All transcripts are written to :\n $output_directory/Bridger.fasta      \n";
    print "############################################################################\n\n";


    if ($clean) {
        &clean(\@raw_graphs);
    }
    system("rm -f $target_fa");

    exit(0);	
}



### prepare sequence
sub prep_seqs {

    my ($initial_file, $seqType, $file_prefix, $SS_lib_type) = @_;

    if ($seqType eq "fq") {

	# make fasta
        my $perlcmd = "$PERL_DIR/fastQ_to_fastA.pl -I $initial_file";
        my $fastool_cmd = "$FASTOOL_DIR/fastool";
        if ($SS_lib_type && $SS_lib_type eq "R") {
            $perlcmd .= " --rev";
            $fastool_cmd .= " --rev";
        }
        $fastool_cmd .= " --to-fasta $initial_file > $file_prefix.fa";
        $perlcmd .= " > $file_prefix.fa";  
	&process_cmd($fastool_cmd) unless (-e "$file_prefix.fa");
    } elsif ($seqType eq "fa") {

	if ($SS_lib_type && $SS_lib_type eq "R") {
  	    my $cmd = "$PERL_DIR/revcomp_fasta.pl $initial_file > $file_prefix.fa";
	    &process_cmd($cmd) unless (-s "$file_prefix.fa");
	} else {
	    ## just symlink it here:
	    my $cmd = "ln -s $initial_file $file_prefix.fa";
	    &process_cmd($cmd) unless (-s "$file_prefix.fa");
	}
    } elsif ($seqType eq "cfa" | $seqType eq "cfq"){

	#make double-encoded fasta
	my $cmd = "$PERL_DIR/csfastX_to_defastA.pl -I $initial_file";
        if ($SS_lib_type && $SS_lib_type eq "R") {
                $cmd .= " --rev ";
        }
        $cmd .= "> $file_prefix.fa";
        &process_cmd($cmd) unless (-e "$file_prefix.fa");
    } else{

	print "Error, illegal argument for '--seqType'!\n";
	exit(1);
    }	

    return;
}


###
sub create_full_path {
    my ($file) = @_;
    my $cwd = cwd();
    if ($file !~ m|^/|) { # must be a relative path
      $file = $cwd . "/$file";
    }
    return ($file);
}

####
sub process_cmd {
    my ($cmd) = @_;
    print "CMD: $cmd\n";
    my $start_time = `date +%s`;
    my $ret = system($cmd);
    my $end_time = `date +%s`;
    if ($ret) {
      die "Error, cmd: $cmd died with ret $ret !\n";
    }
    print "CMD finished (" . ($end_time - $start_time) . " seconds)\n";
    return;
}

###
sub show_lit_citation {

    print "\n";
    print "############################################################\n";
    print "  Tools or Codes Used within This Softwares\n";
    print "############################################################\n\n";

    print "* Fastool (for fast fastQ-to-fastA conversion):\n"
         ."Francesco Strozzi\n"
         ."Code: https://github.com/fstrozzi/Fastool\n\n";

    print "* Trinity\n"
        . "Full-length transcriptome assembly from RNA-Seq data without a reference genome.\n"
        . "Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L,\n"
        . "Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F,\n"
        . "Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A.\n"
        . "Nature Biotechnology 29, 644–652 (2011)\n"
        . "Paper: http://www.nature.com/nbt/journal/v29/n7/full/nbt.1883.html\n"
        . "Code:  http://trinityrnaseq.sf.net\n\n";

    print "* Cufflinks\n"
        . "Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and \n"
        . "isoform switching during cell differentiation.\n"
        . "Cole Trapnell, Brian A Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Marijke J van Baren,\n"
        . "Steven L Salzberg, Barbara J Wold & Lior Pachter\n"
        . "Nature Biotechnology 28, 511–515 (2010)\n"
        . "Paper: http://www.nature.com/nbt/journal/v28/n5/full/nbt.1621.html\n"
        . "Code: http://cufflinks.cbcb.umd.edu/\n\n";

    print "################################################################\n\n";
}


###
sub collect_transcripts {
 
    my @rgs = @{$_[0]};  
    my $output_file = $_[1];

    my $begin = `date +%s`;
    print "\n### Collecting all transcripts ###\n\n...\n\n";  

    open (FASTA, ">$output_file") || die "cannot open $output_file.\n";

    foreach (@rgs) {

        my $rg_basename = $_;
        my $fasta_file = "./transcripts/$rg_basename.transcripts.fasta";

        if (-e $fasta_file) {
            open (my $fh, $fasta_file) or die "Error, cannot open file $fasta_file";
            while (<$fh>) {
                print FASTA $_;
            }
            close $fh;
        } else {
            print STDERR "Error, no fasta file reported as: $fasta_file\n";
        }
    }

    close FASTA;

    my $end = `date +%s`;
    print "Done! (", $end-$begin, " seconds)\n\n";
 
    return ;
}

###
sub clean {

    my @rgs = @{$_[0]};

    my $begin = `date +%s`;
    print "\n### Cleaning all intermediate files ###\n\n...\n\n";
    foreach (@rgs) {
        my $delete_rg = "RawGraphs/$_.rg";
        my $delete_ta = "transcripts/$_.transcripts.fasta";
        system("rm -rf $delete_rg") if (-e $delete_rg);
        system("rm -rf $delete_ta") if (-e $delete_ta);
    }
    system("rm -rf RawGraphs") if (-d "RawGraphs");
    system("rm -rf transcripts") if (-d "transcripts");

    my $end = `date +%s`;
    print "Done! (", $end-$begin, " seconds)\n\n";

    return;
}



