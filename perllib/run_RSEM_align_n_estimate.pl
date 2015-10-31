#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<__EOUSAGE__;

#########################################################################
#
#  --transcripts <string>           transcript fasta file
#  --seqType <string>              fq|fa
# 
#  If Paired-end:
#
#  --left <string>
#  --right <string>
#  
#    or Single-end:
#
#  --single <string>
#
#
#
# Optional:
# 
# --prefix <string>                prefix for RSEM output files (default: 'RSEM')
#
# --SS_lib_type <string>           strand-specific library type:  paired('RF' or 'FR'), single('F' or 'R').
#
# --no_group_by_component          Trinity-mode, using 'components' as 'genes'
#
# --thread_count                   number of threads to use (default = 4)
#
# --debug                  retain intermediate files
#  
#####################
#  Non-Trinity options:
#
#  --gene_trans_map <string>        file containing 'gene(tab)transcript' identifiers per line.
#
#
#########################################################################
#  
#  To pass additional parameters to rsem-calculate-expression, 
#    type ' -- ' followed by additional pass-through params
#
#########################################################################




__EOUSAGE__

    ;


my $help_flag;
my $transcripts;
my $bam_file;
my $paired_flag;
my $DEBUG_flag = 0;
my $SS_lib_type;
my $no_group_by_component = 0;
my $thread_count = 4;
my $seqType;
my $left;
my $right;
my $single;
my $gene_trans_map_file;
my $prefix = "RSEM";

&GetOptions ( 'h' => \$help_flag,
              'transcripts=s' => \$transcripts,
              'name_sorted_bam=s' => \$bam_file,
              'paired' => \$paired_flag,
              'debug' => \$DEBUG_flag,
              'SS_lib_type=s' => \$SS_lib_type,
              'no_group_by_component' => \$no_group_by_component,
              'thread_count=i' => \$thread_count,
              'gene_trans_map=s' => \$gene_trans_map_file,

              'seqType=s' => \$seqType,
              'left=s' => \$left,
              'right=s' => \$right,
              'single=s' => \$single,
              
              'prefix=s' => \$prefix,
              
              );



if ($help_flag) {
    die $usage;
}

unless ($transcripts && $seqType && ($single || ($left && $right))) {
    die $usage;
}




if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(RF|FR|R|F)$/) {
        die "Error, do not recognize SS_lib_type: [$SS_lib_type]\n";
    }
    if ($left && $right && length($SS_lib_type) != 2 ) {
        die "Error, SS_lib_type [$SS_lib_type] is not compatible with paired reads";
    }
}

if ( $thread_count !~ /^\d+$/ ) {
    die "Error, --thread_count value must be an integer";
}

my $RSEM_dir = "$FindBin::Bin/../plugins/rsem";

main: {

    my $cmd = "$RSEM_dir/rsem-prepare-reference";
    
    unless (-s "TRANS.1.ebwt") { ## this step already run

        if ($gene_trans_map_file) {
            $cmd .= " --transcript-to-gene-map $gene_trans_map_file ";
        }
        elsif (! $no_group_by_component) {
            
            # create Trinity component-to-transcript mapping.
            my $trans_to_gene_map_file = &write_gene_to_trans_map_file($transcripts);
            
            $cmd .= " --transcript-to-gene-map $trans_to_gene_map_file";
        
        }
        $cmd .= " $transcripts TRANS";
        
        &process_cmd($cmd);
    }
    
    
    my $keep_intermediate_files_opt = ($DEBUG_flag) ? "--keep-intermediate-files" : "";


    my $SS_opt = "";
    if ($SS_lib_type) {
        if ($SS_lib_type =~ /^F/) {
            $SS_opt = "--forward-prob 1.0";
        }
        else {
            $SS_opt = "--forward-prob 0";
        }
    }
    
    $cmd = "$RSEM_dir/rsem-calculate-expression @ARGV " ## allow for custom params
        . " -p $thread_count"
        . " $keep_intermediate_files_opt"
        . " $SS_opt";
    
    if ($seqType eq "fa")  {
        $cmd .= " --no-qualities";
    }

    if ($left && $right) {
        $cmd .= " --paired-end $left $right";
    }
    else {
        $cmd .= " $single";
    }
            
    $cmd .= " TRANS $prefix";
    
    &process_cmd($cmd);
   
    exit(0);
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "**********************************************************************\n";
    print STDERR "**  Running RSEM Command:\n";
    print STDERR "**  $cmd\n";
    print STDERR "**********************************************************************\n\n\n";

    
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret: $ret";
    }
    
    return;
}


####
sub write_gene_to_trans_map_file {
    my ($transcripts_fasta_file) = @_;    
        
    open (my $fh, $transcripts_fasta_file) or die "Error, cannot open file $transcripts_fasta_file";
    
    my $mapping_file = "$transcripts_fasta_file.component_to_trans_map";
    open (my $ofh, ">$mapping_file") or die "Error, cannot write to file: $mapping_file";
    
    while (<$fh>) {
        if (/>(comp\S+)/) {
            my $acc = $1;
            $acc =~ /^(comp\d+)_seq\d+/ or die "Error, cannot parse the trinity component ID from $acc";
            my $comp_id = $1;
            print $ofh "$comp_id\t$acc\n";
        }
    }
    close $fh;
    close $ofh;

    return($mapping_file);
}

