#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;


## we delete all files we don't need.

chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";


my @files_to_keep = qw (clean.pl
                        README
                        reads.left.fq.gz
                        reads.right.fq.gz
                        run_Me.sh
                        run_Me_as_DS.sh
                        run_abundance_estimation.sh                       
                        );


my %keep = map { + $_ => 1 } @files_to_keep;


`rm -rf bridger_out_dir/` if (-d "bridger_out_dir");
`rm -rf RSEM.stat` if (-d "RSEM.stat");


foreach my $file (<*>) {

        if (! $keep{$file}) {
                print STDERR "-removing file: $file\n";
                unlink($file);
        }
}

exit(0);

