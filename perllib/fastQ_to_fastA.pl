#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
#use lib ("$FindBin::Bin/../PerlLib");
use lib("$FindBin::Bin/modules");
use Nuc_translator;

use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

##########################################################
#
#  -I <string>     input.fq   
#
#  --ignoreDirty   ignores poorly formed entries
#
#  -a <int>        append "/num" to the accession name.
#
#  --rev           reverse complement nucleotide sequence.
#
###########################################################

_EOUSAGE_

	;

my $inputFile;
my $ignore_dirty = 0;
my $append_num;
my $revcomp_flag = 0;

&GetOptions( 'I=s' => \$inputFile,
		 'ignore_dirty' => \$ignore_dirty,
		 'a=i' => \$append_num,
		 'rev' => \$revcomp_flag,
	);


unless ($inputFile) {
	die $usage;
}



open (my $fh, $inputFile) or die "Error, cannot open $inputFile";

my $counter = 0;
my $num_clean = 0;
my $num_dirty = 0;

my @rec;

my $line = <$fh>;

while ($line) {

	if ($line =~ /^\@/) {
		$counter++;
		
		# print STDERR "\r[$counter] [$num_clean clean] [$num_dirty dirty]       " if ($counter % 10000 == 0);
		
		push (@rec, $line);
		
		$line = <$fh>;
		for (1..3) {
			push (@rec, $line);
			$line = <$fh>;
		}
		
		my $record_text = join("", @rec);

		my $header = shift @rec;
		my $seq = shift @rec;
		my $qual_header = shift @rec;
		my $qual_line = shift @rec;
				
		chomp $header;
		chomp $seq if $seq;
		chomp $qual_header if $qual_header;
		chomp $qual_line if $qual_line;
		
		if ($header && $seq && $qual_header && $qual_line && 
			$qual_header =~ /^\+/ && length($seq) == length($qual_line)) {
			
			# can do some more checks here if needed to be sure that the lines are formatted as expected.

			$header =~ s/^\@//;
			

			if (defined $append_num) {
				$header .= "/$append_num";
			}
			
			if ($revcomp_flag) {
				$seq = &reverse_complement($seq);
			}
			

			print ">$header\n$seq\n";
			$num_clean++;
		}
		else {

			$num_dirty++;

			unless ($ignore_dirty) {
				die "Error, improperly formatted entry:\n\n$record_text  ";
			}
			
		}
		@rec = ();
	} else {
		$line = <$fh>;
	}
	
}

exit(0);


		
		
		
		
