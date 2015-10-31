#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/modules");

use Nuc_translator;
use Fasta_reader;

my $usage = "\nUsage: $0 fasta_file\n";

if (@ARGV == 0) {
 die "$usage\n";
}

my $fasta_file = $ARGV[0] || *STDIN{IO};

main: {
	
	my $fasta_reader = new Fasta_reader($fasta_file);
	
	while (my $seq_obj = $fasta_reader->next()) {

		my $header = $seq_obj->get_header();
		my $sequence = $seq_obj->get_sequence();

		$sequence = &reverse_complement($sequence);

		print ">$header (reverse-complemented)\n";

		#$sequence =~ s/(\S{60})/$1\n/g;

		#chomp $sequence;
		print "$sequence\n";
	}

	exit(0);

}

