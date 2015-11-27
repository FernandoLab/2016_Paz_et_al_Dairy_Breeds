#!/usr/bin/perl -w

#remove any seqs less than 100 bp
#anything longer than 467, truncaate at 467
#calc length distribution

use strict;
use Getopt::Long;

#Command line parameters:
my $fasta = "";
my $min = "";
my $max = "";

#Setup the command line options using Getopt:Long
my $commandline = GetOptions("fasta:s", \$fasta,
							 "min:s", \$min,
							 "max:s", \$max);

if (!$commandline || $fasta eq "" || $min eq "" || $max eq "") {
	print STDERR "Usage: $0 -fasta -min -max \n";
	print STDERR "example: ./header_lines_sff_split.pl -fasta=16S.fasta -min=100 -max=467 \n\n";
	exit;
}

open (my $FASTA_FILE, "$fasta") or die "Can't open FASTA file!";
open (my $NEW_FASTA, ">seqs_trimmed.fasta") or die "Can't open output FASTA file";

my $input_read_count = 0;
my $output_read_count = 0;
my $header;
my $index_max = ($max - 1);
my $sequence;

while (my $line = readline($FASTA_FILE)) {
	chomp $line;
	my $check_line = substr ($line, 0, 1);
	
	if ($check_line eq ">") {
		$header = $line;
		$input_read_count ++;
		next;
	}
	
	my $seq_length = length ($line);
    #print "$seq_length\n" ;
	
	if ($seq_length < $min) {
		next;
	}
	
	if ($seq_length > $max) {
		my $max_line = substr ($line, 0, $index_max);
		if ($output_read_count == 0) {
			print $NEW_FASTA "$header\n";
			print $NEW_FASTA "$max_line";
		}
		print $NEW_FASTA "\n$header\n";
		print $NEW_FASTA "$max_line";
		$output_read_count ++;
	} else {
		if ($output_read_count == 0) {
			print $NEW_FASTA "$header\n";
			print $NEW_FASTA "$line";
		}
		print $NEW_FASTA "\n$header\n";
		print $NEW_FASTA "$line";
		$output_read_count ++;
	}
}

print "\nInput read count: $input_read_count \n";
print "Output read count: $output_read_count \n\n";

close ($FASTA_FILE) or die "Can’t close the FASTA file!";
close ($NEW_FASTA) or die "Can’t close the NEW_FASTA file!";