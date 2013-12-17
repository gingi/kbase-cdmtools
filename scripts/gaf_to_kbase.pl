#!/usr/bin/env perl

use strict;
use warnings;

use feature 'state';
use Getopt::Long;
use Pod::Usage;
use Readonly;

Readonly my $ONTOLOGY_TYPE     => "GO";
Readonly my $EMPTYSTR          => q{};

MAIN: {
    my %options = (output => "ontology.tab");
    GetOptions(\%options, 'output|o=s', 'help|h', 'man|m')
        || pod2usage(-verbose => 1);

    pod2usage(-verbose => 1) if $options{help};
    pod2usage(-verbose => 2) if $options{man};

    my ($infh, $outfh);
    my $filename = $ARGV[0];

    if ($filename) {
        open($infh, "<", $filename)
            or die "Cannot read from input file $filename: $!\n";
    } else {
        $infh = \*STDIN;
    }
    open($outfh, ">", $options{output})
        or
        die "Cannot open ontology file $options{output} for writing: $!\n";
    print $outfh join("\t",
        qw/GeneID TranscriptID OntologyID OntologyDescription OntologyDomain
        OntologyEvidenceCode OntologyType/), "\n";
    while (my $line = <$infh>) {
        chomp $line;
		next if ($line =~ m/^\s*\!/);
        my ($gene, $term, $evidence_code)
            = (split(/\t/, $line, 11))[1,2,4];
        print $outfh join("\t", $gene, $EMPTYSTR, $term, $EMPTYSTR, $EMPTYSTR,
            $evidence_code, $ONTOLOGY_TYPE), "\n";
    }
    close($outfh);
}

1;

=head1 NAME

gaf_to_kbase.pl

=head1 SYNOPSIS

gaf_to_kbase.pl [options] ontology.gaf

Options:

    --output, -o : Output ontology file (ontology.tab)
    --help,   -h : Brief usage instructions
    --man,    -m : More verbose usage instructions

=head1 OPTIONS

B<--help>, B<-h>
    Brief usage instructions.

B<--man>, B<-m>
    More verbose usage instructions.   

=head1 DESCRIPTION

Convert GAF data into a KBase-formatted ontology.tab file

=head1 AUTHOR

Shiran Pasternak <shiranpasternak@gmail.com>

=cut    
