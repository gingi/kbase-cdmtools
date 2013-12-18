#!/usr/bin/env perl

use strict;
use warnings;

use feature 'state';
use Getopt::Long;
use Pod::Usage;
use Readonly;

Readonly my $ONTOLOGY_TYPE     => "GO";
Readonly my $ONTOLOGY_EVIDENCE => "IEA";
Readonly my $EMPTYSTR          => q{};

MAIN: {
    my %options = (
        ontology  => "ontology.tab",
        functions => "functions.tab"
    );
    GetOptions(\%options, 'ontology|O=s', 'functions|F=s', 'help|h', 'man|m')
        || pod2usage(-verbose => 1);

    pod2usage(-verbose => 1) if $options{help};
    pod2usage(-verbose => 2) if $options{man};

    my ($infh, $onfh, $fnfh);
    my $filename = $ARGV[0];

    if ($filename) {
        open($infh, "<", $filename)
            or die "Cannot read from input file $filename: $!\n";
    } else {
        $infh = \*STDIN;
    }
    open($onfh, ">", $options{ontology})
        or
        die "Cannot open ontology file $options{ontology} for writing: $!\n";
    open($fnfh, ">", $options{functions})
        or die
        "Cannot open functions file $options{functions} for writing: $!\n";
    print $onfh join("\t",
        qw/GeneID TranscriptID OntologyID OntologyDescription OntologyDomain
        OntologyEvidenceCode OntologyType/), "\n";
    while (my $line = <$infh>) {
        chomp $line;
        my ($gene, $trans, $cds, $pfam, $panther, $kog, $k_ec, $k_orth, $go)
            = (split(/\t/, $line, 11))[1 .. 9];
        for my $term (split(',', $go)) {
            print $onfh join("\t", $gene, $trans, $term, $EMPTYSTR, $EMPTYSTR,
                $ONTOLOGY_EVIDENCE, $ONTOLOGY_TYPE), "\n";
        }
        for my $fngrp ($pfam, $panther, $kog, $k_ec, $k_orth) {
            for my $fn (split(',', $fngrp)) {
                if (defined $fn && $fn ne "") {
                    print $fnfh join("\t", $trans, $fn), "\n";
                }
            }
        }
    }
    close($onfh);
    close($fnfh);
}

1;

=head1 NAME

convert_jgi_info.pl

=head1 SYNOPSIS

convert_jgi_info.pl [options] annotation_info_file.txt

Options:

    --ontology,  -O : Output ontology file (ontology.tab)
    --functions, -F : The functions file (functions.tab)
    --help,      -h : Brief usage instructions
    --man,       -m : More verbose usage instructions

=head1 OPTIONS

B<--help>, B<-h>
    Brief usage instructions.

B<--man>, B<-m>
    More verbose usage instructions.   

=head1 DESCRIPTION

Convert JGI annotation info into requisite KBase Exchange files

=head1 AUTHOR

Shiran Pasternak <shiranpasternak@gmail.com>

=cut    
