#!/usr/bin/env perl

use strict;
use warnings;

use Readonly;
use Getopt::Long;
use Pod::Usage;

Readonly my $GRMURL => "http://gramene.org";

MAIN: {
    my %options = (output => "provenance.json");
    GetOptions(\%options, 'output|o=s', 'help|h', 'man|m')
        || pod2usage(-verbose => 1);

    pod2usage(-verbose => 1) if $options{help};
    pod2usage(-verbose => 2) if $options{man};

    my $filename = $ARGV[0];

	my $meta = load_metadata($filename);
    open(my $fh, ">", $options{output}) or
        die "Cannot open provenance file $options{output} for writing: $!\n";

	my $species = $meta->{name}; $species =~ s/ /_/g;
	my $species_url = "$GRMURL/$species";
	my $common = $meta->{"common-name"};
	my $version = $meta->{"assembly-version"};
	$version =~ s/^$species\.//;

	# Use first if multiple values
	$common = $common->[0] if (ref $common eq "ARRAY");
	
	print $fh <<EOJSON;
{
    "scope" : ["data build"],
    "name"  : "$meta->{name}",
    "version" : {
        "external" : "$version"
    },
    "origin_url": "$species_url",
    "format": "FASTA",
    "type": "NUCLEOTIDE",
    "content": "ORGANISM",
    "source": {
        "name": "$meta->{source}",
        "url": "$GRMURL"
    },
    "description": "Genome data for $meta->{name} ($common)",
    "citation": "$species_url"
}
EOJSON
	close($fh);
}

sub load_metadata {
	my ($filename) = @_;
	my ($fh);
	my $meta = {};
    if ($filename) {
        open($fh, "<", $filename)
            or die "Cannot read from input file $filename: $!\n";
    } else {
        $fh = \*STDIN;
    }
	my $status = "key";
	my $key    = undef;
	while (my $line = <$fh>) {
		chomp $line;
		if ($line =~ m|^//|) {
			$status = "key";
		} elsif ($status eq "key") {
			$key = $line;
			$meta->{$key} = undef;
			$status = "value";
		} elsif ($status eq "value") {
			my $value = $meta->{$key};
			if (not defined $value) {
				$meta->{$key} = $line;
			} elsif (ref $value eq "ARRAY") {
				push @{ $meta->{$key} }, $line;
			} elsif (ref \$value eq "SCALAR") {
				$meta->{$key} = [ $value, $line ];
			}
		}
	}
	close($fh);
	return $meta;
}

1;

=head1 NAME

gramene-provenance.pl

=head1 SYNOPSIS

gramene-provenance.pl [options] metadata.tbl

Options:

    --output, -o : Output ontology file (provenance.json)
    --help,   -h : Brief usage instructions
    --man,    -m : More verbose usage instructions

=head1 OPTIONS

B<--help>, B<-h>
    Brief usage instructions.

B<--man>, B<-m>
    More verbose usage instructions.   

=head1 DESCRIPTION

Creates a KBase provenance file for Gramene genomes.

=head1 AUTHOR

Shiran Pasternak <shiranpasternak@gmail.com>

=cut    

