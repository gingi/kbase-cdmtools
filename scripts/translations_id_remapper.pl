#!/usr/local/bin/perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Text::RecordParser;

my %args = ('mapfield' => 'name');

GetOptions(
    \%args,
    'features=s',
    'translations=s',
    'output=s',
    'mapfield=s',
);

die "No features file"                                           unless $args{'features'};
die "No translation file"                                        unless $args{'translations'};
warn "No output file specified...updating translations in place" unless $args{'output'};

$args{'output'} ||= $args{'translations'};


my $featuresP = Text::RecordParser->new($args{'features'});
$featuresP->field_separator("\t");
$featuresP->bind_fields(qw(id type parent name excess)); 

my $mapping = {};
my $transcript_id_mapping = {};

while (my $feature = $featuresP->fetchrow_hashref) {

    if ($feature->{'type'} eq 'mRNA') {
        $transcript_id_mapping->{ $feature->{'id'} } = $feature->{'name'};
    }

    if ($feature->{'type'} eq 'CDS') {
        if (defined $mapping->{ $feature->{'parent'} }) {
            die "Multiple CDses provides for $feature->{'parent'}";
        }
        $mapping->{ $transcript_id_mapping->{ $feature->{'parent'} } } = $feature->{ $args{'mapfield'} };
    }
}

my @output = ();

open (my $translationsfh, "<", $args{'translations'});
$/ = "\n>";

while (my $record = <$translationsfh>) {
    chomp($record);
    my ($id, $sequence) = split /\n/, $record, 2;
    $id =~ s/^>//;
    
    if ($mapping->{$id}) {
        $id = $mapping->{$id};
    }
    else {
        warn "No mapping provided for $id on line $.";
    }
    
    push @output, ">$id\n$sequence";
}

close $translationsfh;

warn "Overwriting $args{'output'}" if -e $args{'output'};

open(my $outputfh, ">", $args{'output'});
print $outputfh join("\n", @output);
close $outputfh;

__DATA__

=pod

=head1 DESCRIPTION

Given a features file and a translations file, will remap all transcript ids in
the translation file to CDS ids provided in the parent file. Each transcript
must be associated with only one CDS.

Arguments:

 --features      /path/to/features.file
 --translations  /path/to/translations.file
 --output        /path/to/output.file  (optional. Will default to --translations
                                        file and update in place if not set)
--mapfield       field to map transcript ids to. Defaults to "name"

=cut