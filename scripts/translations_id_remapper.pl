#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Text::RecordParser;
use Pod::Usage;

MAIN: {
    my %options = (mapfield => 'feature_id');

    GetOptions(\%options, 'features=s', 'translations=s', 'output=s',
        'mapfield=s', 'help|h')
        || pod2usage(-verbose => 1);

    pod2usage(-verbose => 1) if $options{help};

    die "No features file"    unless $options{features};
    die "No translation file" unless $options{translations};

    $options{output} ||= $options{translations};

    my $featuresP = Text::RecordParser->new($options{features});
    $featuresP->field_separator("\t");
    $featuresP->bind_fields(qw(feature_id type location parent subset name));

    my %mapping;
    my %ancestors = (undef => undef);

    while (my $feature = $featuresP->fetchrow_hashref) {
        my $id     = $feature->{feature_id};
        my $parent = $feature->{parent};
        $parent = undef if ($parent eq ".");
        $mapping{$id}   = $id;
        $ancestors{$id} = $parent;
        while (defined $parent) {
            $mapping{$parent} = $feature->{feature_id};
            $parent = $ancestors{$parent};
        }
    }

    my @output = ();

    open(my $translationsfh, "<", $options{translations});
    $/ = "\n>";

    while (my $record = <$translationsfh>) {
        chomp($record);
        my ($id, $sequence) = split /\n/, $record, 2;
        $id =~ s/^>//;

        if ($mapping{$id}) {
            $id = $mapping{$id};
        }
        else {
            print STDERR "No mapping provided for $id on line $.\n";
        }

        push @output, ">$id\n$sequence";
    }

    close $translationsfh;

    # warn "Overwriting $options{output}" if -e $options{output};

    open(my $outputfh, ">", $options{output});
    print $outputfh join("\n", @output);
    close $outputfh;
}

1;

=pod

=head1 NAME

    translations_id_mapper.pl

=head1 SYNOPSIS

    translations_id_mapper.pl [options]

Given a features file and a translations file, will remap all transcript ids in
the translation file to CDS ids provided in the parent file. Each transcript
must be associated with only one CDS.

Options:

 --help, -h      This help message
 --features      /path/to/features.file
 --translations  /path/to/translations.file
 --output        /path/to/output.file  (optional. Will default to --translations
                                        file and update in place if not set)
 --mapfield      field to map transcript ids to. Defaults to "feature_id"

=cut
