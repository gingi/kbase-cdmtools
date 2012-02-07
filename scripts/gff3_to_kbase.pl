#!/usr/bin/env perl

use strict;
use warnings;

use feature 'state';

use Getopt::Long;
use Pod::Usage;

use Readonly;

Readonly my %GFF3FIELDS => (
    map         => 0,
    source      => 1,
    type        => 2,
    start       => 3,
    end         => 4,
    phase       => 5,
    orientation => 6,
    score       => 7,
    extra       => 8
);
Readonly my %KBASE_FTYPE => (gene => 'locus',);
Readonly my @USABLE_TYPES => qw/gene exon mRNA tRNA rRNA CDS/;

MAIN: {
    my %options;
    GetOptions(\%options, 'gff=s', 'output|o=s', 'help|h', 'man|m')
        || pod2usage(-verbose => 1);

    pod2usage(-verbose => 1) if $options{help};
    pod2usage(-verbose => 2) if $options{man};

    my ($infh, $outfh);

    if (!defined($options{gff})) {
        $options{gff} = $ARGV[0];
    }

    if ($options{gff}) {
        open($infh, '<', $options{gff})
            or die "Cannot open GFF file $options{gff} for reading: $!\n";
    } else {
        $infh = \*STDIN;
    }
    my $features = read_gff3($infh);
    close($infh);

    if ($options{output}) {
        open($outfh, '>', $options{output})
            or die "Cannot open file $options{output} for writing: $!\n";
    } else {
        $outfh = \*STDOUT;
    }
    for my $feature (@$features) {
        print_children($outfh, $feature);
    }
    close($outfh);
}

sub print_children {
    my ($fh, $feature) = @_;
    print $fh $feature->to_string, "\n";
    for my $child (@{ $feature->children }) {
        print_children($fh, $child);
    }
}

sub read_gff3 {
    my ($fh) = @_;
    my (@features, %parents, %orphans);
    my (@fields);
    my $f = sub { $fields[ $GFF3FIELDS{ shift() } ] };

    # Read features from GFF3
    while (my $line = <$fh>) {
        next if ($line =~ m/^\s*#/);
        chomp $line;
        @fields = split(/\s+/, $line, 9);
        next unless (grep { $_ eq $f->('type') } @USABLE_TYPES);
        my %extra = map { m/(\S+)\s*=\s*(\S+)/; lc($1) => $2; }
            split(";", $f->('extra'));

        # Hash by type/location, not ID (which could be undefined)
        my $hashkey = join(':', map { $f->($_) } qw/type map start end/);
        my $id      = $extra{id};
        my $feature = KbaseExchangeFeature->new();
        $feature->id   = $id;
        $feature->type = $f->('type');
        $feature->add_location(
            map         => $f->('map'),
            start       => $f->('start'),
            end         => $f->('end'),
            orientation => $f->('orientation')
        );
        $feature->name = $extra{name};
        my (@parent_ids);

        if ($extra{parent}) {
            @parent_ids = split(',', $extra{parent});
        }

        if (scalar @parent_ids) {
            for my $parent_id (@parent_ids) {
                my $parent = $parents{$parent_id};
                if (!defined $parent) {
                    $orphans{$hashkey} = 1;
                } else {
                    $parent->add_child($feature);
                    $feature->add_parent($parent);
                    if ($feature->id) { $parents{ $feature->id } = $feature; }
                }
            }
        } else {
            $parents{ $feature->id } = $feature;
            push @features, $feature;
        }
    }
    for my $feature (@features) {
        compress_features($feature);
    }
    return \@features;
}

sub print_feature {
    my ($feature, $indent) = @_;
    $indent ||= '';
    print $indent, $feature->to_string, "\n";
    for my $child (@{ $feature->children }) {
        print_feature($child, "   $indent");
    }
}

sub compress_features {
    my ($feature) = @_;

    return unless scalar @{ $feature->children };

    if ($feature->type =~ m/RNA$/) {
        merge_cds($feature);
        assign_mrna_locations($feature);
    } else {
        for my $child (@{ $feature->children }) {
            compress_features($child);
        }
    }
}

sub merge_cds {
    my ($feature)    = @_;
    my %new_children = ();
    my $new_id       = ();
    for my $cds (@{ $feature->children_by_type('CDS') }) {
        if (!defined $cds->id || $cds->id eq '') {
            if (defined $new_id) {
                $cds->id = $new_id;
            } else {
                assign_feature_id($cds);
                $new_id = $cds->id;
            }
        }
        my $new_child = $new_children{ $cds->id };
        if (defined $new_child) {
            $new_child->merge_locations($cds);
        } else {
            $new_children{ $cds->id } = $cds;
        }
    }
    $feature->delete_children_by_type('CDS');
    for my $cds (sort { $a->id cmp $b->id } values %new_children) {
        $feature->add_child($cds);
    }
}

sub assign_mrna_locations {
    my ($mrna) = @_;
    $mrna->clear_locations();
    for my $exon (@{ $mrna->children_by_type('exon') }) {
        $mrna->merge_locations($exon);
    }
    $mrna->delete_children_by_type('exon');
}

sub assign_feature_id {
    my $feature = shift;

    state %Identifiers;
    if (scalar keys %Identifiers == 0) {
        %Identifiers = map { $_ => 1 } @USABLE_TYPES;
    }
    $feature->id = sprintf('%s%.5d', lc $feature->type,
        $Identifiers{ $feature->type }++);
}

package KbaseExchangeFeature;

sub new {
    my $class = shift;
    $class = ref $class if ref $class;
    my $self = bless {}, $class;
    $self->{children}  = [];
    $self->{locations} = [];
    $self->{parents}   = [];
    $self;
}

sub id : lvalue {
    shift->{id};
}

sub type : lvalue {
    shift->{type};
}

sub name : lvalue {
    shift->{name};
}

sub parents : lvalue {
    shift->{parents};
}

sub add_parent {
    my ($self, $parent) = @_;
    push @{ $self->{parents} }, $parent;
}

sub children : lvalue {
    shift->{children};
}

sub children_by_type {
    my ($self, $type) = @_;
    return [ grep { $_->type eq $type } @{ $self->children } ];
}

sub add_child {
    my ($self, $child) = @_;
    push @{ $self->{children} }, $child;
}

sub delete_children_by_type {
    my ($self, $type) = @_;
    for (my $i = $#{ $self->children }; $i >= 0; $i--) {
        if ($self->children->[$i]->type eq $type) {
            splice(@{ $self->children }, $i, 1);
        }
    }
}

sub clear_locations {
    shift->{locations} = [];
}

sub add_location {
    my $self      = shift;
    my %params    = ref($_[0]) eq 'HASH' ? %{ $_[0] } : @_;
    my $found     = 0;
    my $locations = $self->{locations};
    for my $loc (@$locations) {
        if (   $loc->{map} eq $params{map}
            && $loc->{start} == $params{start}
            && $loc->{end} == $params{end}
            && $loc->{orientation} eq $params{orientation})
        {
            $found = 1;
        }
    }
    if (not $found) {
        push @$locations, {%params};
        @$locations = sort { $a->{start} <=> $b->{start} } @$locations;
    }
}

sub location {
    my $self = shift;
    return join(
        ',',
        map {
            sprintf('%s_%d%s%d',
                $_->{map}, $_->{orientation} eq '+' ? $_->{start} : $_->{end},
                $_->{orientation}, $_->{end} - $_->{start} + 1);
            } @{ $self->{locations} }
    );
}

sub nth_location {
    my ($self, $n) = @_;
    return $self->{locations}->[ $n - 1 ];
}

sub merge_locations {
    my ($self, $other) = @_;
    for my $loc (@{ $other->{locations} }) {
        $self->add_location($loc);
    }
}

sub to_string {
    my $self = shift;
    my $parent_id
        = scalar @{ $self->parents } ? $self->parents->[0]->id : '.';
    return join("\t",
        $self->id                   || '.',
        $KBASE_FTYPE{ $self->type } || $self->type,
        $parent_id, $self->name     || '.',
        $self->location);
}

1;

=head1 NAME

gff3_to_kbase.pl - Converts GFF3 files of gene annotations to Kbase Exchange format

=head1 SYNOPSIS

gff3_to_kbase.pl --gff FEATURES.GFF --output FEATURES.TSV
gff3_to_kbase.pl < FEATURES.GFF > FEATURES.TSV

Options:

    --gff,    -g  : The input GFF file
    --output, -o  : The output file
    --help,   -h  : Brief usage instructions
    --man,    -m  : More verbose usage instructions

=head1 OPTIONS

B<--gff>,B<-g>
    The input GFF file. If not specified, the first positional parameter is used. If neither works, the script will try to read GFF output from stdin.

B<--output>,B<-o>
    The output file. Will overwrite the file if already exists. If not specified, the script prints to stdout.
    
B<--help>,B<-h>
    Brief usage instructions.

B<--help>,B<-h>
    More verbose usage instructions.   

=head1 DESCRIPTION

This script converts GFF3 output into the Kbase Exchange format.

=head1 AUTHOR

Shiran Pasternak <shiran@cshl.edu>

=cut    
