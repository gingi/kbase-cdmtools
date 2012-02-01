#!/usr/local/bin/perl

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
Readonly my %CONVERT_FTYPE => (
    gene => 'locus',
    mRNA => 'mRNA',
    CDS  => 'CDS',
);
Readonly my @USABLE_TYPES => qw/gene mRNA CDS/;

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

    my $lines = 0;

    # Read features from GFF3
    while (my $line = <$fh>) {
        next if ($line =~ m/^\s*#/);
        chomp $line;
        @fields = split(/\s+/, $line, 9);
        next unless (grep { $_ eq $f->('type') } @USABLE_TYPES);
        my %extra = map { m/(\S+)\s*=\s*(\S+)/; lc $1 => $2; }
            split(";", $f->('extra'));

        # Hash by type/location, not ID (which could be undefined)
        my $hashkey = join(':', map { $f->($_) } qw/type map start end/);
        my $id      = $extra{id};
        my $feature = Feature->new();
        $feature->id          = $id;
        $feature->type        = $f->('type');
        $feature->map         = $f->('map');
        $feature->start       = $f->('start');
        $feature->end         = $f->('end');
        $feature->orientation = $f->('orientation');
        $feature->name        = $extra{name};
        $feature->parent      = $extra{parent};
        if (defined $feature->parent) {
            $feature->parent
                =~ s/,.*$//;   # Stripping extra parents. TODO: Is this smart?
            my $parent = $parents{ $feature->parent };
            if (!defined $parent) {
                $orphans{$hashkey} = 1;
            } else {
                my $children = $parent->children;
                push @$children, $feature;
                if ($feature->id) { $parents{ $feature->id } = $feature; }
            }
        } else {
            $parents{ $feature->id } = $feature;
            push @features, $feature;
        }
        $lines++;

        # print $feature->to_string, "\n";
    }
    compress_features(\@features);
    return \@features;
}

sub compress_features {
    my ($features) = @_;

    for my $feature (@$features) {
        my @children = @{ $feature->children };
        next unless scalar @children;
        my $first_child = $children[0];
        if (defined $first_child->id && $first_child->id ne '') {
            compress_features(\@children);
        } else {
            assign_feature_id($first_child);
        }

        # Assign location of mRNA based on component CDS features
        if ($feature->type eq 'mRNA' && $first_child->type eq 'CDS') {
            my $mrna_end = $feature->end;
            $feature->end = $first_child->end;
            for (my $i = 1; $i <= $#children; $i++) {
                $feature->merge_location($children[$i]);
            }
            $feature->nth_location($#children)->{end} = $mrna_end;
        }

        for (my $i = 1; $i <= $#children; $i++) {
            $first_child->merge_location($children[$i]);
            delete $feature->children->[$i];
        }
        $feature->set_children([$first_child]);
    }
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

package Feature;

sub new {
    my $class = shift;
    $class = ref $class if ref $class;
    my $self = bless {}, $class;
    $self->{children}  = [];
    $self->{locations} = [];
    $self;
}

sub id : lvalue {
    shift->{id};
}

sub map : lvalue {
    shift->{map};
}

sub start : lvalue {
    shift->{start};
}

sub end : lvalue {
    shift->{end};
}

sub orientation : lvalue {
    shift->{ori};
}

sub type : lvalue {
    shift->{type};
}

sub name : lvalue {
    shift->{name};
}

sub parent : lvalue {
    shift->{parent};
}
sub children { shift->{children} }

sub set_children {
    my ($self, $children) = @_;
    $self->{children} = $children;
}

sub location {
    my $self = shift;
    if (scalar @{ $self->{locations} } > 0) {
        return join(
            ';',
            map {
                sprintf('%s_%d%s%d',
                    $_->{map}, $_->{start}, $_->{orientation},
                    $_->{end} - $_->{start});
                } @{ $self->{locations} }
        );
    }
    return sprintf('%s_%d%s%d',
        $self->map, $self->start, $self->orientation,
        $self->end - $self->start);
}

sub nth_location {
    my ($self, $n) = @_;
    return $self->location_tuple if (scalar @{ $self->{locations} } == 0);
    return $self->{locations}->[ $n - 1 ];
}

sub merge_location {
    my ($self, $other) = @_;
    $self->{locations} ||= [ $self->location_tuple ];

# print "Merging location ", join(" ", map { $other->location_tuple->{$_} } qw/map start end/), "\n";
    push @{ $self->{locations} }, $other->location_tuple;
}

sub location_tuple {
    my $self = shift;
    return +{ map { $_ => $self->$_ } qw/map start end orientation/ };
}

sub to_string {
    my $self = shift;
    return join("\t",
        $self->id,
        $CONVERT_FTYPE{ $self->type },
        $self->parent || '.',
        $self->name   || '.',
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
