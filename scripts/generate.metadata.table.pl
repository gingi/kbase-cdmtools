#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Readonly;
use DBI;
use List::MoreUtils qw(uniq);
use Bio::DB::Taxonomy;
use File::Spec;

my $options = {
    'port'       => 3306,
    'output_dir' => './',
};

GetOptions(
    $options,
    'help',
    'man',
    
    'input_file=s',
    'output_dir=s',
    
    #OR a raw DSN
    'dsn=s',
    #OR a raw database name
    'database=s',
    
    #user and password, if necessary
    'user=s',
    'pass=s',

    #host and port, if necessary
    'host=s',
    'port=i',
    
) or pod2usage(2);


if ( $options->{'help'} || $options->{'man'} ) {
    pod2usage({
        -exitval => 0,
        -verbose => $options->{'man'} ? 2 : 1
    });
}; 


#first case - we were given a file that contains a list of databases
if ($options->{'input_file'}) {

    open (my $inputfh, "<", $options->{'input_file'});
    while (defined (my $dsn = <$inputfh>)) {
        chomp($dsn);
        
        #not actually a dsn? Assume it's a database name and build a DSN.
        if ($dsn !~ /^dbi:/) {
            $dsn = sprintf(
                'dbi:mysql:%s;host=%s;port=%d',
                    $dsn,
                    $options->{'host'},
                    $options->{'port'}
            );
        }
        
        my $dbh = DBI->connect(
            $dsn,
            $options->{'user'},
            $options->{'pass'},
        );
        
        print STDERR "Dumping database : $dsn\n";
        dump_meta_table($dbh, undef, $dsn);
    }

    #we're done. Bow out.
    exit;
}

# if a dsn was NOT passed, but we do have a database, then build
# a dsn from the database, host, and port
if (! defined $options->{'dsn'} && defined $options->{'database'}) {
    $options->{'dsn'} = sprintf(
        'dbi:mysql:%s;host=%s;port=%d',
            $options->{'database'},
            $options->{'host'},
            $options->{'port'}
    );
}

die "No DSN provided. Cannot connect." unless $options->{'dsn'};

#finally, connect to the db
my $dbh = DBI->connect(@$options{qw(dsn user pass)});

#and spit it out to STDOUT
dump_meta_table($dbh, \*STDOUT);


#--------------
# the one and only function necessary to do the work. Accepts a DBI database
# handle to read from and a filehandle to print it out to. Also accepts a 3rd
# argument - a fallback file name if it can't figure one out from the data
# in the database. Yay integrity!

sub dump_meta_table {
    my ($dbh, $fh, $default_file_name) = @_;

    my $meta_sth = $dbh->prepare_cached(<<"eSQL");
        SELECT
            meta_key,
            meta_value
        FROM
            meta
        WHERE
            meta_key like 'species.%'
            OR meta_key like 'assembly.%'
eSQL
    
    $meta_sth->execute();
    
    #we need to map from ugly ensembl keys to pretty metadata keys.
    my $data_key_map = {
        'assembly.name'             => {key => 'assembly-version', order => 1},
        'assembly.date'             => {key => 'assembly-date',    order => 2},
        'species.scientific_name'   => {key => 'name',             order => 3},
        'species.common_name'       => {key => 'common-name',      order => 4},
        'species.alias'             => {key => 'common-name',      order => 5},
        'species.classification'    => {key => 'classification',   order => 6},
        'species.division'          => {key => 'source',           order => 7},
        'species.taxonomy_id'       => {key => 'taxonomy-id',      order => 8}
    };
    
    #store our output info
    my $output = {};
    
    # actually loop through the data, map the keys, and toss it into our output
    # hash. ONLY keep keys that we're aware of and have mapped.
    while (my $row = $meta_sth->fetchrow_hashref('NAME_lc')) {
    
        #don't auto-vivify the key!
        my $output_set = $data_key_map->{ $row->{'meta_key'} } or next;
        my $output_key = $output_set->{'key'};
    
        push @{ $output->{$output_key} }, $row->{'meta_value'};

    }
    
    
    {
        # But we don't -really- trust the ensembl database to properly
        # populate the classification. So we use BioPerl instead.
        my $biodb = Bio::DB::Taxonomy->new('-source' => 'entrez');

        my $node =$biodb->get_Taxonomy_Node(-taxonid=>$output->{'taxonomy-id'});
        if (defined $node) {
            my @classification = ();
            do {
                unshift @classification, $node->name('scientific')->[0];
            } while ($node = $node->ancestor);
            
            $output->{'classification'} = \@classification if @classification;
        }
        else {
            print STDERR "\tNo entrez classification. Using unsorted DB values";
        }
    }
    
    #if we were given no output filehandle, open one now:
    if (! $fh) {
        my $filename = lc (
               $output->{'name'}->[0]
            || $output->{'common-name'}->[0]
            || $default_file_name
            || time);   #if we got absolutely no clue, just timestamp it.
        $filename =~ s/\W+/_/g;
        
        open (
            $fh,
            '>',
            File::Spec->catfile(
                $options->{'output_dir'},
                "$filename.tbl"
            )
        ) or die $!;
    }

    # some input keys map to the same output key. only print 'em out once. Not
    # grepped in the for loop for legibility's sake.
    my %seen = ();
    
    # finally, loop through our $output structure, sort it, and spit out
    # the data in exchange format.
    foreach my $output_key (
        grep { defined $output->{$_} }
        map  { $data_key_map->{$_}->{'key'} }
        sort {$data_key_map->{$a}->{'order'} <=> $data_key_map->{$b}->{'order'}}
        keys %$data_key_map
            ) {
        
        next if $seen{$output_key}++;
        
        my @values;
        
        #classification is a special case. We pulled it out of entrez up above.
        if ($output_key eq 'classification') {
            @values = @{ $output->{$output_key} };
        }
        else {
            @values = sort uniq( @{ $output->{$output_key} });
        }

        print $fh join("\n",
            $output_key,
            @values,
            "//\n",
        );
    }
}

__END__

# ----------------------------------------------------

=pod

=head1 NAME

generate.metadata.table.pl - Generate metadata.tbl files
                             from Ensembl 'meta' table.

=head1 SYNOPSIS

  generate.metadata.table.pl 

Options:

  --help   Show brief help and exit
  --man    Show full documentation
  
  There are several ways to get to a database:
  
  You may pass:
  
  --input_file for a list of databases (or dsns) one per line to connect to
  --output_dir for a directory to output files to. Defaults to $CWD.
               If not using --input_file, will send info to STDOUT.
  
  --dsn a full mysql dsn driver string
  or
  --database the name of the database to connect to, which will use
            default connection settings
  
  
  and additionally, you'll probably need to specify:
  
  --user the user to connect with
  --pass user's password
  --host host to connect to (will default to cabot)
  --port port to connect to (will default to 3306)
  
  
=head1 DESCRIPTION

Will export a file in exchange format of ensembl meta data. Prints to STDOUT,
so direct to where you want to go:

./generate.metadata.table.pl --database some_arabidopsis_db > arabidopsis.tbl

=head1 SEE ALSO

perl.

=head1 AUTHOR

Jim Thomason E<lt>thomason@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2012-2013 Cold Spring Harbor Laboratory

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
