#!/usr/local/bin/perl

use strict;
use warnings;
use autodie;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Readonly;
use DBI;
use List::MoreUtils qw(uniq);
use Gramene::Config;
use Bio::DB::Taxonomy;
use File::Spec;

my $options = {
    'host'       => 'cabot',
    'port'       => 3306,
    'output_dir' => './',
};

GetOptions(
    $options,
    'help',
    'man',
    
    # flag to spit out all ensembl dbs in the gramene conf, and the
    # corresponding output dir, which defaults to ./ if not specified
    'all_gramene_dbs',
    'output_dir=s',
    
    #OR a specific gramene DB
    'gramene_db=s',

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

#first case, the user requested ALL gramene dbs. Grab and dump 'em all out.
if ($options->{'all_gramene_dbs'}) {

    my $gramene_config  = Gramene::Config->new->getall;
    foreach my $key (keys %$gramene_config) {
        next unless $key =~ /^ensembl_/;
    
        my $dbh = DBI->connect(
            $gramene_config->{$key}->{'db_dsn'},
            $gramene_config->{$key}->{'db_user'},
            $gramene_config->{$key}->{'db_pass'},
            
        );
        
        print STDERR "Dumping database : $key\n";
        
        dump_meta_table($dbh, undef, $key);

    }

    #we're done. Bow out.
    exit;
}

# Next case, the user requested just one gramene db. Pull out all info from the
# gramene config file
elsif ($options->{'gramene_db'}) {
    my $gramene_config  = Gramene::Config->new->getall;
    
    if (my $db_params = $gramene_config->{ $options->{'gramene_db'} } ) {
        $options->{'dsn'}  = $db_params->{'db_dsn'},
        $options->{'user'} = $db_params->{'db_user'},
        $options->{'pass'} = $db_params->{'db_pass'},
    }
    else {
        die "Cannot connect to database: $options->{'gramene_db'} : Imaginary.";
    }
}

# otherwise, if a dsn was NOT passed, but we do have a database, then build
# a dsn from the database, host, and port
elsif (! defined $options->{'dsn'} && defined $options->{'database'}) {
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
        'species.scientific_name'   => {key => 'species',          order => 3},
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
        my $taxonid = $biodb->get_taxonid('Homo sapiens');
        my $node = $biodb->get_Taxonomy_Node(-taxonid => $output->{'taxonomy-id'});
        if (defined $node) {
            my @classification = ();
            do {
                unshift @classification, $node->name('scientific')->[0];
            } while ($node = $node->ancestor);
            
            $output->{'classification'} = \@classification if @classification;
        }
        else {
            print STDERR "\tNo entrez classification. Using unsorted Ensembl DB values";
        }
    }
    
    #if we were given no output filehandle, open one now:
    if (! $fh) {
        my $filename = lc (
               $output->{'species'}->[0]
            || $output->{'common-name'}->[0]
            || $default_file_name
            || time);   #if we got absolutely no clue, just timestamp it.
        $filename =~ s/\s+/_/g;
        
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
        sort { $data_key_map->{$a}->{'order'} <=> $data_key_map->{$b}->{'order'} }
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

generate.metadata.table.pl - Generate metadata.tbl files from Ensembl 'meta' table.

=head1 SYNOPSIS

  meta2.table.pl 

Options:

  --help   Show brief help and exit
  --man    Show full documentation
  
  There are several ways to get to a database:
  
  --gramene_db specifies a key in the gramene.conf file. All connection
               info will be read from there.
  
  
  Alternatively, you may pass:
  
  --dsn a full mysql dsn driver string
  or
  --database the name of the database to connect to, which will use
            default connection settings
  
  
  and additional to dsn or database:
  
  --user the user to connect with
  --pass user's password
  --host host to connect to (will default to cabot)
  --port port to connect to (will default to 3306)
  
  
  Finally, you an just specify all known gramene ensembl databases:
  
  --all_gramene_dbs
  
  In that case, it'll name each file based off of the species name, or failing that
  the first common name, or failing that the conf key, or failing that just a timestamp.
  Will dump all
  files to the current directory, or the output_dir option:
  
  --all_gramene_dbs --output_dir=/tmp/kbase-metadata

=head1 DESCRIPTION

Will export a file in exchange format of ensembl meta data. Prints to STDOUT,
so direct to where you want to go:

./generate.metadata.table.pl --gramene_db ensembl_arabidopsis > arabidopsis.meta.tbl

Or, just dump 'em all:

./generate.metadata.table.pl --all_gramene_dbs --output_dir=/tmp/kbase-meta

=head1 SEE ALSO

perl.

=head1 AUTHOR

Jim Thomason E<lt>thomason@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2012 Cold Spring Harbor Laboratory

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
