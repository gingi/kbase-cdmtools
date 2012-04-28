#!/usr/local/bin/perl

use lib '/usr/local/gramene/lib/perl';


use lib map { "/usr/local/ensembl-live/$_" } qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-draw/modules ensembl-compara/modules ensembl-variation/modules/ gramene-live/ensembl-plugins/gramene/modules gramene-live/ensembl-plugins/maize/modules conf);

use strict;
use warnings;

use Getopt::Long;
use File::Spec;
use Bio::EnsEMBL::Registry;

my %args = (
    'slimtype' => 'GO-slim',
    'registry' => '/usr/local/gramene/conf/ensembl.registry'
);

GetOptions(
    \%args,
    'species=s',
    'transcript=s',
    'registry=s',
    'output_dir=s',
    'slimtype=s',
    'term=s@',
);

$args{'term'} = { map { $_ => 1} @{ $args{'term'} || [] } };

my @all = Bio::EnsEMBL::Registry->load_all($args{'registry'});
my @db_adaptors = @{ Bio::EnsEMBL::Registry->get_all_DBAdaptors() };

my $goA         = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'ontology', 'goterm');
my %transcript_adaptors = ();

if ($args{'species'}) {
    $transcript_adaptors{ $args{'species'} } = (Bio::EnsEMBL::Registry->get_adaptor($args{'species'}, 'Core', 'Transcript'));
}
else {

    die "Must specify output_dir if no species!" unless $args{'output_dir'};

    foreach my $db (@db_adaptors) {
        eval {
            my $ta = $db->get_TranscriptAdaptor();
            $transcript_adaptors{ $db->species } = $ta;
        }
    }
}

warn "Dumping all terms" unless keys %{$args{'term'}};

foreach my $species (keys %transcript_adaptors) {

    warn "Dumping $species";

    my $transcriptA = $transcript_adaptors{$species};

    my $outputfh = \*STDOUT;

    if ($args{'output_dir'}) {
        open $outputfh, ">", File::Spec->catfile($args{'output_dir'}, "$species.exf");
    }

    my $transcripts = $args{'transcript'}
        ? [$transcriptA->fetch_by_stable_id($args{'transcript'})]
        : $transcriptA->fetch_all;
    
    my $delim = "\t";
    my $altdelim = ",";
    
    print $outputfh join($delim, qw( TranscriptID OntologyID OntologyDescription OntologyDomain OntologyEvidenceCode OntologyType )), "\n";  
    
    foreach my $transcript (@$transcripts) {

        my $links = $transcript->get_all_DBLinks();
        foreach my $link (@$links) {
        
            my $ontology_term = $goA->fetch_by_accession($link->display_id);
            my $ancestors = $ontology_term->ancestors || [];
            my $domain = @$ancestors
                ? $ancestors->[-1]->name
                : '';
        
            next unless (ref $link) =~ /OntologyXref/;
            print $outputfh 
                join(
                    $delim,
                        $transcript->display_id,
                        $link->display_id,
                        $link->description,
                        $domain,
                        join(
                            $altdelim,
                                @{$link->get_all_linkage_types}
                        ),
                        $link->dbname,
                ),
                "\n" unless scalar keys %{$args{'term'}} && ! $args{'term'}{ $link->dbname };

            foreach my $slimterm ( @{ $ancestors } ) {

                my $s = join(',', @{$slimterm->subsets});
                next unless $s =~/goslim_generic/;

                print $outputfh 
                    join(
                        $delim,
                            $transcript->display_id,
                            $slimterm->accession,
                            $slimterm->name,
                            $domain,
                            '', #use term evidence code?
                            $args{'slimtype'}, # hardwired to go-slim instead of $slimterm->ontology?,
                    ),
                    "\n" unless scalar keys %{$args{'term'}} && ! $args{'term'}{ $args{'slimtype'} };
            }

        }

    }
}

=pod

=head1 DESCRIPTION

spits out transcript ontology links (accession, term, linkage type) to STDOUT. 
give it options :

 --registry (registry file to use (defaults to Gramene registry file in known location)
 --species (species you're hitting) (optional leave blank to dump all in registry)
 --transcript (the transcript ID to dump) (optional. If omitted will dump all transcripts)
 --output_dir (optional directory to dump to. Required if no species given, otherwise goes to STDOUT)
 --term (specify multiple times). Terms to include. If omitted, then all terms
 --slimtype The ontology type to print out for go-slim terms. Optional. Defaults to 'GO-slim'. 
            Be careful mixing this with --term up above! Must specify the --slimtype term!

=cut
