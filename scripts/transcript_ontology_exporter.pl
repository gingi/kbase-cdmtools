#!/usr/local/bin/perl

use lib '/usr/local/gramene/lib/perl';


use lib map { "/usr/local/ensembl-live/$_" } qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-draw/modules ensembl-compara/modules ensembl-variation/modules/ gramene-live/ensembl-plugins/gramene/modules gramene-live/ensembl-plugins/maize/modules conf);

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::Registry;

my %args = (
    'registry' => '/usr/local/gramene/conf/ensembl.registry'
);

GetOptions(
    \%args,
    'species=s',
    'transcript=s',
    'registry=s',
);

Bio::EnsEMBL::Registry->load_all($args{'registry'});

my $transcriptA = Bio::EnsEMBL::Registry->get_adaptor($args{'species'}, 'Core', 'Transcript');
my $goA         = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'ontology', 'goterm');

my $transcripts = $args{'transcript'}
    ? [$transcriptA->fetch_by_stable_id($args{'transcript'})]
    : $transcriptA->fetch_all;

my $delim = "\t";
my $altdelim = ",";

print join($delim, qw( TranscriptID OntologyID OntologyDescription OntologyEvidenceCode OntologyType )), "\n";  

foreach my $transcript (@$transcripts) {

    my $links = $transcript->get_all_DBLinks();
    foreach my $link (@$links) {
        next unless (ref $link) =~ /OntologyXref/;
        print 
            join(
                $delim,
                    $transcript->display_id,
                    $link->display_id,
                    $link->description,
                    join(
                        $altdelim,
                            @{$link->get_all_linkage_types}
                    ),
                    $link->dbname,
            ),
            "\n";
            
        my $ontology_term = $goA->fetch_by_accession($link->display_id);
        foreach my $slimterm ( @{ $ontology_term->ancestors } ) {
            my $s = join(',', @{$slimterm->subsets});
            next unless $s =~/goslim_generic/;
            #print Dumper($slimterm); last; use Data::Dumper;
            print 
                join(
                    $delim,
                        $transcript->display_id,
                        $slimterm->accession,
                        $slimterm->name,
                        '', #use term evidence code?
                        $slimterm->ontology,
                ),
                "\n";
        }
    }
}

=pod

=head1 DESCRIPTION

spits out transcript ontology links (accession, term, linkage type) to STDOUT. 
give it options :

 --registry (registry file to use (defaults to Gramene registry file in known location)
 --species (species you're hitting)
 --transcript (the transcript ID to dump) (optional. If omitted will dump all transcripts)

=cut

