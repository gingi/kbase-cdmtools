#!/usr/local/bin/perl

# spits out transcript ontology links (accession, term, linkage type) to STDOUT. 
# give it options :
#
# --registry (registry file to use (defaults to Gramene registry file in known location)
# --species (species you're hitting)
# --transcript (the transcript ID to dump)


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


my $transcript = $transcriptA->fetch_by_stable_id($args{'transcript'});

my $delim = "\t";
my $altdelim = ",";

print join($delim, qw(Accession Term Linkage)), "\n";  

my $links = $transcript->get_all_DBLinks();
foreach my $link (@$links) {
    next unless (ref $link) =~ /OntologyXref/;
    print 
        join(
            $delim,
                $link->display_id,
                $link->description,
                join(
                    $altdelim,
                        @{$link->get_all_linkage_types}
                )
        ),
        "\n";
}

