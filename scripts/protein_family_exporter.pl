#!/usr/local/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Spec;

use lib '/usr/local/gramene/lib/perl';


use lib map { "/usr/local/ensembl-live/$_" } qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-draw/modules ensembl-compara/modules ensembl-variation/modules/ gramene-live/ensembl-plugins/gramene/modules gramene-live/ensembl-plugins/maize/modules conf);

use Bio::EnsEMBL::Registry;

my %args = (
    'registry' => '/usr/local/gramene/conf/ensembl.registry',
    'pair_id' => 1,
);

GetOptions(
    \%args,
    'species=s@',
    'registry=s',
    'output_dir=s',
    'pair_id=s',
);

#print Dumper(\%args); use Data::Dumper; exit;

Bio::EnsEMBL::Registry->load_all($args{'registry'});

#my $pair_id = $args{'pair_id'};
my $orth_id=1;
my $para_id=1;

my $seen_genes = {};

my $family_ids = {};

my $relationships = {};

my $throttle = 0;

#my %valid_ids = map {$_ => 1} ('BRADI0007S00200');
my %valid_ids = ();

open my $families_fh, ">>", File::Spec->catfile($args{'output_dir'}, "plants.families.2c");          #$args{'species'}.families.2c");
open my $functions_fh, ">>", File::Spec->catfile($args{'output_dir'}, "plants.families.functions");  #$args{'species'}.families.functions");

foreach my $species (@{$args{'species'}}) {

    print $families_fh "# homologues for $species\n";

    my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');

    my $ga = $dba->get_GeneAdaptor();

    my $genes = $ga->fetch_all();

    while (my $gene = shift @$genes) {
        my $seen_homologues = {};

        #last if $throttle++ > 100;
        print "@ ", $throttle, "\n" unless ++$throttle % 500;
        #print "$throttle\n";

    #    next unless $valid_ids{ $gene->stable_id };

        next if $seen_genes->{$gene->stable_id}++;

        my $orthologues = $gene->get_all_homologous_Genes();

        while (my $o = shift @$orthologues) {

            next if $seen_homologues->{$o->[0]->stable_id . $o->[2]}++;

            my ($gene_id_1, $gene_id_2) = ($gene->stable_id, $o->[0]->stable_id);

            $valid_ids{$gene_id_1}++;
            $valid_ids{$gene_id_2}++;

            my $type = $o->[1]->description;# =~ /ortholog/ ? 'ortholog' : 'paralog';

            my $prefix;

            if ($type =~ /ortholog/) {
                $prefix = 'EPCO_';
    #            print $gene->species, " vs " , $o->[0]->species, "<-O\n";
            }
            elsif ($type =~ /paralog/) {
                $prefix = 'EPCP_';
    #            print $gene->species, " vs " , $o->[0]->species, "\n";
            }
            else {
    #            print $gene->species, " vs " , $o->[0]->species, "<-P\n";
                warn "UNKNOWN RELATIONSHIP - $type for $gene_id_1 -> $gene_id_2";
                next;
            }

            my $local_id = $family_ids->{$gene_id_1} || $family_ids->{$gene_id_2};
            if (! defined $local_id) {
		$local_id = $prefix . $orth_id++ if($prefix eq "EPCO_";
		$local_id = $prefix . $para_id++ if($prefix eq "EPCP_";
                my $description = join(';', $type, 'EnsemblPlantCompara');
                print $functions_fh join("\t", $local_id, $description), "\n";
            }

            if (! defined $family_ids->{$gene_id_1}) {
                $family_ids->{$gene_id_1} = $local_id;
                print $families_fh join("\t", $local_id, $gene_id_1), "\n";
            }

            if (! defined $family_ids->{$gene_id_2}) {
                $family_ids->{$gene_id_2} = $local_id;
                print $families_fh join("\t", $local_id, $gene_id_2), "\n";
            }

            $relationships->{$local_id}->{$gene_id_1}->{$gene_id_2}++;
            $relationships->{$local_id}->{$gene_id_2}->{$gene_id_1}++;

        }

    }
}

close $families_fh;
close $functions_fh;

open my $test_fh, ">", File::Spec->catfile($args{'output_dir'}, "plants.families.test");
foreach my $rel (keys %$relationships) {
    print $test_fh $rel,"\t",join("|",sort keys %{$relationships{$rel}}),"\n";
}
close $test_fh;

__DATA__
foreach my $rel (keys %$relationships) {
    my $rel_count = 0;
    my $local_rel = $relationships->{$rel};
    print Dumper($local_rel); use Data::Dumper;
    foreach my $gene_id ( keys %$local_rel ) {
        if (! $rel_count) {
            $rel_count = keys %{ $local_rel->{$gene_id} };
            next;
        }
        elsif ($rel_count != keys %{ $local_rel->{$gene_id} }) {
            die "OMFG - mismatched count on $gene_id ($rel_count != ", scalar keys %{ $local_rel->{$gene_id} }, " for $rel";
        }
    }
}
