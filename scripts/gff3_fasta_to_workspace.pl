#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();
my $header=1;


################################################################
##
## Input is 1 contig fasta file, 1 protein fasta file, 
## and 1 annotation file from SEED
##
################################################################

my $contigs = $ARGV[0];
exit if !$contigs || !-f $contigs;
exit if $contigs !~ /^[\w\|.-]+$/;

my $proteins = $ARGV[1];
exit if !$proteins || !-f $proteins;
exit if $proteins !~ /^[\w\|.-]+$/;

my $annotation = $ARGV[1];
exit if !$annotation || ! -f $annotation;


################################################################
##
## Read files
##
################################################################

open(FH, "< $contigs");
my $Contigs="";
while(<FH>){
    $Contigs.=$_;
}
close(FH);

my @contigss = map { $_->[2] =~ tr/ \n\r\t//d; $_ }
               map { /^(\S+)([ \t]+([^\n\r]+)?)?[\n\r]+(.*)$/s ? [ $1, $3 || '', $4 || '' ] : () }
               split /[\n\r]+>[ \t]*/m, $Contigs;

open(FH, "< $proteins");
my $Proteins="";
while(<FH>){
    $Proteins.=$_;
}
close(FH);

my @proteinss = map { $_->[2] =~ tr/ \n\r\t//d; $_ }
               map { /^(\S+)([ \t]+([^\n\r]+)?)?[\n\r]+(.*)$/s ? [ $1, $3 || '', $4 || '' ] : () }
               split /[\n\r]+>[ \t]*/m, $Proteins;

open(FH, "< ".$annotation);
my %Functions=();
while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    my ($id,$function)=($temp[1],$temp[2]);
    $Functions{$id}=$function;
}

################################################################
##
## Build Contigs and Features JSON objects
##
################################################################

my @Contigs=();
foreach my $seq (@contigs){
    $seq->[0] =~ s/^>//;

    my $contigHash = {id => $seq->[0],
		      seq => $seq->[2]};

    push(@Contigs,$contigHash);
}    

my @Features=();
foreach my $seq (@proteins){

#Location arrayref is structured as
#location => [[$seq->[0],$start,[+-],$length]]

    my $featureHash = {id => $seq->[0],
		       location => [],
		       type => "CDS",
		       function => "",
		       alternative_functions => [],
		       protein_translation => $seq->[2],
		       aliases => [$seq->[0]],
		       annotations => []};

    if(exists($Functions{$seq->[0]})){
	$featureHash->{function}=$Functions{$seq->[0]};
    }

    push(@Features,$featureHash);
}

################################################################
##
## Build Entire Genome JSON Object and load into workspace "Plants"
##
################################################################

my %genomeData = (id=>"Test",
                  scientific_name=>"unknown",
                  domain=>"Plant",
                  genetic_code=>11,
                  source=>"Private-User",
                  source_id=>"unknown",
                  contigs=>\@Contigs,
                  features=>\@Features);

use Bio::KBase::workspaceService::Helpers qw(auth);
my $AToken = auth();

use Bio::KBase::fbaModelServices::Client;
my $FBAClient = Bio::KBase::fbaModelServices::Client->new('http://www.kbase.us/services/KBaseFBAModeling/');
$FBAClient->genome_object_to_workspace({workspace=>"Plants",auth=>$AToken,genomeobj=>\%genomeData,mapping=>"PlantSEED",mapping_workspace=>"PlantSEED"});
