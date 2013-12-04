#!/usr/bin/env perl
use warnings;
use strict;
my @temp=();

my $Dir=$ARGV[0];
exit if !$Dir || !-d $Dir || !-f $Dir."/features.tab";
$Dir.="/" if $Dir !~ /\/$/;

print "Testing $Dir\n";
open(FH, "< ".$Dir."/features.tab");
my %Loci=();
my %mRNAs=();
my %CDSs=();
my %Untranslated=();
while(<FH>){
    chomp;
    @temp=split(/\s+/,$_,-1);
    $Loci{$temp[0]}=0 if $temp[1] eq "locus";
    $mRNAs{$temp[0]}={defined=>0,parent=>$temp[3]} if $temp[1] eq "mRNA";
    $CDSs{$temp[0]}=$temp[3] if $temp[1] eq "CDS";
}
close(FH);

foreach my $mRNA ( keys %mRNAs ){
    if(exists($Loci{$mRNAs{$mRNA}{parent}})){
	$Loci{$mRNAs{$mRNA}{parent}}=$mRNA;
    }
}

foreach my $CDS ( keys %CDSs ){
    if(exists($mRNAs{$CDSs{$CDS}})){
	$mRNAs{$CDSs{$CDS}}{defined}=$CDS;
    }
}

my %Missing=();
open(FH, "< ".$Dir."/proteins.fa");
while(<FH>){
    chomp;
    next unless $_ =~ s/^>//;

    $Missing{$_}=1 if !exists($CDSs{$_});
}
close(FH);

mkdir $Dir."../Logs/" if !-d $Dir."../Logs/";
open(OUT, "> ".$Dir."../Logs/exchange_format_identifiers.out");
open(ERR, "> ".$Dir."../Logs/exchange_format_identifiers.err");
print OUT scalar(keys %Loci)," genes\n";
print OUT scalar(keys %mRNAs)," mRNAs\n";
print OUT scalar(keys %CDSs)," proteins\n";
print scalar(keys %Loci)," genes\n";
print scalar(keys %mRNAs)," mRNAs\n";
print scalar(keys %CDSs)," proteins\n";
print ERR scalar( grep { $Loci{$_} eq "0" } keys %Loci)," loci not transcribed into mRNA\n";
print ERR scalar( grep { $mRNAs{$_}{defined} eq "0" } keys %mRNAs)," mRNAs not translated into proteins\n";
print scalar( grep { $Loci{$_} eq "0" } keys %Loci)," loci not transcribed into mRNA\n";
print scalar( grep { $mRNAs{$_}{defined} eq "0" } keys %mRNAs)," mRNAs not translated into proteins\n";
if(scalar( grep { $Loci{$_} eq "0" } keys %Loci )>0){
    print ERR join("\n", map { "Untranscribed locus $_" } grep { $Loci{$_} eq "0" } keys %Loci),"\n";
}
if(scalar( grep { $mRNAs{$_}{defined} eq "0" } keys %mRNAs )>0){
    print ERR join("\n", map { "Untranslated mRNA ".$mRNAs{$_}{parent} } grep { $mRNAs{$_}{defined} eq "0" } keys %mRNAs),"\n";
}
if(scalar(keys %Missing)>0){
    print ERR join("\n",map {"Missing protein sequence $_"} keys %Missing),"\n";
}
