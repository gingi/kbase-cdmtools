#!/usr/bin/perl

#*************************************************
#name           art_to_gff.pl
#
#description    Convert the exchange format lines to .gff
#
#author         Miriam Land - ORNL
#
#history        Feb. 2012 - Initial code
#
#*************************************************

my $root_dir = '/home/ml3/dev/KBase';

my $gff  = "$root_dir/Ath/TAIR10_GFF3_genes_gff";
my $fun  = "$root_dir/Ath/functions.Athaliana.TAIR10.tsv";
my $prot = "$root_dir/Ath/proteins.Athaliana.TAIR10.fa";

open (GFF,$gff) || die "Could not open $gff";   # Input GFF

my %List;
my %Prot;

while (my $buf = <GFF>)
{
	chomp $buf;
	my @ary = split(/	/,$buf);
	next unless ($ary[2] eq 'gene' || $ary[2] eq 'CDS');
	@id_str = split(/[;=,]/,$ary[8]);
	$id = $id_str[1];
	$type = $id_str[3];
	$List{$id} = $type if ($ary[2] eq 'gene' );
	$Prot{$id} = 'protein' if ($ary[2] eq 'CDS'  );
#	print "DEBUG: ID=$id PROT=$Prot{$id} LIST=$List{$id} ARY2=$ary[2]\n" if ($id =~ /AT1G0103/);
}

close GFF;


open (FH,$prot) || die "Could not open $prot";   # proteins

while (my $id = <FH>)
{
	next unless ($id =~ /^\>/);
	chomp $id;
	$id =~ s/^\>//;
	($parent,$dummy) = split(/\./,$id);

	if (exists $Prot{$id}) { 
		$Prot{$id} = 'Found';
		print "Inconsistent type for $id\n" if ($List{$parent} ne 'protein_coding_gene');
	}
	else {
		print "CDS for $id not found in the GFF\n";
	}
}
close FH;

#
#	For every CDS in the GFF, did we find a protein?
#
foreach $key (keys(%Prot))
{
#	print "DEBUG: ID=$id PROT=$Prot{$key} \n" if ($id =~ /AT1G0103/);
	next if ($Prot{$key} eq 'Found' );
	print "Did not find a Protein for $key\n";
}



open (FH,$fun) || die "Could not open $fun";   # Input functions

while (my $buf = <FH>)
{
	chomp $buf;
	($id,$function) = split(/	/,$buf);

#
#	Was the ID found in the GFF above
#	If it was, is the function consistent with the type of gene
#
	if (exists $List{$id}) { 
		$type = $List{$id};
		$List{$id} = 'found';

		next if ($type eq 'protein_coding_gene' || $function =~ /$type/);
		next if ( ($type eq 'other_RNA' || $type eq 'miRNA' ) && ($function =~ /other RNA/ || $function =~ /otherRNA/) );
		next if ($type eq 'rRNA' && $function =~ /ribosomal/);
		next if ( ($type eq 'other_RNA' ) && ($function =~ /Potential natural antisense gene/ ) );
		next if ( ($type eq 'other_RNA' ) && ($function =~ /Unknown gene/ ) );

		
		print "INCONSISTENT: $buf -- TYPE=$type\n" ;
	}
	else  {
		print "Did not find a GFF gene for $buf\n";
	}
}
close FH;

#
#	For every gene in the GFF, did we find a function?
#
foreach $key (keys(%List))
{
	next if ($List{$key} eq 'found' );
	print "Did not find a function for $key $List{$key}\n";
}


