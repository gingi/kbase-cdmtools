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

my $root_dir = '/homes/oakridge/landml/dev/KBase/plants';
my $in_name = shift || die "Usage: $0 <file>";

unless (-e $in_name)
{
	die "Did not find $in_name\n";
}

my $gff = "$in_name.test";
my $source = 'TAIR10';

open (FH,$in_name) || die "Could not find $in_name"; # Export Formatted file (input)
open (GFF,">$gff") || die "Could not create $gff";   # Convert back to GFF (output)

print GFF "##gff-version 3\n";

#
#  Hashes for the range of features - Used to get the range of a parent
#	mRNA when processing a CDS.  Assumes that features such as mRNA and
#	CDS might span two nucleotide sequences.
#
#	Key = feature_id:sequence_id  e.g. AT1G01020:Chr1
#

undef my %First ;  # Left most base
undef my %Last  ;  # Right most base
undef my %Strand;  # Strand
undef my %Use_loc; # List of locations to use in GFF

# Input file format
#	feature_id
#	type - e.g., mRNA, locus, CDS
#	parent_feature_id
#	name - name be identical to feature_id
#	location - commas separated list of ranges.  Each range consists of:
#		-	an ID for the nucleotide sequence (e.g. Chr1)
#		-	an underscore
#		-	the start site
#		-	strand (plus or minus)
#		-	length

while (my $buf = <FH>)
{
	chomp $buf;
	my ($id,$type,$parent,$name,$location) = split(/	/,$buf);
#next unless ($buf =~ /AT1G01020/);
#last if ($buf =~ /AT1G01040/);

#
#	Split the location into ranges
#
	@loc_ary = split(/,/,$location);
	$Use_loc{$id} = {}; # List of locations to use in GFF

#
#	Parse the locations.  Find that First, Last, and Strand for each
#	nucleotide sequence where the feature is found
#	Convert the start/length/strand information into the
#	left/right/score/strand used by GFF3 and make a list for %Use_loc
#
	foreach my $loc (@loc_ary)
	{
		my ($seqid,$start,$length) = split(/[_+-]/,$loc);
		my $left   = 0;
		my $right  = 0;
		my $strand = '';
		my $key = "$id:$seqid";

#print "ID=$id  TYPE=$type PARENT=$parent SEQID=$seqid KEY=$key START=$start LENGTH=$length\n" ;

		unless (exists $Use_loc{$id}{$seqid} )
		{
			$Use_loc{$id}{$seqid} = [] ;
#print "DEBUG: ID=$id SEQID=$seqid Initialize Use_loc{id}{seqid} $Use_loc{$id}{$seqid}\n";
		}

		if ($loc =~ /\+/)
		{
			$left   = $start;
			$right  = $start + $length - 1;
			$strand = '+';
			push ($Use_loc{$id}{$seqid},"$left\t$right\t.\t$strand");
		}
		else
		{
			$right  = $start;
			$left   = $start - $length + 1;
			$strand = '-';
			unshift ($Use_loc{$id}{$seqid},"$left\t$right\t.\t$strand"); # Reverse the order
		}

		if ( (exists $First{$key} && $left < $First{$key}) || 
             ( ! exists $First{$key} ))  { 
			$First{$key} = $left;
		}
		if ( (exists $Last{$key} && $right > $Last{$key}) || 
             ( ! exists $Last{$key} ))  { 
             $Last{$key} = $right;
		}
		$Strand{$key} = $strand;
	}

	if ($type eq 'locus')
	{
#		Change display type from locus to gene
#		Comment = No parent here
		my $dis_type = 'gene';
		foreach my $seqid (keys(%{ $Use_loc{$id} }))
		{
			foreach my $loc (@{ $Use_loc{$id}{$seqid} })
			{
				print GFF "$seqid\t$source\t$dis_type\t$loc\t.\tID=$id;Note=protein_coding_gene;Name=$name\n";
			}
		}
	}
	elsif ($type eq 'mRNA')
	{
#		Comment = Has parent 
#		mRNA is First/Last/Strand
#		exons are reach of the individual regions
		foreach my $seqid (keys(%{ $Use_loc{$id} }))
		{
			$key = "$id:$seqid";
#
#			mRNA
#
			my $loc = "$First{$key}\t$Last{$key}\t.\t$Strand{$key}";
			print GFF "$seqid\t$source\t$type\t$loc\t.\tID=$id;Parent=$parent;Name=$name;Index=1\n";
#
#			exons
#
			my $dis_type = 'exon';
			foreach my $loc (@{ $Use_loc{$id}{$seqid} })
			{
				print GFF "$seqid\t$source\t$dis_type\t$loc\t.\tParent=$id\n";
			}
		}
	}

	elsif ($type eq 'CDS')
	{
#		Comment = Has parent 
#		protein is the First/Last of the CDS
#		CDS is each of the ranges individually
#		UTRs - for each exon, the segment not part of the CDS
#		5' UTR - The UTRs before the CDS start
#		3' UTR - The UTRs after the CDS stop

		foreach my $seqid (keys(%{ $Use_loc{$id} }))
		{
			$key = "$id:$seqid";
			my $first      = $First{$key};
			my $last       = $Last{$key};

#
#			Protein
#
			my $dis_type = 'protein';
			my $loc = "$First{$key}\t$Last{$key}\t.\t$Strand{$key}";
			print GFF "$seqid\t$source\t$dis_type\t$loc\t.\tID=$parent-Protein;Name=$parent;Derives_from=$parent\n";

#
#			5' UTR
#
			my $mrna = "$parent:$seqid";
#print "DEBUG:PARENT=$parent SEQID=$seqid MRNA=$mrna \n";
			foreach my $loc (@{ $Use_loc{$parent}{$seqid} })
			{
				my ($left,$right,$score,$strand) = split(/	/,$loc);
#print "DEBUG:LOC=$loc KEY=$key FIRST=$First{$key} LAST=$Last{$key} LEFT=$left RIGHT=$right\n";
#					Contained within CDS
				next if ($left >= $First{$key} && $right <= $Last{$key} );  
#print "DEBUG:\t\tFIRST=$First{$key} LAST=$Last{$key} LEFT=$left RIGHT=$right\n";
#					Exon Left of the CDS
				if ($right < $First{$key})
				{
					$dis_type = ($strand eq '+') ? 'five_prime_UTR' : 'three_prime_UTR';
					my $loc = "$left\t$right\t.\t$Strand{$mrna}";
					print GFF "$seqid\t$source\t$dis_type\t$loc\t.\tParent=$parent\n";
				}
#					Exon Overlaps Left end of the CDS
				elsif ($left <= $First{$key})
				{
					my $right = $First{$key} - 1;
					$dis_type = ($strand eq '+') ? 'five_prime_UTR' : 'three_prime_UTR';
					my $loc = "$left\t$right\t.\t$Strand{$mrna}";
					print GFF "$seqid\t$source\t$dis_type\t$loc\t.\tParent=$parent\n";
				}
#					Exon Right of the CDS
				if ($left > $Last{$key})
				{
					$dis_type = ($strand eq '+') ? 'three_prime_UTR' : 'five_prime_UTR';
					my $loc = "$left\t$right\t.\t$Strand{$mrna}";
					print GFF "$seqid\t$source\t$dis_type\t$loc\t.\tParent=$parent\n";
				}
#					Exon Overlaps right end of the CDS
				elsif ($right >= $Last{$key})
				{
					my $left = $Last{$key} + 1;
					$dis_type = ($strand eq '+') ? 'three_prime_UTR' : 'five_prime_UTR';
					my $loc = "$left\t$right\t.\t$Strand{$mrna}";
					print GFF "$seqid\t$source\t$dis_type\t$loc\t.\tParent=$parent\n";
				}
			}

#
#			CDS
#
			my $dis_type = 'CDS';
			my $phase = 0;
			foreach my $loc (@{ $Use_loc{$id}{$seqid} })
			{
				print GFF "$seqid\t$source\t$dis_type\t$loc\t$phase\tParent=$parent,$parent-Protein;\n";
				my ($left,$right,$score,$strand) = split(/	/,$loc);
#$length = $right - $left + 1;
#print "DEBUG:\tPHASE=$phase LEFT=$left RIGHT=$right LENGTH=$length\n";
				$phase = ($phase + 2 + $left - $right) % 3 ;
#print "DEBUG:\t\tPHASE=$phase LEFT=$left RIGHT=$right LENGTH=$length\n";
			}

		}
	}
	else
	{
		print "Did not find $buf\n";
	}
}

close FH;
close GFF;

