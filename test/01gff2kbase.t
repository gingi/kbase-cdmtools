use strict;
use Test::Simple tests => 1;
use Text::Diff;
use File::Temp;

my $SCRIPTDIR  = './scripts';
my $FIXTUREDIR = './fixtures';


for my $fixnum (1..1) {
    my $expfile = "$FIXTUREDIR/features$fixnum.tsv";
    my $gfffile = "$FIXTUREDIR/features$fixnum.gff";
    my $fh = File::Temp->new();
    my $outfile = $fh->filename;

    system "$SCRIPTDIR/gff3_to_kbase.pl --gff $gfffile --output $outfile";
    my $diff = diff $expfile, $outfile;
    ok($diff eq '', "Output does not match expected:\n$diff");
}