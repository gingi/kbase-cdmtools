use strict;
use Test::Simple tests => 1;
use Text::Diff;
use File::Temp;

my $SCRIPTDIR  = './scripts';
my $FIXTUREDIR = './fixtures';


my $expfile = "$FIXTUREDIR/features1.tsv";
my $fh = File::Temp->new();
my $outfile = $fh->filename;

system "$SCRIPTDIR/gff3_to_kbase.pl --gff $FIXTUREDIR/features1.gff --output $outfile";
my $diff = diff $expfile, $outfile;
ok($diff eq '', "Output does not match expected: $diff");