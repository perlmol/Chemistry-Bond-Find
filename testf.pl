#!/home/ivan/bin/perl -ws
$|=1;
use strict;
use Chemistry::Mol;
use Chemistry::File::Mopac;
use Chemistry::File::MDLMol;
use Chemistry::File::PDB;
use blib;
use Chemistry::Bond::Find;
use Benchmark ':all';

if ($a) {
    my $mol = Chemistry::Mol->read($ARGV[0] || "test1717.pdb");
    Chemistry::Bond::Find::find_bonds($mol, margin => $b || 1.5);
    printf "found %d bonds\n", scalar($mol->bonds);
    exit;
}

my @files = glob "test*.pdb";

for my $file (@files) {
    print "FILENAME: $file\n";
    my $mol = Chemistry::Mol->read($file);
    my $n = $mol->atoms;
    my $times = int(5000/$n) || 1;    
    cmpthese($times,
        {
            clone => sub { $mol->clone },
            rec => sub { Chemistry::Bond::Find::find_bonds_rec($mol->clone) },
            grid => sub { Chemistry::Bond::Find::find_bonds($mol->clone) },
        },
    );
    print "\n\n";
}
