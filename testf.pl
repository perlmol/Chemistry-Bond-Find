#!/home/ivan/bin/perl -w

use strict;
use Chemistry::Mol;
use Chemistry::File::Mopac;
use blib;
use Chemistry::Bond::Find 'find_bonds';

my $mol = Chemistry::Mol->read($ARGV[0] || "02.zt");
print $mol->formula, "\n";

find_bonds($mol);
for my $bond ($mol->bonds) {
    print $bond, ": ", $bond->atoms, "\n";
}
