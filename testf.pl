#!/home/ivan/bin/perl -w

use strict;
use Chemistry::Mol;
use Chemistry::File::Mopac;
use Chemistry::File::MDLMol;
use Chemistry::File::PDB;
use blib;
use Chemistry::Bond::Find 'find_bonds';

my $mol = Chemistry::Mol->read($ARGV[0] || "test750.pdb");
#print $mol->formula, "\n";

#find_bonds($mol);
#for my $bond ($mol->bonds) { print $bond, ": ", $bond->atoms, "\n"; }
#print $mol->print(format => "mdl");

Chemistry::Bond::Find::find_bonds($mol);
$mol->write("rec.mol", format => 'mdl');
