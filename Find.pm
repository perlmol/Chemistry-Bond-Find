package Chemistry::Bond::Find;
$VERSION = '0.05';

=head1 NAME

Chemistry::Bond::Find - Detect bonds in a molecule from atomic 3D coordinates

=head1 SYNOPSIS

    use Chemistry::Bond::Find 'find_bonds';

    find_bonds($mol, order => 1);

=cut


use strict;
use Chemistry::Mol;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(find_bonds);

=head1 DESCRIPTION

Detect the bonds in a molecule from its 3D coordinates by using
simple cutoffs. Currently uses a silly brute-force, N^2 algorithm, so
it's probably not a good idea to run it on a protein.

=cut


# table taken from
# http://environmentalchemistry.com/yogi/periodic/covalentradius.html
my %Covalent_Radius = (
    Ag => 1.34, Al => 1.18, Ar => 0.98, As => 1.20, At => 1.45, Au => 1.34,
    B  => 0.82, Ba => 1.98, Be => 0.90, Bi => 1.46, Br => 1.14, C  => 0.77,
    Ca => 1.74, Cd => 1.48, Ce => 1.65, Cl => 0.99, Co => 1.16, Cr => 1.18,
    Cs => 2.35, Cu => 1.17, Dy => 1.59, Er => 1.57, Eu => 1.85, F  => 0.72,
    Fe => 1.17, Ga => 1.26, Gd => 1.61, Ge => 1.22, H  => 0.32, He => 0.93,
    Hf => 1.44, Hg => 1.49, Ho => 1.58, I  => 1.33, In => 1.44, Ir => 1.27,
    K  => 2.03, Kr => 1.12, La => 1.69, Li => 1.23, Lu => 1.56, Mg => 1.36,
    Mn => 1.17, Mo => 1.30, N  => 0.75, Na => 1.54, Nb => 1.34, Nd => 1.64,
    Ne => 0.71, Ni => 1.15, O  => 0.73, Os => 1.26, P  => 1.06, Pb => 1.47,
    Pd => 1.28, Pm => 1.63, Po => 1.46, Pr => 1.65, Pt => 1.30, Rb => 2.16,
    Re => 1.28, Rh => 1.25, Ru => 1.25, S  => 1.02, Sb => 1.40, Sc => 1.44,
    Se => 1.16, Si => 1.11, Sm => 1.62, Sn => 1.41, Sr => 1.91, Ta => 1.34,
    Tb => 1.59, Tc => 1.27, Te => 1.36, Th => 1.65, Ti => 1.32, Tl => 1.48,
    Tm => 1.56, U  => 1.42, V  => 1.22, W  => 1.30, Xe => 1.31, Y  => 1.62,
    Yb => 1.74, Zn => 1.25, Zr => 1.45,
);

my $Tolerance = 0.3; # as a fraction
my $Default_Radius = 1.5;

my %Cutoff; # Table of cutoffs. Example: key = "H C", value = [0.76, 1.42]

sub are_bonded {
    my ($a, $b, $r) = @_;
    my ($min, $max);
    my $pair = "$a $b";
    if($Cutoff{$pair}) {
        ($min, $max) = @{$Cutoff{$pair}};
    } else {
        my $bond_r = ($Covalent_Radius{$a} || $Default_Radius) 
                + ($Covalent_Radius{$b} || $Default_Radius);
        $min = $bond_r * (1-$Tolerance);
        $max = $bond_r * (1+$Tolerance);
        $Cutoff{$pair} = [$min, $max];
    }
    return ($r > $min and $r < $max);
}


sub find_bonds {
    my ($mol, %opts) = @_;
    my @a = $mol->atoms;
    for (my $i = 0; $i < @a; ++$i) {
	for (my $j = $i + 1; $j < @a; ++$j) {
	    my ($a1, $a2) = ($a[$i], $a[$j]);
	    if (are_bonded($a1->symbol, $a2->symbol, $a1->distance($a2))) {
		$mol->add_bond(Chemistry::Bond->new(atoms=>[$a1, $a2]));
	    }
	}
    }
}

# Future work:
# assign hybridization
# assign bond order
# let's think about how to deal with aromaticity...

my %Groups = (
    H => 1, B => 3, C => 4, N => 5, O => 6, F => 7, Cl => 7, => Br => 7, 
    I => 7, Si => 4, P => 5, S => 6, 
);

my %Valencies = (
    H => 1, B => 3, C => 4, N => 3, O => 2, F => 1, Cl => 1, => Br => 1, 
    I => 1, Si => 4, P => 3, S => 2, 
);


# under construction...
sub assign_orders {
    my ($mol, %opts) = @_;

    #assign_unambiguous($mol);
    #assign_ambiguous($mol);
    for my $atom ($mol->atoms) {
        my $charge = 0;
        my $free = $Valencies{$atom->symbol} - $atom->bonds;
        $free = 0 if $free < 0;
        $atom->attr("bond-find/free", $free);
    }

}

1;

=head1 SEE ALSO

Chemistry::Mol, Chemistry::Atom, Chemistry::Bond

=head1 AUTHOR

Ivan Tubert-Brohman <ivan@tubert.org>

=head1 VERSION

$Id$

=cut
