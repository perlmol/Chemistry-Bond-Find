package Chemistry::Bond::Find;

$VERSION = '0.05';
# $Id$

=head1 NAME

Chemistry::Bond::Find - Detect bonds in a molecule from atomic 3D coordinates

=head1 SYNOPSIS

    use Chemistry::Bond::Find 'find_bonds';

    find_bonds($mol, %options);

=head1 DESCRIPTION

Detects the bonds in a molecule from its 3D coordinates by using simple
cutoffs.  The current version does not guess the bond orders; all bonds will
have a bond order of 1.

This module is part of the PerlMol project, L<http://www.perlmol.org/>.

=head1 FUNCTIONS

=over

=cut

use strict;
use Chemistry::Mol;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(find_bonds);

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

my $Default_Radius = 1.5; # radius for unknown elements

# I considered inlining this function, but the performance gain was minimal
# ( < 5 % ), so it's probably better to leave it here
# $opts->{cuttof} Hash Table of cutoffs. 
# Example: key = "H C", value = [0.76, 1.42]
sub are_bonded {
    my ($a, $b, $r, $opts) = @_;
    return $r < ($opts->{cuttoffs}{"$a $b"} ||= 
        (($Covalent_Radius{$a} || $Default_Radius) 
         + ($Covalent_Radius{$b} || $Default_Radius)) * $opts->{tolerance});
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



=item find_bonds($mol, %options)

Find and add the bonds in a molecule. Only use it in molecules that have no 
explicit bonds; for example, after reading a file with 3D coordinates but no
bond orders.

Available options:

=over

=item tolerance

Defaults to 1.1. Two atoms are considered to be bound if the distance between
them is less than the sum of their covalent radii multiplied by the tolerance.

=item margin

NOTE: in general setting this option is not recommended, unless you know what
you are doing. Used by the recursive partitioning algorithm when "stitching"
the interfaces between partitions. It defaults to 2 * Rmax * tolerance, where
Rmax is the largest covalent radius among the elements found in the molecule.
For example, if a molecule has C, H, N, O, and I, Rmax = R(I) = 1.33, so the
margin defaults to 2 * 1.33 * 1.1 = 2.926. This margin ensures that no bonds
are missed by the partitioning algorithm. Using a smaller value gives faster
results, but at the risk of missing some bonds. In this exaple, if you are
certain that your molecule doesn't contain I-I bonds (but it has C-I bonds),
you can set margin to (0.77 + 1.33) * 1.1 = 2.31 and you still won't miss any
bonds. All this only has a real impact for molecules with a thousand atoms or
more.

=item min_atoms

Defaults to 20. Partitions with fewer that min_atoms will not be partitioned
further, but their bonds will be found by a simple n^2 method. The default
value seems to be optimal, but you can play with it if you want.

=back

=cut

sub find_bonds {
    my ($mol, %opts) = @_;
    %opts = (min_atoms => 20, tolerance => 1.1, %opts,
        cutoffs => {});  # set defaults
    my $margin = guess_margin($mol, \%opts);    
    %opts = (margin => $margin, %opts);
    _partition($mol, [$mol->atoms], 0, \%opts);
}

sub guess_margin {
    my ($mol, $opts) = @_;
    my $formula = $mol->formula_hash;
    my $max = 0;
    for my $elem (keys %$formula) {
        $max = $Covalent_Radius{$elem} if $Covalent_Radius{$elem} > $max;
    }
    $max *= 2 * $opts->{tolerance};
    #printf "MARGIN guessed at (%.2f)\n", $max;
    $max;
}

sub _partition {
    my ($mol, $atoms, $dir, $opts) = @_;

    #printf "_partition(%s, $dir)\n", scalar(@$atoms);

    if (@$atoms < $opts->{min_atoms}) {
        #print "BOTTOM!\n";
        find_bonds_n2_one_set($mol, $atoms, $opts);
        return;
    }

    my $min = ($atoms->[0]->coords->array)[$dir];
    my $max = $min;
    my $center;
    for my $atom (@$atoms) {
        my $coord = ($atom->coords->array)[$dir]; 
        $center += $coord;
        $min = $coord if $coord < $min;
        $max = $coord if $coord > $min;
    }
    $center /= @$atoms;
    #printf "center($dir)=%.2f; range=(%.2f, %.2f)\n", $center, $min, $max;

    my @left = grep { ($_->coords->array)[$dir] < $center } @$atoms;
    my @right = grep { ($_->coords->array)[$dir] >= $center } @$atoms;
    _partition($mol, \@left, ($dir + 1) % 3, $opts);
    _partition($mol, \@right, ($dir + 1) % 3, $opts);

    # merge the interface between the two halves
    my $margin = $opts->{margin};
    my @left_margin = 
        grep { ($_->coords->array)[$dir] > $center - $margin } @left;
    my @right_margin = 
        grep { ($_->coords->array)[$dir] < $center + $margin } @right;
    find_bonds_n2_two_sets($mol, \@left_margin, \@right_margin, $opts);
}

sub find_bonds_n2_one_set {
    my ($mol, $atoms, $opts) = @_;
    for (my $i = 0; $i < @$atoms; ++$i) {
	for (my $j = $i + 1; $j < @$atoms; ++$j) {
	    my ($a1, $a2) = ($atoms->[$i], $atoms->[$j]);
	    if (are_bonded($a1->symbol, $a2->symbol, scalar $a1->distance($a2), $opts)) {
		$mol->add_bond(Chemistry::Bond->new(atoms=>[$a1, $a2]));
	    }
	}
    }
}

sub find_bonds_n2_two_sets {
    my ($mol, $atoms1, $atoms2, $opts) = @_;
    for my $a1 (@$atoms1) {
	for my $a2 (@$atoms2) {
	    if (are_bonded($a1->symbol, $a2->symbol, scalar $a1->distance($a2), $opts)) {
		$mol->add_bond(Chemistry::Bond->new(atoms=>[$a1, $a2]));
	    }
	}
    }
}

1;

=back

=head1 VERSION

0.05

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::Atom>, L<Chemistry::Bond>,
L<http://www.perlmol.org/>.

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

