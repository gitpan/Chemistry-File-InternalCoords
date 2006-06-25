package Chemistry::File::InternalCoords;

use warnings;
use strict;

our $VERSION = '0.01';

use base qw(Chemistry::File);

use Chemistry::Mol;
use Chemistry::InternalCoords::Builder 'build_zmat';

my  $EXT = 'zmat';
Chemistry::Mol->register_format( zmat => __PACKAGE__ );

sub name_is {
    my ($class, $fname) = @_;
    $fname =~ /\.$EXT$/i;
}

sub file_is {
    my ($class, $fname) = @_;
    $fname =~ /\.$EXT$/i;
}

sub parse_string {
    my ($class, $s, %opts) = @_;

    my $mol_class  = $opts{mol_class}  || 'Chemistry::Mol';
    my ($s_atoms, $s_vars) = split /^\s*$/m, $s;

    my $mol = $mol_class->new();

    my @lines = split /(?:\n|\r\n?)/, $s_vars;
    my %vars;
    foreach (@lines){
      next unless /(\w+)\s+(-?\d+\.?\d*)/;
      $vars{$1} = $2;
    }

    @lines = split /(?:\n|\r\n?)/, $s_atoms;
    foreach my $i (0..$#lines){
      my ($elem, @internal_coords) = split ' ', $lines[$i];
      $_ = $vars{$_} for @internal_coords[ grep { $_ <= $#internal_coords } 1,3,5 ];
      my $atom = $mol->new_atom(
        parent => $mol,
        id => $i+1,
        ($elem =~ /^\d+$/ ? "Z" : "symbol") => $elem,
        internal_coords => \@internal_coords,
      );
      $atom->internal_coords->add_cartesians;
      $mol->new_bond(
        id => "b$i",
        atoms => [ $atom, $mol->atoms($internal_coords[0]) ],
        length => $internal_coords[1],
      ) if $internal_coords[0];
    }

    return $mol;
}

sub write_string {
    my ($class, $mol, %opts) = @_;

    %opts = (symbol => 1, vars => 1, bfs => 0, sort => 1, %opts);

    build_zmat($mol, %opts);

    my $s = '';
    my ( @B, @A, @D );
    my %index;
    foreach my $i ( 1 .. scalar $mol->atoms ){
        my $atom = $mol->atoms($i);
	$index{ $atom->id } = $i;
        my $ic = $atom->internal_coords;
        my @ic = ($ic->distance, $ic->angle, $ic->dihedral);
        if( $opts{vars} ){
          if(defined $ic[1]){
            push @B, $ic[1];
            $ic[1] = 'B'.scalar(@B);
          }
          if(defined $ic[3]){
            push @A, $ic[3];
            $ic[3] = 'A'.scalar(@A);
          }
          if(defined $ic[5]){
            push @D, $ic[5];
            $ic[5] = 'D'.scalar(@D);
          }
        }else{
          $_ = sprintf "%.8f", $_ for @ic[ grep {defined $ic[$_]} 1,3,5 ];
        }
        $_ = $index{ $_->id } for @ic[ grep {defined $ic[$_]} 0,2,4 ];
        pop @ic while( @ic && !defined $ic[-1] );
        $s .= sprintf "%-2s" . " %5d %15s" x int(@ic/2) . "\n",
                $opts{symbol} ? $atom->symbol : $atom->Z,
                @ic,
        ;
    }
    if( $opts{vars} ){
      $s .= "\n";
      $s .= join "", map { sprintf "  B%-4d %25.8f\n", $_+1, $B[$_] } 0..$#B;
      $s .= join "", map { sprintf "  A%-4d %25.8f\n", $_+1, $A[$_] } 0..$#A;
      $s .= join "", map { sprintf "  D%-4d %25.8f\n", $_+1, $D[$_] } 0..$#D;
    }

    return $s;
}

1;

=pod

=head1 NAME

Chemistry::File::InternalCoords - Internal coordinates (z-matrix) molecule format reader/writer

=head1 VERSION

Version 0.01

=head1 SYNOPSIS

This module is not intended for direct use -- it is intended for use via L<Chemistry::Mol>.

  use Chemistry::File qw/ InternalCoords XYZ /;
  my $mol = Chemistry::Mol->read("foo.zmat");
  warn $mol->print;
  $mol->write(\*STDOUT, format => 'zmat');
  $mol->write('foo.xyz', format => 'xyz');

=head1 DESCRIPTION

This is a subclass of L<Chemistry::File> for reading and writing the zmatriz (aka Z-matrix aka InternalCoords) molecule data format. It registers the 'zmat' file extension with L<Chemistry::Mol>.

For example, here is hydrogen:

  H
  H      1              B1
  
  B1                   0.70000000

and water:

  O
  H      1              B1
  H      1              B2     2              A1
  
  B1                   0.96659987
  B2                   0.96659987
  A1                 107.67114171


=head1 METHODS

This class inherits from L<Chemistry::File>.  The following methods are overloaded:

=over 2

=item name_is

Checks if the filename extension is 'zmat'.

=item file_is

Checks if the filename extension is 'zmat'.

=item parse_string

Expects a plain zmatrix format. Variables are support. No special options.

=item write_string

Creates a plain zmatrix formatted string. Any options are also passed to L<Chemistry::InternalCoords::Builder>'s I<build_zmat> function (defaults to bfs off and sort on).  Also recognizes I<symbol> (default on; if on uses the element instead of the atomic number) and I<vars> (default on; if on uses variables for the bond lengths & angles) options, which affect the output.

=back

=head1 PREREQUISITES

=over 4

=item *

L<Chemistry::Mol>

=item *

L<Chemistry::File>

=item *

L<Chemistry::InternalCoords>

=item *

L<Chemistry::InternalCoords::Builder>

=back

=head1 SEE ALSO

L<http://www.perlmol.org/>

=head1 AUTHOR

David Westbrook, C<< <dwestbrook at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-chemistry-file-internalcoords at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Chemistry-File-InternalCoords>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

I'm also available by email or via '/msg davidrw' on L<http://perlmonks.org>.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Chemistry::File::InternalCoords

You can also look for information at:

=over 4

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Chemistry-File-InternalCoords>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Chemistry-File-InternalCoords>

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Chemistry-File-InternalCoords>

=item * Search CPAN

L<http://search.cpan.org/dist/Chemistry-File-InternalCoords>

=back

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2006 David Westbrook, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

