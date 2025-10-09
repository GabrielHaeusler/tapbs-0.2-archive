#!/usr/bin/perl -w

# Perl script that creates PQR file from CHARMM files.
# 2006 Gernot Kieseritzky

use strict;
use PDB;

# parse command line arguments

die "Usage: create-pqr -psf <psf> -crd <crd> -rtf <mass.rtf> -prm <parm.prm> [-renumber <offset>]\n" if( @ARGV<4 );

my %args = 
(
 -renumber=>-1,
 -prm=>'',
 -rtf=>'',
 -psf=>'',
 -crd=>'',
 @ARGV
);

my $pqr = PQR->new;
$pqr->parseFromCharmmFiles($args{-psf}, $args{-crd}, $args{-rtf}, $args{-prm});
$pqr->renumber( $args{-renumber} ) if ( $args{-renumber}>=0 );
$pqr->fold( sub{my $atom = shift; $atom->chain_name = 'A' } );
print $pqr->toString;

exit 1;
