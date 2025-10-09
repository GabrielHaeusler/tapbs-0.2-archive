#!/usr/bin/perl -w

use strict;
use PDB;
use ST;

unless(@ARGV)
{
 print "Usage: pqr2sites.pl <PQR filename>\n";
 exit 0;
}

my $file = shift @ARGV;

# SITES file name
my $sitesfile = 'prot.sites';
   $sitesfile = "$1.sites" if ( $file =~ /(.*)\.+.*/ );

# parse pqr file
my $pqr  = PQR->parseFromFile($file);
#my $pqr  = PQR->parseFromFileNonStrict($file);
buildST($pqr, $sitesfile);

exit(1);

sub buildST
{
 my ($pqr, $sitesfile) = @_;
 my @titratables = SingleConfTitratable->getTitratableResidues($pqr);

 # write SITES
 open(SITES, "> $sitesfile") or die "buildST: Could not open file $sitesfile";
 foreach my $t (@titratables)
 {
  my $segid = $t->segid;
  my $resid = $t->resid;
  my $name  = $t->name;
  printf(SITES "%s %i %s %s.st\n", $segid, $resid, $name, $name);
 }
 close SITES;

 # write ST files
 my %hash = map {$_->name, $_} @titratables;
 foreach my $name (keys %hash)
 {
  my $stfile = "$name.st";
  open(ST, "> $stfile") or die "buildST: Could not open file $stfile";
   print ST $hash{$name}->stringify; 
  close ST;
 }
}
