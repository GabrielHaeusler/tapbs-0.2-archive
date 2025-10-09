use strict;
use Cartesian;

package PDB;

use overload '""' => 'stringify',
             '==' => 'equals';

our %CODE = (
            'ALA'=>'A',
            'ARG'=>'R',
            'RPP'=>'R',
            'ASP'=>'D',
            'DPP'=>'D',
            'ASN'=>'N',
            'CYS'=>'C',
            'GLU'=>'E',
            'EPP'=>'E',
            'GLN'=>'Q',
            'GLY'=>'G',
            'HIS'=>'H',
            'ILE'=>'I',
            'LEU'=>'L',
            'LYS'=>'K',
            'MET'=>'M',
            'PHE'=>'F',
            'PRO'=>'P',
            'SER'=>'S',
            'THR'=>'T',
            'TRP'=>'W',
            'TYR'=>'Y',
            'VAL'=>'V',
            'HSP'=>'H'
           );

sub new
{
 my $CLASS = shift;
 my $struc;
 
 $struc->{'atoms'}  = [];
 $struc->{'resids'} = {};
 $struc->{'icodes'} = {};
 $struc->{'chains'} = {};
 $struc->{'segnames'} = {};
 $struc->{'confs'}  = {};
 $struc->{'coors'}  = {};

 bless $struc, $CLASS;
}

sub copy
{
 my $CLASS = shift;
 my $pdb   = shift;
 my $copy  = PDB->new;

 foreach my $atom (@{$pdb->all_atoms})
 {
  my $atomCopy = PDB_atom->copy($atom);
  $copy->add_atom($atomCopy);
 }
 
 $copy->update;
 
 return $copy;
}

sub parseFromFile
{
 my $self     = shift;
 my $filename = shift;
 my $noupdate = shift;

 $self = $self->new unless( ref($self) ); # if class provided behave as constructor

 open(IN, "< $filename") or die "Could not open file $filename!";
 
 foreach my $line (<IN>)
 {
  my $atom = PDB_atom->new;
  if ( $atom->parseFromString($line) )
  {
   $self->add_atom($atom);
  }
 }
 
 close IN;
 
 # update hash tables for resids, chains and conformers
 $self->update unless($noupdate);
 
 return $self;
}

sub parseFromNMRFile
{
 my $class    = shift;
 my $filename = shift;

 open(IN, "< $filename") or die "Could not open file $filename!";
 
 my @pdbs;
 my $pdb;
 
 foreach my $line (<IN>)
 {
  next unless($line =~ /^(ATOM|HETATM|MODEL|ENDMDL)/);
  if ($line =~ /MODEL/)
  {
   $pdb = $class->new;
   next;
  }
  if ($line =~ /ENDMDL/)
  {
   $pdb->update;
   push @pdbs, $pdb;
	next;
  }
  my $atom = PDB_atom->new;
  if ( $atom->parseFromString($line) )
  {
   $pdb->add_atom($atom);
  }
 }
 
 close IN;
 
 return @pdbs;
}

sub parseFromCRDFile
{
 my $self     = shift;
 my $filename = shift;

 open(IN, "< $filename") or die "Could not open file $filename!";
 
 foreach my $line (<IN>)
 {
  my $atom = PDB_atom->new();
  if ( $atom->parseFromCRDString($line) )
  {
   $self->add_atom($atom);
  }
 }
 
 close IN;
 
 # update hash tables for resids, chains and conformers
 $self->update();
}

sub empty
{
 my $self = shift;
 (@{$self->{'atoms'}}>0)?0:1;
}

sub add_atom
{
 my $self     = shift;
 my $pdb_atom = shift;
 
 push @{$self->{'atoms'}}, $pdb_atom;
}

sub add_atoms
{
 my $self      = shift;
 my @pdb_atoms = @_;
 
 foreach my $atom (@pdb_atoms)
 {
  push @{$self->{'atoms'}}, $atom;
 }
}

sub update
{
# adds hash tables for resids, chains and conformers
 my $self = shift;

 # sort atoms
 my @atoms = sort byAtomID @{$self->{'atoms'}};
 $self->{'atoms'} = \@atoms;

 # remove old hash tables
 $self->{'resids'} = {};
 $self->{'icodes'} = {};
 $self->{'chains'} = {};
 $self->{'segnames'} = {};
 $self->{'confs'}  = {};
 $self->{'coors'}  = {};

 foreach my $atom ( @{$self->{'atoms'}} )
 {
  my $resid       = $atom->resid;
  my $icode       = $atom->icode;
  my $chain       = $atom->chain_name;
  my $conf        = $atom->conformer_id;
  my $segid       = $atom->segment_name;
  my ($x, $y, $z) = @{$atom->coordinates};
  push @{$self->{'resids'}->{$resid}},    $atom if (defined($resid));
  push @{$self->{'icodes'}->{$icode}},    $atom unless($icode eq '');
  push @{$self->{'chains'}->{$chain}},    $atom unless($chain eq '');
  push @{$self->{'segnames'}->{$segid}},  $atom unless($segid eq '');
  push @{$self->{'confs'}->{$conf}},      $atom unless($conf eq '');
  push @{$self->{'coors'}->{"$x,$y,$z"}}, $atom;
 }
 
 return $self;
}

sub renumber
{
 my ($self, $roffset, $aoffset) = @_;
 my ($atomid, $resid) = ( defined($aoffset)?$aoffset:0, defined($roffset)?$roffset:0 );
 my ($oldchain, $oldresid) = (0, 0);

 foreach my $atom (@{$self->{'atoms'}})
 {
  $atomid++;
  if ( ($atom->resid != $oldresid) || ($atom->chain_name ne $oldchain) )
  {
   $resid++;
	$oldresid = $atom->resid;
	$oldchain = $atom->chain_name;
  }
  $atom->atomid = $atomid;
  $atom->resid = $resid;
 }
 $self->update;
}

sub renumberBySegid
{
 my ($self) = @_;
 my ($atomid, $resid)    = (0, 0);
 my ($oldseg, $oldresid) = (0, 0);

 foreach my $segid (sort keys %{$self->{'segnames'}})
 {
  foreach my $atom (@{$self->{'segnames'}->{$segid}})
  {
   $atomid++;
   if ( ($atom->resid != $oldresid) )
   {
    $resid++;
 	 $oldresid = $atom->resid;
   }
	if ( $atom->segment_name ne $oldseg )
	{
	 $resid = 1;
	 $oldseg = $atom->segment_name;
	}
   $atom->atomid = $atomid;
   $atom->resid = $resid;
  }
 }
 $self->update;
}

sub renumberByChain
{
 my ($self, $roffset, $aoffset) = @_;
 my ($atomid, $resid) = ( defined($aoffset)?$aoffset:0, defined($roffset)?$roffset:0 );
 my ($oldchain, $oldresid) = (0, 0);

 foreach my $c (sort keys %{$self->{'chains'}})
 {
  foreach my $atom (@{$self->{'segnames'}->{$c}})
  {
   $atomid++;
   if ( $atom->resid != $oldresid )
   {
    $resid++;
 	 $oldresid = $atom->resid;
   }
	if ( $atom->chain_name ne $oldchain )
	{
	 $resid = 1;
	 $oldchain = $atom->chain_name;
	}
   $atom->atomid = $atomid;
   $atom->resid = $resid;
  }
 }
 $self->update;
}

sub toString
{
 my $self = shift;
 my $str = '';
 
 foreach my $atom (@{$self->{'atoms'}})
 {
  $str .= sprintf( "%s\n", $atom);
 }
 #$str .= sprintf("TER\n");
 #$str .= sprintf("END");
 
 return $str;
}

sub stringify
{
 $_[0]->toString;
}

sub toStringCRD
{
 my $self = shift;
 my $str = '';
 my $acnt = 0;
 my $rcnt = 0;
 my $oldresid = 0;
 
 # title
 $str .= "* PERL generated CRD file\n";
 $str .= "* by conversion of PDB file.\n";
 $str .= "* \n";
 
 # nr of atoms
 $str .= sprintf("%5i \n", $self->nrOfAtoms);
 
 # c ATOMNO RESNO RES TYPE X Y Z SEGID RESID Weighting
 # c I5 I5 1X A4 1X A4 F10.5 F10.5 F10.5 1X A4 1X A4 F10.5
 foreach my $atom ( @{$self->all_atoms} )
 {
  if($oldresid ne $atom->resid)
  {
   $rcnt++;
   $oldresid = $atom->resid;
  }
  $str .= sprintf
          (
           "%5i%5i %-4s %-4s%10.5f%10.5f%10.5f %-4s%5i%10.5f\n",
			  ++$acnt, $rcnt, $atom->residue_name, $atom->atom_name, 
			  $atom->coordinates->[0], $atom->coordinates->[1], $atom->coordinates->[2],
			  $atom->segment_name, $atom->resid, 0
          );
 }
 
 return $str;
}

sub all_atoms
{
 my $self = shift;
 
 return $self->{'atoms'};
}

sub coordinates
{
 my $self = shift;
 
 return $self->{'coors'};
}

sub append
{
 my $self = shift;
 my $other= shift;
 
 foreach my $atom (@{$other->all_atoms})
 {
  $self->add_atom($atom);
 }
 
 $self->update;
}

sub removeAtomName
{
 my $self = shift;
 my %names = map {$_, 1} @_;
 my @atoms;
 
 foreach my $atom (@{$self->all_atoms})
 {
  push @atoms, $atom unless( exists $names{$atom->atom_name} );
 }

 $self->{'atoms'} = \@atoms;
 $self->update(); 
}

sub removeAtomPositions
{
 my $self      = shift;
 my $positions = shift;
 my @atoms     = ();

 foreach my $atom (@{$self->all_atoms})
 {
  my ($x, $y, $z) = @{$atom->coordinates};
  push @atoms, $atom unless( exists($positions->{"$x,$y,$z"}) );
 }
 
 $self->{'atoms'} = \@atoms;
 
 $self->update();
}

sub residue_numbers
{
 my $self  = shift;
 my $resids= [];
 foreach my $resid ( sort {$a <=> $b} keys %{$self->{'resids'}} )
 {
  push @$resids, $resid;
 }
 
 return $resids;
}

sub chains
{
 my $self  = shift;
 my $chains= [];
 
 foreach my $chain ( sort keys %{$self->{'chains'}} )
 {
  push @$chains, $chain;
 }
 
 return $chains;
}

sub segment_names
{
 my $self  = shift;
 my $segids = [];
 
 foreach my $segment ( sort keys %{$self->{'segnames'}} )
 {
  push @$segids, $segment;
 }
 
 return $segids;
}

sub conformers
{
 my $self  = shift;
 my $confs = [];
 
 foreach my $conf ( sort keys %{$self->{'confs'}} )
 {
  push @$confs, $conf;
 }
 
 return $confs;
}

sub icodes
{
 my $self  = shift;
 my $icodes = [];
 
 foreach my $icode ( sort keys %{$self->{'icodes'}} )
 {
  push @$icodes, $icode;
 }
 
 return $icodes;
}

sub get_resid
{
 my $self  = shift;
 my $resid = shift;

 my $atoms = $self->{'resids'}->{$resid} or die "Unknown residue number: $resid!";
 my $pdb   = PDB->new();

 foreach my $atom ( @{$atoms} )
 {
  $pdb->add_atom( $atom );
 }

 $pdb->update();

 return $pdb;
}

# returns residue atoms in incomplete PDB instance
sub get_resid_fast
{
 my $self  = shift;
 my $resid = shift;
 my $pdb   = PDB->new;
 
 $pdb->{'atoms'} = $self->{'resids'}->{$resid};
 
 return $pdb;
}

sub get_chain
{
 my $self  = shift;
 my $chain = shift;

 my $atoms = $self->{'chains'}->{$chain} or die "Unknown chain identifier: $chain!";
 my $pdb   = PDB->new();

 foreach my $atom ( @{$atoms} )
 {
  $pdb->add_atom( $atom );
 }

 $pdb->update();

 return $pdb;
}

sub get_segment
{
 my $self  = shift;
 my $segid = shift;

 my $atoms = $self->{'segnames'}->{$segid} or die "Unknown segment name: $segid!";
 my $pdb   = PDB->new();

 foreach my $atom ( @{$atoms} )
 {
  $pdb->add_atom( $atom );
 }

 $pdb->update();

 return $pdb;
}

sub get_conformer
{
 my $self = shift;
 my $conf = shift;
 
 my $atoms = $self->{'confs'}->{$conf} or die "Unknown conformer identifier: $conf!";
 my $pdb   = PDB->new();

 foreach my $atom ( @{$atoms} )
 {
  $pdb->add_atom( $atom );
 }

 $pdb->update();
 
 return $pdb;
}

sub get_icode
{
 my $self = shift;
 my $icode = shift;
 
 my $atoms = $self->{'icodes'}->{$icode} or die "Unknown insertion code: $icode!";
 my $pdb   = PDB->new();

 foreach my $atom ( @{$atoms} )
 {
  $pdb->add_atom( $atom );
 }

 $pdb->update();
 
 return $pdb;
}

sub size
{
 my $self  = shift;
 my @min  = @{$self->all_atoms->[0]->coordinates};
 my @max  = @{$self->all_atoms->[0]->coordinates};
 my @diff = (0,0,0);
 
 foreach my $atom (@{$self->all_atoms})
 {
  my ($x, $y, $z) = @{$atom->coordinates};
  $min[0] = $x if ($x<$min[0]);
  $min[1] = $y if ($y<$min[1]);
  $min[2] = $z if ($z<$min[2]);
  $max[0] = $x if ($x>$max[0]);
  $max[1] = $y if ($y>$max[1]);
  $max[2] = $z if ($z>$max[2]);
 }
 
 @diff = ($max[0]-$min[0], $max[1]-$min[1], $max[2]-$min[2]);
 
 return (\@min, \@max, \@diff);
}

# special functions

sub equals
{
 my $self = shift;
 my $other = shift;
 
 #print "self != other because nr of atoms differ!\n" if( @{$self->all_atoms}!=@{$other->all_atoms} );
 return 0 if( @{$self->all_atoms}!=@{$other->all_atoms} );

 for( my $i=0; $i<@{$self->all_atoms}; $i++ )
 {
  #print "self != other because coordinates differ!\n" unless($self->all_atoms->[$i]->coordinates->equal($other->all_atoms->[$i]->coordinates));
  return 0 unless($self->all_atoms->[$i]->coordinates->equal($other->all_atoms->[$i]->coordinates));
 }
 
 return 1;
}

sub rmsd
{
 my $self   = shift;
 my $other  = shift;
 my $rmsd   = 0.00;
 
 my @atoms1 = @{ $self->all_atoms };
 my @atoms2 = @{ $other->all_atoms };
 
 die "PDB::rmsd: Arguments should have the same number of atoms!" unless( @atoms1==@atoms2 );

 for(my $i=0; $i<@atoms1; $i++)
 {
  my $diff = $atoms1[$i]->coordinates-$atoms2[$i]->coordinates;
  $rmsd   += $diff*$diff;
 }
 
 $rmsd /= @atoms1;
 
 return sqrt($rmsd);
}

# @brief Like map, applies callback function to all atoms
sub fold
{
 my $self     = shift;
 my $callback = shift;
 
 foreach my $atom (@{$self->all_atoms})
 {
  &$callback($atom);
 }
}

# @brief Calculates the common set of atoms with same coordinates from structures.
# assumes that there no identical atoms within structures.
sub intersection
{
 my @pdbs  = @_;
 my %coors;
 my $intersection = PDB->new;

 # useless if only one pdb given
 return $intersection unless(@pdbs>1);
 
 # take all coordinates and store atoms
 foreach my $pdb (@pdbs)
 {
  $pdb->fold( sub {
               my $atom = shift;
               my ($x, $y, $z) = @{$atom->coordinates};
               push @{$coors{"$x,$y,$z"}}, $atom;
              }
            );
 }

 # take only hash positions that contain as many atoms as structures in @pdbs
 while( my ($coor, $array) = each %coors )
 {
  if (@$array==@pdbs)
  {
   foreach my $atom (@$array)
	{
    $intersection->add_atom($atom);
	}
  }
 }
 
 $intersection->update();
 
 return $intersection;
}

# useful macros

sub writeToFile
{
 my $self     = shift;
 my $filename = shift;
 open(OUT, "> $filename") or die "PDB::writeToFile: Could not open file $filename!";
 print OUT $self->toString;
 print OUT "END\n";
 close OUT;
}

sub writeToCRDFile
{
 my $self     = shift;
 my $filename = shift;
 open(OUT, "> $filename") or die "PDB::writeToCRDFile: Could not open file $filename!";
 print OUT $self->toStringCRD;
 close OUT;
}

sub writeChainsToFiles
{
 my $self     = shift;
 my $filestem = shift;
 
 foreach my $chainID ( keys %{$self->{'chains'}} )
 {
  my $chain = $self->get_chain($chainID);
  open(OUT, "> chain$chainID.$filestem");
  print OUT $chain->toString();
  close OUT;
 }
}

sub writeConformersToFiles
{
 my $self     = shift;
 my $filestem = shift;

 foreach my $chainID ( keys %{$self->{'chains'}} )
 {
  my $chain = $self->get_chain($chainID);
  foreach my $confID ( @{$chain->conformers()} )
  {
   next if $confID eq '_';
   my $conf = $chain->get_conformer($confID);
   foreach my $resID ( @{$conf->residue_numbers()} )
   {
    open(OUT, "> conf$confID.res$resID.chain$chainID.$filestem");
    print OUT $conf->get_resid($resID)->toString();
    close OUT;
   }
  }
 }
}

sub getSequence
{
 my $self = shift;
 my $chain = shift;
    $chain = $self->chains->[0] unless($chain);

 my @residues;
 my @sequence;
 my $seq;
 
 # collect sequence
 $self->get_chain($chain)->fold
 ( 
  sub
  {
   $seq->{ $_[0]->resid } = $_[0]->residue_name; #if ( exists($CODE{$_[0]->residue_name}) );
  }
 );

 # fill sequence gaps
 my @temp = sort {$a<=>$b} keys %$seq;
 @residues = (1 .. pop @temp);
 foreach my $resid (@residues)
 {
  if ( exists($seq->{$resid}) )
  {
   push @sequence, $seq->{$resid};
  }
  else
  {
   push @sequence, 'DUM';
  }
 }
 
 return \@sequence;
}

sub removeNonAAResidues
{
 my $self = shift;
 my $pdb = PDB->new;

 $self->fold
 (
  sub
  {
   $pdb->add_atom($_[0]) if( exists($CODE{$_[0]->residue_name}) );
  }
 );
 
 $self->{'atoms'} = $pdb->{'atoms'};
 $self->update;
}

sub filter_multiple_coordinates
{
 my $self = shift;
 my $hash = {};
 my $pdb  = PDB->new;

 foreach my $atom (@{$self->all_atoms})
 {
  my ($x, $y, $z) = @{$atom->coordinates};
  push @{$hash->{"$x,$y,$z"}}, $atom;
 }
 
 while(my ($coor, $atoms) = each %$hash)
 {
  $pdb->add_atom($atoms->[0]);
 }
 
 $pdb->update();
 
 return $pdb;
}

sub undefCoordinates
{
 my $self = shift;

 foreach my $atom (@{$self->all_atoms})
 {
  $atom->coordinates = Cartesian->new(9999.999, 9999.999, 9999.999);
 }
 
 $self->update;
}

sub getDISU
{
 my $self = shift;
 my %args = 
 (
  -cutoff=> 2.5,
  -chain=>($self->chains)->[0],
  @_
 );
 my $chain = $args{-chain};
 my $cutoff= $args{-cutoff};
 
 my $disu = {};

 # first filter cys sulfurs
 my $sulfurs = PDB->new;
 foreach my $atom ( @{$self->get_chain($chain)->all_atoms} )
 {
  $sulfurs->add_atom($atom) if ( ($atom->residue_name eq 'CYS') and ($atom->atom_name eq 'SG') );
 }
 $sulfurs->update;

 # check pair distances
 for(my $i=0; $i<@{$sulfurs->all_atoms}; $i++ )
 {
  for(my $j=$i+1; $j<@{$sulfurs->all_atoms}; $j++)
  {
   my $s1 = $sulfurs->all_atoms->[$i];
	my $s2 = $sulfurs->all_atoms->[$j];
	my $resid1 = $s1->resid;
	my $resid2 = $s2->resid;
	my ($x1, $y1, $z1) = @{ $s1->coordinates };
	my ($x2, $y2, $z2) = @{ $s2->coordinates };
	my ($d1, $d2, $d3) = ( $x1-$x2, $y1-$y2, $z1-$z2 );
	my $r = sqrt($d1*$d1 + $d2*$d2 + $d3*$d3);
	$disu->{$resid1} = $resid2 if($r<=$cutoff);
  }
 }
 
 return $disu;
}

sub getTermini
{

 # returns resids for N- and C-terminus.
 # C-terminus is only assumed to be present
 # if OT1 and OT2 atoms exist (O/OXT).

 my $pdb   = shift;
 my $chain = shift;

 # get chain
 $chain = $pdb->chains->[0] unless($chain);
 $pdb   = $pdb->get_chain($chain);

 # get N-Terminus
 my $nter = $pdb->all_atoms->[0]->resid;

 # check C-terminus
 my $cnt = 0;
 my $cter = ($pdb->all_atoms->[$#{$pdb->all_atoms}])->resid;
 foreach my $atom ( @{$pdb->get_resid($cter)->all_atoms} )
 {
  $cnt++ if ( $atom->atom_name eq 'OT1' or $atom->atom_name eq 'O' );
  $cnt++ if ( $atom->atom_name eq 'OT2' or $atom->atom_name eq 'OXT' );
 }
 if( $cnt<2 )
 {
  $cter = -1;
 }

 return ($nter, $cter);
}

sub nrOfAtoms
{
 my $self = shift;
 
 return scalar @{$self->all_atoms};
}

# --------------------
# private methods
# --------------------

sub byAtomID
{
 $a->atomid <=> $b->atomid;
}

package PQR;

use base qw(PDB);

sub copy
{
 my $CLASS = shift;
 my $pqr   = shift;
 my $copy  = PQR->new;

 foreach my $atom (@{$pqr->all_atoms})
 {
  my $atomCopy = PQR_atom->copy($atom);
  $copy->add_atom($atomCopy);
 }
 
 $copy->update;
 
 return $copy;
}

sub parseFromFile
{
 my $self     = shift;
 my $filename = shift;
 my $noupdate = shift;
 
 $self = $self->new unless( ref($self) ); # if class provided behave as constructor
 my $lineNr = 0;

 open(IN, "< $filename") or die "Could not open file $filename!";
 
 foreach my $line (<IN>)
 {
  my $atom = PQR_atom->new;
  $lineNr++;
  if ( $atom->parseFromString($line) )
  {
   $self->add_atom($atom);
  }
  #else
  #{
  # die "PDB::parseFromFile: Parse error in $filename at line $lineNr";
  #}
 }

 close IN;
 
 # update hash tables for resids, chains and conformers
 $self->update unless($noupdate);
 
 return $self;
}

sub parseFromFileNonStrict
{
 my $self     = shift;
 my $filename = shift;

 open(IN, "< $filename") or die "Could not open file $filename!";
 
 foreach my $line (<IN>)
 {
  my $atom = PQR_atom->new();
  my @temp = split(' ', $line);
  $atom->atomid        = $temp[1];
  $atom->atom_name     = $temp[2];
  $atom->residue_name  = $temp[3];
  $atom->resid         = $temp[4];
  $atom->coordinates   = Cartesian->new($temp[5], $temp[6], $temp[7]);
  $atom->charge        = $temp[8];
  $atom->radius        = $temp[9];
  $atom->chain_name    = ' ';
  $atom->segment_name  = ' ';
  $atom->element_name  = ' ';
  $atom->conformer_id  = '_';
  $self->add_atom($atom);
 }

 close IN;
 
 # update hash tables for resids, chains and conformers
 $self->update();
}

sub fromCharmmData
{
 my ($self, $psf, $crd, $rtf, $prm) = @_;

 # build PQR structure
 for(my $i=0; $i<@$psf; $i++)
 {
  my $psfatom  = $psf->[$i];
  my $coord    = $crd->[$i] or die "Numbers of atoms in PSF and CRD do not match!";

  # find radius
  my $radius = find_radius($psfatom, $rtf, $prm);

  # define PQR atom 
  my $pqr_atom = PQR_atom->new();
  $pqr_atom->atomid             = $psfatom->{'atomid'};
  $pqr_atom->atom_name          = $psfatom->{'name'};
  $pqr_atom->residue_name       = $psfatom->{'residue'};
  $pqr_atom->resid              = $psfatom->{'resid'};
  $pqr_atom->chain_name         = $psfatom->{'segid'};
  $pqr_atom->segment_name       = $psfatom->{'segid'};
  $pqr_atom->occupancy          = $psfatom->{'charge'};
  $pqr_atom->temperature_factor = $radius;
  $pqr_atom->coordinates        = Cartesian->new(@$coord);
  $pqr_atom->element_name       = ' ';
  $pqr_atom->conformer_id       = ' ';
 
  # add atom to structure
  $self->add_atom($pqr_atom);
 }

 $self->update();
}

sub parseFromCharmmFiles
{
 my $self = shift;

 my ($atomid, $oldresid, $resid) = (0, 0, 0);

 # read & parse CHARMM files
 my $psf  = read_psf(shift);
 my $crd  = read_crd(shift);
 my $rtf  = read_rtf(shift);
 my $prm  = read_prm(shift);

 $self->fromCharmmData($psf, $crd, $rtf, $prm);
}

sub parseFromPDBCharmmFiles
{
 my $self = shift;
 my $pdb  = shift;

 my ($atomid, $oldresid, $resid) = (0, 0, 0);

 # get coordinates
 my $crd = [];
 foreach my $pdb_atom (@{$pdb->all_atoms})
 {
  push @$crd, $pdb_atom->coordinates;
 }

 # read & parse CHARMM files
 my $psf  = read_psf(shift);
 my $rtf  = read_rtf(shift);
 my $prm  = read_prm(shift);

 $self->fromCharmmData($psf, $crd, $rtf, $prm);
}

sub totalCharge
{
 my $self   = shift;
 my $charge = 0.0;
 
 foreach my $atom ( @{$self->all_atoms} )
 {
  $charge += $atom->charge;
 }
 
 return $charge;
}

sub coulomb
{
 my $self = shift;
 my $eps  = shift;

 # electrostatic constants
 my $pi    = 3.141592653589793238;
 my $ec    = 1.6021773e-19; # As
 my $e0    = 8.8541878176E-12; # As^2/Jm
 my $Angstrom2Meter = 1E-10; # m/A
 my $const = $ec*$ec/(4*$pi*$e0*$eps)/$Angstrom2Meter;
 my $NA    = 6.0221367e+23; # mol^-1

 my ($i, $j, $atomi, $atomj);
 my ($r1, $r2, $r);
 my ($q1, $q2);
 my $E = 0;

 for($i=0; $i<@{$self->{'atoms'}}; $i++)
 {
  $atomi = $self->{'atoms'}->[$i];
  #$r1 = $atomi->coordinates;
  #$q1 = $atomi->charge;
  $r1 = $atomi->{'coord'};
  $q1 = $atomi->{'occup'};
  for($j=$i+1; $j<@{$self->{'atoms'}}; $j++)
  {
   $atomj = $self->{'atoms'}->[$j];
   #$r2 = $atomj->coordinates;
	#$q2 = $atomj->charge;
	$r2  = $atomj->{'coord'};
   $q2  = $atomj->{'occup'};
	#$r  = ($r1-$r2)->value;
	$r = sqrt( ($r1->[0]-$r2->[0])*($r1->[0]-$r2->[0]) + ($r1->[1]-$r2->[1])*($r1->[1]-$r2->[1]) + ($r1->[2]-$r2->[2])*($r1->[2]-$r2->[2]) );
	$E += $q1*$q2/$r;
  }
 }
 
 return $const*$NA*$E/1000; # kJ/mol;
}

sub undefCharges
{
 my $self = shift;

 foreach my $atom (@{$self->all_atoms})
 {
  $atom->charge = 99.999;
 }
}

sub undefRadii
{
 my $self = shift;

 foreach my $atom (@{$self->all_atoms})
 {
  $atom->radius = 99.999;
 }
}

sub equals
{
 my $self = shift;
 my $other = shift;
 
 #print "PQR: self != other because nr of atoms differ!\n" if( @{$self->all_atoms}!=@{$other->all_atoms} );
 return 0 if( @{$self->all_atoms}!=@{$other->all_atoms} );

 for( my $i=0; $i<@{$self->all_atoms}; $i++ )
 {
  #print "PQR: self != other because coordinates differ!\n" unless($self->all_atoms->[$i]->coordinates->equal($other->all_atoms->[$i]->coordinates));
  #print "PQR: self != other because charges differ!\n" unless($self->all_atoms->[$i]->charge==$other->all_atoms->[$i]->charge);
  return 0 unless($self->all_atoms->[$i]->charge==$other->all_atoms->[$i]->charge);
  return 0 unless($self->all_atoms->[$i]->coordinates->equal($other->all_atoms->[$i]->coordinates));
 }
 
 return 1;
}

sub toStringAPBS
{
 my $self = shift;
 my $str = '';
 
 foreach my $pdb_atom (@{$self->{'atoms'}})
 {
  $str .= sprintf( "%s\n", $pdb_atom->stringify_APBS);
 }
 
 return $str;
}

sub writeToFileAPBS
{
 my $self = shift;
 my $name = shift;
 
 open(OUT, "> $name") or die "Could not write to file $name";
  print OUT $self->toStringAPBS;
 close OUT;
}

# private helper methods

sub read_psf
{
 my $filename = shift;
 my $flag = 0;
 my $psf;
 
 open(PSF, "< $filename") or die "Could not open file $filename!";
 
 while(my $line = <PSF>)
 {
  last if ( $line =~ /!NBOND/ );
  if ($line =~ /!NATOM/)
  {
   $flag = 1;
	next;
  }
  my @temp = split(' ', $line);
  if ( @temp > 0 && $flag )
  {
   my $atom = {'atomid'  => $temp[0],
               'segid'   => $temp[1],
               'resid'   => $temp[2],
               'residue' => $temp[3],
					'name'    => $temp[4],
					'massid'  => $temp[5],
					'charge'  => sprintf("%6.3f", $temp[6])};
   push @{$psf}, $atom;
  }
 }
 
 close PSF;
 
 return $psf;
}

sub read_crd
{
 my $filename = shift;
 my $crd;
 my ($x, $y, $z);
 
 open(CRD, "< $filename") or die "Could not open file $filename!";
 
 while(my $line = <CRD>)
 {
  next if ($line =~ /^\*/); # ignore title
  if ( length($line)>60 )
  {
   $x = substr($line, 20, 10);
   $y = substr($line, 30, 10);
   $z = substr($line, 40, 10);
   # convert to 3 digits after comma
   push @{$crd}, [sprintf("%8.3f", $x), sprintf("%8.3f", $y), sprintf("%8.3f", $z)]; # x, y, z coordinate
  }
 }
 
 close CRD;
 
 return $crd;
}

sub read_rtf
{
 my $filename = shift;
 my $atomtype;
 
 open(RTF, "< $filename") or die "Could not open file $filename!";
 
 while(my $line = <RTF>)
 {
  if ( $line =~ /^MASS/ )
  {
   my @temp = split(' ', $line);
	$atomtype->{$temp[1]} = $temp[2];
  }
 }
 
 close RTF;
 
 return $atomtype;
}

sub read_prm
{
 my $filename = shift;
 my $flag = 0;
 my $radius;
 
 open(PRM, "< $filename") or die "Could not open file $filename!";
 
 while(my $line = <PRM>)
 {
  $flag = 1 if ($line =~ /^NONBONDED/);
  last if ($line =~ /^HBOND/);
  
  if ( $flag && $line =~ /^[A-Z0-9]{1,4}\s*\d\.\d*/ )
  {
   my @temp = split(' ', $line);
	$radius->{$temp[0]} = $temp[3];
	#print STDERR "Radius for atom type ".$temp[0]." is ".$temp[3]."\n";
  }
 }
 
 close PRM;
 
 return $radius;
}

sub find_radius
{
 my ($psfatom, $rtf, $prm) = @_;

 my ($massname, $radius) = ('', 0.0);

 # get massid of current atom
 my $massid   = $psfatom->{'massid'};

 # translate massid into massname
 if ( exists($rtf->{$massid}) )
 {
  $massname = $rtf->{$massid};
 }
 else
 {
  die "Could not find atom type with mass id $massid!";
 }

 #  use massname to get radius from prm-file
 if ( exists($prm->{$massname}) )
 {
  $radius = $prm->{$massname};
 }
 else
 {
  print STDERR "Could not find radius for atom type $massname! Using 0.00 A\n";
 }
 
 return $radius;
}

package PDB_atom;
use overload '""' => \&stringify;

sub new
{
 my $CLASS = shift;
 my $self  = {};
 
 $self = {'atmid' => 0,
          'aname' => '',
          'rname' => '',
          'resid' => 0,
			 'icode' => '',
          'chain' => '',
          'segna' => '',
          'coord' => [],
          'occup' => 0.0,
          'bfact' => 0.0,
          'elemn' => '',
          'confm' => ''
         };
 
 bless $self, $CLASS;
}

sub copy
{
 my $CLASS = shift;
 my $atom  = shift;
 my $self  = {};

 $self->{'atmid'} = $atom->atomid;
 $self->{'aname'} = $atom->atom_name;
 $self->{'rname'} = $atom->residue_name;
 $self->{'resid'} = $atom->resid;
 $self->{'icode'} = $atom->icode;
 $self->{'chain'} = $atom->chain_name;
 $self->{'segna'} = $atom->segment_name;
 #$self->{'coord'} = [@{$atom->coordinates}];
 $self->{'coord'} = Cartesian->new(@{$atom->coordinates});
 $self->{'occup'} = $atom->occupancy;
 $self->{'bfact'} = $atom->temperature_factor;
 $self->{'elemn'} = $atom->element_name;
 $self->{'confm'} = $atom->conformer_id;

 bless $self, $CLASS;
}

sub parseFromString
{
 my $self   = shift;
 my $string = shift;
 my $temp   = '';

 if ($string =~ /^(ATOM|HETATM)/)
 {
  my $atomno   = sprintf( "%u", substr($string, 6, 5) );
  my ($atomna) = split( ' ', substr($string, 12, 4) );
  my $resnam   = substr($string, 17, 3);
  my $resid    = sprintf( "%u", substr($string, 22, 4) );
  my $icode    = substr($string, 26, 1);
  my ($chain)  = split( ' ', substr($string, 21, 1) ); $chain = defined($chain)?$chain:'A';
  my ($segna)  = split( ' ', substr($string, 72, 4) ); $segna = defined($segna)?$segna:' ';

  $temp        = substr($string, 30, 8);
  my $x        = sprintf( "%8.3f", ($temp ne '********')?$temp:'9999.999' );

  $temp        = substr($string, 38, 8);
  my $y        = sprintf( "%8.3f", ($temp ne '********')?$temp:'9999.999' );

  $temp        = substr($string, 46, 8);
  my $z        = sprintf( "%8.3f", ($temp ne '********')?$temp:'9999.999' );

  my $occup    = sprintf( "%6.2f", substr($string, 54, 6) );
  my $bfac     = sprintf( "%6.2f", substr($string, 60, 6) );
  my ($elem)   = split( ' ', substr($string, 76, 2) ); $elem = defined($elem)?$elem:substr($atomna,0,1);
  my ($conf)   = split( ' ', substr($string, 16, 1) ); $conf = defined($conf)?$conf:' ';
  $self->{'atmid'} = $atomno;
  $self->{'aname'} = $atomna;
  $self->{'rname'} = $resnam;
  $self->{'resid'} = $resid;
  $self->{'icode'} = $icode;
  $self->{'chain'} = $chain;
  $self->{'segna'} = $segna;
  #$self->{'coord'} = [$x, $y, $z];
  $self->{'coord'} = Cartesian->new($x, $y, $z);
  $self->{'occup'} = $occup;
  $self->{'bfact'} = $bfac;
  $self->{'elemn'} = $elem;
  $self->{'confm'} = $conf;
  return 1;
 }
 else
 {
  return 0;
 }
}

sub parseFromCRDString
{
 my $self   = shift;
 my $string = shift;
 my @temp;

 return 0 if ($string =~ /^\*/); # ignore title
 
 # sloppy parsing with split
 @temp = split(' ', $string);
 
 return 0 unless (@temp>1);
 
 $self->{'atmid'} = $temp[0];
 $self->{'aname'} = $temp[3];
 $self->{'rname'} = $temp[2];
 $self->{'resid'} = $temp[1];
 $self->{'chain'} = 'A';
 $self->{'segna'} = $temp[7];
 #$self->{'coord'} = [$temp[4], $temp[5], $temp[6]];
 $self->{'coord'} = Cartesian->new($temp[4], $temp[5], $temp[6]);
 $self->{'occup'} = 1.00;
 $self->{'bfact'} = 0.00;
 $self->{'elemn'} = ' ';
 $self->{'confm'} = ' ';
 
 return 1;
}

sub atomid : lvalue
{
 my $self   = shift;
 
 $self->{'atmid'};
}

sub resid : lvalue
{
 my $self   = shift;

 $self->{'resid'};
}

sub icode : lvalue
{
 my $self   = shift;

 $self->{'icode'};
}

sub conformer_id : lvalue
{
 my $self   = shift;

 $self->{'confm'};
}

sub coordinates : lvalue
{
 my $self   = shift;
 
 $self->{'coord'};
}

sub atom_name : lvalue
{
 my $self   = shift;

 $self->{'aname'};
}

sub element_name : lvalue
{
 my $self   = shift;
 
 $self->{'elemn'};
}

sub residue_name : lvalue
{
 my $self   = shift;
 
 $self->{'rname'};
}

sub chain_name : lvalue
{
 my $self   = shift;
 
 $self->{'chain'};
}

sub segment_name : lvalue
{
 my $self   = shift;
 
 $self->{'segna'};
}

sub occupancy : lvalue
{
 my $self = shift;
 
 $self->{'occup'};
}

sub temperature_factor : lvalue
{
 my $self = shift;
 
 $self->{'bfact'};
}

sub stringify
{
 my $self  = shift;
 my $str   = '';

 if($self->{'aname'} !~ /^\dH/ && length($self->{'aname'})<4) # NMR format issue
 {
  $str = sprintf("ATOM  %5u %-4s%1.1s%3.3s %1.1s%4u%1.1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4.4s%2.2s", $self->{'atmid'}, " ".($self->{'aname'}), $self->{'confm'}, $self->{'rname'}, $self->{'chain'}, $self->{'resid'}, $self->{'icode'}, $self->{'coord'}->[0], $self->{'coord'}->[1], $self->{'coord'}->[2], $self->{'occup'}, $self->{'bfact'}, $self->{'segna'}, $self->{'elemn'});
 }
 else
 {
  $str = sprintf("ATOM  %5u %-4s%1.1s%3.3s %1.1s%4u%1.1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4.4s%2.2s", $self->{'atmid'}, $self->{'aname'}, $self->{'confm'}, $self->{'rname'}, $self->{'chain'}, $self->{'resid'}, $self->{'icode'}, $self->{'coord'}->[0], $self->{'coord'}->[1], $self->{'coord'}->[2], $self->{'occup'}, $self->{'bfact'}, $self->{'segna'}, $self->{'elemn'});
 }

 return $str;
}

package PQR_atom;

use base qw(PDB_atom);
use overload '""' => \&stringify;

sub new
{
 my $CLASS = shift;
 my $self  = {};
 
 $self = {'atmid' => 0,
          'aname' => '',
          'rname' => '',
          'resid' => 0,
			 'icode' => '',
          'chain' => '',
          'segna' => '',
          'coord' => [],
          'occup' => 0.0,
          'bfact' => 0.0,
          'elemn' => '',
          'confm' => ''
         };

 bless $self, $CLASS;
}

# everything the same as in PDB_atom except
# that charges and radii are parsed in format
# %6.3f

sub parseFromString
{
 my $self   = shift;
 my $string = shift;
 my $temp   = '';

 if ($string =~ /^(ATOM|HETATM)/)
 {
  my $atomno   = sprintf( "%u", substr($string, 6, 5) );
  my ($atomna) = split( ' ', substr($string, 12, 4) );
  my $resnam   = substr($string, 17, 3);
  my $resid    = sprintf( "%u", substr($string, 22, 4) );
  my $icode    = substr($string, 26, 1);
  my ($chain)  = split( ' ', substr($string, 21, 1) ); $chain = defined($chain)?$chain:'A';
  my ($segna)  = split( ' ', substr($string, 72, 4) ); $segna = defined($segna)?$segna:' ';

  $temp        = substr($string, 30, 8);
  my $x        = sprintf( "%8.3f", ($temp ne '********')?$temp:'9999.999' );

  $temp        = substr($string, 38, 8);
  my $y        = sprintf( "%8.3f", ($temp ne '********')?$temp:'9999.999' );

  $temp        = substr($string, 46, 8);
  my $z        = sprintf( "%8.3f", ($temp ne '********')?$temp:'9999.999' );

  $temp        = substr($string, 54, 6);
  #my ($occup)  = split( ' ', substr($string, 54, 6) );
  my ($occup)  = sprintf("%6.3f", $temp);
  $temp        = substr($string, 60, 6);
  #my ($bfac)   = split( ' ', substr($string, 60, 6) );
  my ($bfac)   = sprintf("%6.3f", $temp);
  my ($elem)   = split( ' ', substr($string, 76, 2) ); $elem = defined($elem)?$elem:substr($atomna,0,1);
  my ($conf)   = split( ' ', substr($string, 16, 1) ); $conf = defined($conf)?$conf:' ';
  $self->{'atmid'} = $atomno;
  $self->{'aname'} = $atomna;
  $self->{'rname'} = $resnam;
  $self->{'resid'} = $resid;
  $self->{'icode'} = $icode;
  $self->{'chain'} = $chain;
  $self->{'segna'} = $segna;
  #$self->{'coord'} = [$x, $y, $z];
  $self->{'coord'} = Cartesian->new($x, $y, $z);
  $self->{'occup'} = $occup;
  $self->{'bfact'} = $bfac;
  $self->{'elemn'} = $elem;
  $self->{'confm'} = $conf;

  return 1;
 }
 else
 {
  return 0;
 }
}

sub charge : lvalue
{
 my $self = shift;
 $self->{'occup'};
}

sub radius : lvalue
{
 my $self = shift;
 $self->{'bfact'};
}

# everything the same as in PDB_atom except
# that charges and radii are parsed in format
# %5.3f instead of %6.2f

sub stringify
{
 my $self  = shift;
 my $str   = '';
 if($self->{'aname'} !~ /^\dH/ && length($self->{'aname'})<4) # NMR format issue
 {
  $str = sprintf("ATOM  %5u %-4s%1.1s%3.3s %1.1s%4u%1.1s   %8.3f%8.3f%8.3f%6.3f%6.3f      %-4.4s%2.2s", $self->{'atmid'}, " ".($self->{'aname'}), $self->{'confm'}, $self->{'rname'}, $self->{'chain'}, $self->{'resid'}, $self->{'icode'}, $self->{'coord'}->[0], $self->{'coord'}->[1], $self->{'coord'}->[2], $self->{'occup'}, $self->{'bfact'}, $self->{'segna'}, $self->{'elemn'});
 }
 else
 {
  $str = sprintf("ATOM  %5u %-4s%1.1s%3.3s %1.1s%4u%1.1s   %8.3f%8.3f%8.3f%6.3f%6.3f      %-4.4s%2.2s", $self->{'atmid'}, $self->{'aname'}, $self->{'confm'}, $self->{'rname'}, $self->{'chain'}, $self->{'resid'}, $self->{'icode'}, $self->{'coord'}->[0], $self->{'coord'}->[1], $self->{'coord'}->[2], $self->{'occup'}, $self->{'bfact'}, $self->{'segna'}, $self->{'elemn'});
 }
 #$str = sprintf("ATOM  %5u %-4s%1.1s%3.3s %1.1s%4u    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4.4s%2.2s", $self->{'atmid'}, $self->{'aname'}, $self->{'confm'}, $self->{'rname'}, $self->{'chain'}, $self->{'resid'}, $self->{'coord'}->[0], $self->{'coord'}->[1], $self->{'coord'}->[2], $self->{'occup'}, $self->{'bfact'}, $self->{'segna'}, $self->{'elemn'});

 return $str;
}

sub stringify_APBS
{
 my $self  = shift;
 my $str   = '';
 $str = sprintf("ATOM  %5u %-4s %3.3s %4u    %8.3f %8.3f %8.3f %6.3f %6.3f      %-4.4s %2.2s", $self->{'atmid'}, $self->{'aname'}, $self->{'rname'}, $self->{'resid'}, $self->{'coord'}->[0], $self->{'coord'}->[1], $self->{'coord'}->[2], $self->{'occup'}, $self->{'bfact'}, $self->{'segna'}, $self->{'elemn'});

 return $str;
}

1;
