BEGIN
{
 @INC = (@INC, '/user/gernotf/knappstuff/exe');
}

use strict;
use PDB;

{ 
package State;

sub titratableResidues
{

my $id = 
{
 'NTE'  => ['R','P'],
 'LYS'  => ['R','P'],
 'ARG'  => ['R','P'],
 'HSP'  => ['R','P','0'],
 'GLU'  => ['R','D'],
 'ASP'  => ['R','D'],
 'CTE'  => ['R','D'],
 'PXA'  => ['R','D'],
 'PXD'  => ['R','D'],
 'CYS'  => ['R','D'],
 'TYR'  => ['R','D'],
 'SEP'  => ['R','D'],
 'DPP'  => ['R','D','0'],
 'EPP'  => ['R','D','0'],
 'PAP'  => ['R','D','0'],
 'PDP'  => ['R','D','0']
};

my $pka = 
{
 'NTE'  => [0.00, -7.50],
 'LYS'  => [0.00, -10.40],
 'ARG'  => [0.00, -12.00],
 'HSP'  => [0.00, -7.00, -0.40],
 'GLU'  => [0.00, 4.40],
 'ASP'  => [0.00, 4.00],
 'PXA'  => [0.00, 4.80],
 'PXD'  => [0.00, 4.80],
 'CTE'  => [0.00, 3.80],
 'CYS'  => [0.00, 9.50],
 'TYR'  => [0.00, 9.60],
 'SEP'  => [0.00, 15.70],
 'DPP'  => [0.00, 4.00, 0.00],
 'EPP'  => [0.00, 4.40, 0.00],
 'PAP'  => [0.00, 4.80, 0.00],
 'PDP'  => [0.00, 4.80, 0.00]
};

my $charge = 
{
 'NTE' =>
 {
  'R' => {'N'=>-0.97, 'HT1'=>0.22, 'HT2'=>0.22, 'HT3'=>0.22},
  'P' => {'N'=>-0.30, 'HT1'=>0.33, 'HT2'=>0.33, 'HT3'=>0.33}
 },
 'LYS' =>
 {
  'R' => {'NZ'=>-0.97, 'HZ1'=>0.22, 'HZ2'=>0.22, 'HZ3'=>0.22},
  'P' => {'NZ'=>-0.30, 'HZ1'=>0.33, 'HZ2'=>0.33, 'HZ3'=>0.33}
 },
 'ARG' =>
 {
  'R' => {'NE'=>-0.81, 'HE'=>0.44, 'CZ'=>0.71, 'NH1'=>-0.90, 'HH11'=>0.27, 'HH12'=>0.27, 'NH2'=>-0.90, 'HH21'=>0.27, 'HH22'=>0.27},
  'P' => {'NE'=>-0.70, 'HE'=>0.44, 'CZ'=>0.64, 'NH1'=>-0.80, 'HH11'=>0.46, 'HH12'=>0.46, 'NH2'=>-0.80, 'HH21'=>0.46, 'HH22'=>0.46}
 },
 'HSP' =>
 {
  'R' => {'NE2'=>-0.70, 'CG'=>-0.05, 'ND1'=>-0.36, 'HD1'=>0.32, 'CE1'=>0.25, 'HE2'=>0.00, 'CD2'=>0.22, 'CB'=>-0.09, 'HD2'=>0.10, 'HE1'=>0.13}, # delta
  'P' => {'NE2'=>-0.51, 'CG'=>0.19, 'ND1'=>-0.51, 'HD1'=>0.44, 'CE1'=>0.32, 'HE2'=>0.44, 'CD2'=>0.19, 'CB'=>-0.05, 'HD2'=>0.13, 'HE1'=>0.18},  # prot
  '0' => {'NE2'=>-0.36, 'CG'=>0.22, 'ND1'=>-0.70, 'HD1'=>0.00, 'CE1'=>0.25, 'HE2'=>0.32, 'CD2'=>-0.05, 'CB'=>-0.08, 'HD2'=>0.09, 'HE1'=>0.13}  # epsilon 
 },
 'GLU' => 
 {
  'R' => {'CG'=>-0.21, 'CD'=>0.75, 'OE1'=>-0.36, 'OE2'=>-0.36}, 
  'D' => {'CG'=>-0.28, 'CD'=>0.62, 'OE1'=>-0.76, 'OE2'=>-0.76}
 },
 'ASP' => 
 {
  'R' => {'CB'=>-0.21, 'CG'=>0.75, 'OD1'=>-0.36, 'OD2'=>-0.36}, 
  'D' => {'CB'=>-0.28, 'CG'=>0.62, 'OD1'=>-0.76, 'OD2'=>-0.76}
 },
 'CTE' =>
 {
  'R' => {'C'=>0.34, 'OT1'=>-0.17, 'OT2'=>-0.17},
  'D' => {'C'=>0.34, 'OT1'=>-0.67, 'OT2'=>-0.67}
  },
 'PXA' =>
 {
 'R' => {'CBA'=>-0.21, 'CGA'=>0.75, 'O1A'=>-0.36, 'O2A'=>-0.36},
 'D' => {'CBA'=>-0.28, 'CGA'=>0.62, 'O1A'=>-0.76, 'O2A'=>-0.76}
 },
 'PXD' =>
 {
  'R' => {'CBD'=>-0.21, 'CGD'=>0.75, 'O1D'=>-0.36, 'O2D'=>-0.36},
  'D' => {'CBD'=>-0.28, 'CGD'=>0.62, 'O1D'=>-0.76, 'O2D'=>-0.76}
 },
 'CYS' =>
 {
  'R' => {'CB'=>-0.11, 'SG'=>-0.23, 'HG1'=>0.16},
  'D' => {'CB'=>-0.25, 'SG'=>-0.93, 'HG1'=>0.00}
 },
#  'TYR' =>
#  {
#   'R' => {'CZ'=>0.11, 'OH'=>-0.54, 'HH'=>0.43},
#	'D' => {'CZ'=>-0.18, 'OH'=>-0.82, 'HH'=>0.00}
#  }
 'TYR' =>
 {
  'R' => {'CG'=>0.000, 'CD1'=>-0.115, 'HD1'=>0.115, 'CE1'=>-0.115, 'HE1'=>0.115, 'CZ'=>0.110, 'OH'=>-0.540, 'HH'=>0.430, 'CE2'=>-0.115, 'HE2'=>0.115, 'CD2'=>-0.115, 'HD2'=>0.115},
  'D' => {'CG'=>-0.341, 'CD1'=>0.028, 'HD1'=>0.072, 'CE1'=>-0.525, 'HE1'=>0.124, 'CZ'=>0.769, 'OH'=>-0.826, 'HH'=>0.00, 'CE2'=>-0.525, 'HE2'=>0.124, 'CD2'=>0.028, 'HD2'=>0.072}
 },
 'SEP' =>
 {
  'R' => {'CB'=> 0.050, 'OG'=>-0.660, 'HG1'=>0.430},
  'D' => {'CB'=>-0.480, 'OG'=>-0.700, 'HG1'=>0.000}
 },
 'DPP' => 
 {
  'R' => {'CB'=>-0.21, 'CG'=>0.75, 'OD1'=>-0.61, 'OD2'=>-0.55, 'HD1'=>0.44, 'HD2'=>0.00},
  'D' => {'CB'=>-0.28, 'CG'=>0.62, 'OD1'=>-0.76, 'OD2'=>-0.76, 'HD1'=>0.00, 'HD2'=>0.00},
  '0' => {'CB'=>-0.21, 'CG'=>0.75, 'OD1'=>-0.55, 'OD2'=>-0.61, 'HD1'=>0.00, 'HD2'=>0.44}
 },
 'EPP' => 
 {
  'R' => {'CG'=>-0.21, 'CD'=>0.75, 'OE1'=>-0.61, 'OE2'=>-0.55, 'HE1'=>0.44, 'HE2'=>0.00},
  'D' => {'CG'=>-0.28, 'CD'=>0.62, 'OE1'=>-0.76, 'OE2'=>-0.76, 'HE1'=>0.00, 'HE2'=>0.00},
  '0' => {'CG'=>-0.21, 'CD'=>0.75, 'OE1'=>-0.55, 'OE2'=>-0.61, 'HE1'=>0.00, 'HE2'=>0.44}
 },
 'PAP' =>
 {
  'R' => {'CBA'=>-0.21, 'CGA'=>0.75, 'O1A'=>-0.61, 'O2A'=>-0.55, 'HO1A'=>0.44, 'HO2A'=>0.00},
  'D' => {'CBA'=>-0.28, 'CGA'=>0.62, 'O1A'=>-0.76, 'O2A'=>-0.76, 'HO1A'=>0.00, 'HO2A'=>0.00},
  '0' => {'CBA'=>-0.21, 'CGA'=>0.75, 'O1A'=>-0.55, 'O2A'=>-0.61, 'HO1A'=>0.00, 'HO2A'=>0.44}
 },
 'PDP' =>
 {
  'R' => {'CBD'=>-0.21, 'CGD'=>0.75, 'O1D'=>-0.61, 'O2D'=>-0.55, 'HO1D'=>0.44, 'HO2D'=>0.00},
  'D' => {'CBD'=>-0.28, 'CGD'=>0.62, 'O1D'=>-0.76, 'O2D'=>-0.76, 'HO1D'=>0.00, 'HO2D'=>0.00},
  '0' => {'CBD'=>-0.21, 'CGD'=>0.75, 'O1D'=>-0.55, 'O2D'=>-0.61, 'HO1D'=>0.00, 'HO2D'=>0.44}
 }
};

return ($id, $pka, $charge)
}

sub getChargedStates
{
 my ($class, $name, $pqr) = @_;
 my ($id, $pka, $charge)  = $class->titratableResidues;

 return () unless( exists($id->{$name}) );
 
 my @ids = @{$id->{$name}};
 my @pks = @{$pka->{$name}};
 my @states;
 
 for(my $i=0; $i<@ids; $i++)
 {
  my $state;
  my $alist = PQR->new;
  foreach my $atom ( @{$pqr->all_atoms} )
  {
   my $c = $charge->{$name}->{$ids[$i]}->{$atom->atom_name};
	if ( defined($c) )
	{
	 my $acopy = PQR_atom->copy($atom);
    $acopy->charge = $c;
	 $alist->add_atom($acopy);
	}
  }
  $state = {'id'=>$ids[$i], 'pka'=>$pks[$i], 'alist'=>$alist}; bless $state, $class;
  push @states, $state;
 }
 
 return @states;
}

sub id : lvalue
{
 $_[0]->{'id'}
}

sub pka : lvalue
{
 $_[0]->{'pka'}
}

sub alist
{
 $_[0]->{'alist'};
}

} # class State

{
package SingleConfTitratable;

sub new
{
 my ($class, $name, $resid, $segid, $pqr) = @_;
 my @states;
 if($segid)
 {
  @states = State->getChargedStates($name, $pqr->get_segment($segid)->get_resid_fast($resid));
 }
 else
 {
  @states = State->getChargedStates($name, $pqr->get_resid_fast($resid));
 }
 map {$_->alist->undefCoordinates} @states;
 map {$_->alist->undefRadii} @states;
 my $self = {'name'=>$name, 'resid'=>$resid, 'segid'=>$segid, '@states'=>\@states};
 
 bless $self, $class;
}

sub getTitratableResidues
{
 my ($class, $pqr) = @_;
 my @titratables;

 foreach my $segid ( @{$pqr->segment_names} )
 {
  my $seg     = $pqr->get_segment($segid);
  my ($n, $c) = $class->NorCTerminus($seg);
  my %disu    = %{$seg->getDISU}; %disu = (%disu, reverse %disu);
  foreach my $resid ( @{$seg->residue_numbers} )
  {
   next if( exists($disu{$resid}) ); # do not titrate disulfide bridges
   my $res    = $seg->get_resid_fast($resid);
   my $name   = $res->all_atoms->[0]->residue_name;
   my @states = State->getChargedStates($name, $res);
	if(@states)
	{
    map {$_->alist->undefCoordinates} @states;
    map {$_->alist->undefRadii} @states;
    my $titr   = {'name'=>$name, 'resid'=>$resid, 'segid'=>$segid, '@states'=>\@states};
    bless $titr, $class;
	 push @titratables, $titr;
	}
  }
  # N and C terminus
  if($n!=-1)
  {
   my $res  = $seg->get_resid_fast($n);
	my $name = 'NTE';
	my @states = State->getChargedStates($name, $res);
   map {$_->alist->undefCoordinates} @states;
   map {$_->alist->undefRadii} @states;
   my $titr   = {'name'=>$name, 'resid'=>$n, 'segid'=>$segid, '@states'=>\@states};
   bless $titr, $class;
   push @titratables, $titr;
  }
  if($c!=-1)
  {
   my $res  = $seg->get_resid_fast($c);
	my $name = 'CTE';
	my @states = State->getChargedStates($name, $res);
   map {$_->alist->undefCoordinates} @states;
   map {$_->alist->undefRadii} @states;
   my $titr   = {'name'=>$name, 'resid'=>$c, 'segid'=>$segid, '@states'=>\@states};
   bless $titr, $class;
   push @titratables, $titr;
  }
 }

 return @titratables;
}

sub states
{
 @{$_[0]->{'@states'}};
}

sub name
{
 $_[0]->{'name'}
}

sub resid 
{
 $_[0]->{'resid'}
}

sub segid
{
 $_[0]->{'segid'}
}

sub NorCTerminus
{
 my ($class, $seg) = @_;

 my $n     = $seg->all_atoms->[0]->resid;
 my $nname = $seg->all_atoms->[0]->residue_name;
 my $c     = ($seg->all_atoms->[$#{$seg->all_atoms}])->resid;
 my $cnt   = 0;
 
 # detect N-terminus
 foreach my $atom ( @{$seg->get_resid($n)->all_atoms} )
 {
  $cnt++ if ( $atom->atom_name eq 'HT1' );
  $cnt++ if ( $atom->atom_name eq 'HT2' );
  $cnt++ if ( $atom->atom_name eq 'HT3' );
  $cnt++ if ( $atom->atom_name eq 'N' );
 }
 if($cnt<4)
 {
  $n = -1;
 }
 elsif($nname eq 'PRO')
 {
  $n = -1;
 }

 # detect C-terminus
 $cnt = 0;
 foreach my $atom ( @{$seg->get_resid($c)->all_atoms} )
 {
  $cnt++ if ( $atom->atom_name eq 'OT1' );
  $cnt++ if ( $atom->atom_name eq 'OT2' );
 }
 $c = -1 if($cnt<2);

 return ($n, $c);
}

sub switchToState
{
 my ($self, $idx, $pqr) = @_;

 my @states = $self->states;
 my $segid  = $self->segid;
 my $resid  = $self->resid;

 die "SingleConfTitratable::switchToState: Unknown state $idx" if( ($idx<0) || ($idx>=@states) );

 my %hash = map {$_->atom_name, $_} @{$states[$idx]->alist->all_atoms};
 my $res;
 if($segid)
 {
  $res = $pqr->get_segment($segid)->get_resid_fast($resid)->all_atoms;
 }
 else
 {
  $res = $pqr->get_resid_fast($resid)->all_atoms;
 }
 foreach my $atom ( @{$res} )
 {
  $atom->charge = $hash{$atom->atom_name}->charge if( exists($hash{$atom->atom_name}) );
 }
}

sub weightedAverage
{
 my $self   = shift;
 my $pqr    = shift;
 my @probs  = @_;

 my @states = $self->states;
 my $segid  = $self->segid;
 my $resid  = $self->resid;
 
 my $res;
 if($segid)
 {
  $res = $pqr->get_segment($segid)->get_resid_fast($resid);
 }
 else
 {
  $res = $pqr->get_resid_fast($resid);
 }
 
 # linear interpolation between state charges according to
 # user supplied occupancies
 my %hash;
 for(my $i=0; $i<@states; $i++)
 {
  foreach my $atom ( @{$states[$i]->alist->all_atoms} )
  {
   $hash{$atom->atom_name} += $probs[$i] * $atom->charge;
  }
 }
 
 # set interpolated charges
 foreach my $atom ( @{$res->all_atoms} )
 {
  $atom->charge = $hash{$atom->atom_name} if( exists($hash{$atom->atom_name}) );
 }
}

sub mostProbable
{
 my $self   = shift;
 my $pqr    = shift;
 my @probs  = @_;

 my @states = $self->states;
 my $segid  = $self->segid;
 my $resid  = $self->resid;
 
 die "SingleConfTitratable::mostProbable: #states @states != #probs @probs" unless( @states==@probs );
 
 # take largest probability
 my $max = [-1,0];
 for(my $i=0; $i<@states; $i++)
 {
  $max = [$i,$probs[$i]] if( $probs[$i] > $$max[1] );
 }
 
 # switch to state with highest probability
 $self->switchToState($$max[0], $pqr);
}

sub stringify
{
 my $self  = shift;
 my $str = '';

 foreach my $state ($self->states)
 {
  $str .= sprintf "%.2f pK %s\n%s", $state->pka, $state->id, $state->alist->toString;
 }
 
 return $str;
}


} # class SingleConfTitratable
1;
