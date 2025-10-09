#!/usr/bin/perl -w

use strict;

my $prefix = shift;

die "usage: make_curves2.pl 'file prefix'" unless ( defined($prefix) );

my $files     = get_filelist($prefix);
my $conffiles = get_conf_filelist($prefix);

# print files
print "I found following files:\n";
print array2string([@$files, @$conffiles]);
print "\n";

my $sites = get_sites($files->[0]);

# print sites
print "I found following sites from ".$files->[0].":\n";
print array2string([sort sites_by_res keys %$sites]);
print "\n";

# read data
print "Reading data from files ... ";
my ($biased, $error, $rebiased)  = get_data($sites, $files);
my ($confoccup, $conferror)      = get_conf_occup($conffiles);
print "finished.\n";

# print curves
foreach my $t (sort numerically keys %{$biased})
{
 # site files
 foreach my $site (sort sites_by_res keys %$sites)
 {
  foreach my $state (sort numerically keys %{$biased->{$t}->{$site}})
  {
   my $filename = $prefix."_$site"."_$t"."K_s$state.curve";
   print "Creating file $filename ... ";
   open(OUT, "> $filename");
   print OUT combined_curve($site, $biased->{$t}->{$site}->{$state}, $error->{$t}->{$site}->{$state}, $rebiased->{$t}->{$site}->{$state});
   close OUT;
   print "finished.\n";
  }
 }
 # conf occupancies
 if (@$conffiles != 0)
 {
  my $filename = $prefix."_conf_$t"."K.curve";
  print "Creating file $filename ... ";
  open(OUT, "> $filename");
  print OUT conf_plot($confoccup->{$t}, $conferror->{$t});
  close OUT;
  print "finished.\n";
 }
}

# print protonation table
foreach my $t (sort numerically keys %{$biased})
{
 my $filename = $prefix."_$t"."K.table.html";
 print "Writing protonation table into file $filename ... ";
 open(OUT, "> $filename");
 print OUT protonation_table_html($biased->{$t});
 close OUT;
 print "finished.\n";
}

exit 1;

sub combined_curve
{
 my $site     = shift;
 my $biased   = shift;
 my $error    = shift;
 my $rebiased = shift;
 my $string   = '';
 my @phs      = sort numerically keys %$biased;
 my @redoxs   = sort numerically keys %{$biased->{$phs[0]}};

 # print occupancy of current state
 
 foreach my $ph (@phs)
 {
  foreach my $redox (@redoxs)
  {
   my $occup = $biased->{$ph}->{$redox};
	my $err   = $error->{$ph}->{$redox};
   $string .= sprintf("%3.2f %4.1f %1.6f %1.6f\n", $ph, $redox, $occup, $err);
  }
 }

 # pK1/2
 foreach my $redox (@redoxs)
 {
  my ($oldpH, $currpH, $oldOcc, $currOcc);
  foreach my $ph (@phs)
  {
   $currpH  = $ph;
   $currOcc = $biased->{$currpH}->{$redox};
   $oldpH   = $currpH unless ( defined($oldpH) );
   $oldOcc  = $currOcc unless ( defined($oldOcc) );
   if ( (($currOcc<=0.5) && ($oldOcc>0.5)) || (($currOcc>=0.5) && ($oldOcc<0.5)) )
   {
    last;
   }
   else
   {
    $oldpH  = $currpH;
	 $oldOcc = $currOcc;
   }
  }
  my ($slope, $pka);
  if ($oldpH==$currpH)
  {
   $string .= sprintf("# pKa1/2 %s%3.1f @ %3.2f mV\n", (($currOcc-$oldOcc)>=0.0)?'>':'<', $oldpH, $redox);
  }
  else
  {
   $slope = ($currOcc - $oldOcc)/($currpH - $oldpH);
   $pka   = (0.5 - $oldOcc)/$slope + $oldpH;
   $string .= sprintf("# pKa1/2 = %3.1f @ %3.2f mV\n", $pka, $redox);
  }
 }
 
  #Em1/2
 foreach my $ph (@phs)
 {
  my ($oldmV, $currmV, $oldOcc, $currOcc);
  foreach my $mV (@redoxs)
  {
   $currmV  = $mV;
   $currOcc = $biased->{$ph}->{$currmV};
   $oldmV   = $currmV unless ( defined($oldmV) );
   $oldOcc  = $currOcc unless ( defined($oldOcc) );
   if ( (($currOcc<=0.5) && ($oldOcc>0.5)) || (($currOcc>=0.5) && ($oldOcc<0.5)) )
   {
    last;
   }
   else
   {
    $oldmV  = $currmV;
	 $oldOcc = $currOcc;
   }
  }
  my ($slope, $Em);
  if ($oldmV==$currmV)
  {
   $string .= sprintf("# E1/2 = %s%3.1f mV @ pH %3.2f\n", (($currOcc-$oldOcc)>=0.0)?'>':'<',$oldmV, $ph);
  }
  else
  {
   $slope = ($currOcc - $oldOcc)/($currmV - $oldmV);
   $Em    = (0.5 - $oldOcc)/$slope + $oldmV;
   $string .= sprintf("# E1/2 = %3.1f mV @ pH %3.2f\n", $Em, $ph);
  }
 }

 return $string;
}

sub protonation_table_html
{
 my $prob = shift;

 my $string = '';
 my @color  = ('#FFFFFF', '#DDDDDD');
 my $index  = 0;

 # get various key lists to get list of redox potentials
 my @sites  = sort sites_by_res keys %$prob;
 my @states = sort numerically keys %{$prob->{$sites[0]}};
 my @phs    = sort numerically keys %{$prob->{$sites[0]}->{$states[0]}};
 my @redoxs = sort numerically keys %{$prob->{$sites[0]}->{$states[0]}->{$phs[0]}};
 
 # print header
 $string .= "<HTML>\n<HEAD>\n<TITLE>$prefix</TITLE>\n</HEAD>\n<BODY>\n";
 $string .= "<H1>$prefix</H1>\n";
 
 # foreach redox potential
 foreach my $redox (@redoxs)
 {
  # print table ------------->
  $string .= sprintf("<P>Solution redox potential: %3.2f</P>\n", $redox); 
  $string .= "<TABLE BORDER='1'>\n";
  $string .= "<TR>\n";
  $string .= "<TH>&nbsp;</TH>";
  # header line
  foreach my $ph (@phs)
  {
   $string .= "<TH>$ph</TH>";
  }
  $string .= "\n</TR>\n";
  # data
  foreach my $site (@sites)
  {
   $index = ($index==0)?1:0;
   foreach my $state ( sort numerically keys %{$prob->{$site}} )
   {
    $string .= "<TR BGCOLOR='".$color[$index]."'>\n";
    $string .= "<TH ALIGN='LEFT'>$site ($state):</TH>";
    foreach my $ph (@phs)
    {
     my $occup = $prob->{$site}->{$state}->{$ph}->{$redox};
     $string .= sprintf("<TD>%1.2f</TD>", $occup);
    }
    $string .= "\n</TR>\n";
   }
  }
  $string .= "</TABLE>\n";
  # <------------- print table
 }

 # print footer
 $string .= "</BODY>\n";
 $string .= "</HTML>\n";
 
 return $string;
}

sub conf_plot
{
 my ($occup, $error) = (shift, shift);
 my @phs    = sort numerically keys %{$occup->{1}};
 my @redoxs = sort numerically keys %{$occup->{1}->{$phs[0]}};
 
 my $string = '';
 
 foreach my $ph ( @phs )
 {
  foreach my $redox ( @redoxs )
  {
   my $sum = 0;
	foreach my $conf (sort numerically keys %$occup)
	{
	 $sum += $occup->{$conf}->{$ph}->{$redox};
	}
   $string .= sprintf("% 3.2f % 3.2f %3.2f ", $ph, $redox, $sum);
   foreach my $conf (sort numerically keys %$occup)
	{
    $string .= sprintf("%3.2f ", $occup->{$conf}->{$ph}->{$redox});
	}
	$string .= "\n";
  }
 }
 
 return $string;
}

sub get_data
{
 # reads all data in files
 # and returns 3 complex hashes:
 #
 # biased->temp->site->state->pH->redox
 # error->temp->site->state->pH->redox
 # rebiased->temp->site->state->pH->redox

 my $biased;
 my $error;
 my $rebiased;

 my $sites = shift;
 my $files = shift;

 foreach my $file (@$files)
 {
  # get pH and redox potential from file name
  $file =~ /_pH(.*)_(.*)mV_(.*)K/; # $1 => pH, $2 => redox, $3 -> temperature
  my $ph    = $1;
  my $redox = $2;
  my $t     = $3;
  open(IN, "< $file") or die "Could not open file $file!";
  foreach my $line (<IN>)
  {
   my @temp = split(' ', $line);
	next if (@temp < 8); # 8 columns minimum
	my $site   = $temp[1];
	if ( exists($sites->{$site}) )
	{
    # loop through states
    my $states = 0;
    for(my $i=2; $i<(@temp-2); $i+=3)
    {
     $biased->{$t}->{$site}->{$states}->{$ph}->{$redox} = $temp[$i];
     $error->{$t}->{$site}->{$states}->{$ph}->{$redox}  = $temp[$i+1];
     $rebiased->{$t}->{$site}->{$states}{$ph}->{$redox} = $temp[$i+2];
     $states++;
	 }
	}
   else
   {
    printf(STDERR "Found unknown site %s in file %s!\n", $site, $file);
    exit 0;
   }
  }
  close IN;
 }

 return ($biased, $error, $rebiased);
}

sub get_conf_occup
{
 my $occup;
 my $error;
 
 my $files = shift;

 foreach my $file (@$files)
 {
  # get pH, mV, temp and conformer number
  $file =~ /_pH(.*)_(.*)mV_(.*)K_conf(\d*)$/;
  #print STDERR "$file, $1, $2, $3, $4\n";
  # open file & read first line
  open(IN, "< $file") or die "Could not open file $file!";
  my $firstline = <IN>;
  close IN;
  my @temp = split(' ', $firstline);
  # save occupancy
  $occup->{$3}->{$4}->{$1}->{$2} = $temp[0];
  $error->{$3}->{$4}->{$1}->{$2} = $temp[1];
 }
 
 return ($occup, $error);
}

sub get_filelist
{
 my $prefix = shift;
 my @temp = ();
 my $filelist = [];
 
 opendir DIR, "." or die "ERROR: Can't read working directory!";
 @temp = readdir DIR;
 closedir DIR;
 
 foreach (@temp)
 {
  if (/^($prefix)(_pH\d*\.\d*|\d*\.\d*|\.?\d*)_.+?mV_.+?K$/) # ignores multiconf
  {
   push @$filelist, $_;
  }
 }
 
 return $filelist;
}

sub get_conf_filelist
{
 my $prefix = shift;
 my @temp = ();
 my $filelist = [];
 
 opendir DIR, "." or die "ERROR: Can't read working directory!";
 @temp = readdir DIR;
 closedir DIR;
 
 foreach (@temp)
 {
  if (/^($prefix)(_pH\d*\.\d*|\d*\.\d*|\.?\d*)_.+?mV_.+?K_conf\d*$/) # ignores multiconf
  {
   push @$filelist, $_;
  }
 }
 
 return $filelist;
}

sub get_sites
{
 # reads one file and extracts site names
 
 my $filename = shift;
 my $names    = {};
 
 die "Not a valid file name!" unless (defined($filename));
 
 open(IN, "< $filename") or die "Could not open file $filename!";
 foreach my $line (<IN>)
 {
  my @temp = split(' ', $line);
  next if (@temp < 8); # 8 columns minimum
  $names->{$temp[1]} = 1;
 }
 close IN;
 
 return $names;
}

sub array2string
{
 my $array_ref = shift;
 my $string = '';
 
 foreach my $element (@$array_ref)
 {
  $string .= sprintf("%s\n", $element);
 }
 
 # remove last new line
 chop $string;
 
 return $string;
}

sub numerically {$a <=> $b}

sub sites_by_res
{
 $a =~ /([A-Z0-9]*)-(\d*)/;
 my $res1 = $1;
 my $idx1 = $2;
 $b =~ /([A-Z0-9]*)-(\d*)/;
 my $res2 = $1;
 my $idx2 = $2;
 
 my $cmp1 = $res1 cmp $res2;
 my $cmp2 = $idx1 <=> $idx2;
 
 if ( !$cmp1 )
 {
  return $cmp2;
 }

 if ( !$cmp2 )
 {
  return $cmp1;
 }
 
 return $cmp1;
}
