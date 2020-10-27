#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use File::Basename;

my $file_pdb = shift;

use Cwd 'abs_path';
my $program_dssp = dirname(abs_path($0))."/dssp-2.0.4-linux-amd64";

confess "ERROR! $file_pdb does not exist!" if not -f $file_pdb;

my %RESIDUE = ();
my %SS  = ();
my %PHI = ();
my %PSI = ();

my $command = "$program_dssp $file_pdb | grep -C 0 -A 1500 \'  #  RESIDUE\' | tail -n +2";
my @dssp_rows = `$command`;

foreach(@dssp_rows){
	my $rnum = substr($_,  5, 6);
	my $res  = substr($_, 13, 1);
	my $sstr = substr($_, 16, 1);
	my $phia = substr($_,103, 6);
	my $psia = substr($_,109, 6);
	$rnum =~ s/\s+//g;
	$res  =~ s/\s+//g;
	$sstr =~ s/\s+//g;
	$phia =~ s/\s+//g;
	$psia =~ s/\s+//g;
	# alternate residue locations may have alphabets
	#$rnum =~ s/[^0-9]//g;
	#$res  =~ s/[^A-Z]//g;
	$sstr =~ s/[^A-Z]//g;
	next if length($rnum) < 1;
	#confess " :( residue not defined for $rnum" if length($res) < 1;
	confess " :( phi not defined for $rnum" if length($phia) < 1;
	confess " :( psi not defined for $rnum" if length($psia) < 1;
	$sstr = "C" if length($sstr) < 1;
	$sstr =~ s/\./C/g;
	$sstr =~ s/I/C/g;
	$sstr =~ s/S/C/g;
	$sstr =~ s/T/C/g;
	$sstr =~ s/B/C/g;
	$sstr =~ s/G/C/g;
	$RESIDUE{$rnum} = $res;
	$SS{$rnum}      = $sstr;
	$PHI{$rnum}     = $phia;
	$PSI{$rnum}     = $psia;
}

foreach (sort { $a <=> $b} keys %RESIDUE){
	print $_." ".$RESIDUE{$_}." ".$SS{$_}." ".$PHI{$_}." ".$PSI{$_}."\n";
}
