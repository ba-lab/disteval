#!/usr/bin/perl -w
# CONFOLD updated to accept distance constraints
# version 1.1, Badri Adhikari, 10/28/2020

use strict;
use warnings;
use Carp;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

# Location of CNSsuite and DSSP
my $program_dssp   = dirname(abs_path($0))."/dssp-2.0.4-linux-amd64";
my $cns_suite      = dirname(abs_path($0))."/cns_solve_1.3";
my $cns_executable = "$cns_suite/intel-x86_64bit-linux/bin/cns_solve";

# User inputs 
my ($help, $dir_out, $file_fasta, $file_pair, $file_rr, $file_ss);
my ($selectrr, $rrtype, $lambda, $stage2, $contwt, $sswt, $mcount, $rep2, $pthres);

GetOptions(
	"h"			=> \$help,
	"o=s"		=> \$dir_out,
	"seq=s"		=> \$file_fasta,
	"pair=s"	=> \$file_pair,
	"rr=s"		=> \$file_rr,
	"ss=s"		=> \$file_ss,
	"selectrr=s"=> \$selectrr,
	"rrtype=s"	=> \$rrtype,
	"lambda=s"	=> \$lambda,
	"stage2=i"	=> \$stage2, 
	"contwt=i"	=> \$contwt, 
	"sswt=s" 	=> \$sswt,
	"mcount=i"	=> \$mcount,
	"rep2=s"	=> \$rep2,
	"pthres=s"	=> \$pthres)
or confess "ERROR! Error in command line arguments!";

print_usage() if defined $help;
print_usage() if not defined $dir_out;

if (! -e $program_dssp) {
    confess "ERROR! program dssp not found at ".$program_dssp."\n";
}
if (! -e $cns_executable) {
    confess "ERROR! cns_solve not found at ".$cns_executable."\n";
}

# Default parameters
$selectrr = "all" if !$selectrr;
$rrtype   = "cb"  if !$rrtype;
$lambda   = 0.4   if !$lambda;
$stage2   = 0     if !$stage2;
$contwt   = 10    if !$contwt;
$sswt     = 5     if !$sswt;
$mcount   = 20    if !$mcount;
$rep2     = 0.85  if !$rep2;
$pthres   = 7.0   if !$pthres;

# Other parameters
my $rep1         = 1.0;    # initial van der Waals repel radius (md.cool.init.rad)
my $mini;                  # number of minimization steps
my $pthres_hbond = $pthres + 3.0; 
my $atomselect   = 2;      # See below
my $mode         = "trial";

# Standard Residues http://proteopedia.org/wiki/index.php/Standard_Residues
my %AA3TO1 = qw/ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V/;
my %CBATOM = qw/A cb N cb C cb Q cb H cb L cb M cb P cb T cb Y cb R cb D cb E cb G ca I cb K cb F cb S cb W cb V cb/;
my %AA1TO3 = reverse %AA3TO1;

# Global variables for overall script
my $id;
my $L;
my %residues = ();

# Global variables for generating secondary structure restraints
my %ATOMTYPE     = qw( CA 1 N 1 C 1 O 1 );
my %SHIFT        = qw( 0 1 +1 1 -1 1 );
my %res_dihe     = ();
my %res_strnd_OO = ();
my %res_dist     = ();
my %res_hbnd     = ();

check_programs();
process_parameters();

chdir "$dir_out/stage1" or confess $!;
system_cmd("cp $dir_out/input/* $dir_out/stage1/");
load_ss_restraints($lambda, "ssrestraints.log");

print "\nStart [$0]: ".(localtime)."\n";

print "\nBuild extended mtf and pdb..\n";
build_extended() if not -f "extended.pdb";
my %defined_atoms = xyz_pdb("extended.pdb", "all");

# Hydrogen or Nitrogen must be present for all atoms, in order to apply hydrogen bond restraints
foreach (keys %residues){
	next if defined $defined_atoms{$_." H"};
	next if defined $defined_atoms{$_." N"};
	confess "ERROR!! Something went wrong with residue $_ extended must have H or N for each residue!";
}

my %stage_list = qw/stage1 1/;
%stage_list    = qw/stage1 1 stage2 2/ if $stage2 != 0;
foreach my $stage (sort keys %stage_list){
	my $dir_stage = "$dir_out/$stage";
	mkdir $dir_stage or confess $! if $stage eq "stage2"; 
	chdir $dir_stage or confess $!;
	print "\nStart $stage job..";
	if ($stage eq "stage2"){
		confess "ERROR! Stage 1 top model not found! Something went wrong!" if not -f "$dir_out/stage1/${id}_model1.pdb";
		system_cmd("cp $dir_out/stage1/extended.* ./");
		system_cmd("cp $dir_out/stage1/$id.fasta  ./");
        if (-e "$dir_out/stage1/$id.ss"){ 
            system_cmd("cp $dir_out/stage1/$id.ss ./");
        }
	}
	contact_restraints($stage);
	sec_restraints($stage);
	build_models($stage);
	assess_dgsa($stage);
}

# Log script finish time
print "\nFinished [$0]: ".(localtime)."\n";

sub check_programs{
	confess "ERROR! Cannot find DSSP program at $program_dssp" if not -f $program_dssp;
	confess "ERROR! Cannot find CNS suite folder at $cns_suite! Check CNS installation!" if not -d $cns_suite;
	confess "ERROR! cns_solve_env.sh not found inside $cns_suite/! Check CNS installation!" if not -f "$cns_suite/cns_solve_env.sh";
	confess "ERROR! CNS executable $cns_executable not found! Check CNS installation!" if not -f $cns_executable;
	confess "ERROR! protein.param not found inside $cns_suite/libraries/toppar/! Check CNS installation!" if not -f "$cns_suite/libraries/toppar/protein.param";
}

sub process_parameters{
	confess "ERROR! Output directory not supplied" if !$dir_out;
	confess "ERROR! Fasta file not supplied" if !$file_fasta;
	confess "ERROR! Either RR file or SS need to be supplied" if (!$file_ss and !$file_rr);
	confess "ERROR! Fasta file $file_fasta does not exist!" if not -f $file_fasta;
	confess "ERROR! Fasta file $file_rr does not exist!" if ($file_rr and not -f $file_rr);
	confess "ERROR! Fasta file $file_ss does not exist!" if ($file_ss and not -f $file_ss);
	confess "ERROR! Pairing file $file_pair does not exist!" if ($file_pair and not -f $file_pair);
	$L = length(seq_fasta($file_fasta));
	$mini = 15 * $L;
	$id = basename($file_fasta, ".fasta");
	$selectrr =~ s/L//;
	$selectrr = int($selectrr * $L + 0.5) if $selectrr ne "all";
	$selectrr = 100000 if $selectrr eq "all";
	confess "ERROR! Select rr option does not look right!" if !looks_like_number($selectrr);
	confess "ERROR! RR type must be ca or cb!" if ($rrtype ne "ca" and $rrtype ne "cb");
	confess "ERROR! Lambda must be between 0.1 and 10.0!" if ($lambda > 10.0 or $lambda < 0.1);
	confess "ERROR! Stage2 selection must be 1,2,or 3!" if ($stage2 > 3 or $stage2 < 0);
	confess "ERROR! Contact restraints weight must be between 0.1 and 10000!" if ($contwt > 1000.0 or $contwt < 0.1);
	confess "ERROR! Secondary structure restraints weight must be between 0.1 and 100!" if ($contwt > 100.0 or $contwt < 0.1);
	confess "ERROR! Model-count must be between 5 and 50!" if ($mcount > 50 or $mcount < 5);
	confess "ERROR! Rep2 must be between 0.6 and 1.5!" if ($rep2 > 1.5 or $rep2 < 0.5);
	confess "ERROR! Pair detection threshold must be between 5.0 and 10.0!" if ($pthres > 10.0 or $pthres < 5.0);
	$dir_out = abs_path($dir_out);
	mkdir $dir_out or confess $! if not -d $dir_out;
	system_cmd("rm -r $dir_out/input") if -d "$dir_out/input";
	system_cmd("rm -r $dir_out/stage1") if -d "$dir_out/stage1";
	system_cmd("rm -r $dir_out/stage2") if (-d "$dir_out/stage2" and $stage2 > 0);
	mkdir "$dir_out/input" or confess $!;
	mkdir "$dir_out/stage1" or confess $!;
	system_cmd("cp $file_fasta $dir_out/input/$id.fasta");
	system_cmd("cp $file_ss    $dir_out/input/$id.ss") if $file_ss;
	system_cmd("cp $file_rr    $dir_out/input/$id.rr") if $file_rr;
	system_cmd("cp $file_pair  $dir_out/input/$id.pair") if $file_pair;
	$file_fasta = "$id.fasta";
	$file_ss    = "$id.ss" if $file_ss;
	$file_rr    = "$id.rr" if $file_rr;
	$file_pair  = "$id.pair" if $file_pair;
	chdir "$dir_out/input" or confess $!;
	flatten_fasta($file_fasta);
	%residues = fasta2residues_hash($file_fasta);
	# Check whether or not the sequences match
	my $seq = seq_fasta($file_fasta);
	chk_errors_seq($seq);
	if ($file_rr){
		my $a = seq_rr($file_rr);
		my $b = seq_fasta($file_fasta);
		confess "ERROR! Fasta and rr sequence do not match!\nFasta:$b\nRR   :$a" if $a ne $b;
	}
	if ($file_ss){
		my $a = seq_fasta($file_ss);
		my $b = seq_fasta($file_fasta);
		confess "ERROR! Fasta and ss sequence lengths do not match!\nFasta:$b\nSS   :$a" if length($a) != length($b);
	}
	if ($file_rr){
		# Sort rr file by confidence
		system_cmd("rm -f sorted.rr");
		# Keep contact rows only
		system_cmd("sed -i '/^[A-Z]/d' $file_rr");
		# Some files have leading white spaces
		system_cmd("sed -i 's/^ *//' $file_rr");
		# Stable sort with -s option, i.e. maintain order in case confidence are equal
		system_cmd("sort -nr -s -k5 $file_rr > sorted.rr");
		system_cmd("rm -f $file_rr");
		print2file($file_rr, $seq) if defined $seq;
		system_cmd("cat sorted.rr >> $file_rr");
		system_cmd("rm -f sorted.rr");
		# Check if any of the contacts are out of sequence
		my %rr = rr2contacts_hash($file_rr, 1);
		foreach (keys %rr){
			my @C = split /\s+/, $_;
			confess "ERROR! Residue ".$C[0]." in contact $_ is out of sequence!" if not defined $residues{$C[0]};
			confess "ERROR! Residue ".$C[1]." in contact $_ is out of sequence!" if not defined $residues{$C[1]};
		}
	}
	# Check SS characters for non-standard characters
	if($file_ss){
		my $sec = seq_fasta($file_ss);
		for(my $i = 1; $i <= length($sec); $i++){
			my $char = substr $sec, $i-1, 1;
			if (not ($char eq "H" or $char eq "C" or $char eq "E")){
				confess "ERROR! undefined secondary structure unit $char in $sec";
			}
		}
	}
	# Check sequence for non-standard amino acids
	for(my $i = 1; $i <= $L; $i++){
		my $char = substr $seq, $i-1, 1;
		if (not defined $AA1TO3{$char}){
			confess "ERROR! undefined amino acid $char in $seq";
		}
	}
	# Check pairing file contents
	if ($file_pair){
		open PAIR, $file_pair or confess $!;
		while (<PAIR>){
			my $row = $_;
			chomp $row;
			$row =~ s/^\s+//;
			next unless $row;
			confess "ERROR! Pairing file Error! Numbers expected at beginning in $row" if $row !~ /^[0-9]/; 
			my @C = split /\s+/, $row;
			for( my $i=0; $i <= 4; $i++){
				confess "ERROR! Pairing file Error! Column $i is missing in $row" if not defined $C[$i];
			}
			confess "ERROR! Pairing file Error! a must be less than b in $row" if ($C[0] >= $C[1]);
			confess "ERROR! Pairing file Error! Pairing must be A or P in $row" if ($C[4] ne "A" and $C[4] ne "P");
			confess "ERROR! Pairing file Error! c must be less than d for parallel pairing in $row" if ($C[2] >= $C[3] and $C[4] eq "P");
			confess "ERROR! Pairing file Error! c must be greater than d for anti-parallel pairing in $row" if ($C[2] <= $C[3] and $C[4] eq "A");
		}
		close PAIR;
	}
	print "dir_out     $dir_out     "; 
	print "\nfile_fasta  $file_fasta  ";
	print "\nfile_pair   $file_pair   " if $file_pair;
	print "\nfile_rr     $file_rr     " if $file_rr;
	print "\nfile_ss     $file_ss     " if $file_ss;
	print "\nselectrr    $selectrr    ";
	print "\nrrtype      $rrtype      ";
	print "\nlambda      $lambda      ";
	print "\nstage2      $stage2      ";
	print "\ncontwt      $contwt      ";
	print "\nsswt        $sswt        ";
	print "\nmcount      $mcount      ";
	print "\nrep2        $rep2        ";
	print "\npthres      $pthres      ";
	print "\nrep1        $rep1        ";
	print "\nmini        $mini        ";
	print "\natomselect  $atomselect  ";
	print "\nmode        $mode        ";
	print "\nid          $id          ";
	print "\nL           $L           ";
	print "\nsequence    ".seq_fasta($file_fasta);
	print "\nss          ".seq_fasta($file_ss) if $file_ss;
	print "\n";
	chk_errors_rr($file_rr) if $file_rr;
}

sub build_extended{
	my $file_seq = "$id.fasta";
	confess "ERROR! expected fasta file not found in the extended directory!" if not -f $file_seq;
	flatten_fasta($file_seq);
	write_cns_seq($file_seq, "input.seq");
	write_cns_generate_seq_file();
	write_cns_generate_extended_file();
	open  JOB, ">job.sh" or confess "ERROR! Could not open job.sh $!";
	print JOB "#!/bin/bash                                     \n";
	print JOB "# CNS-CONFIGURATION                             \n";
	print JOB "source $cns_suite/cns_solve_env.sh        \n";
	print JOB "export KMP_AFFINITY=none                        \n";
	print JOB "$cns_executable < gseq.inp > gseq.log \n";
	print JOB "$cns_executable < extn.inp > extn.log \n";
	close JOB;
	system_cmd("chmod +x job.sh");
	`./job.sh > job.log`;
	confess "FAILED! extended.mtf not found!" if not -f "extended.pdb";
	system_cmd("rm -f gseq.log");
	system_cmd("rm -f extn.log");
}

sub contact_restraints{
	my $stage = shift;
	if($stage eq "stage1"){
		return if not $file_rr;
		my $xL = $selectrr + 1; # +1 to account for header line
		system_cmd("head -n $xL $file_rr > temp.rr");
		system_cmd("rm $file_rr");
		system_cmd("mv temp.rr $file_rr");
		rr2tbl($file_rr, "contact.tbl", $rrtype);
		return;
	}
	my $stg1_rr = "$dir_out/stage1/$id.rr";
	return if not -f $stg1_rr;
	my $stg1_model = "$dir_out/stage1/${id}_model1.pdb";
	if($stage2 == 1 or $stage2 == 2){
		rr_remove_unsatisfied_top_model($stg1_rr, $stg1_model, "$id.rr", "rr_filter.log");
	}
	else{
		system_cmd("cp $stg1_rr ./$id.rr");
	}
	rr2tbl("$id.rr", "contact.tbl", $rrtype);
}

sub build_models{
	my $stage = shift;
	# load restraints and print to logs
	my %tbl_list = ();
	$tbl_list{"contact.tbl"}  = count_lines("contact.tbl")  if count_lines("contact.tbl");
	$tbl_list{"ssnoe.tbl"}    = count_lines("ssnoe.tbl")    if count_lines("ssnoe.tbl");
	$tbl_list{"hbond.tbl"}    = count_lines("hbond.tbl")    if count_lines("hbond.tbl");
	$tbl_list{"dihedral.tbl"} = count_lines("dihedral.tbl") if count_lines("dihedral.tbl");
	print "\nWARNING! contact/rr restraints not found or empty!" if not defined $tbl_list{"contact.tbl"};
	print "\nWARNING! secondary structure restraints not found or empty!" if not defined $tbl_list{"dihedral.tbl"};
	print "\n".seq_fasta("$id.fasta");
	print "\n".seq_fasta("$id.ss") if -f "$id.ss";
	foreach my $tbl (sort keys %tbl_list){
		my $flag_dihe = 0;
		$flag_dihe = 1 if $tbl eq "dihedral.tbl";
		print "\n".coverage_tbl("${id}.fasta", $tbl, $flag_dihe);
	}
	# prepare CNS task files
	system_cmd("rm -f iam.*");
	write_cns_customized_modules();
	write_cns_dgsa_file();
	system_cmd("sed -i s/contact.tbl//g dgsa.inp")  if not defined $tbl_list{"contact.tbl"};
	system_cmd("sed -i s/hbond.tbl//g dgsa.inp")    if not defined $tbl_list{"hbond.tbl"};
	system_cmd("sed -i s/dihedral.tbl//g dgsa.inp") if not defined $tbl_list{"dihedral.tbl"};
	system_cmd("sed -i s/ssnoe.tbl//g dgsa.inp")    if not defined $tbl_list{"ssnoe.tbl"};
	open  JOB, ">job.sh" or confess "ERROR! Could not open job.sh $!";
	print JOB "#!/bin/bash                                       \n";
	print JOB "echo \"starting cns..\"                           \n";
	print JOB "touch iam.running                                 \n";
	print JOB "# CNS-CONFIGURATION                               \n";
	print JOB "source $cns_suite/cns_solve_env.sh                \n";
	print JOB "export KMP_AFFINITY=none                          \n";
	print JOB "export CNS_CUSTOMMODULE=$dir_out/$stage           \n";
	print JOB "$cns_suite/intel-x86_64bit-linux/bin/cns_solve < dgsa.inp > dgsa.log \n";
	print JOB "if [ -f \"${id}_${mcount}.pdb\" ]; then      \n";
	print JOB "   rm iam.running                                 \n";
	print JOB "   echo \"trial structures written.\"             \n";
	print JOB "   rm *embed*                                     \n";
	print JOB "   exit                                           \n";
	print JOB "fi                                                \n";
	print JOB "if [ -f \"${id}a_${mcount}.pdb\" ]; then \n";
	print JOB "   rm iam.running                                 \n";
	print JOB "   echo \"accepted structures written.\"          \n";
	print JOB "   rm *embed*                                     \n";
	print JOB "   exit                                           \n";
	print JOB "fi                                                \n";
	print JOB "tail -n 30 dgsa.log                              \n";
	print JOB "echo \"ERROR! Final structures not found!\"       \n";
	print JOB "echo \"CNS FAILED!\"                              \n";
	print JOB "mv iam.running iam.failed                         \n";
	close JOB;
	print "\n"."Starting job [$dir_out/$stage/job.sh > job.log]";
	system_cmd("chmod +x $dir_out/$stage/job.sh");
	system_cmd("$dir_out/$stage/job.sh > job.log");
	confess "ERROR! Something went wrong while running CNS! Check job.log and dgsa.log!" if -f "iam.failed";
}

sub sec_restraints{
	my $stage = shift;
	return if not $file_ss;
	system_cmd("rm -f ssrestraints.log");
	# helix restraints; start writing to the files hbond.tbl, ssnoe.tbl and dihedral.tbl
	print_helix_restraints();
	# strand and/or sheet restraints
	if (strand_count($file_ss) == 0){
		print2file("ssrestraints.log", "no strand restraints");
		return;
	}
	if (strand_count($file_ss) == 1){
		warn "\nWARNING! only 1 strand! please check your secondary structure file!";
		print2file("ssrestraints.log", "WARNING! only 1 strand! please check your secondary structure file!");
		return;
	}
	my $pair_hbonds;
	if ($stage eq "stage1"){
		if ($file_pair){
			pairing2hbonds($file_pair, "pairing_hbonds.txt");
			$pair_hbonds = "pairing_hbonds.txt";
			confess "ERROR! Pairing file is defined but does not exist!" if not -f $file_pair;
			confess "ERROR! Pairing hbond file is defined but does not exist!" if not -f $pair_hbonds;
			confess "ERROR! Pairing file is empty!" if not -s $file_pair;
			confess "ERROR! Pairing hbond file is empty!" if not -s $pair_hbonds;
			update_sec_using_pairing_info($file_ss, $file_pair, $file_ss);
			# continue to write to hbond.tbl, ssnoe.tbl and dihedral.tbl
			strand_and_sheet_tbl($file_ss, $file_pair, $pair_hbonds);
		}
		else{
			# continue to write to hbond.tbl, ssnoe.tbl and dihedral.tbl
			strand_and_sheet_tbl($file_ss, undef, undef);
		}
	}
	else{
		my $stg1_model = "$dir_out/stage1/${id}_model1.pdb";
		if ($stage2 == 1 or $stage2 == 3){
			detect_hbonds_from_model($stg1_model, "pairing.txt", "pairing_hbonds.txt");
			if (-f "pairing_hbonds.txt"){
				$file_pair = "pairing.txt";
				$pair_hbonds = "pairing_hbonds.txt";
				confess "ERROR! Pairing file is defined but does not exist!" if not -f $file_pair;
				confess "ERROR! Pairing hbond file is defined but does not exist!" if not -f $pair_hbonds;
				confess "ERROR! Pairing file is empty!" if not -s $file_pair;
				confess "ERROR! Pairing hbond file is empty!" if not -s $pair_hbonds;
				update_sec_using_pairing_info($file_ss, $file_pair, $file_ss);
			}
			else{
				print "\n\nWARNING! Could not detect any sheet pairings from the top model in stage1!\n";
			}
		}
		# continue to write to hbond.tbl, ssnoe.tbl and dihedral.tbl
		strand_and_sheet_tbl($file_ss, $file_pair, $pair_hbonds);
	}
}

sub print_helix_restraints{
	my %res_sec = fasta2residues_hash($file_ss);
	foreach (keys %res_sec){
		delete $res_sec{$_} if $res_sec{$_} ne "H";
	}
	if (scalar keys %res_sec){
		print2file("ssrestraints.log", "write helix tbl restraints");
		foreach my $i (sort {$a <=> $b} keys %res_sec){
			my @PHI = split /\s+/, $res_dihe{"H PHI"};
			my @PSI = split /\s+/, $res_dihe{"H PSI"};
			# 9/1/2014: If phi angle is removed from the pre-first residue, the first residue cannot form helix most of the times. Checked in 4/5 proteins and 2 dozen helices.
			# So, adding the pre-first phi angle as well (by removing the condition for  res_sec{$i-1} to exist)
			print2file("dihedral.tbl", (sprintf "assign (resid %3d and name c) (resid %3d and name  n) (resid %3d and name ca) (resid %3d and name c) 5.0 %7s %7s 2 !helix phi", $i-1, $i, $i, $i, $PHI[0], $PHI[1])) if defined $residues{$i-1};
			# 9/1/2014: Even if we don't have PSI for the last residue, the last residue form helix, almost in all cases.
			# 9/2/2014: In some cases like target Tc767, only 106 helix residues were formed while the input is 111. And all of it is because of the psi angle at the end.
			# So, adding the post-last psi angle as well (by removing the condition for res_sec{$i+1} to exist)
			print2file("dihedral.tbl", (sprintf "assign (resid %3d and name n) (resid %3d and name ca) (resid %3d and name  c) (resid %3d and name n) 5.0 %7s %7s 2 !helix psi", $i, $i, $i, $i+1, $PSI[0], $PSI[1])) if defined $residues{$i+1};
		}
		foreach my $i (sort {$a <=> $b} keys %res_sec){
			next if not defined $res_sec{$i+1};
			next if not defined $res_sec{$i+2};
			next if not defined $res_sec{$i+3};
			next if not defined $res_sec{$i+4};
			my @HR = split /\s+/, $res_hbnd{"H"};
			# In case of Proline, use N atom instead of H
			my $HATOM = "H";
			$HR[0] -= 1.0;
			if($residues{$i+4} eq "P"){
				$HR[0] += 1.0;
				$HATOM  = "N";
			}
			print2file("hbond.tbl", (sprintf "assign (resid %3d and name O) (resid %3d and name $HATOM) %4s %4s %4s !helix", $i, $i+4, $HR[0], $HR[1], $HR[2]));
		}
		foreach my $i (sort {$a <=> $b} keys %res_sec){
			next if not defined $res_sec{$i+1};
			next if not defined $res_sec{$i+2};
			next if not defined $res_sec{$i+3};
			next if not defined $res_sec{$i+4};
			foreach my $A (sort keys %ATOMTYPE){
				foreach my $SH (sort keys %SHIFT){
					my @AR = split /\s+/, $res_dist{"H $A-$A O $SH"};
					# found restraints even for non-helix SS
					next if not defined $res_sec{$i+4+$SH};
					confess ":(" if not defined ($i and $A and ($i+4+$SH) and $AR[0] and $AR[1] and $AR[2]);
					print2file("ssnoe.tbl", (sprintf "assign (resid %3d and name %2s) (resid %3d and name %2s) %.2f %.2f %.2f !helix", $i, $A, ($i+4+$SH), $A, $AR[0], $AR[1], $AR[2]));
				}
			}
		}
	}
	else{
		print2file("ssrestraints.log", "no helix predictions!");
	}
}

sub pairing2hbonds{
	my $file_pair = shift;
	my $file_hbonds = shift;
	my %res_ss = fasta2residues_hash($file_ss);
	my %pairs = pairing_info2hash($file_pair);
	print2file("ssrestraints.log", "initial pairing");
	print2file("ssrestraints.log", "A/P  R1-R2   R1-R2  Confidence");
	foreach (sort {$pairs{$a} <=> $pairs{$b}} keys %pairs){
		my @A = split /\s+/, $_;
		print2file("ssrestraints.log", sprintf "%3s %3d-%-3d %3d-%-3d", $A[4], $A[0], $A[1], $A[2], $A[3]);
	}
	#build hydrogen bonding pattern, and shift strands if necessary, for pairing
	my %hbond_connectors = ();
	my %hbonds = ();
	my %slided_pairs = ();
	my %possible_slides = qw( 0 1 L1 2 R1 3 L2 4 R2 5 ); # total 10 possible configurations
	my %possible_configs = qw( a 1 b 2 c 3 d 4 ); # parallel can have c and d hydrogen bonding pattern as well
	print2file("ssrestraints.log", "Selected Configurations:");
	foreach my $pair (sort {$pairs{$a} <=> $pairs{$b}} keys %pairs){
		my @A = split /\s+/, $pair;
		if ((abs($A[0] - $A[1]) < 2) or abs($A[2] - $A[3]) < 2){
			delete $pairs{$pair};
			next;
		}
		my $selected_shift_config = undef;
		foreach my $slide (sort {$possible_slides{$a} <=> $possible_slides{$b}} keys %possible_slides){
			last if defined $selected_shift_config;
			foreach my $config (sort {$possible_configs{$a} <=> $possible_configs{$b}} keys %possible_configs){
				last if defined $selected_shift_config;
				my @A = split /\s+/, $pair;
				next if ($A[4] eq "A" and not ($config eq "a" or $config eq "b")); # anti-parallel do not have hbond-patterns c and d
				my $no_of_residues = $slide;
				$no_of_residues =~ s/[A-Z]//g;
				# first strand is sliding forward 
				if ($slide =~ /^R/){
					$A[1] = $A[1] - $no_of_residues if $A[4] eq "A";
					$A[2] = $A[2] - $no_of_residues if $A[4] eq "A";
					$A[1] = $A[1] - $no_of_residues if $A[4] eq "P";
					$A[2] = $A[2] + $no_of_residues if $A[4] eq "P";
				}
				# first strand is sliding backward
				elsif ($slide =~ /^L/){
					$A[0] = $A[0] + $no_of_residues if $A[4] eq "A";
					$A[3] = $A[3] + $no_of_residues if $A[4] eq "A";
					$A[0] = $A[0] + $no_of_residues if $A[4] eq "P";
					$A[3] = $A[3] - $no_of_residues if $A[4] eq "P";
				}
				next if ((abs($A[0] - $A[1]) < 2) or abs($A[2] - $A[3]) < 2);
				warn "\nWARNING! 2 residue shifting is needed for pair $pair slide $slide config $config!" if $no_of_residues == 2;
				my %ideal_hbonds = make_ideal_hbonds($A[0], $A[1], $A[2], $A[3], $A[4], $config);
				my %connectors = hbonds2connectors(\%ideal_hbonds);
				my $flag_reject = 0;
				foreach (keys %connectors){
					$flag_reject = 1 if (defined $hbond_connectors{$_});
				}
				if (not $flag_reject){
					$selected_shift_config = $A[0]." ".$A[1]." ".$A[2]." ".$A[3]." ".$A[4]." ".$config;
					print2file("ssrestraints.log", "pair:$pair slide:$slide selected:".$A[0]." ".$A[1]." ".$A[2]." ".$A[3]." ".$A[4]." ".$config);
				}
			}
		}
		confess ":( $_" if not defined $selected_shift_config;
		my @S = split /\s+/, $selected_shift_config;
		my %ideal_hbonds = make_ideal_hbonds($S[0], $S[1], $S[2], $S[3], $S[4], $S[5]);
		my %connectors = hbonds2connectors(\%ideal_hbonds);
		foreach (keys %ideal_hbonds){
			$hbonds{$_} = $S[4];
		}
		foreach (keys %connectors){
			$hbond_connectors{$_} = 1;
		}
		$slided_pairs{$S[0]." ".$S[1]." ".$S[2]." ".$S[3]." ".$S[4]." ".$S[5]} = $S[0];
	}
	foreach (sort keys %hbonds){
		print2file($file_hbonds, $_." ".$hbonds{$_});
	}
}

sub strand_and_sheet_tbl{
	my $file_ss = shift;
	my $file_pair = shift;
	my $file_pairing_hbond = shift;
	my %res_ss = fasta2residues_hash($file_ss);
	my %unpaired_residues = %res_ss;
	my %paired_residues = ();
	my %res_ssE = %res_ss;
	foreach (sort keys %res_ssE){
		delete $res_ssE{$_} if $res_ssE{$_} ne "E";
	}
	my %hbonds = ();
	if ($file_pair){
		%paired_residues = paired_residues($file_pair);
		foreach (sort keys %unpaired_residues){
			delete $unpaired_residues{$_} if $unpaired_residues{$_} ne "E";
		}
		foreach (sort keys %unpaired_residues){
			delete $unpaired_residues{$_} if defined $paired_residues{$_};
		}
		%hbonds = load_hbonds($file_pairing_hbond);
		confess ":(" if not scalar %paired_residues;
		confess ":(" if not scalar %hbonds;
		foreach (sort keys %hbonds){
			my @H = split /\s+/, $_;
			my @HR = split /\s+/, $res_hbnd{$hbonds{$_}};
			confess ":( distance not defined ".$hbonds{$_} if not (defined $HR[0] and defined $HR[1] and defined $HR[2]);
			# In case of Proline, use N atom instead of H
			$HR[0] -= 1.0;
			if(not defined $defined_atoms{$H[0]." ".$H[1]}){
				confess "ERROR! N must be defined when H is not for ".$H[0] if not defined $defined_atoms{$H[0]." N"};
				$HR[0] += 1.0;
				$H[1]  = "N";
			}
			if(not defined $defined_atoms{$H[2]." ".$H[3]}){
				confess "ERROR! N must be defined when H is not for ".$H[2] if not defined $defined_atoms{$H[2]." N"};
				$HR[0] += 1.0;
				$H[3]  = "N";
			}
			print2file("hbond.tbl", (sprintf "assign (resid %3d and name %1s) (resid %3d and name %1s) %4s %4s %4s !sheet", $H[0], $H[1], $H[2], $H[3], $HR[0], $HR[1], $HR[2]));
		}
		foreach my $hb (sort keys %hbonds){
			my @H = split /\s+/, $hb;
			my @hbondConnectors = ( $H[0]." ".$H[1]." ".$H[2], $H[2]." ".$H[3]." ".$H[0]); 
			foreach (@hbondConnectors){
				my @HBC = split /\s+/, $_;
				foreach my $A (sort keys %ATOMTYPE){
					foreach my $S (sort keys %SHIFT){
						# A   O-O   O -1 
						my $distances = $res_dist{$hbonds{$hb}." ".$A."-".$A." ".$HBC[1]." ".$S};
						my @DIST = split /\s+/, $distances;
						confess ":( distance not defined ".$hbonds{$hb}." ".$A."-".$A." ".$HBC[1]." ".$S if not defined $DIST[2];
						next if not defined $paired_residues{($HBC[2]+$S)};
						print2file("ssnoe.tbl", (sprintf "assign (resid %3d and name %2s) (resid %3d and name %2s) %2.2f %2.2f %2.2f !sheet", $HBC[0], $A, ($HBC[2]+$S), $A, $DIST[0], $DIST[1], $DIST[2]));
					}
				}
			}
		}
	}
	# Identify the strands that are not used for pairing and generate generic dihedral restraints for them
	foreach my $i (sort {$a <=> $b} keys %res_ssE){
		my @SPHI = ();
		my @SPSI = ();
		my $strand_type = "unpaired E residue";
		if (defined $paired_residues{$i}){
			@SPHI = split /\s+/, $res_dihe{$paired_residues{$i}." PHI"};
			@SPSI = split /\s+/, $res_dihe{$paired_residues{$i}." PSI"};
			$strand_type = "paired E residue";
		}
		else{
			@SPHI = split /\s+/, $res_dihe{"U PHI"};
			@SPSI = split /\s+/, $res_dihe{"U PSI"};
		}
		if (defined $res_ss{$i-1} and $res_ss{$i-1} eq "E"){
			print2file("dihedral.tbl", (sprintf "assign (resid %3d and name c) (resid %3d and name  n) (resid %3d and name ca) (resid %3d and name c) 5.0 %7s %7s 2 !$strand_type phi", $i-1, $i, $i, $i, $SPHI[0], $SPHI[1]));
		}
		if (defined $res_ss{$i+1} and $res_ss{$i+1} eq "E"){
			print2file("dihedral.tbl", (sprintf "assign (resid %3d and name n) (resid %3d and name ca) (resid %3d and name  c) (resid %3d and name n) 5.0 %7s %7s 2 !$strand_type psi", $i, $i, $i, $i+1, $SPSI[0], $SPSI[1]));
		}
	}
	# Identify the strands that are not used for pairing and generate generic restraints for them
	foreach my $i (sort {$a <=> $b} keys %res_ssE){
		my @SD = ();
		my $strand_type = "unpaired E residue";
		if (defined $paired_residues{$i}){
			@SD = split /\s+/, $res_strnd_OO{$paired_residues{$i}};
			confess ":(" if (!$SD[0] or !$SD[1] or !$SD[2]);
			$strand_type = "paired E residue";
		}
		else{
			@SD = split /\s+/, $res_strnd_OO{"U"};
			confess ":(" if (!$SD[0] or !$SD[1] or !$SD[2]);
		}
		next if not defined $res_ssE{$i+1};
		next if $res_ssE{$i+1} ne "E";
		print2file("ssnoe.tbl", (sprintf "assign (resid %3d and name %2s) (resid %3d and name %2s) %.2f %.2f %.2f !$strand_type", $i, "O", $i+1, "O", $SD[0], $SD[1], $SD[2]));
	}
}

sub load_ss_restraints{
	my $lamda = shift;
	my $log_reference = shift;
	# T      Helix or Parallel or anti-parallel or Unknown Strand Type
	# A1_A2  Atom1-Atom2 pair
	# Ref    Hydrogen bonding connector atom (reference hbond connector)
	# N      Neighborhood residue shifting on the hbond connector of R2 side. For example, If R1:N and R2:O have hbond and S = +1, the restraint A1-A2 are for R1 and (R2+1)
	# Note: hbond distances are the distances between Nitrogen and Oxygen
	# Places to verify:
	# http://www.beta-sheet.org/page29/page51/page53/index.html
	# In this model, the HP sheet is composed of identical straight helical chains with phi = -122 degrees, psi = 135 degrees, and a slightly non-linear interchain H-bond angle delta of 170 degrees.
	# http://en.wikipedia.org/wiki/Alpha_helix
	# Residues in α-helices typically adopt backbone (φ, ψ) dihedral angles around (-60°, -45°), as shown in the image at right
	# the H to O distance is about 2 Å (0.20 nm

	# T A Mean Standard_deviation
	$res_dihe{"A PSI"} = "136.91 ".($lamda * 17.39);
	$res_dihe{"A PHI"} = "-120.89 ".($lamda * 21.98);
	$res_dihe{"P PSI"} = "130.96 ".($lamda * 16.66);
	$res_dihe{"P PHI"} = "-115.00 ".($lamda * 20.31);
	$res_dihe{"U PSI"} = "134.95 ".($lamda * 17.65);
	$res_dihe{"U PHI"} = "-118.91 ".($lamda * 21.73);
	$res_dihe{"H PSI"} = "-41.51 ".($lamda * 9.84);
	$res_dihe{"H PHI"} = "-63.47 ".($lamda * 9.20);

	#T A1_A2 Ref N Mean Standard_deviation
	$res_dist{"A O-O O +1"}   = get_dist_neg_pos(7.73, 0.59, $lamda);
	$res_dist{"A O-O O -1"}   = get_dist_neg_pos(4.84, 0.16, $lamda);
	$res_dist{"A O-O O 0"}    = get_dist_neg_pos(3.57, 0.28, $lamda);
	$res_dist{"A O-O H +1"}   = get_dist_neg_pos(7.76, 0.60, $lamda);
	$res_dist{"A O-O H -1"}   = get_dist_neg_pos(4.90, 0.45, $lamda);
	$res_dist{"A O-O H 0"}    = get_dist_neg_pos(3.58, 0.31, $lamda);
	$res_dist{"A C-C O +1"}   = get_dist_neg_pos(7.66, 0.52, $lamda);
	$res_dist{"A C-C O -1"}   = get_dist_neg_pos(4.80, 0.17, $lamda);
	$res_dist{"A C-C O 0"}    = get_dist_neg_pos(4.96, 0.21, $lamda);
	$res_dist{"A C-C H +1"}   = get_dist_neg_pos(7.65, 0.51, $lamda);
	$res_dist{"A C-C H -1"}   = get_dist_neg_pos(4.85, 0.34, $lamda);
	$res_dist{"A C-C H 0"}    = get_dist_neg_pos(4.96, 0.21, $lamda);
	$res_dist{"A N-N O +1"}   = get_dist_neg_pos(5.09, 0.34, $lamda);
	$res_dist{"A N-N O -1"}   = get_dist_neg_pos(6.86, 0.40, $lamda);
	$res_dist{"A N-N O 0"}    = get_dist_neg_pos(4.42, 0.24, $lamda);
	$res_dist{"A N-N H +1"}   = get_dist_neg_pos(5.04, 0.21, $lamda);
	$res_dist{"A N-N H -1"}   = get_dist_neg_pos(6.85, 0.45, $lamda);
	$res_dist{"A N-N H 0"}    = get_dist_neg_pos(4.43, 0.25, $lamda);
	$res_dist{"A CA-CA O +1"} = get_dist_neg_pos(6.43, 0.41, $lamda);
	$res_dist{"A CA-CA O -1"} = get_dist_neg_pos(5.67, 0.28, $lamda);
	$res_dist{"A CA-CA O 0"}  = get_dist_neg_pos(5.26, 0.24, $lamda);
	$res_dist{"A CA-CA H +1"} = get_dist_neg_pos(6.38, 0.36, $lamda);
	$res_dist{"A CA-CA H -1"} = get_dist_neg_pos(5.71, 0.40, $lamda);
	$res_dist{"A CA-CA H 0"}  = get_dist_neg_pos(5.27, 0.25, $lamda);
	$res_dist{"P O-O O +1"}   = get_dist_neg_pos(7.90, 0.61, $lamda);
	$res_dist{"P O-O O -1"}   = get_dist_neg_pos(4.86, 0.16, $lamda);
	$res_dist{"P O-O O 0"}    = get_dist_neg_pos(3.78, 0.34, $lamda);
	$res_dist{"P O-O H +1"}   = get_dist_neg_pos(4.92, 0.40, $lamda);
	$res_dist{"P O-O H -1"}   = get_dist_neg_pos(8.02, 0.60, $lamda);
	$res_dist{"P O-O H 0"}    = get_dist_neg_pos(3.78, 0.32, $lamda);
	$res_dist{"P C-C O +1"}   = get_dist_neg_pos(8.03, 0.51, $lamda);
	$res_dist{"P C-C O -1"}   = get_dist_neg_pos(4.82, 0.17, $lamda);
	$res_dist{"P C-C O 0"}    = get_dist_neg_pos(5.21, 0.25, $lamda);
	$res_dist{"P C-C H +1"}   = get_dist_neg_pos(4.88, 0.34, $lamda);
	$res_dist{"P C-C H -1"}   = get_dist_neg_pos(7.87, 0.44, $lamda);
	$res_dist{"P C-C H 0"}    = get_dist_neg_pos(5.22, 0.22, $lamda);
	$res_dist{"P N-N O +1"}   = get_dist_neg_pos(8.14, 0.35, $lamda);
	$res_dist{"P N-N O -1"}   = get_dist_neg_pos(4.86, 0.40, $lamda);
	$res_dist{"P N-N O 0"}    = get_dist_neg_pos(5.13, 0.32, $lamda);
	$res_dist{"P N-N H +1"}   = get_dist_neg_pos(4.80, 0.18, $lamda);
	$res_dist{"P N-N H -1"}   = get_dist_neg_pos(7.54, 0.69, $lamda);
	$res_dist{"P N-N H 0"}    = get_dist_neg_pos(5.10, 0.28, $lamda);
	$res_dist{"P CA-CA O +1"} = get_dist_neg_pos(8.55, 0.37, $lamda);
	$res_dist{"P CA-CA O -1"} = get_dist_neg_pos(4.90, 0.29, $lamda);
	$res_dist{"P CA-CA O 0"}  = get_dist_neg_pos(6.21, 0.26, $lamda);
	$res_dist{"P CA-CA H +1"} = get_dist_neg_pos(4.90, 0.28, $lamda);
	$res_dist{"P CA-CA H -1"} = get_dist_neg_pos(7.49, 0.60, $lamda);
	$res_dist{"P CA-CA H 0"}  = get_dist_neg_pos(6.24, 0.24, $lamda);
	$res_dist{"H O-O O +1"}   = get_dist_neg_pos(8.40, 0.27, $lamda);
	$res_dist{"H O-O O -1"}   = get_dist_neg_pos(4.99, 0.16, $lamda);
	$res_dist{"H O-O O 0"}    = get_dist_neg_pos(6.12, 0.26, $lamda);
	$res_dist{"H O-O H +1"}   = get_dist_neg_pos(5.03, 0.31, $lamda);
	$res_dist{"H O-O H -1"}   = get_dist_neg_pos(8.43, 0.32, $lamda);
	$res_dist{"H O-O H 0"}    = get_dist_neg_pos(6.12, 0.26, $lamda);
	$res_dist{"H C-C O +1"}   = get_dist_neg_pos(8.16, 0.24, $lamda);
	$res_dist{"H C-C O -1"}   = get_dist_neg_pos(4.87, 0.13, $lamda);
	$res_dist{"H C-C O 0"}    = get_dist_neg_pos(6.09, 0.23, $lamda);
	$res_dist{"H C-C H +1"}   = get_dist_neg_pos(4.89, 0.23, $lamda);
	$res_dist{"H C-C H -1"}   = get_dist_neg_pos(8.17, 0.25, $lamda);
	$res_dist{"H C-C H 0"}    = get_dist_neg_pos(6.09, 0.23, $lamda);
	$res_dist{"H N-N O +1"}   = get_dist_neg_pos(8.07, 0.23, $lamda);
	$res_dist{"H N-N O -1"}   = get_dist_neg_pos(4.84, 0.19, $lamda);
	$res_dist{"H N-N O 0"}    = get_dist_neg_pos(6.10, 0.20, $lamda);
	$res_dist{"H N-N H +1"}   = get_dist_neg_pos(4.81, 0.13, $lamda);
	$res_dist{"H N-N H -1"}   = get_dist_neg_pos(8.08, 0.21, $lamda);
	$res_dist{"H N-N H 0"}    = get_dist_neg_pos(6.10, 0.20, $lamda);
	$res_dist{"H CA-CA O +1"} = get_dist_neg_pos(8.63, 0.28, $lamda);
	$res_dist{"H CA-CA O -1"} = get_dist_neg_pos(5.13, 0.20, $lamda);
	$res_dist{"H CA-CA O 0"}  = get_dist_neg_pos(6.16, 0.26, $lamda);
	$res_dist{"H CA-CA H +1"} = get_dist_neg_pos(5.14, 0.21, $lamda);
	$res_dist{"H CA-CA H -1"} = get_dist_neg_pos(8.64, 0.26, $lamda);
	$res_dist{"H CA-CA H 0"}  = get_dist_neg_pos(6.16, 0.26, $lamda);

	# T Mean Standard_deviation
	$res_strnd_OO{"A"} = get_dist_neg_pos(4.57, 0.30, $lamda);
	$res_strnd_OO{"P"} = get_dist_neg_pos(4.57, 0.29, $lamda);
	$res_strnd_OO{"U"} = get_dist_neg_pos(4.57, 0.30, $lamda);

	# T Mean Standard_deviation
	$res_hbnd{"A"} = get_dist_neg_pos(2.92, 0.16, $lamda);
	$res_hbnd{"P"} = get_dist_neg_pos(2.93, 0.16, $lamda);
	$res_hbnd{"H"} = get_dist_neg_pos(2.99, 0.17, $lamda);

	confess ":( dihe restraints could not be loaded" if not scalar keys %res_dihe;
	confess ":( dstr restraints could not be loaded" if not scalar keys %res_strnd_OO;
	confess ":( dist restraints could not be loaded" if not scalar keys %res_dist;
	confess ":( hbnd restraints could not be loaded" if not scalar keys %res_hbnd;

	system_cmd("rm -f $log_reference");
	foreach (keys %res_dihe){ print2file($log_reference, $_."\t".$res_dihe{$_});}
	foreach (keys %res_strnd_OO){ print2file($log_reference, $_."\t".$res_strnd_OO{$_});}
	foreach (keys %res_dist){ print2file($log_reference, $_."\t".$res_dist{$_});}
	foreach (keys %res_hbnd){ print2file($log_reference, $_."\t".$res_hbnd{$_});}
}

sub get_dist_neg_pos{
	my $mean = shift;
	my $devi = shift;
	my $lamda = shift;
	my $lamda_devi = sprintf "%.1f", $lamda * $devi;
	return $mean." ".$lamda_devi." ".$lamda_devi;
}

sub detect_hbonds_from_model{
	my $file_model = shift;
	my $file_pair  = shift;
	my $file_hbond = shift;
	my %strands = strand_list($file_ss);
	my %final_pairs = ();
	my %final_aligned_pairs = ();
	confess ":( model does not exist $file_model " if not -f $file_model;
	my %pdb_xyz = xyz_pdb($file_model, "all");
	# list all the possible pairings
	my %pairs_subsets_info = ();
	my %pairs_subsets_dist = ();
	foreach my $s1 (keys %strands){
		foreach my $s2 (keys %strands){
			next if $s2 eq $s1;
			my @A = split /\s+/, $s1;
			my @B = split /\s+/, $s2;
			next if $A[0] > $B[0];
			next if abs($A[0]-$A[1]) < 2; # implies at least 3 residues
			next if abs($B[0]-$B[1]) < 2; # implies at least 3 residues
			my %strand1_trims = ();
			my %strand2_trims = ();
			# (a) they are already of equal length
			if ($strands{$s1} == $strands{$s2}){  
				$strand1_trims{$s1} = 1;
				$strand2_trims{$s2} = 1;
			}
			# (b) find subsets of s1 to match the length as s2
			elsif ($strands{$s1} >= $strands{$s2}){
				$strand2_trims{$s2} = 1;
				for (my $i = 0; $i <= ($strands{$s1} - $strands{$s2}); $i++){
					confess ":(" if ($A[0]+$i+$strands{$s2}-1) > $A[1];
					$strand1_trims{($A[0]+$i)." ".($A[0]+$i+$strands{$s2}-1)} = 1;
				}
			}
			# (c) find subsets of s2 of same length as s1
			else{
				$strand1_trims{$s1} = 1;
				for (my $i = 0; $i <= ($strands{$s2} - $strands{$s1}); $i++){
					confess ":(" if ($B[0]+$i+$strands{$s1}-1) > $B[1];
					$strand2_trims{($B[0]+$i)." ".($B[0]+$i+$strands{$s1}-1)} = 1;
				}
			}
			my %subsets = ();
			my %possible_slides = qw( 0 1 L1 2 R1 3 L2 4 R2 5 ); # total 10 possible configurations
			my %possible_configs = qw( a 1 b 2 c 3 d 4);
			my %possible_types = qw( P 1 A 2 );
			foreach my $subset1 (keys %strand1_trims){
				foreach my $subset2 (keys %strand2_trims){
					my @X = split /\s+/, $subset1;
					my @Y = split /\s+/, $subset2;
					confess ":( strand size must be equal" if abs($X[0]-$X[1]) != abs($Y[0]-$Y[1]);
					my $i = 1;
					foreach my $slide (sort {$possible_slides{$a} <=> $possible_slides{$b}} keys %possible_slides){
						foreach my $config (sort {$possible_configs{$a} <=> $possible_configs{$b}} keys %possible_configs){
							foreach my $type (keys %possible_types){
								my @A = ($X[0], $X[1], $Y[0], $Y[1]);
								next if ($type eq "A" and not ($config eq "a" or $config eq "b")); # anti-parallel do not have hbond-patterns c and d
								my $no_of_residues = $slide;
								$no_of_residues =~ s/[A-Z]//g;
								# first strand is sliding forward
								if ($slide =~ /^R/){
									$A[1] = $A[1] - $no_of_residues;
									$A[3] = $A[3] - $no_of_residues if $type eq "A";
									$A[2] = $A[2] + $no_of_residues if $type eq "P";
								}
								# first strand is sliding backward
								elsif ($slide =~ /^L/){
									$A[0] = $A[0] + $no_of_residues;
									$A[2] = $A[2] + $no_of_residues if $type eq "A";
									$A[3] = $A[3] - $no_of_residues if $type eq "P";
								}
								confess "ERROR! strand size must be equal!" if abs($A[0]-$A[1]) != abs($A[2]-$A[3]);
								next if abs($A[0]-$A[1]) < 2; # must have at least 3 residues
								next if abs($A[2]-$A[3]) < 2; # must have at least 3 residues
								$subsets{$A[0]." ".$A[1]." ".$A[2]." ".$A[3]." $type $config"} = "[pattern $i] slide $slide type $type config $config " if $type eq "P";
								$subsets{$A[0]." ".$A[1]." ".$A[3]." ".$A[2]." $type $config"} = "[pattern $i] slide $slide type $type config $config " if $type eq "A";
								$i++;
							}
						}
					}
				}
			}
			$pairs_subsets_info{$s1." ".$s2} = \%subsets;
		}
	}
	system_cmd("rm -f pairs_and_distances.log");
	print2file("pairs_and_distances.log", "original_pair subset avg_ca_ca_dist avg_hbond_dist overall_farness pair_information");
	# for each strand pair, for each subset pair, compute average of distances in pdb files, for each A and P
	foreach my $pairing (keys %pairs_subsets_info){
		my %subsets_info = %{$pairs_subsets_info{$pairing}};
		my %subsets_dist = ();
		foreach my $subset (keys %subsets_info){
			my @A = split /\s+/, $subset;
			# Step (1) measure the average backbone Ca-Ca distance
			my $avg_ca_ca = avg_ca_ca_dist(\%pdb_xyz, $A[0], $A[1], $A[2], $A[3], $A[4]);
			print2file("pairs_and_distances.log", sprintf "%15s : %15s : %.3f [%s]", $pairing, $subset, $avg_ca_ca, $subsets_info{$subset});
			next if $avg_ca_ca > $pthres;
			# Step (2) measure the average hbond distance for this configuration
			my $avg_hb_dist = avg_hbond_dist(\%pdb_xyz, $A[0], $A[1], $A[2], $A[3], $A[4], $A[5]);
			print2file("pairs_and_distances.log", sprintf "%15s : %15s : %.3f : %.3f [%s]", $pairing, $subset, $avg_ca_ca, $avg_hb_dist, $subsets_info{$subset});
			next if $avg_hb_dist > $pthres_hbond;
			# Step (3) average them.
			my $farness = ($avg_ca_ca + $avg_hb_dist)/2;
			print2file("pairs_and_distances.log", sprintf "%15s : %15s : %.3f : %.3f : %.3f : %s (selected)", $pairing, $subset, $avg_ca_ca, $avg_hb_dist, $farness, $subsets_info{$subset});
			$subsets_dist{$subset} = $farness if not defined $subsets_dist{$subset};
			$subsets_dist{$subset} = $farness if $subsets_dist{$subset} > $farness;
		}
		$pairs_subsets_dist{$pairing} = \%subsets_dist if scalar keys %subsets_dist;
	}
	if (not scalar keys %pairs_subsets_dist){
		print2file("pairs_and_distances.log", ":( not any eligible pairs!");
		return;
	}
	# best looking pairs
	print2file("pairs_and_distances.log", "\nBest looking pairs:");
	foreach my $pairing (sort keys %pairs_subsets_dist){
		my %subsets_info = %{$pairs_subsets_info{$pairing}};
		my %subsets_dist = %{$pairs_subsets_dist{$pairing}};
		foreach (sort {$subsets_dist{$a} <=> $subsets_dist{$b}} keys %subsets_dist){
			print2file("pairs_and_distances.log", sprintf "%15s : %15s : %s : %s", $pairing, $_, $subsets_dist{$_}, $subsets_info{$_});
			last;
		}
	}
	# rank pairs
	my %pairs_rank = ();
	foreach my $pairing (keys %pairs_subsets_dist){
		my %subsets_dist = %{$pairs_subsets_dist{$pairing}};
		my $best_subset = (sort {$subsets_dist{$a} <=> $subsets_dist{$b}} keys %subsets_dist)[0];
		$pairs_rank{$pairing} = $subsets_dist{$best_subset};
	}
	# keep the best subset having minimum distance AND not violating the hbond pattern
	my %hbonds = ();
	my %final_pairing_info = ();
	my %final_pairing_dist = ();
	print2file("pairs_and_distances.log", "\nSelecting pairs:");
	foreach my $pairing (sort {$pairs_rank{$a} <=> $pairs_rank{$b}} keys %pairs_rank){
		print2file("pairs_and_distances.log", "\t$pairing");
		my %subsets_info = %{$pairs_subsets_info{$pairing}};
		my %subsets_dist = %{$pairs_subsets_dist{$pairing}};
		# find the direction of top 1
		my $top1 = (sort {$subsets_dist{$a} <=> $subsets_dist{$b}} keys %subsets_dist)[0];
		my @T = split /\s+/, $top1;
		my $direction_of_top1 = $T[4];
		foreach my $best_subset (sort {$subsets_dist{$a} <=> $subsets_dist{$b}} keys %subsets_dist){
			print2file("pairs_and_distances.log", "\t\t$best_subset", "");
			my @T = split /\s+/, $best_subset;
			if($T[4] ne $direction_of_top1){
				print2file("pairs_and_distances.log", " -> direction is different.. ignoring rest!", "");
				last;
			}
			my $addable = 0;
			$addable = is_addable_to_hbond_list($best_subset, %hbonds);
			if ($addable){
				my @S = split /\s+/, $best_subset;
				my %this_hbonds = make_ideal_hbonds($S[0], $S[1], $S[2], $S[3], $S[4], $S[5]);
				foreach (keys %this_hbonds){
					$hbonds{$_} = $S[4];
				}
				$final_pairing_info{$best_subset} = $subsets_info{$best_subset};
				$final_pairing_dist{$best_subset} = $subsets_dist{$best_subset};
				last;
			}
			else{
				print2file("pairs_and_distances.log", " -> could not add because of hbond issues!", "");
			}
			print2file("pairs_and_distances.log", "");
		}
		print2file("pairs_and_distances.log", "");
	}
	if (not scalar keys %hbonds){
		print2file("pairs_and_distances.log", "Quitting because no strands could be paired!");
		return;
	}
	print2file("pairs_and_distances.log", "\nHbonds:");
	foreach (keys %hbonds){
		print2file("pairs_and_distances.log", $_." ".$hbonds{$_});
	}
	# log final pairs
	print2file("pairs_and_distances.log", "\nFinal pairs:");
	foreach (sort {$final_pairing_dist{$a} <=> $final_pairing_dist{$b}} keys %final_pairing_dist){
		print2file("pairs_and_distances.log", sprintf "%15s : %.3f : %s", $_, $final_pairing_dist{$_}, $final_pairing_info{$_});
	}
	# print to the output files, finally
	system_cmd("rm -f $file_pair");
	foreach (sort {$final_pairing_dist{$a} <=> $final_pairing_dist{$b}} keys %final_pairing_dist){
		my @A = split /\s+/, $_;
		# converting distance to confidence
		print2file("pairing.txt", sprintf $A[0]." ".$A[1]." ".$A[2]." ".$A[3]." ".$A[4]." %.2f ".$A[5]."  # average distance ".$final_pairing_dist{$_}." [".$final_pairing_info{$_}."]", (1/$final_pairing_dist{$_}));
	}
	system_cmd("rm -f $file_hbond");
	foreach (sort keys %hbonds){
		print2file($file_hbond, $_." ".$hbonds{$_});
	}
}

sub avg_ca_ca_dist{
	my %pdb_xyz = %{shift @_};
	my $a = shift;
	my $b = shift;
	my $c = shift;
	my $d = shift;
	my $T = shift;
	confess "ERROR! Empty hash!" if not scalar keys %pdb_xyz;
	my $sum_dist = 0.0;
	my $count_dist = 0;
	confess "ERROR! type must be A or P in $a $b $c $d $T!" if not ($T eq "A" or $T eq "P");
	if ($T eq "P"){
		foreach (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i++, $j++){
			my @row1 = split /\s+/, $pdb_xyz{$i." CA"};
			my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
			confess "$i not defined in ".$pdb_xyz{$i." CA"} if not defined $x1;
			my @row2 = split /\s+/, $pdb_xyz{$j." CA"};
			my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
			confess "$j not defined in ".$pdb_xyz{$j." CA"} if not defined $x2;
			$sum_dist += sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
			$count_dist++;
		}
	}
	else{
		confess "ERROR! anti-parallel config must have c > d in $a $b $c $d $T!" if $c < $d;
		foreach (my $i = $a, my $j = $c; $i <= $b and $j >= $d; $i++, $j--){
			my @row1 = split /\s+/, $pdb_xyz{$i." CA"};
			my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
			my @row2 = split /\s+/, $pdb_xyz{$j." CA"};
			my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
			$sum_dist += sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
			$count_dist++;
		}
	}
	return sprintf "%.3f", ($sum_dist/$count_dist);
}

sub avg_hbond_dist{
	my %pdb_xyz = %{shift @_};
	my $a = shift;
	my $b = shift;
	my $c = shift;
	my $d = shift;
	my $T = shift;
	my $C = shift;
	my $sum_dist = 0.0;
	my $count_dist = 0;
	my %ideal_hbonds = make_ideal_hbonds($a, $b, $c, $d, $T, $C);
	foreach (keys %ideal_hbonds){
		# some complications here because (a) Proline, and (b) First residue, do not have H atom.
		my @H = split /\s+/, $_;
		my @row1;
		my @row2;
		if(defined $defined_atoms{$H[0]." ".$H[1]}){
			@row1 = split /\s+/, $pdb_xyz{$H[0]." ".$H[1]};
		}
		else{
			confess "ERROR! N must be defined when H is not for ".$H[0] if not defined $defined_atoms{$H[0]." N"};
			@row1 = split /\s+/, $pdb_xyz{$H[0]." N"};
		}
		if(defined $defined_atoms{$H[2]." ".$H[3]}){
			@row2 = split /\s+/, $pdb_xyz{$H[2]." ".$H[3]};
		}
		else{
			confess "ERROR! N must be defined when H is not for ".$H[2] if not defined $defined_atoms{$H[2]." N"};
			@row2 = split /\s+/, $pdb_xyz{$H[2]." N"};
		}
		my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
		confess "x1 not defined in $_" if not defined $x1;
		my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
		confess "x2 not defined in $_" if not defined $x2;
		$sum_dist += sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
		$count_dist++;
	}
	return sprintf "%.3f", ($sum_dist/$count_dist);
}

sub update_sec_using_pairing_info{
	my $file_sec = shift;
	my $file_pairing_info = shift;
	my $file_updated_ss = shift;
	confess "ERROR! file sec $file_sec does not exist!" if not -f $file_sec;
	confess "ERROR! file pairing info $file_pairing_info does not exist!" if not -f $file_pairing_info;
	confess "ERROR! file updated ss not defined!" if not defined $file_updated_ss;
	my %pairing_info = pairing_info2hash($file_pairing_info);
	my %strands = strand_list($file_sec);
	my %paired_residues = ();
	# If the pairing information is out of the strand residues, extend the secondary structure information
	foreach my $pair (keys %pairing_info){
		my @A = split /\s+/, $pair;
		confess "ERROR! $pair pair not defined!" if not (defined $A[0] and defined $A[1]);
		for(my $i = $A[0]; $i <= $A[1]; $i++){
			$paired_residues{$i} = 1;
		}
		for(my $i = $A[2]; $i <= $A[3]; $i++){
			$paired_residues{$i} = 1;
		}
		for(my $i = $A[3]; $i <= $A[2]; $i++){
			$paired_residues{$i} = 1;
		}
	}
	my %res_ss = fasta2residues_hash($file_sec);
	foreach (sort {$a <=> $b} keys %res_ss){
		$res_ss{$_} = "E" if defined $paired_residues{$_};
	}
	my $line_sec = "";
	foreach (sort {$a <=> $b} keys %res_ss){
		$line_sec .= $res_ss{$_};
	}
	my $line_updated_E = "";
	foreach (sort {$a <=> $b} keys %res_ss){
		if (defined $paired_residues{$_}){
			$line_updated_E .= "E";
			next;
		}
		$line_updated_E .= "-";
	}
	confess "\n $line_sec :( " if length(seq_fasta($file_sec)) ne length($line_sec);
	print "\nStage 2 update strands using pairing information..";
	print "\n".seq_fasta($file_sec)." [".basename($file_sec)."]";
	print "\n".$line_sec." [updated]";
	print "\n".$line_updated_E." [paired strands]\n";
	system_cmd("rm -f $file_updated_ss");
	print2file($file_updated_ss, ">$id");
	print2file($file_updated_ss, $line_sec);
}

sub pairing_info2hash{
	my $file_pairing_info = shift;
	confess "ERROR! file pairing info $file_pairing_info does not exist!" if not -f $file_pairing_info;
	my %pairing_rows = ();
	open PAIR, $file_pairing_info or confess $!;
	while (<PAIR>){
		chomp $_;
		next if not defined $_;
		$_ =~ s/^\s+//;
		my @A = split /\s+/, $_;
		confess "ERROR! all columns not defined in row $_ in pairing file!" if not (defined $A[0] and defined $A[1] and defined $A[2] and defined $A[3] and defined $A[4]); 
		$pairing_rows{$_} = $A[0];
	}
	close PAIR;
	return %pairing_rows;
}

sub paired_residues{
	my $file_pairing_info = shift;
	confess "ERROR! file pairing info $file_pairing_info does not exist!" if not -f $file_pairing_info;
	my %pairing_info = pairing_info2hash($file_pairing_info);
	my %paired_residues = ();
	foreach my $pair (keys %pairing_info){
		my @A = split /\s+/, $pair;
		for(my $i = $A[0]; $i <= $A[1]; $i++){
			$paired_residues{$i} = $A[4];
		}
		for(my $i = $A[2]; $i <= $A[3]; $i++){
			$paired_residues{$i} = $A[4];
		}
		for(my $i = $A[3]; $i <= $A[2]; $i++){
			$paired_residues{$i} = $A[4];
		}
	}
	return %paired_residues;
}

sub load_hbonds{
	my $file_hbonds = shift;
	confess "ERROR! file hbonds $file_hbonds does not exist!" if not -f $file_hbonds;
	my %hbond_rows = ();
	open HBOND, $file_hbonds or confess $!;
	while (<HBOND>){
		my @A = split /\s+/, $_;
		$hbond_rows{$A[0]." ".$A[1]." ".$A[2]." ".$A[3]} = $A[4];
	}
	close HBOND;
	return %hbond_rows;
}

sub hbonds2connectors{
	my $ref_hbonds = shift;
	my %hbonds = %{$ref_hbonds};
	my %connectors = ();
	foreach (sort keys %hbonds){
		my @H = split /\s+/, $_;
		$connectors{$H[0]." ".$H[1]} = 1;
		$connectors{$H[2]." ".$H[3]} = 1;
	}
	return %connectors;
}

sub strand_list{
	my $file_sec = shift;
	confess "ERROR! file ss $file_sec does not exist!" if not -f $file_sec;
	my %each_residue_ss = fasta2residues_hash($file_sec);
	my %strand_list = ();
	my $start = 1;
	foreach my $r (sort {$a <=> $b} keys %each_residue_ss){
		next if $each_residue_ss{$r} ne "E";
		$start = $r if (defined $each_residue_ss{$r-1} && $each_residue_ss{$r-1} ne "E");
		$strand_list{$start." ".$r} = $r - $start + 1 if (defined $each_residue_ss{$r+1} && $each_residue_ss{$r+1} ne "E");
	}
	return %strand_list;
}

sub make_ideal_hbonds{
	my $a = shift;
	my $b = shift;
	my $c = shift;
	my $d = shift;
	my $T = shift;
	my $C = shift; # configuration a or b (look at the pictures folder) 
	confess "ERROR! unknown type! type must be P or A" if not ($T eq "P" or $T eq "A"); 
	confess "ERROR! unknown configuration! configuration must be a or b" if not ($C eq "a" or $C eq "b" or $C eq "c" or $C eq "d"); 
	my %hbond_list = ();
	if ($T eq "A"){
		confess "ERROR! in anti-parallel strand pairing in [$T $a-$b:$c-$d]. In [a-b:c-d] c must be greater than d" if $c < $d;
		if ($C eq "a"){
			for (my $i = $a, my $j = $c; $i <= $b and $j >= $d; $i = $i+2,  $j = $j-2){
				$hbond_list{"$i H $j O"} = $T;
				$hbond_list{"$i O $j H"} = $T;
			}
		}
		else{
			for (my $i = $a+1, my $j = $c-1; $i <= $b and $j >= $d; $i = $i+2,  $j = $j-2){
				$hbond_list{"$i H $j O"} = $T;
				$hbond_list{"$i O $j H"} = $T;
			}
		}
	}
	else{
		confess "ERROR! in parallel strand pairing in [$T $a-$b:$c-$d]. In [a-b:c-d] c must be less than d" if $d < $c;
		# 9/2/2014:
		# I noticed that my output number of EEE in parallel pairs is less than input, in my "test1" results.
		# Looking at how Chimera shows beta sheets and how DSSP calculates E for a PDB, it seems this:
		# In a native PDB (I used 1a05a's first parallel), if there are two strands 3-6 and 36-39, forming a parallel beta sheet, then the actual residues involved in hydrogen bonding are 2-6 and 36-40.
		# But, somehow, DSSP and Chimera both ignore the first hydrogen-bonded residue and last hydrogen-bonded residue.
		# This makes me start hydrogen bonding from $a-1 and go up to $d+1
		# 12/24/2014:
		# Seeing 1yipA's hydrogen bonding 122N-158O, in reindexed pdb, and looking at DSSP's interpretations, it 
		# is clear that pairs starting with N-O hydrogen bonds are ignored, and so my assumption of ideal h-bond starting with O-H bonding may
		# be all right.
		# 12/26/2014:
		# I know that this is not appropriate in terms of reconstruction.
		# I must go up to 1 more and so should j. But, that way, it will work great for reconstruction but not for prediction (my intuition)
		if ($C eq "a"){
			for (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i = $i+2, $j = $j+2){
				$hbond_list{"$i O ".($j+1)." H"} = $T if $j+1 <= $d;
				$hbond_list{"".($i+2)." H ".($j+1)." O"} = $T if ($i+2 <= $b and $j+1 <= $d);
			}
		}
		elsif ($C eq "b"){
			for (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i = $i+2, $j = $j+2){
				$hbond_list{($i+1)." H $j O"} = $T if $i+1 <= $b;
				$hbond_list{($i+1)." O ".($j+2)." H"} = $T if ($i+1 <= $b and $j+2 <= $d);
			}
		}
		elsif ($C eq "c"){
			for (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i = $i+2, $j = $j+2){
				$hbond_list{"$i H ".($j+1)." O"} = $T if $j+1 <= $d;
				$hbond_list{"".($i+2)." O ".($j+1)." H"} = $T if ($i+2 <= $b and $j+1 <= $d);
			}
		}
		else{
			for (my $i = $a, my $j = $c; $i <= $b and $j <= $d; $i = $i+2, $j = $j+2){
				$hbond_list{($i+1)." O $j H"} = $T if $i+1 <= $b;
				$hbond_list{($i+1)." H ".($j+2)." O"} = $T if ($i+1 <= $b and $j+2 <= $d);
			}
		}
	}
	return %hbond_list;
}

sub is_addable_to_hbond_list{
	my $pairing = shift;
	my %hbonds = @_;
	confess "ERROR! pairing not defined!" if not defined $pairing;
	my %connectors = ();
	my @P = split /\s+/, $pairing;
	my %ideal_hbonds = make_ideal_hbonds($P[0], $P[1], $P[2], $P[3], $P[4], $P[5]);
	foreach (keys %hbonds){
		my @H = split /\s+/, $_;
		$connectors{$H[0]." ".$H[1]} = 1;
		$connectors{$H[2]." ".$H[3]} = 1;
	}
	my $flag_problem = 0;
	foreach (keys %ideal_hbonds){
		my @H = split /\s+/, $_;
		$flag_problem = 1 if defined $connectors{$H[0]." ".$H[1]};
		$flag_problem = 1 if defined $connectors{$H[2]." ".$H[3]};
	}
	return 1 if not $flag_problem;
	return 0;
}

sub system_cmd{
	my $command = shift;
	my $log = shift;
	confess "EXECUTE [$command]?\n" if (length($command) < 5  and $command =~ m/^rm/);
	if(defined $log){
		system("$command &> $log");
	}
	else{
		system($command);
	}
	if($? != 0){
		my $exit_code  = $? >> 8;
		confess "ERROR!! Could not execute [$command]! \nError message: [$!]";
	}
}

sub count_lines{
	my $file = shift;
	my $lines = 0;
	return 0 if not -f $file;
	open FILE, $file or confess "ERROR! Could not open $file! $!";
	while (<FILE>){
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		next if not defined $_;
		next if length($_) < 1;
		$lines ++;
	}
	close FILE;
	return $lines;
}

sub print2file{
	my $file = shift;
	my $message = shift;
	my $newline = shift;
	$newline = "\n" if not defined $newline;
	if (-f $file){
		open  FILE, ">>$file" or confess $!;
		print FILE $message.$newline;
		close FILE;
	}
	else{
		open  FILE, ">$file" or confess $!;
		print FILE $message.$newline;
		close FILE;
	}
}

sub tbl2rows_hash{
	my $file_tbl = shift;
	confess "ERROR! tbl file $file_tbl does not exist!" if not -f $file_tbl;
	my %noe = ();
	open NOETBL, $file_tbl or confess $!;
	while (<NOETBL>){
		$_ =~ s/^\s+//;
		next if $_ !~ m/^assign/;
		$_ =~ s/^\s+//;
		$_ =~ s/\(/ /g;
		$_ =~ s/\)/ /g;
		$noe{$_} = 1;
	}
	close NOETBL;
	confess "$file_tbl seems empty!" if not %noe;
	return %noe;
}

sub ssnoe_tbl_min_pdb_dist{
	my $file_tbl = shift;
	my $file_pdb = shift;
	confess "ERROR! tbl file $file_tbl does not exist!" if not -f $file_tbl;
	confess "ERROR! tbl pdb $file_pdb does not exist!" if not -f $file_pdb;
	my %noe_hash = ();
	open NOETBL, $file_tbl or confess $!;
	while (<NOETBL>){
		chomp $_;
		$_ =~ s/^\s+//;
		$_ =~ s/\)/ /g;
		$_ =~ s/\(/ /g;
		confess "$_" if $_ !~ m/^assign/;
		my @C = split /\s+/, $_;
		if ($C[6] eq "or" and $C[17] eq "or"){
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5]." ".$C[8]." ".$C[11];
			$noe_hash{$_}{"right"}    = $C[13]." ".$C[16]." ".$C[19]." ".$C[22];
			$noe_hash{$_}{"distance"} = $C[23]." ".$C[24]." ".$C[25];
		}
		elsif ($C[6] eq "or" and $C[17] ne "or"){
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5]." ".$C[8]." ".$C[11];
			$noe_hash{$_}{"right"}    = $C[13]." ".$C[16];
			$noe_hash{$_}{"distance"} = $C[17]." ".$C[18]." ".$C[19];
		}
		elsif ($C[6] ne "or" and $C[11] eq "or"){
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5];
			$noe_hash{$_}{"right"}    = $C[7]." ".$C[10]." ".$C[13]." ".$C[16];
			$noe_hash{$_}{"distance"} = $C[17]." ".$C[18]." ".$C[19];
		}
		else{
			$noe_hash{$_}{"left"}     = $C[2]." ".$C[5];
			$noe_hash{$_}{"right"}    = $C[7]." ".$C[10];
			$noe_hash{$_}{"distance"} = $C[11]." ".$C[12]." ".$C[13];
		}
	}
	close NOETBL;
	confess "$file_tbl seems empty!" if not %noe_hash;
	my %xyzPDB = xyz_pdb($file_pdb, "all");
	foreach (sort keys %noe_hash){
		my $left = $noe_hash{$_}{"left"};
		my $right = $noe_hash{$_}{"right"};
		my $distance = $noe_hash{$_}{"distance"};
		my @L = split /\s+/, $left;
		my @R = split /\s+/, $right;
		my @D = split /\s+/, $distance;
		my %left_list = ();
		my %right_list = ();
		for(my $i = 0; $i <= $#L; $i = $i+2){
			$left_list{$L[$i]." ".$L[$i+1]} = 1;
		}
		for(my $i = 0; $i <= $#R; $i = $i+2){
			$right_list{$R[$i]." ".$R[$i+1]} = 1;
		}
		my $distance_pdb = 1000.0;
		foreach my $le (keys %left_list){
			foreach my $ri (keys %right_list){
				my @L = split /\s+/, $le;
				my @R = split /\s+/, $ri;
				confess "$file_pdb does not have ".$L[0]." ".uc($L[1])."\n" if not defined $xyzPDB{$L[0]." ".uc($L[1])};
				confess "$file_pdb does not have ".$R[0]." ".uc($R[1])."\n" if not defined $xyzPDB{$R[0]." ".uc($R[1])};
				my $d = calc_dist($xyzPDB{$L[0]." ".uc($L[1])}, $xyzPDB{$R[0]." ".uc($R[1])});
				$distance_pdb = $d if $d < $distance_pdb;
			}
		}
		$noe_hash{$_}{"pdb_distance"} = $distance_pdb;
	}
	return %noe_hash;
}

sub coverage_tbl{
	my $file_fasta = shift;
	my $file_tbl = shift;
	my $flag_dihe = shift;
	confess "ERROR! fasta file $file_fasta does not exist!" if not -f $file_fasta;
	confess "ERROR! tbl file $file_tbl does not exist!" if not -f $file_tbl;
	$flag_dihe = 0 if not defined $flag_dihe;
	my $seq = seq_fasta($file_fasta);
	my $L = length $seq;
	my $cov = $seq;
	$cov =~ s/[A-Z]/-/g;
	my %noe = tbl2rows_hash($file_tbl);
	foreach (keys %noe){
		#assign (resid 123 and name CA) (resid  58 and name CA) 3.60 0.10 3.40
		my @C = split /\s+/, $_;
		my $r1 = $C[2]; my $r2 = $C[7];
		# the case when "ca or cb" restraints are provided
		#assign ((resid 123 and name ca) or (resid 123 and name cb)) ((resid 58 and name ca) or (resid 58 and name cb)) 3.60 0.10 3.40
		$r2 = $C[13] if $C[6] eq "or";
		$r2 = $C[17] if $flag_dihe;
		my $c1 = substr $cov, ($r1 - 1), 1;
		my $c2 = substr $cov, ($r2 - 1), 1;
		if ($c1 eq "-" ){
			$c1 = 1;
		}
		elsif ($c1 eq "*" ){
			$c1 = "*";
		}
		else{
			$c1++;
			$c1 = "*" if ($c1 > 9);
		}
		if ($c2 eq "-" ){
			$c2 = 1;
		}
		elsif ($c2 eq "*" ){
			$c2 = "*";
		}
		else{
			$c2++;
			$c2 = "*" if ($c2 > 9);
		}
		substr $cov, ($r1 - 1), 1, $c1;
		substr $cov, ($r2 - 1), 1, $c2;
	}
	my $cov2 = $cov;
	$cov2 =~ s/-//g;
	return sprintf "$cov [%12s : %3s restraints touching %s residues]", basename($file_tbl), (scalar keys %noe), length($cov2);
}

sub count_ss_match{
	my $file_pdb = shift;
	my $file_sec = shift;
	my $char = shift; # secondary structure element
	confess "ERROR! file pdb $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! file sec $file_sec does not exist!" if not -f $file_sec;
	confess "Invalid character!" if not ($char eq "H" or $char eq "E" or $char eq "C");
	my %residue_ss1 = fasta2residues_hash($file_sec);
	my %residue_ss2 = dssp_result($file_pdb, "ss");
	my $count = 0;
	foreach my $r (keys %residue_ss1){
		next if $residue_ss1{$r} ne $char;
		$count++ if $residue_ss1{$r} eq $residue_ss2{$r};
	}
	return $count;
}

sub count_satisfied_tbl_rows{
	my $file_pdb = shift;
	my $file_tbl = shift;
	my $type    = shift;
	confess "ERROR! file pdb $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! file tbl $file_tbl does not exist!" if not -f $file_tbl;
	confess "ERROR! Invalid type!" if not ($type eq "noe" or $type eq "dihedral");
	my $count = 0;
	my $total = 0;
	my %log_rows = ();
	if ($type eq "dihedral"){
		my %noe = tbl2rows_hash($file_tbl);;
		my %phi = dssp_result($file_pdb, "phi");
		my %psi = dssp_result($file_pdb, "psi");
		foreach (sort keys %noe){
			chomp $_;
			#assign resid   1 and name ca resid   2 and name ca resid   3 and name ca resid   4 and name ca  5.0 -164.6 5.0 2 ! alpha for 2
			my @C = split /\s+/, $_;
			my $angle_true = 0.0;
			if (uc($C[5]) eq "C" and uc($C[10]) eq "N" and uc($C[15]) eq "CA" and uc($C[20]) eq "C"){
				$angle_true = $phi{$C[2]};
				confess "Undefined dssp phi angle for ".$C[2]." in $file_pdb!" if not defined $angle_true;
			}
			elsif (uc($C[5]) eq "N" and uc($C[10]) eq "CA" and uc($C[15]) eq "C" and uc($C[20]) eq "N"){
				$angle_true = $psi{$C[2]};
				confess "Undefined dssp phi angle for ".$C[2]." in $file_pdb!" if not defined $angle_true;
			}
			else{
				confess "Undefined dihedral angle $_";
			}
			confess "Undefined dihedral angle $_" if not defined $C[22];
			my $viol_flag = 1;
			my $d = abs($angle_true - $C[22]);
			$d = abs(360.0 - $d) if $d > 180.0;
			if ($d < ($C[23] + 2.0)){
				$count ++;
				$viol_flag = 0;
			}
			$total++;
			$log_rows{(sprintf "%3s\t%.2f\t%.2f # $_", $viol_flag, $d, $angle_true)} =  $viol_flag;
		}
	}
	else{
		my %tbl_hash = ssnoe_tbl_min_pdb_dist($file_tbl, $file_pdb);
		foreach (sort keys %tbl_hash){
			my $viol_flag = 1;
			my $distance = $tbl_hash{$_}{"distance"};
			my @D = split /\s+/, $distance;
			my $pdb_distance = $tbl_hash{$_}{"pdb_distance"};
			my $deviation = $pdb_distance - ( $D[0] + $D[2] );
			if ($pdb_distance < ( $D[0] + $D[2] + 0.2) ){
				$count ++;
				$viol_flag = 0;
				$deviation = 0.0;
			}
			if ($pdb_distance < ( $D[0] - $D[1] - 0.2) ){
				$count--;
				$viol_flag = 1;
				$deviation = -($D[0] - $D[1] - $pdb_distance);
			}
			
			$log_rows{(sprintf "%3s\t%.2f\t%.2f # $_", $viol_flag, $deviation, $pdb_distance)} = $viol_flag;
			$total++;
		}
	}
	my $file_viol = "".basename($file_tbl,".tbl")."_violation.txt";
	print2file($file_viol, "#NOE violation check; $file_pdb against $file_tbl");
	print2file($file_viol, "#violation-flag, deviation, actual-measurement, Input-NOE-restraint");
	foreach (sort {$log_rows{$b} <=> $log_rows{$a}} keys %log_rows){
		print2file($file_viol, $_);
	}
	return $count."/".$total;
}

sub noe_tbl_violation_coverage{
	my $file_pdb = shift;
	my $file_tbl = shift;
	confess "ERROR! file pdb $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! file tbl $file_tbl does not exist!" if not -f $file_tbl;
	my $cov = seq_chain($file_pdb);
	$cov =~ s/[A-Z]/-/g;
	my %tbl_hash = ssnoe_tbl_min_pdb_dist($file_tbl, $file_pdb);
	my %xyz = xyz_pdb($file_pdb, "all");
	foreach (sort keys %tbl_hash){
		my $left = $tbl_hash{$_}{"left"};
		my $right = $tbl_hash{$_}{"right"};
		my $distance = $tbl_hash{$_}{"distance"};
		my @L = split /\s+/, $left;
		my @R = split /\s+/, $right;
		my @D = split /\s+/, $distance;
		my $pdb_distance = $tbl_hash{$_}{"pdb_distance"};
		substr $cov, $L[0] - 1, 1, "x" if $pdb_distance > $D[0] + $D[2] + 0.2;
		substr $cov, $R[0] - 1, 1, "x" if $pdb_distance > $D[0] + $D[2] + 0.2;
		substr $cov, $L[0] - 1, 1, "x" if $pdb_distance < $D[0] - $D[1] - 0.2;
		substr $cov, $R[0] - 1, 1, "x" if $pdb_distance < $D[0] - $D[1] - 0.2;
	}
	return $cov;
}

sub sum_noe_dev{
	my $file_pdb = shift;
	my $file_tbl = shift;
	confess "ERROR! file pdb $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! file tbl $file_tbl does not exist!" if not -f $file_tbl;
	my $sum_dev = 0.0;
	my %tbl_hash = ssnoe_tbl_min_pdb_dist($file_tbl, $file_pdb);
	foreach (sort keys %tbl_hash){
		my $viol_flag = 1;
		my @D = split /\s+/, $tbl_hash{$_}{"distance"};
		my $pdb_distance = $tbl_hash{$_}{"pdb_distance"};
		if ($pdb_distance > ( $D[0] + $D[2] + 0.2) ){
			$sum_dev += $pdb_distance - ( $D[0] + $D[2] );
		}
		if ($pdb_distance < ( $D[0] - $D[1] - 0.2) ){
			$sum_dev += ( $D[0] - $D[1] ) - $pdb_distance;
		}
	}
	return sprintf "%.2f", $sum_dev;
}

sub get_cns_energy{
	my $cns_pdb = shift;
	my $energy_term = shift; # "overall", "vdw", "bon", "noe"
	confess "ERROR! file $cns_pdb does not exist!" if not -f $cns_pdb;
	confess "ERROR! energy term must be overall vdw bon or noe" if not ($energy_term eq "overall" or $energy_term eq "bon" or $energy_term eq "vdw" or $energy_term eq "noe");
	my $value = "X";
	open CHAIN, $cns_pdb or confess $!;
	while(<CHAIN>){
		chomp $_;
		next if $_ !~ m/^REMARK\ $energy_term/;
		$_ =~ s/\s+//g;
		my @C = split /=/, $_;
		$value = $C[1];
	}
	close CHAIN;
	return int($value);
}

sub pdb_each_residue_atoms{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %rnum_atoms = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_atoms{parse_pdb_row($_,"rnum")} = "";
	}
	close CHAIN;
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_atoms{parse_pdb_row($_,"rnum")} .= parse_pdb_row($_,"aname")." ";
	}
	close CHAIN;
	confess "ERROR! no atoms in file $chain!" if not scalar keys %rnum_atoms;
	return %rnum_atoms;
}

sub load_pdb{
	my $dir_chains = shift;
	confess ":( directory $dir_chains does not exist!" if not -d $dir_chains;
	my @pdb_list = <$dir_chains/*.pdb>;
	if(not (@pdb_list)){
		@pdb_list = <$dir_chains/*.ent>;
	}
	confess "ERROR! Directory $dir_chains has no pdb files!\n" unless(@pdb_list);
	return @pdb_list;
}

sub pdb2rnum_rname{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my %rnum_rname = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$rnum_rname{parse_pdb_row($_,"rnum")} = parse_pdb_row($_,"rname");
	}
	close CHAIN;
	confess ":(" if not scalar keys %rnum_rname;
	return %rnum_rname;
}

sub xyz_pdb{
	my $chain = shift;
	my $atom_selection = shift; # ca or cb or all
	$atom_selection = uc($atom_selection);
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	confess "ERROR! Selection must be ca or cb or all" if not ($atom_selection eq "CA" or $atom_selection eq "ALL" or $atom_selection eq "CB");
	my %xyz_pdb = ();
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		$xyz_pdb{"".parse_pdb_row($_,"rnum")." ".parse_pdb_row($_,"aname")} = "".parse_pdb_row($_,"x")." ".parse_pdb_row($_,"y")." ".parse_pdb_row($_,"z");
	}
	close CHAIN;
	confess "ERROR!: xyz_pdb is empty\n" if (not scalar keys %xyz_pdb);
	return %xyz_pdb if $atom_selection eq "ALL";

	my %rnum_rname = pdb2rnum_rname($chain);
	my %selected_xyz = ();
	foreach (sort keys %xyz_pdb){
		my @C = split /\s+/, $_;
		my $this_atom = $atom_selection;
		$this_atom = "CA" if ($atom_selection eq "CB" and $rnum_rname{$C[0]} eq "GLY");
		next if $C[1] ne $this_atom;
		$selected_xyz{$C[0]} = $xyz_pdb{$_};
	}
	confess ":(" if not scalar keys %selected_xyz;
	return %selected_xyz;
}

sub pdb2rr{
	my $chain = shift;
	my $rr = shift;
	my $seq_separation = shift;
	my $dist_threshold = shift;
	my $atom = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	confess "ERROR! file rr not defined!" if not defined $rr;
	$seq_separation = 6 if not defined $seq_separation;
	$dist_threshold = 8 if not defined $dist_threshold;
	$atom = "cb" if not defined $atom;
	my %contacts_pdb = ();
	my %xyz;
	%xyz = xyz_pdb($chain, $atom);
	confess "ERROR! atom type must be ca or cb!" if not scalar keys %xyz;
	foreach my $r1 (sort keys %xyz){
		foreach my $r2 (sort keys %xyz){
			next if ($r1 >= $r2);
			next if (($r2 - $r1) < $seq_separation);
			my $d = calc_dist($xyz{$r1}, $xyz{$r2});
			$contacts_pdb{"$r1 $r2"} = sprintf "%.3f", $d if ($d <= $dist_threshold);
		}
	}
	open RR, ">$rr" or confess $!;
	print RR "".seq_chain_with_gaps($chain)."\n";
	foreach (sort keys %contacts_pdb){
		print RR "$_ 0 $dist_threshold ".$contacts_pdb{$_}."\n";
	}
	close RR;
}

sub seq_chain_with_gaps{
	my $chain = shift;
	my $flag = shift; # flag 1 if the left side dashes of the sequence are not wanted
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $start = 1;
	# if flagged, keep start for trimming
	if (defined $flag){
		open CHAIN, $chain or confess $!;
		while(<CHAIN>){
			next if $_ !~ m/^ATOM/;
			next unless (parse_pdb_row($_,"aname") eq "CA");
			confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
			$start = parse_pdb_row($_,"rnum");
			last;
		}
		close CHAIN;
	}
	# 1.find end residue number
	my $end;
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		next unless (parse_pdb_row($_,"aname") eq "CA");
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		$end = parse_pdb_row($_,"rnum");
	}
	close CHAIN;
	# 2.initialize
	my $seq = "";
	for (my $i = 1; $i <= $end; $i++){
		$seq .= "-";
	}
	# 3.replace with residues
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		next unless (parse_pdb_row($_,"aname") eq "CA");
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $rnum = parse_pdb_row($_,"rnum");
		$rnum =~ s/[A-G]//g;
		substr $seq, ($rnum - 1), 1, $AA3TO1{parse_pdb_row($_,"rname")}; 
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return (substr $seq, $start - 1);
}

sub seq_chain{
	my $chain = shift;
	confess "ERROR! file $chain does not exist!" if not -f $chain;
	my $seq = "";
	open CHAIN, $chain or confess $!;
	while(<CHAIN>){
		next if $_ !~ m/^ATOM/;
		next unless (parse_pdb_row($_,"aname") eq "CA");
		confess "ERROR!: ".parse_pdb_row($_,"rname")." residue not defined! \nFile: $chain! \nLine : $_" if (not defined $AA3TO1{parse_pdb_row($_,"rname")});
		my $res = $AA3TO1{parse_pdb_row($_,"rname")};
		$seq .= $res;
	}
	close CHAIN;
	confess "$chain has less than 1 residue!" if (length($seq) < 1);
	return $seq;
}

sub parse_pdb_row{
	my $row = shift;
	my $param = shift;
	my $result;
	$result = substr($row,6,5) if ($param eq "anum");
	$result = substr($row,12,4) if ($param eq "aname");
	$result = substr($row,16,1) if ($param eq "altloc");
	$result = substr($row,17,3) if ($param eq "rname");
	$result = substr($row,22,5) if ($param eq "rnum");
	$result = substr($row,26,1) if ($param eq "insertion");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	confess "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}

sub clash_count{
	my $file_pdb = shift;
	my $threshold = shift;
	confess "ERROR! pdb file $file_pdb not defined!" if not -f $file_pdb;
	confess "ERROR! threshold not defined!" if not defined $threshold;
	my $count = 0;
	my %ca_xyz = xyz_pdb($file_pdb, "ca");
	foreach my $r1 (sort keys %ca_xyz){
		my @R1 = split /\s+/, $ca_xyz{$r1};
		my $x1 = $R1[0]; my $y1 = $R1[1]; my $z1 = $R1[2];
		foreach my $r2 (sort keys %ca_xyz){
			next if $r1 >= $r2;
			my @R2 = split /\s+/, $ca_xyz{$r2};
			my $x2 = $R2[0]; my $y2 = $R2[1]; my $z2 = $R2[2];
			my $d = sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
			if ($d <= $threshold){
				$count++;
			}
		}
	}
	return $count;
}

sub calc_dist{
	my $x1y1z1 = shift;
	my $x2y2z2 = shift;
	confess "ERROR! x1y1z1 not defined!" if !$x1y1z1;
	confess "ERROR! x2y2z2 not defined!" if !$x2y2z2;
	my @row1 = split(/\s+/, $x1y1z1);
	my $x1 = $row1[0]; my $y1 = $row1[1]; my $z1 = $row1[2];
	confess "ERROR! One of the xyz in $x1y1z1 is not a number!" if (not looks_like_number($x1) or not looks_like_number($y1) or not looks_like_number($z1));
	my @row2 = split(/\s+/, $x2y2z2);
	my $x2 = $row2[0]; my $y2 = $row2[1]; my $z2 = $row2[2];
	confess "ERROR! One of the xyz in $x2y2z2 is not a number!" if (not looks_like_number($x2) or not looks_like_number($y2) or not looks_like_number($z2));
	my $d = sprintf "%.3f", sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
	return $d;
}

sub rr2contacts_hash{
    my $file_rr = shift;
    my $seq_sep = shift;
    my $count = shift;
    my $hashvalue = shift;
	confess "ERROR! file_rr $file_rr does not exist!" if not -f $file_rr;
	confess "ERROR! seq_sep not defined!" if !$seq_sep;
	$count = 100000 if not defined $count;
	$hashvalue = "threshold" if not defined $hashvalue;
	my %contacts = ();
	open RR, $file_rr or confess $!;
	while(<RR>){
		next unless ($_ =~ /[0-9]/);
		$_ =~ s/^\s+//;
		next unless ($_ =~ /^[0-9]/);
		my @C = split(/\s+/, $_);
		confess "ERROR! Expecting a pair in row [".$_."]!\n" if (not defined $C[0] || not defined $C[1]);
		next if (abs($C[1] - $C[0]) < $seq_sep);
		if(defined $C[3]){
			$contacts{$C[0]." ".$C[1]} = $C[3];
			$contacts{$C[0]." ".$C[1]} = $C[4] if $hashvalue eq "confidence";
			$contacts{$C[0]." ".$C[1]} = $C[2] if $hashvalue eq "lowerbound";
		}
		else{
			confess "ERROR! Confidence column not defined in row [".$_."] in file $file_rr!\n";
		}
		last if (scalar keys %contacts) == $count;
	}
	close RR;
	return %contacts;
}

sub seq_rr{
	my $file_rr = shift;
	confess ":(" if not -f $file_rr;
	my $seq;
	open RR, $file_rr or confess "ERROR! Could not open $file_rr! $!";
	while(<RR>){
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$_ =~ s/^\s+//;
		next if ($_ =~ /^PFRMAT/);
		next if ($_ =~ /^TARGET/);
		next if ($_ =~ /^AUTHOR/);
		next if ($_ =~ /^SCORE/); 
		next if ($_ =~ /^REMARK/);
		next if ($_ =~ /^METHOD/);
		next if ($_ =~ /^MODEL/); 
		next if ($_ =~ /^PARENT/);
		last if ($_ =~ /^TER/);   
		last if ($_ =~ /^END/);
		# Now, I can directly merge to RR files with sequences on top
		last if ($_ =~ /^[0-9]/);
		$seq .= $_;
	}
	close RR;
	confess ":( no sequence header in $file_rr" if not defined $seq;
	return $seq;
}

sub rr2tbl{
	my $file_rr  = shift;
	my $file_tbl = shift;
	my $rrtype   = shift;
	confess ":(" if not -f $file_rr;
	confess ":(" if not $file_tbl;
	confess ":(" if not ($rrtype eq "ca" or $rrtype eq "cb");
	my %r1a1r2a2 = rr2r1a1r2a2($file_rr, $rrtype);
	my %lowerbound = rr2contacts_hash($file_rr, 1, 100000, "lowerbound");
	my %upperbound = rr2contacts_hash($file_rr, 1, 100000);
	my %rows_and_weights = ();
	foreach (keys %r1a1r2a2){
		my @C = split /\s+/, $_;
		my $lbound   = $lowerbound{$C[0]." ".$C[2]};
		my $ubound   = $upperbound{$C[0]." ".$C[2]};
		my $distance = sprintf("%.2f", 3.6);
		my $negdev   = sprintf("%.2f", 0.1);
		my $posdev   = sprintf("%.2f", ($r1a1r2a2{$_} - 3.6));
		# This is probably a non-contact information
		if ($lbound > 4){
			$distance = sprintf("%.2f", ($lbound + $r1a1r2a2{$_})/2);
			$negdev   = sprintf("%.2f", $distance - $lbound);
			$posdev   = sprintf("%.2f", $distance - $lbound);
		}
		$rows_and_weights{(sprintf "assign (resid %3d and name %2s) (resid %3d and name %2s) %.2f %.2f %.2f", $C[0], $C[1], $C[2], $C[3], $distance, $negdev, $posdev)} = $C[0];
	}
	system_cmd("rm -f $file_tbl");
	foreach (sort {$rows_and_weights{$a} <=> $rows_and_weights{$b}} keys %rows_and_weights){
		print2file($file_tbl, $_);
	}
}

sub rr2r1a1r2a2{
	my $file_rr  = shift;
	my $rrtype = shift;
	confess "ERROR! rr file $file_rr does not exist!" if not -f $file_rr;
	confess "ERROR! rrtype must be ca or cb" if not ($rrtype eq "ca" or $rrtype eq "cb");
	my %contacts = rr2contacts_hash($file_rr, 1);
	my %r1a1r2a2 = ();
	my $seq = seq_rr($file_rr);
	foreach (keys %contacts) {
		my @pair = split(/\s+/,$_);
		my $r1 = substr($seq, ($pair[0] -1), 1);
		my $r2 = substr($seq, ($pair[1] -1), 1);
		my $ca1 = $rrtype;
		my $ca2 = $rrtype;
		if ($rrtype eq "cb"){
			$ca1 = $CBATOM{$r1};
			$ca2 = $CBATOM{$r2};
		}
		$r1a1r2a2{$pair[0]." ".$ca1." ".$pair[1]." ".$ca2} = $contacts{$_};
	}
	return %r1a1r2a2;
}

sub rr_remove_unsatisfied_top_model{
	my $file_rr  = shift;
	my $file_pdb = shift;
	my $file_new_rr  = shift;
	my $file_log = shift;
	confess "ERROR! rr file $file_rr does not exist!" if not -f $file_rr;
	confess "ERROR! pdb file $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! file_new_rr not defined!" if not defined $file_new_rr;
	confess "ERROR! log file not defined!" if not defined $file_log;
	my %contacts = min_all_atom_dist_in_top_model($file_rr, $file_pdb, "log_min_all_atom_dist_in_top_model.txt");
	my %confidence = rr2contacts_hash($file_rr, 1, 10000, "confidence");
	my %rr_thres = rr2contacts_hash($file_rr, 1);
	system_cmd("rm -f $file_new_rr");
	system_cmd("rm -f $file_log");
	print2file($file_new_rr, seq_rr($file_rr));
	my $count = 0;
	foreach (keys %contacts){
		print2file($file_new_rr, $_." 0 ".$rr_thres{$_}." ".$confidence{$_}) if $contacts{$_} <= $rr_thres{$_};
		if ($contacts{$_} > $rr_thres{$_}){
			print2file($file_log, "$_ ".$contacts{$_}." ignored");
			$count ++;
		}
		else{
			print2file($file_log, "$_ ".$contacts{$_}." accepted");
		}
	}
	print "\n\nStage 2 contact filtering.. removed $count contacts!\n";
}

sub min_all_atom_dist_in_top_model{
	# find the minimum of all atom to all atom distance between two residues, in model1.pdb
	my $file_rr = shift;
	my $file_pdb = shift;
	my $file_log  = shift;
	confess "ERROR! rr file $file_rr does not exist!" if not -f $file_rr;
	confess "ERROR! pdb file $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! log file not defined!" if not defined $file_log;
	my %contacts = rr2contacts_hash($file_rr, 1);
	my %contacts_atom_list = ();
	confess ":(" if not %contacts;
	my %pdb_atom_list = pdb_each_residue_atoms($file_pdb);
	foreach (keys %contacts){
		my @C = split /\s+/, $_;
		my %left_list = ();
		my %right_list = ();
		my @A1 = split /\s+/, $pdb_atom_list{$C[0]};
		my @A2 = split /\s+/, $pdb_atom_list{$C[1]};
		foreach (@A1){
			$left_list{$C[0]." ".$_} = 1;
		}
		foreach (@A2){
			$right_list{$C[1]." ".$_} = 1;
		}
		$contacts_atom_list{$_." left"} = \%left_list;
		$contacts_atom_list{$_." right"} = \%right_list;
	}
	my %pdb_xyz = xyz_pdb($file_pdb, "all");
	print2file($file_log, "Minimum all atom distance in $file_pdb:");
	print2file($file_log, (sprintf "\n%-10s %-10s %-20s", "contact", "min_dist", "atom-pair"));
	foreach (keys %contacts){
		my %left_list = %{$contacts_atom_list{$_." left"}};
		my %right_list = %{$contacts_atom_list{$_." right"}};
		my $min_dist = 1000;
		my $atom_pair = "";
		foreach my $le (keys %left_list){
			foreach my $ri (keys %right_list){
				my @L = split /\s+/, $le;
				my @R = split /\s+/, $ri;
				confess "$file_pdb does not have ".$L[0]." ".uc($L[1])."\n" if not defined $pdb_xyz{$L[0]." ".uc($L[1])};
				confess "$file_pdb does not have ".$R[0]." ".uc($R[1])."\n" if not defined $pdb_xyz{$R[0]." ".uc($R[1])};
				my @c1 = split(/\s+/,$pdb_xyz{uc($L[0])." ".uc($L[1])});
				my $x1 = $c1[0]; my $y1 = $c1[1]; my $z1 = $c1[2];
				my @c2 = split(/\s+/,$pdb_xyz{uc($R[0])." ".uc($R[1])});
				my $x2 = $c2[0]; my $y2 = $c2[1]; my $z2 = $c2[2];
				my $dist = sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
				if($min_dist > $dist){
					$min_dist = sprintf "%.2f", $dist;
					$atom_pair = $le." - ".$ri;
				}
			}
		}
		$contacts{$_} = $min_dist;
		print2file($file_log, (sprintf "\n%-10s %-10s %-20s", $_, $min_dist, $atom_pair));
	}
	return %contacts;
}

sub strand_count{
	my $file_ss = shift;
	confess "ERROR! SS file $file_ss does not exist!" if not -f $file_ss;
	my $ss = seq_fasta($file_ss);
	my $count = 0;
	for(my $i = 0; $i <= length($ss); $i++){
		next if (substr $ss, $i, 1) ne "E";
		next if not defined (substr $ss, $i-1, 1);
		next if not defined (substr $ss, $i, 1);
		next if not defined (substr $ss, $i+1, 1);
		$count++ if ((substr $ss, $i-1, 1) eq "E" and (substr $ss, $i, 1) eq "E" and (substr $ss, $i+1, 1) ne "E");
	}
	return $count;
}

sub dssp_result{
	my $file_pdb = shift;
	my $selection = shift;
	confess "ERROR! $file_pdb does not exist!" if not -f $file_pdb;
	confess "ERROR! selection not defined!" if not defined $selection;
	my %RESIDUE = ();
	my %SS  = ();
	my %PHI = ();
	my %PSI = ();
	my $command = "$program_dssp $file_pdb | grep -C 0 -A 1000 \'  #  RESIDUE\' | tail -n +2";
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
		$res  =~ s/[^A-Z]//g;
		$sstr =~ s/[^A-Z]//g;
		next if length($rnum) < 1;
		confess ":( residue not defined for $rnum" if length($res) < 1;
		confess ":( phi not defined for $rnum" if length($phia) < 1;
		confess ":( psi not defined for $rnum" if length($psia) < 1;
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
	# special temporary change
	#return if not scalar keys %PHI;
	confess "$file_pdb seems empty!" if not scalar keys %PHI;
	return %SS if ($selection eq "ss");
	return %PHI if ($selection eq "phi");
	return %PSI if ($selection eq "psi");
	confess "ERROR! Invalid selection string!";
}

sub pdb2ss{
	my $file_pdb = shift;
	confess "ERROR! Fasta file $file_pdb does not exist!" if not -f $file_pdb;
	my %ss = dssp_result($file_pdb, "ss");
	my $ssrow = "";
	foreach(sort {$a <=> $b} keys %ss){
		$ssrow .= $ss{$_};
	}
	confess "looks like DSSP failed!" if length($ssrow) < 2;
	return $ssrow;
}

sub seq_fasta{
	my $file_fasta = shift;
	confess "ERROR! Fasta file $file_fasta does not exist!" if not -f $file_fasta;
	my $seq = "";
	open FASTA, $file_fasta or confess $!;
	while (<FASTA>){
		next if (substr($_,0,1) eq ">"); 
		chomp $_;
		$_ =~ tr/\r//d; # chomp does not remove \r
		$seq .= $_;
	}
	close FASTA;
	return $seq;
}

sub write_cns_seq{
	my $file_fasta = shift;
	my $file_cns_seq = shift;
	confess "ERROR! Fasta file $file_fasta does not exist!" if not -f $file_fasta;
	confess "ERROR! Output not defined!" if not defined $file_cns_seq;
	my @seq = split //, seq_fasta($file_fasta);
	my $three_letter_seq;
	foreach (@seq) {
		$three_letter_seq .= $AA1TO3{$_}." ";
	}
	system_cmd("rm -f $file_cns_seq");
	while($three_letter_seq){
		if(length ($three_letter_seq) <= 64 ){
			print2file($file_cns_seq, $three_letter_seq);
			$three_letter_seq = ();
		}
		else{
			print2file($file_cns_seq, substr($three_letter_seq, 0, 64));
			$three_letter_seq = substr($three_letter_seq, 64);
		}
	}
}

sub flatten_fasta{
	my $file_fasta = shift;
	confess "ERROR! Fasta file $file_fasta does not exist!" if not -f $file_fasta;
	my $seq = seq_fasta($file_fasta);
	open FASTA, $file_fasta or confess $!;
	my $header = <FASTA>;
	chomp $header;
	close FASTA;
	system_cmd("rm -f $file_fasta");
	print2file($file_fasta, $header);
	print2file($file_fasta, $seq);
}

sub fasta2residues_hash{
	my $file_fasta = shift;
	confess "ERROR! Fasta file $file_fasta does not exist!" if not -f $file_fasta;
	my $seq = seq_fasta($file_fasta);
	my %res = ();
	foreach (my $i = 0; $i < length($seq); $i++){
		$res{$i+1} = substr $seq, $i, 1;
	}
	return %res;
}

sub assess_dgsa{
	my $stage = shift;
	my $seq = seq_fasta($file_fasta);
	my @pdb_list = load_pdb("$dir_out/$stage");
	if($#pdb_list < 2){
		confess "ERROR!! Something went wrong while running CNS in $stage! Try a different atom selection scheme!";
	}
	if($#pdb_list < ($mcount - 1)){
		warn "\n\nWARNING! There were some issues! $mcount models were not genenerated in $stage! Try a different atom selection scheme!\n";
	}
	my %tbl_list = ();
	$tbl_list{"contact.tbl"}  = count_lines("contact.tbl")  if count_lines("contact.tbl");
	$tbl_list{"ssnoe.tbl"}    = count_lines("ssnoe.tbl")    if count_lines("ssnoe.tbl");
	$tbl_list{"hbond.tbl"}    = count_lines("hbond.tbl")    if count_lines("hbond.tbl");
	$tbl_list{"dihedral.tbl"} = count_lines("dihedral.tbl") if count_lines("dihedral.tbl");
	foreach my $tbl (sort keys %tbl_list){
		# also verify if CNS actually accepted all the distance restraints provided
		next if $tbl =~ m/dihedral.tbl/;
		my $search_string = "N1";
		$search_string = "N2" if $tbl eq "ssnoe.tbl";
		$search_string = "HBND" if $tbl eq "hbond.tbl";
		if(not -f "dgsa.log"){
			system_cmd("touch assess.failed");
			confess "ERROR! Something went wrong! Log file not found (dgsa.log)!";
		}
		my $result = `grep NOEPRI dgsa.log | grep $search_string | head -n 1`;
		my @C = split /\s+/, $result;
		my $count = $C[($#C)-1];
		if ($count != count_lines($tbl)){
			system_cmd("touch assess.failed");
			confess ":( CNS did not accept all restraints of $tbl! Something wrong somewhere! Only $count accepted";
		}
	}
	# remove "trial" structure of corresponding "accepted" structure because they are same
	for(my $i = 0; $i <= 1000; $i++){
		next if not -f "${id}a_$i.pdb";
		print "\ndeleting ${id}_$i.pdb because ${id}a_$i.pdb exists!";
		system_cmd("rm ${id}_$i.pdb");
	}
	my %energy_noe = ();
	foreach my $pdb (@pdb_list) {
		next if $pdb =~ m/sub_embed/;
		next if $pdb =~ m/extended/;
		$energy_noe{$pdb} = get_cns_energy($pdb, "noe");
	}
	my $top_pdb = (sort {$energy_noe{$a} <=> $energy_noe{$b}} keys %energy_noe)[0];
	print "\n";
	print "\n";
	print "           ENERGY            CLASH     SS              NOE SATISFIED(±0.2A)            SUM OF DEVIATIONS >= 0.2     PDB\n";
	print "--------------------------  -------  -------  ---------------------------------------  -------------------------  --------\n";
	print "TOTAL  VDW    BOND   NOE    2.5 3.5  H   E    CONTACTS  SS-NOE    HBONDS    DIHEDRAL   CONTACTS SS-NOE   HBONDS\n"; 
	foreach my $pdb (sort {$energy_noe{$a} <=> $energy_noe{$b}} keys %energy_noe){
		my $e1 = "-";
		my $e2 = "-";
		my $e3 = "-";
		my $e4 = "-";
		my $c1 = "-";
		my $c2 = "-";
		my $h  = "-";
		my $e  = "-";
		my $n1 = "-";
		my $n2 = "-";
		my $n3 = "-";
		my $n4 = "-";
		my $s1 = "-";
		my $s2 = "-";
		my $s3 = "-";
		$e1 = get_cns_energy($pdb, "overall");
		$e2 = get_cns_energy($pdb, "vdw");
		$e3 = get_cns_energy($pdb, "bon");
		$e4 = $energy_noe{$pdb};
		$c1 = clash_count($pdb, 2.5);
		$c2 = clash_count($pdb, 3.5);
		$h  = count_ss_match($pdb, $file_ss, "H") if $file_ss;
		$e  = count_ss_match($pdb, $file_ss, "E") if $file_ss;
		$n1 = count_satisfied_tbl_rows($pdb, "contact.tbl", "noe") if defined $tbl_list{"contact.tbl"};
		$n2 = count_satisfied_tbl_rows($pdb, "ssnoe.tbl", "noe")   if defined $tbl_list{"ssnoe.tbl"};
		$n3 = count_satisfied_tbl_rows($pdb, "hbond.tbl", "noe")   if defined $tbl_list{"hbond.tbl"};
		$n4 = count_satisfied_tbl_rows($pdb, "dihedral.tbl", "dihedral")   if defined $tbl_list{"dihedral.tbl"};
		$s1 = sum_noe_dev($pdb, "contact.tbl") if defined $tbl_list{"contact.tbl"};
		$s2 = sum_noe_dev($pdb, "ssnoe.tbl")   if defined $tbl_list{"ssnoe.tbl"};
		$s3 = sum_noe_dev($pdb, "hbond.tbl")   if defined $tbl_list{"hbond.tbl"};
		printf "%-6s %-6s %-6s %-6s %-3s %-3s  %-3s %-3s  ",  $e1, $e2, $e3, $e4, $c1, $c2, $h, $e;
		printf "%-9s %-9s %-9s %-9s  %-8s %-8s %-8s %-25s\n", $n1, $n2, $n3, $n4, $s1, $s2, $s3, basename($pdb); 
	}
	print "\n";
	foreach my $pdb (sort {$energy_noe{$a} <=> $energy_noe{$b}} keys %energy_noe){
		my $ss = pdb2ss($pdb);
		$ss =~ s/C/-/g;
		printf $ss." [".basename($pdb)."]\n";
	}
	foreach my $tbl (sort keys %tbl_list){
		next if $tbl =~ m/dihedral/;
		print "\n";
		foreach my $pdb (sort {$energy_noe{$a} <=> $energy_noe{$b}} keys %energy_noe){
			print "".noe_tbl_violation_coverage($pdb, $tbl)." [ violation of ".basename($tbl)." in ".basename($pdb)." ]\n";
		}
	}
	print "\n";
	my $i = 1;
	foreach( sort {$energy_noe{$a} <=> $energy_noe{$b}} keys %energy_noe){
		print "model$i.pdb <= $_\n";
		system_cmd("mv $_ ${id}_model$i.pdb");
		$i++;
		last if $i > 5;
		last if $i > $mcount;
	}
	system_cmd("rm dgsa.log");
}

sub chk_errors_seq{
	my $seq = shift;
	for(my $i = 1; $i <= length($seq); $i++){
		my $char = substr $seq, $i-1, 1;
		if (not defined $AA1TO3{$char}){
			confess "Undefined amino acid $char in $seq";
		}
	}
}

sub chk_errors_rr{
	my $input_rr = shift;
	open RR, $input_rr or error_exit($!);
	my $no_rows = 1;
	while(<RR>){
		my $line = $_;
		my @RR = split /\n/, $line;
		next if length($line) < 2;
		next if $line =~ /^[A-Z]/;
		my @C = split /\s+/, $line;
		confess "column 1 not defined in $line in rr input" if not defined $C[0];
		confess "column 2 not defined in $line in rr input" if not defined $C[1];
		confess "column 3 not defined in $line in rr input" if not defined $C[2];
		confess "column 4 not defined in $line in rr input" if not defined $C[3];
		confess "column 5 not defined in $line in rr input" if not defined $C[4];
		$no_rows = 0;
	}
	close RR;
	confess "Input RR file has 0 rows! Something went wrong!" if $no_rows;
}

sub write_cns_customized_modules{
	# The scalecoolsetup module edited for scaling implementation
	print2file("scalecoolsetupedited", "! module file: scalecoolsetup");
	print2file("scalecoolsetupedited", "!");
	print2file("scalecoolsetupedited", "! CNS MODULE");
	print2file("scalecoolsetupedited", "! **********");
	print2file("scalecoolsetupedited", "!");
	print2file("scalecoolsetupedited", "! Authors: Gregory L. Warren and Axel T. Brunger");
	print2file("scalecoolsetupedited", "!");
	print2file("scalecoolsetupedited", "! copyright Yale University");
	print2file("scalecoolsetupedited", "!");
	print2file("scalecoolsetupedited", "! version 2/25/98");
	print2file("scalecoolsetupedited", "!");
	print2file("scalecoolsetupedited", "! Function:");
	print2file("scalecoolsetupedited", "!    This module sets up the NMR restraint data weighting ");
	print2file("scalecoolsetupedited", "!    scheme for either the first or second slow-cooling ");
	print2file("scalecoolsetupedited", "!    dynamics sections");
	print2file("scalecoolsetupedited", "!");
	print2file("scalecoolsetupedited", "! Requirements:");
	print2file("scalecoolsetupedited", "!    See the top level script file for the ");
	print2file("scalecoolsetupedited", "!    defined parameter(s). ");
	print2file("scalecoolsetupedited", "! ");
	print2file("scalecoolsetupedited", "!");
	print2file("scalecoolsetupedited", "");
	print2file("scalecoolsetupedited", "module { scalecoolsetup }");
	print2file("scalecoolsetupedited", "( ");
	print2file("scalecoolsetupedited", "   &kdani=kdani;             {INPUT: dani force and class numbers}");
	print2file("scalecoolsetupedited", "   &ksani=ksani;             {INPUT: sani force and class numbers}");
	print2file("scalecoolsetupedited", "   &nmr=nmr;                 {INPUT: nmr restraints parameters}");
	print2file("scalecoolsetupedited", "   &input.noe.scale=input.noe.scale;    {INPUT: noe weight}");
	print2file("scalecoolsetupedited", "   &input.cdih.scale=input.cdih.scale;  {INPUT: cdih weight}");
	print2file("scalecoolsetupedited", "   &input.ncycle=input.ncycle;          {INPUT: number of dynamics cycles}");
	print2file("scalecoolsetupedited", ")");
	print2file("scalecoolsetupedited", "");
	print2file("scalecoolsetupedited", "checkversion 1.3");
	print2file("scalecoolsetupedited", "");
	print2file("scalecoolsetupedited", "set message ? end");
	print2file("scalecoolsetupedited", "evaluate (\$message_old=\$result)");
	print2file("scalecoolsetupedited", "set echo ? end");
	print2file("scalecoolsetupedited", "evaluate (\$echo_old=\$result)");
	print2file("scalecoolsetupedited", "if ( \$log_level = verbose ) then");
	print2file("scalecoolsetupedited", "  set echo=on message=normal end");
	print2file("scalecoolsetupedited", "else");
	print2file("scalecoolsetupedited", "  set echo=off message=off end");
	print2file("scalecoolsetupedited", "end if");
	print2file("scalecoolsetupedited", "");
	print2file("scalecoolsetupedited", "   noe ");
	print2file("scalecoolsetupedited", "      scale * &input.noe.scale");
	print2file("scalecoolsetupedited", "   end");
	print2file("scalecoolsetupedited", "");
	print2file("scalecoolsetupedited", "! start: edited by Badri - 4/22/2014     ");
	print2file("scalecoolsetupedited", "    evaluate(\$newScale = &input.noe.scale * $sswt)");
	print2file("scalecoolsetupedited", "    noe");
	print2file("scalecoolsetupedited", "       scale N2 \$newScale");
	print2file("scalecoolsetupedited", "    end");
	print2file("scalecoolsetupedited", "! end: edited by Badri - 4/22/2014     ");
	print2file("scalecoolsetupedited", "");
	print2file("scalecoolsetupedited", "   restraints dihedral ");
	print2file("scalecoolsetupedited", "      scale = &input.cdih.scale ");
	print2file("scalecoolsetupedited", "   end");
	print2file("scalecoolsetupedited", "");
	print2file("scalecoolsetupedited", "   evaluate (\$count = 1)");
	print2file("scalecoolsetupedited", "   while (&exist%nmr.dani.file.\$count=true) loop nloop");
	print2file("scalecoolsetupedited", "      if (&nmr.dani.file.\$count # \"\" ) then");
	print2file("scalecoolsetupedited", "         if (&kdani.cart.flag=true) then");
	print2file("scalecoolsetupedited", "            evaluate (&kdani.start=1.0)");
	print2file("scalecoolsetupedited", "            evaluate (&kdani.finish=&nmr.dani.force.finl.\$count)");
	print2file("scalecoolsetupedited", "         elseif (&kdani.inter.flag=true) then");
	print2file("scalecoolsetupedited", "            evaluate (&kdani.start=&nmr.dani.force.init.\$count)");
	print2file("scalecoolsetupedited", "            evaluate (&kdani.finish=1.0)");
	print2file("scalecoolsetupedited", "         else");
	print2file("scalecoolsetupedited", "            evaluate (&kdani.start=&nmr.dani.force.init.\$count)");
	print2file("scalecoolsetupedited", "            evaluate (&kdani.finish=&nmr.dani.force.finl.\$count)");
	print2file("scalecoolsetupedited", "         end if");
	print2file("scalecoolsetupedited", "         evaluate (&kdani.\$count  = &kdani.start)");
	print2file("scalecoolsetupedited", "         evaluate (&kdani.fac.\$count = (&kdani.finish/");
	print2file("scalecoolsetupedited", "                                        &kdani.start)^(1/&input.ncycle))");
	print2file("scalecoolsetupedited", "      end if");
	print2file("scalecoolsetupedited", "      evaluate (\$count = \$count + 1)");
	print2file("scalecoolsetupedited", "   end loop nloop");
	print2file("scalecoolsetupedited", "");
	print2file("scalecoolsetupedited", "   evaluate (\$count = 1)");
	print2file("scalecoolsetupedited", "   while (&exist%nmr.sani.file.\$count=true) loop nloop");
	print2file("scalecoolsetupedited", "      if (&nmr.sani.file.\$count # \"\" ) then");
	print2file("scalecoolsetupedited", "         if (&ksani.cart.flag=true) then");
	print2file("scalecoolsetupedited", "            evaluate (&ksani.start=1.0)");
	print2file("scalecoolsetupedited", "            evaluate (&ksani.finish=&nmr.sani.force.finl.\$count)");
	print2file("scalecoolsetupedited", "         elseif (&ksani.inter.flag=true) then");
	print2file("scalecoolsetupedited", "            evaluate (&ksani.start=&nmr.sani.force.init.\$count)");
	print2file("scalecoolsetupedited", "            evaluate (&ksani.finish=1.0)");
	print2file("scalecoolsetupedited", "         else");
	print2file("scalecoolsetupedited", "            evaluate (&ksani.start=&nmr.sani.force.init.\$count)");
	print2file("scalecoolsetupedited", "            evaluate (&ksani.finish=&nmr.sani.force.finl.\$count)");
	print2file("scalecoolsetupedited", "         end if");
	print2file("scalecoolsetupedited", "         evaluate (&ksani.\$count  = &ksani.start)");
	print2file("scalecoolsetupedited", "         evaluate (&ksani.fac.\$count = (&ksani.finish/");
	print2file("scalecoolsetupedited", "                                        &ksani.start)^(1/&input.ncycle))");
	print2file("scalecoolsetupedited", "      end if");
	print2file("scalecoolsetupedited", "      evaluate (\$count = \$count + 1)");
	print2file("scalecoolsetupedited", "   end loop nloop");
	print2file("scalecoolsetupedited", "");
	print2file("scalecoolsetupedited", "set message=\$message_old echo=\$echo_old end");
	# The scalehot module edited for scaling implementation
	print2file("scalehotedited", "! module file: scalehot");
	print2file("scalehotedited", "!");
	print2file("scalehotedited", "! CNS MODULE");
	print2file("scalehotedited", "! **********");
	print2file("scalehotedited", "!");
	print2file("scalehotedited", "! Authors: Gregory L. Warren and Axel T. Brunger");
	print2file("scalehotedited", "!");
	print2file("scalehotedited", "! copyright Yale University");
	print2file("scalehotedited", "!");
	print2file("scalehotedited", "! version 2/25/98");
	print2file("scalehotedited", "!");
	print2file("scalehotedited", "! Function:");
	print2file("scalehotedited", "!    This module sets the NMR restraint data weight ");
	print2file("scalehotedited", "!    for high temperature dynamics");
	print2file("scalehotedited", "!");
	print2file("scalehotedited", "! Requirements:");
	print2file("scalehotedited", "!    See the top level script file for the ");
	print2file("scalehotedited", "!    defined parameter(s) nmr and md. ");
	print2file("scalehotedited", "! ");
	print2file("scalehotedited", "!");
	print2file("scalehotedited", "");
	print2file("scalehotedited", "module { scalehot }");
	print2file("scalehotedited", "( ");
	print2file("scalehotedited", "   &md=md;                   {INPUT: md parameters}");
	print2file("scalehotedited", "   &nmr=nmr;                 {INPUT: nmr restraints parameters}");
	print2file("scalehotedited", "   &input.noe.scale=input.noe.scale;    {INPUT: noe weight}");
	print2file("scalehotedited", "   &input.cdih.scale=input.cdih.scale;  {INPUT: cdih weight}");
	print2file("scalehotedited", ")");
	print2file("scalehotedited", "");
	print2file("scalehotedited", "checkversion 1.3");
	print2file("scalehotedited", "");
	print2file("scalehotedited", "set message ? end");
	print2file("scalehotedited", "evaluate (\$message_old=\$result)");
	print2file("scalehotedited", "set echo ? end");
	print2file("scalehotedited", "evaluate (\$echo_old=\$result)");
	print2file("scalehotedited", "if ( \$log_level = verbose ) then");
	print2file("scalehotedited", "  set echo=on message=normal end");
	print2file("scalehotedited", "else");
	print2file("scalehotedited", "  set echo=off message=off end");
	print2file("scalehotedited", "end if");
	print2file("scalehotedited", "");
	print2file("scalehotedited", "   noe ");
	print2file("scalehotedited", "      scale * &input.noe.scale");
	print2file("scalehotedited", "   end");
	print2file("scalehotedited", "! start: edited by Badri - 4/22/2014     ");
	print2file("scalehotedited", "    evaluate(\$newScale = &input.noe.scale * $sswt)");
	print2file("scalehotedited", "    noe");
	print2file("scalehotedited", "        scale N2 \$newScale");
	print2file("scalehotedited", "    end");
	print2file("scalehotedited", "! end: edited by Badri - 4/22/2014  ");
	print2file("scalehotedited", "   restraints dihedral");
	print2file("scalehotedited", "      scale = &input.cdih.scale");
	print2file("scalehotedited", "   end");
	print2file("scalehotedited", "   evaluate (\$count = 1)");
	print2file("scalehotedited", "   while (&exist%nmr.dani.file.\$count=true) loop nloop");
	print2file("scalehotedited", "      evaluate (\$clsname = \"D\"+encode(\$count))");
	print2file("scalehotedited", "      if (&nmr.dani.file.\$count # \"\" ) then");
	print2file("scalehotedited", "         dani");
	print2file("scalehotedited", "            class \$\$clsname force &nmr.dani.force.init.\$count");
	print2file("scalehotedited", "         end");
	print2file("scalehotedited", "      end if");
	print2file("scalehotedited", "      evaluate (\$count = \$count + 1)");
	print2file("scalehotedited", "   end loop nloop");
	print2file("scalehotedited", "");
	print2file("scalehotedited", "   evaluate (\$count = 1)");
	print2file("scalehotedited", "   while (&exist%nmr.sani.file.\$count=true) loop nloop");
	print2file("scalehotedited", "      evaluate (\$clsname = \"S\"+encode(\$count))");
	print2file("scalehotedited", "      if (&nmr.sani.file.\$count # \"\" ) then");
	print2file("scalehotedited", "         sani");
	print2file("scalehotedited", "            class \$\$clsname force &nmr.sani.force.init.\$count");
	print2file("scalehotedited", "         end");
	print2file("scalehotedited", "      end if");
	print2file("scalehotedited", "      evaluate (\$count = \$count + 1)");
	print2file("scalehotedited", "   end loop nloop");
	print2file("scalehotedited", "");
	print2file("scalehotedited", "set message=\$message_old echo=\$echo_old end");
}

sub write_cns_dgsa_file{
	my $dihed_wt1 = $contwt * $sswt;
	my $dihed_wt2 = $contwt * $sswt;
	print2file("dgsa.inp", "{+ file: dgsa.inp +}");
	print2file("dgsa.inp", "{+ directory: nmr_calc +}");
	print2file("dgsa.inp", "{+ description: distance geometry, full or substructure, with ");
	print2file("dgsa.inp", "                simulated annealing regularization starting from ");
	print2file("dgsa.inp", "                extended strand or pre-folded structures. +}");
	print2file("dgsa.inp", "{+ authors: Gregory Warren, Michael Nilges, John Kuszewski, ");
	print2file("dgsa.inp", "	    Marius Clore and Axel Brunger +}");
	print2file("dgsa.inp", "{+ copyright: Yale University +}");
	print2file("dgsa.inp", "{+ reference: Clore GM, Gronenborn AM, Tjandra N, Direct structure refinement ");
	print2file("dgsa.inp", "              against residual dipolar couplings in the presence of rhombicity");
	print2file("dgsa.inp", "              of unknown magnitude., J. Magn. Reson., 131, In press, (1998) +}");
	print2file("dgsa.inp", "{+ reference: Clore GM, Gronenborn AM, Bax A, A robust method for determining ");
	print2file("dgsa.inp", "              the magnitude of the fully asymmetric alignment tensor of");
	print2file("dgsa.inp", "              oriented macromolecules in the absence of structural");
	print2file("dgsa.inp", "              information., J. Magn. Reson., In press (1998) +}");
	print2file("dgsa.inp", "{+ reference: Garrett DS, Kuszewski J, Hancock TJ, Lodi PJ, Vuister GW,");
	print2file("dgsa.inp", "              Gronenborn AM, Clore GM, The impact of direct refinement against ");
	print2file("dgsa.inp", "              three-bond HN-C alpha H coupling constants on protein structure");
	print2file("dgsa.inp", "              determination by NMR., J. Magn. Reson. Ser. B, 104(1), ");
	print2file("dgsa.inp", "              99-103, (1994) May +}");
	print2file("dgsa.inp", "{+ reference: Kuszewski J, Nilges M, Brunger AT,   Sampling and efficiency ");
	print2file("dgsa.inp", "              of metric matrix distance geometry:  A novel partial metrization ");
	print2file("dgsa.inp", "              algorithm.  J. Biomol. NMR 2, 33-56, (1992). +} ");
	print2file("dgsa.inp", "{+ reference: Kuszewski J, Qin J, Gronenborn AM, Clore GM, The impact of direct");
	print2file("dgsa.inp", "              refinement against 13C alpha and 13C beta chemical shifts on ");
	print2file("dgsa.inp", "              protein structure determination by NMR., J. Magn. Reson. Ser. B,");
	print2file("dgsa.inp", "              106(1), 92-6, (1995) Jan +}");
	print2file("dgsa.inp", "{+ reference: Kuszewski J, Gronenborn AM, Clore GM, The impact of direct");
	print2file("dgsa.inp", "              refinement against proton chemical shifts on protein structure ");
	print2file("dgsa.inp", "              determination by NMR., J. Magn. Reson. Ser. B, 107(3), 293-7, ");
	print2file("dgsa.inp", "              (1995) Jun +}");
	print2file("dgsa.inp", "{+ reference: Kuszewski J, Gronenborn AM, Clore GM, A potential involving ");
	print2file("dgsa.inp", "              multiple proton chemical-shift restraints for ");
	print2file("dgsa.inp", "              nonstereospecifically assigned methyl and methylene protons.");
	print2file("dgsa.inp", "              J. Magn. Reson. Ser. B, 112(1), 79-81, (1996) Jul. +}");
	print2file("dgsa.inp", "{+ reference: Nilges M, Clore GM, Gronenborn AM, Determination of ");
	print2file("dgsa.inp", "              three-dimensional structures of proteins from interproton ");
	print2file("dgsa.inp", "              distance data by hybrid distance geometry-dynamical simulated ");
	print2file("dgsa.inp", "              annealing calculations. FEBS Lett. 229, 317-324 (1988). +}");
	print2file("dgsa.inp", "{+ reference: Nilges M, Clore GM, Gronenborn AM,  Determination of ");
	print2file("dgsa.inp", "              three-dimensional structures of proteins from interproton ");
	print2file("dgsa.inp", "              distance data by dynamical simulated annealing from a random ");
	print2file("dgsa.inp", "              array of atoms. FEBS LEtt. 239, 129-136 (1988). +}");
	print2file("dgsa.inp", "{+ reference: Nilges M, Kuszewski J, Brunger AT, In: Computational Aspects ");
	print2file("dgsa.inp", "              of the Study of Biological Macromolecules by NMR. ");
	print2file("dgsa.inp", "              (J.C. Hoch, ed.),  New York: Plenum Press, (1991). +}");
	print2file("dgsa.inp", "{+ reference: Tjandra N, Garrett DS, Gronenborn AM, Bax A, Clore GM, Defining");
	print2file("dgsa.inp", "              long range order in NMR structure determination from the ");
	print2file("dgsa.inp", "              dependence of heteronuclear relaxation times on rotational ");
	print2file("dgsa.inp", "              diffusion anisotropy. Nature Struct. Biol., 4(6), 443-9,");
	print2file("dgsa.inp", "              (1997) June +}");
	print2file("dgsa.inp", "{+ reference: Tjandra N, Omichinski JG, Gronenborn AM, Clore GM, Bax A, Use of");
	print2file("dgsa.inp", "              dipolar 1H-15N and 1H-13C couplings in the structure");
	print2file("dgsa.inp", "              determination of magnetically oriented macromolecules in");
	print2file("dgsa.inp", "              solution. Nature Struct. Biol., 4(9), 732-8, (1997) Sept +} ");
	print2file("dgsa.inp", "              ");
	print2file("dgsa.inp", "{- Guidelines for using this file:");
	print2file("dgsa.inp", "   - all strings must be quoted by double-quotes");
	print2file("dgsa.inp", "   - logical variables (true/false) are not quoted");
	print2file("dgsa.inp", "   - do not remove any evaluate statements from the file -}");
	print2file("dgsa.inp", "{- begin block parameter definition -} define(");
	print2file("dgsa.inp", "{======================= molecular structure =========================}");
	print2file("dgsa.inp", "{* parameter file(s) *}");
	# CNS_TOPPAR:protein-allhdg5-4.param
	print2file("dgsa.inp", "{===>} par.1=\"CNS_TOPPAR:protein.param\";");
	print2file("dgsa.inp", "{===>} par.2=\"\";");
	print2file("dgsa.inp", "{===>} par.3=\"\";");
	print2file("dgsa.inp", "{===>} par.4=\"\";");
	print2file("dgsa.inp", "{===>} par.5=\"\";");
	print2file("dgsa.inp", "{* structure file(s) *}");
	print2file("dgsa.inp", "{===>} struct.1=\"extended.mtf\";");
	print2file("dgsa.inp", "{===>} struct.2=\"\";");
	print2file("dgsa.inp", "{===>} struct.3=\"\";");
	print2file("dgsa.inp", "{===>} struct.4=\"\";");
	print2file("dgsa.inp", "{===>} struct.5=\"\";");
	print2file("dgsa.inp", "{* input coordinate file(s) *}");
	print2file("dgsa.inp", "{===>} pdb.in.file.1=\"extended.pdb\";");
	print2file("dgsa.inp", "{===>} pdb.in.file.2=\"\";");
	print2file("dgsa.inp", "{===>} pdb.in.file.3=\"\";");
	print2file("dgsa.inp", "{========================== atom selection ===========================}");
	print2file("dgsa.inp", "{* input \"backbone\" selection criteria for average structure generation *}");
	print2file("dgsa.inp", "{* for protein      (name n or name ca or name c)");
	print2file("dgsa.inp", "   for nucleic acid (name O5' or name C5' or name C4' or name C3' ");
	print2file("dgsa.inp", "                     or name O3' or name P) *}");
	print2file("dgsa.inp", "{===>} pdb.atom.select=(name n or name ca or name c);");
	print2file("dgsa.inp", "{======================= refinement parameters ========================}");
	print2file("dgsa.inp", "{* distance geometry *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.dg.flag=true;");
	print2file("dgsa.inp", "{* distance geometry/simualted annealing regularization (DGSA) *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.dgsa.flag=true;");
	print2file("dgsa.inp", "{* if only regularizing coordinate files (no DG) then ");
	print2file("dgsa.inp", "   enter the number of coordinate files to be regularized (DGSA) *}");
	print2file("dgsa.inp", "{===>} pdb.dg.count=$mcount;");
	print2file("dgsa.inp", "{* seed for random number generator *}");
	print2file("dgsa.inp", "{* change to get different initial velocities *}");
	print2file("dgsa.inp", "{===>} md.seed=82364;");
	print2file("dgsa.inp", "{* select whether the number of structures will be either trial or 	");
	print2file("dgsa.inp", "   accepted structures and whether to print only the trial, accepted, 	");
	print2file("dgsa.inp", "   both sets of structures. The printing format is as follows:");
	print2file("dgsa.inp", "   trial = pdb.out.name + _#.pdb , accepted = pdb.out.name + a_#.pdb *} ");
	print2file("dgsa.inp", "{* are the number of structures to be trials or accepted? *}");
	print2file("dgsa.inp", "{+ choice: \"trial\" \"accept\" +}");
	print2file("dgsa.inp", "{===>} flg.trial.struc=\"$mode\";");
	print2file("dgsa.inp", "{* number of trial or accepted structures *}");
	print2file("dgsa.inp", "{===>} pdb.end.count=$mcount;");
	print2file("dgsa.inp", "{* print accepted structures *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.print.accept=true;");
	print2file("dgsa.inp", "{* print trial structures *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.print.trial=true;");
	print2file("dgsa.inp", "{* calculate an average structure for either the trial or 	");
	print2file("dgsa.inp", "   accepted structure.  If calculate accepted average is false then ");
	print2file("dgsa.inp", "   an average for the trial structures will be calculated. *}");
	print2file("dgsa.inp", "{* calculate an average structure? *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.calc.ave.struct=false;");
	print2file("dgsa.inp", "{* calculate an average structure for the accepted structures? *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.calc.ave.accpt=false;");
	print2file("dgsa.inp", "{* minimize average coordinates? *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} flg.min.ave.coor=false;");
	print2file("dgsa.inp", "{============ parameters for the distance geometry stage ==============}");
	print2file("dgsa.inp", "{* shortest path algorithm *}");
	print2file("dgsa.inp", "{+ choice: \"auto\" \"full\" \"sparse\" +}");
	print2file("dgsa.inp", "{===>} md.dg.algo=\"auto\";");
	print2file("dgsa.inp", "{* distance geometry on substructure or complete structure? *}");
	print2file("dgsa.inp", "{* proteins: \"sub\" or \"complete\"; dna/rna: \"complete\" *}");
	print2file("dgsa.inp", "{+ choice: \"sub\" \"complete\" +}");
	print2file("dgsa.inp", "{===>} md.dg.type=\"sub\";");
	print2file("dgsa.inp", "{* input atom selection for substructure  *}");
	# Atom selection scheme for distance geometry
	# (1) as-is / ca + ha + n + hn + c + cb + cg (2) as-is + o (3) as-is + o + h 
	# (2) as-is + o (3) as-is + o + h (4) c + cα + n + o (backbone atoms)
	# (5) bkbone + cβ (6) bkbone + cβ + h (7)  bkbone + cβ + h + cg                              
	# Notes: 
	# (a) CNS recommended is 6, but with option 6, error "%MATROT error encountered: Error in internal consistency check" was observed for some pdbs
	# (b) with option 7, chirality issues were faced on target 1QJP
	# (c) with option 2 (best), maximum reconstruction accuracy and sheet_detection results for 150 fragfold proteins.
	if ($atomselect == 1){
		# as-is
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn");
		print2file("dgsa.inp", "		        or name c or name cb* or name cg*);");
	}
	if ($atomselect == 2){
		# add oxygen
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn");
		print2file("dgsa.inp", "		        or name c or name cb* or name cg* or name o);");
	}
	if ($atomselect == 3){
		# add oxygen and hydrogen
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn or name h");
		print2file("dgsa.inp", "		        or name c or name cb* or name cg* or name o);");
	}
	if ($atomselect == 4){
		# backbone atoms only
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o);");
	}
	if ($atomselect == 5){
		# backbone atoms with cb
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o");
		print2file("dgsa.inp", "		        or name cb);");
	}
	if ($atomselect == 6){
		# backbone atoms with cb and hydrogen
		# atom selection according to instructions in the NIH-XPLORE manual
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o");
		print2file("dgsa.inp", "		        or name cb or name h);");
	}
	if ($atomselect == 7){
		# backbone atoms with cb, hydrogen and cg
		print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name c or name n or name o");
		print2file("dgsa.inp", "		        or name cb or name h or name cg*);");
	}
	print2file("dgsa.inp", "{===>} md.dg.select=(name ca or name ha or name n or name hn");
	print2file("dgsa.inp", "		        or name c or name cb* or name cg*);");
	print2file("dgsa.inp", "{* when using \"complete\" input the rigid group atom selection *}");
	print2file("dgsa.inp", "{===>} md.dg.group.slct=(known);");
	print2file("dgsa.inp", "{* group interatomic error value in angstroms *}");
	print2file("dgsa.inp", "{===>} md.dg.group.err=0.5;");
	print2file("dgsa.inp", "{* use metrization for complete distance geometry? *}");
	print2file("dgsa.inp", "{+ choice: true false +}");
	print2file("dgsa.inp", "{===>} md.dg.metr.flag=false;");
	print2file("dgsa.inp", "{* ordered or random metrization *}");
	print2file("dgsa.inp", "{+ choice: \"ordered\" \"random\" +}");
	print2file("dgsa.inp", "{===>} md.dg.ord=\"random\";");
	print2file("dgsa.inp", "{* input metrization atom selection *}");
	print2file("dgsa.inp", "{===>} md.dg.metr.atom=(all);");
	print2file("dgsa.inp", "{* input number of atoms from selection used during retightening *}");
	print2file("dgsa.inp", "{===>} md.dg.metr.num=4;");
	print2file("dgsa.inp", "{* reference for building the reference data base *}");
	print2file("dgsa.inp", "{+ choice: \"parameter\" \"coordinate\" +}");
	print2file("dgsa.inp", "{===>} md.dg.ref=\"parameter\";");
	print2file("dgsa.inp", "{* scale factor for distance geometry restraint term *}");
	print2file("dgsa.inp", "{===>} md.dg.scale=100.;");
	print2file("dgsa.inp", "{* exponent for distance geometry restraint term *}");
	print2file("dgsa.inp", "{===>} md.dg.expo=2;");
	print2file("dgsa.inp", "{* bond length (in angstroms) error value *}");
	print2file("dgsa.inp", "{===>} md.dg.bacc=0.01;");
	print2file("dgsa.inp", "{* angle (in degrees) error value *}");
	print2file("dgsa.inp", "{===>} md.dg.tacc=2.;");
	print2file("dgsa.inp", "{* improper (in degrees) error value *}");
	print2file("dgsa.inp", "{===>} md.dg.iacc=2.;");
	print2file("dgsa.inp", "{* dihedral (in degrees) error value *}");
	print2file("dgsa.inp", "{===>} md.dg.pacc=2.;");
	print2file("dgsa.inp", "{* number of steps of minimization *}");
	print2file("dgsa.inp", "{===>} md.dg.step=200;");
	print2file("dgsa.inp", "{=== parameters for the distance geometry/simulated annealing stage ===}");
	print2file("dgsa.inp", "{* starting temperature *}");
	print2file("dgsa.inp", "{===>} md.hot.temp=2000;");
	print2file("dgsa.inp", "{* number of steps for high temperature dyanmics *}");
	print2file("dgsa.inp", "{===>} md.hot.step=1000;");
	print2file("dgsa.inp", "{* number of steps for slow-cool annealing *}");
	print2file("dgsa.inp", "{===>} md.cool.step=1000;");
	print2file("dgsa.inp", "{* hot molecular dynamics timestep *}");
	print2file("dgsa.inp", "{===>} md.hot.ss=0.003;");
	print2file("dgsa.inp", "{* slow-cool molecular dynamics timestep *}");
	print2file("dgsa.inp", "{===>} md.cool.ss=0.005;");
	print2file("dgsa.inp", "{* initial scale factor for van der Waals (repel) energy term *}");
	print2file("dgsa.inp", "{===>} md.cool.vdw.init=0.003;");
	print2file("dgsa.inp", "{* final scale factor for van der Waals (repel) energy term *}");
	print2file("dgsa.inp", "{===>} md.cool.vdw.finl=4.0;");
	print2file("dgsa.inp", "{* initial van der Waals repel radius *}");
	print2file("dgsa.inp", "{===>} md.cool.init.rad=$rep1;"); # 0.9 originally
	print2file("dgsa.inp", "{* final van der Waals repel radius *}");
	print2file("dgsa.inp", "{===>} md.cool.fina.rad=$rep2;"); # 0.8 originally
	print2file("dgsa.inp", "{* scale factor for NOE energy term *}");
	print2file("dgsa.inp", "{===>} md.cool.noe=$contwt;"); # 50 originally
	print2file("dgsa.inp", "{* high temperature scale factor for dihedral angle energy term *}");
	print2file("dgsa.inp", "{===>} md.hot.cdih=5;");
	print2file("dgsa.inp", "{* slow-cooling scale factor for dihedral angle energy term *}");
	print2file("dgsa.inp", "{===>} md.cool.cdih=$dihed_wt1;"); # 200 originally
	print2file("dgsa.inp", "{* slow-cool annealing temperature step *}");
	print2file("dgsa.inp", "{===>} md.cool.tmpstp=25.;");
	print2file("dgsa.inp", "{=============== parameters for final minimization stage ==============}");
	print2file("dgsa.inp", "{* scale factor for NOE energy term *}");
	print2file("dgsa.inp", "{===>} md.pow.noe=$contwt;"); # 50 originally
	print2file("dgsa.inp", "{* scale factor for dihedral angle energy term *}");
	print2file("dgsa.inp", "{===>} md.pow.cdih=$dihed_wt2;"); # 400 originally
	print2file("dgsa.inp", "{* number of minimization steps *}");
	print2file("dgsa.inp", "{===>} md.pow.step=$mini;"); # 200 originally
	print2file("dgsa.inp", "{* number of cycles of minimization *}");
	print2file("dgsa.inp", "{===>} md.pow.cycl=10;");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "{============================= noe data ===============================}");
	print2file("dgsa.inp", "{- Important - if you do not have a particular data set then");
	print2file("dgsa.inp", "   set the file name to null (\"\") -}");
	print2file("dgsa.inp", "{* NOE distance restraints files. *}");
	print2file("dgsa.inp", "{* restraint set 1 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.1=\"contact.tbl\";");
	print2file("dgsa.inp", "{* restraint set 2 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.2=\"ssnoe.tbl\";");
	print2file("dgsa.inp", "{* restraint set 3 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.3=\"\";");
	print2file("dgsa.inp", "{* restraint set 4 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.4=\"\";");
	print2file("dgsa.inp", "{* restraint set 5 file *}");
	print2file("dgsa.inp", "{===>} nmr.noe.file.5=\"\";");
	print2file("dgsa.inp", "{* NOE averaging modes *}");
	print2file("dgsa.inp", "{* restraint set 1 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.1=\"cent\";");
	print2file("dgsa.inp", "{* restraint set 2 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.2=\"sum\";");
	print2file("dgsa.inp", "{* restraint set 3 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.3=\"R-6\";");
	print2file("dgsa.inp", "{* restraint set 4 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.4=\"\";");
	print2file("dgsa.inp", "{* restraint set 5 *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.5=\"\";");
	print2file("dgsa.inp", "{======================== hydrogen bond data ==========================}");
	print2file("dgsa.inp", "{* hydrogen-bond distance restraints file. *}");
	print2file("dgsa.inp", "{===>} nmr.noe.hbnd.file=\"hbond.tbl\";");
	print2file("dgsa.inp", "{* enter hydrogen-bond distance averaging mode *}");
	print2file("dgsa.inp", "{+ choice: \"sum\" \"cent\" \"R-6\" \"R-3\" \"symm\" +}");
	print2file("dgsa.inp", "{===>} nmr.noe.ave.mode.hbnd=\"cent\";");
	print2file("dgsa.inp", "{======================= 3-bond J-coupling data =======================}");
	print2file("dgsa.inp", "{* the default setup is for the phi dihedral *}");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* 3-bond J-coupling non-glycine restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.file.1=\"\";");
	print2file("dgsa.inp", "{* 3-bond J-coupling non-glycine potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* 3-bond J-coupling non-glycine force value *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.force.1.1=1;");
	print2file("dgsa.inp", "{* 3-bond j-coupling multiple class force second value *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.force.2.1=0;");
	print2file("dgsa.inp", "{* 3-bond j-coupling Karplus coefficients *}");
	print2file("dgsa.inp", "{* the default values are for phi *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.1.1=6.98;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.2.1=-1.38;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.3.1=1.72;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.4.1=-60.0;");
	print2file("dgsa.inp", "{* Class 2 *}");
	print2file("dgsa.inp", "{* 3-bond j-coupling glycine restraints files *}");
	print2file("dgsa.inp", "{* The potential for the glycine class must be multiple *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.file.2=\"\";");
	print2file("dgsa.inp", "{* 3-bond J-coupling non-glycine potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.pot.2=\"multiple\";");
	print2file("dgsa.inp", "{* 3-bond J-coupling first force value *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.force.1.2=1;");
	print2file("dgsa.inp", "{* 3-bond j-coupling glycine or multiple force second value *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.force.2.2=0;");
	print2file("dgsa.inp", "{* 3-bond j-coupling Karplus coefficients *}");
	print2file("dgsa.inp", "{* the default values are for glycine phi *}");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.1.2=6.98;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.2.2=-1.38;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.3.2=1.72;");
	print2file("dgsa.inp", "{===>} nmr.jcoup.coef.4.2=0.0;");
	print2file("dgsa.inp", "{================ 1-bond heteronuclear J-coupling data ================}");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* 1-bond heteronuclear j-coupling file *}");
	print2file("dgsa.inp", "{===>} nmr.oneb.file.1=\"\";");
	print2file("dgsa.inp", "{* 1-bond heteronuclear j-coupling potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" +}");
	print2file("dgsa.inp", "{===>} nmr.oneb.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* 1-bond heteronuclear j-coupling force value *}");
	print2file("dgsa.inp", "{===>} nmr.oneb.force.1=1.0;");
	print2file("dgsa.inp", "{=============== alpha/beta carbon chemical shift data ================}");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* carbon, alpha and beta, chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.carb.file.1=\"\";");
	print2file("dgsa.inp", "{* carbon, alpha and beta, chemical shift restraint potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" +}");
	print2file("dgsa.inp", "{===>} nmr.carb.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* carbon, alpha and beta, chemical shift restraint force value *}");
	print2file("dgsa.inp", "{===>} nmr.carb.force.1=0.5;");
	print2file("dgsa.inp", "{===================== proton chemical shift data =====================}");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* class 1 proton chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.prot.file.1=\"\";");
	print2file("dgsa.inp", "{* class 1 proton chemical shift potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.prot.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* class 1 proton chemical shift force value *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.1.1=7.5;");
	print2file("dgsa.inp", "{* 2nd class 1 proton chemical shift force value for multi *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.2.1=0;");
	print2file("dgsa.inp", "{* class 1 proton chemical shift violation cutoff threshold *}");
	print2file("dgsa.inp", "{===>} nmr.prot.thresh.1=0.3;");
	print2file("dgsa.inp", "{* Class 2 *}");
	print2file("dgsa.inp", "{* class 2 proton chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.prot.file.2=\"\";");
	print2file("dgsa.inp", "{* class 2 proton chemical shift potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.prot.pot.2=\"harmonic\";");
	print2file("dgsa.inp", "{* class 2 proton chemical shift force value *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.1.2=7.5;");
	print2file("dgsa.inp", "{* 2nd class 2 proton chemical shift force value for multi *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.2.2=0;");
	print2file("dgsa.inp", "{* class 2 proton chemical shift violation cutoff threshold *}");
	print2file("dgsa.inp", "{===>} nmr.prot.thresh.2=0.3;");
	print2file("dgsa.inp", "{* Class 3 *}");
	print2file("dgsa.inp", "{* class 3 proton chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.prot.file.3=\"\";");
	print2file("dgsa.inp", "{* class 3 proton chemical shift potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.prot.pot.3=\"harmonic\";");
	print2file("dgsa.inp", "{* class 3 proton chemical shift force value *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.1.3=7.5;");
	print2file("dgsa.inp", "{* 2nd class 3 proton chemical shift force value for multi *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.2.3=0;");
	print2file("dgsa.inp", "{* class 3 proton chemical shift violation cutoff threshold *}");
	print2file("dgsa.inp", "{===>} nmr.prot.thresh.3=0.3;");
	print2file("dgsa.inp", "{* Class 4 *}");
	print2file("dgsa.inp", "{* class 4 proton chemical shift restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.prot.file.4=\"\";");
	print2file("dgsa.inp", "{* class 4 proton chemical shift potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" \"multiple\" +}");
	print2file("dgsa.inp", "{===>} nmr.prot.pot.4=\"multiple\";");
	print2file("dgsa.inp", "{* class 4 proton chemical shift force value *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.1.4=7.5;");
	print2file("dgsa.inp", "{* 2nd class 4 proton chemical shift force value for multi *}");
	print2file("dgsa.inp", "{===>} nmr.prot.force.2.4=0;");
	print2file("dgsa.inp", "{* class 4 proton chemical shift violation cutoff threshold *}");
	print2file("dgsa.inp", "{===>} nmr.prot.thresh.4=0.3;");
	print2file("dgsa.inp", "{================ diffusion anisotropy restraint data =================}");
	print2file("dgsa.inp", "{* fixed or harmonically restrained external axis *}");
	print2file("dgsa.inp", "{+ choice: \"fixed\" \"harm\" +}");
	print2file("dgsa.inp", "{===>} nmr.dani.axis=\"harm\";");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* diffusion anisotropy restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.dani.file.1=\"\";");
	print2file("dgsa.inp", "{* diffusion anisotropy potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" +}");
	print2file("dgsa.inp", "{===>} nmr.dani.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* diffusion anisotropy initial force value *}");
	print2file("dgsa.inp", "{===>} nmr.dani.force.init.1=0.01;");
	print2file("dgsa.inp", "{* diffusion anisotropy final force value *}");
	print2file("dgsa.inp", "{===>} nmr.dani.force.finl.1=1.0;");
	print2file("dgsa.inp", "{* diffusion anisotropy coefficients *}");
	print2file("dgsa.inp", "{* coef: <Tc> <anis> <rhombicity> <wh> <wn> *}");
	print2file("dgsa.inp", "{* Tc = 1/2(Dx+Dy+Dz) in <ns> *} ");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.1.1=13.1;");
	print2file("dgsa.inp", "{* anis = Dz/0.5*(Dx+Dy) *} ");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.2.1=2.1;");
	print2file("dgsa.inp", "{* rhombicity = 1.5*(Dy-Dx)/(Dz-0.5*(Dy+Dx)) *} ");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.3.1=0.0;");
	print2file("dgsa.inp", "{* wH in <MHz> *} ");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.4.1=600.13;");
	print2file("dgsa.inp", "{* wN in <MHz> *}");
	print2file("dgsa.inp", "{===>} nmr.dani.coef.5.1=60.82;");
	print2file("dgsa.inp", "{============= susceptability anisotropy restraint data ===============}");
	print2file("dgsa.inp", "{* fixed or harmonically restrained external axis *}");
	print2file("dgsa.inp", "{+ choice: \"fixed\" \"harm\" +}");
	print2file("dgsa.inp", "{===>} nmr.sani.axis=\"harm\";");
	print2file("dgsa.inp", "{* Class 1 *}");
	print2file("dgsa.inp", "{* susceptability anisotropy restraints file *}");
	print2file("dgsa.inp", "{===>} nmr.sani.file.1=\"\";");
	print2file("dgsa.inp", "{* susceptability anisotropy potential *}");
	print2file("dgsa.inp", "{+ choice: \"harmonic\" \"square\" +}");
	print2file("dgsa.inp", "{===>} nmr.sani.pot.1=\"harmonic\";");
	print2file("dgsa.inp", "{* susceptability anisotropy initial force value *}");
	print2file("dgsa.inp", "{===>} nmr.sani.force.init.1=0.01;");
	print2file("dgsa.inp", "{* susceptability anisotropy final force value *}");
	print2file("dgsa.inp", "{===>} nmr.sani.force.finl.1=50.0;");
	print2file("dgsa.inp", "{* susceptability anisotropy coefficients *}");
	print2file("dgsa.inp", "{* coef: <DFS> <axial > <rhombicity>;");
	print2file("dgsa.inp", "   a0+a1*(3*cos(theta)^2-1)+a2*(3/2)*sin(theta)^2*cos(2*phi) *}");
	print2file("dgsa.inp", "{* DFS = a0 *}");
	print2file("dgsa.inp", "{===>} nmr.sani.coef.1.1=-0.0601;");
	print2file("dgsa.inp", "{* axial = a0-a1-3/2*a2 *}");
	print2file("dgsa.inp", "{===>} nmr.sani.coef.2.1=-8.02;");
	print2file("dgsa.inp", "{* rhombicity = a2/a1 *}");
	print2file("dgsa.inp", "{===>} nmr.sani.coef.3.1=0.4;");
	print2file("dgsa.inp", "{======================== other restraint data ========================}");
	print2file("dgsa.inp", "{* dihedral angle restraints file *}");
	print2file("dgsa.inp", "{* Note: the restraint file MUST NOT contain restraints ");
	print2file("dgsa.inp", "         dihedral or end *}");
	print2file("dgsa.inp", "{===>} nmr.cdih.file=\"dihedral.tbl\";");
	print2file("dgsa.inp", "{* DNA-RNA base planarity restraints file *}");
	print2file("dgsa.inp", "{* Note: include weights as \$pscale in the restraint file *}");
	print2file("dgsa.inp", "{===>} nmr.plan.file=\"\";");
	print2file("dgsa.inp", "{* input planarity scale factor - this will be written into \$pscale *}");
	print2file("dgsa.inp", "{===>} nmr.plan.scale=150;");
	print2file("dgsa.inp", "{* NCS-restraints file *}");
	print2file("dgsa.inp", "{* example is in inputs/xtal_data/eg1_ncs_restrain.dat *}");
	print2file("dgsa.inp", "{===>} nmr.ncs.file=\"\";");
	print2file("dgsa.inp", "{======================== input/output files ==========================}");
	print2file("dgsa.inp", "{* base name for input coordinate files *}");
	print2file("dgsa.inp", "{* used for simulated annealing when distance geometry is not used *}");
	print2file("dgsa.inp", "{===>} pdb.in.name=\"dg_sub_embed\";");
	print2file("dgsa.inp", "{* base name for output coordinate files *}");
	print2file("dgsa.inp", "{===>} pdb.out.name=\"$id\";");
	print2file("dgsa.inp", "{===========================================================================}");
	print2file("dgsa.inp", "{         things below this line do not normally need to be changed         }");
	print2file("dgsa.inp", "{         except for the torsion angle topology setup if you have           }");
	print2file("dgsa.inp", "{         molecules other than protein or nucleic acid                      }");
	print2file("dgsa.inp", "{===========================================================================}");
	print2file("dgsa.inp", "flg.cv.flag=false;");
	print2file("dgsa.inp", "flg.cv.noe=false;");
	print2file("dgsa.inp", "flg.cv.coup=false;");
	print2file("dgsa.inp", "flg.cv.cdih=false;");
	print2file("dgsa.inp", "nmr.cv.numpart=10;");
	print2file("dgsa.inp", " ) {- end block parameter definition -}");
	print2file("dgsa.inp", "checkversion 1.3");
	print2file("dgsa.inp", "evaluate (\$log_level=quiet)");
	print2file("dgsa.inp", "structure ");
	print2file("dgsa.inp", "   if  (&struct.1 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.1 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if  (&struct.2 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.2 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if  (&struct.3 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.3 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if  (&struct.4 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.4 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if  (&struct.5 # \"\") then");
	print2file("dgsa.inp", "      \@\@&struct.5 ");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "if ( &BLANK%pdb.in.file.1 = false ) then");
	print2file("dgsa.inp", "   coor \@\@&pdb.in.file.1");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "if ( &BLANK%pdb.in.file.2 = false ) then");
	print2file("dgsa.inp", "   coor \@\@&pdb.in.file.2");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "if ( &BLANK%pdb.in.file.3 = false ) then");
	print2file("dgsa.inp", "   coor \@\@&pdb.in.file.3");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "parameter");
	print2file("dgsa.inp", "   if (&par.1 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.1");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if (&par.2 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.2");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if (&par.3 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.3");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if (&par.4 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.4");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "   if (&par.5 # \"\") then");
	print2file("dgsa.inp", "      \@\@&par.5");
	print2file("dgsa.inp", "   end if");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "if ( \$log_level = verbose ) then");
	print2file("dgsa.inp", "  set message=normal echo=on end");
	print2file("dgsa.inp", "else");
	print2file("dgsa.inp", "  set message=off echo=off end");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "parameter");
	print2file("dgsa.inp", "   nbonds");
	print2file("dgsa.inp", "      repel=0.5");
	print2file("dgsa.inp", "      rexp=2 irexp=2 rcon=1.");
	print2file("dgsa.inp", "      nbxmod=-2");
	print2file("dgsa.inp", "      wmin=0.01");
	print2file("dgsa.inp", "      cutnb=4.5 ctonnb=2.99 ctofnb=3.");
	print2file("dgsa.inp", "      tolerance=0.5");
	print2file("dgsa.inp", "   end");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "set seed=&md.seed end");
	print2file("dgsa.inp", "{- Read experimental data -}");
	print2file("dgsa.inp", "   \@CNS_NMRMODULE:readdata ( nmr=&nmr;");
	print2file("dgsa.inp", "                             flag=&flg;");
	print2file("dgsa.inp", "                             output=\$nmr; )");
	print2file("dgsa.inp", "{- Read and store the number of NMR restraints -}");
	print2file("dgsa.inp", "   \@CNS_NMRMODULE:restraintnumber ( num=\$num; )");
	print2file("dgsa.inp", "   ");
	print2file("dgsa.inp", "{- Set mass and parameter values -}");
	print2file("dgsa.inp", "   ");
	print2file("dgsa.inp", "do (fbeta=10) (all)");
	print2file("dgsa.inp", "do (mass=100) (all)");
	print2file("dgsa.inp", "parameter                  ");
	print2file("dgsa.inp", "   nbonds  ");
	print2file("dgsa.inp", "      repel=0.80  ");
	print2file("dgsa.inp", "      rexp=2 irexp=2 rcon=1. ");
	print2file("dgsa.inp", "      nbxmod=3  ");
	print2file("dgsa.inp", "      wmin=0.01  ");
	print2file("dgsa.inp", "      cutnb=6.0 ctonnb=2.99 ctofnb=3.  ");
	print2file("dgsa.inp", "      tolerance=1.5  ");
	print2file("dgsa.inp", "   end  ");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "evaluate (\$nmr.trial.count = 0)    {- Initialize current structure number   -}");
	print2file("dgsa.inp", "evaluate (\$nmr.accept.count = 0)    {- Initialize number accepted            -}");
	print2file("dgsa.inp", "evaluate (\$nmr.counter 	= 0)");
	print2file("dgsa.inp", "evaluate (\$nmr.prev.counter = -1)");
	print2file("dgsa.inp", "\@CNS_NMRMODULE:initave  ( ave=\$ave;");
	print2file("dgsa.inp", "                          ave2=\$ave2;");
	print2file("dgsa.inp", "                          cv=\$cv;");
	print2file("dgsa.inp", "                          ener1=\$ener1;");
	print2file("dgsa.inp", "                          ener2=\$ener2;");
	print2file("dgsa.inp", "                          flag=&flg;");
	print2file("dgsa.inp", "                          nmr.prot=&nmr.prot; )");
	print2file("dgsa.inp", "        ");
	print2file("dgsa.inp", "{- Zero the force constant of disulfide bonds. -}");
	print2file("dgsa.inp", "parameter");
	print2file("dgsa.inp", "   bonds ( name SG ) ( name SG ) 0. TOKEN ");
	print2file("dgsa.inp", "end");
	print2file("dgsa.inp", "{- define a distance restraints for each disulfide bond, i.e., ");
	print2file("dgsa.inp", "   treat it as if it were an NOE and break the bond. -}");
	print2file("dgsa.inp", "for \$ss_rm_id_1 in id ( name SG ) loop STRM");
	print2file("dgsa.inp", "  for \$ss_rm_id_2 in id ( name SG and ");
	print2file("dgsa.inp", "			  bondedto ( id \$ss_rm_id_1 )  ) loop STR2");
	print2file("dgsa.inp", "    if (\$ss_rm_id_1 > \$ss_rm_id_2) then");
	print2file("dgsa.inp", "      pick bond ( id \$ss_rm_id_1 ) ( id \$ss_rm_id_2 ) equil");
	print2file("dgsa.inp", "      evaluate (\$ss_bond=\$result) ");
	print2file("dgsa.inp", "      noe ");
	print2file("dgsa.inp", "         assign ( id \$ss_rm_id_1 ) ( id \$ss_rm_id_2 ) \$ss_bond 0.1 0.1");
	print2file("dgsa.inp", "      end ");
	print2file("dgsa.inp", "    end if");
	print2file("dgsa.inp", "  end loop STR2");
	print2file("dgsa.inp", "end loop STRM");
	print2file("dgsa.inp", "{- Count the number of residues and determine molecule type -}");
	print2file("dgsa.inp", "identify (store9) (tag)");
	print2file("dgsa.inp", "evaluate (\$nmr.rsn.num = \$SELECT)");
	print2file("dgsa.inp", "identify (store9) ( tag and ( resn THY or resn CYT or resn GUA or");
	print2file("dgsa.inp", "                              resn ADE or resn URI ))");
	print2file("dgsa.inp", "evaluate (\$nmr.nucl.num = \$SELECT)    ");
	print2file("dgsa.inp", "if ( &md.dg.ref = \"coordinate\" ) then");
	print2file("dgsa.inp", "   flag exclude * include bond angl impr vdw end ");
	print2file("dgsa.inp", "   minimize lbfgs nstep=2000 drop=10.  nprint=100 end");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "do (refx=x) ( all )");
	print2file("dgsa.inp", "do (refy=y) ( all )");
	print2file("dgsa.inp", "do (refz=z) ( all )");
	print2file("dgsa.inp", "{- generate and store a bounds matrix -}");
	print2file("dgsa.inp", "if (&flg.dg.flag=true) then");
	print2file("dgsa.inp", "   flags exclude * include bond angle dihedral improper vdw noe cdih end");
	print2file("dgsa.inp", "   mmdg");
	print2file("dgsa.inp", "      shortest-path-algorithm=&&md.dg.algo");
	print2file("dgsa.inp", "      scale=&md.dg.scale");
	print2file("dgsa.inp", "      exponent=&md.dg.expo");
	print2file("dgsa.inp", "      baccuracy=&md.dg.bacc");
	print2file("dgsa.inp", "      taccuracy=&md.dg.tacc");
	print2file("dgsa.inp", "      iaccuracy=&md.dg.iacc");
	print2file("dgsa.inp", "      paccuracy=&md.dg.pacc");
	print2file("dgsa.inp", "      if (&md.dg.type=\"sub\") then");
	print2file("dgsa.inp", "   	 reference=&&md.dg.ref");
	print2file("dgsa.inp", "   	 storebounds");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "   	 reference=&&md.dg.ref");
	print2file("dgsa.inp", "   	 group &md.dg.group.slct &md.dg.group.err");
	print2file("dgsa.inp", "   	 storebounds");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "   end");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "   {- Begin protocol to generate structures distance geometry structures -}");
	print2file("dgsa.inp", "   while (&pdb.end.count > \$nmr.counter) loop dg");
	print2file("dgsa.inp", "      evaluate (\$nmr.counter=\$nmr.counter + 1)");
	print2file("dgsa.inp", "      evaluate (\$embedded=false)");
	print2file("dgsa.inp", "      flags exclude * include dg end");
	print2file("dgsa.inp", "      if (&md.dg.type=\"sub\") then");
	print2file("dgsa.inp", "   	 igroup interaction=(&md.dg.select) (&md.dg.select) end");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "      coor init end");
	print2file("dgsa.inp", "      while (\$embedded = false) loop embed");
	print2file("dgsa.inp", "   	 mmdg");
	print2file("dgsa.inp", "   	    if (&md.dg.type=\"sub\") then");
	print2file("dgsa.inp", "   	       recallbounds");
	print2file("dgsa.inp", "   	       substructure=(&md.dg.select)");
	print2file("dgsa.inp", "   	       selection=(&md.dg.select)");
	print2file("dgsa.inp", "   	    else");
	print2file("dgsa.inp", "   	       recallbounds");
	print2file("dgsa.inp", "   	       selection=(all)");
	print2file("dgsa.inp", "   	       if (&md.dg.metr.flag=true) then");
	print2file("dgsa.inp", "   		  &&md.dg.ord");
	print2file("dgsa.inp", "   		  metrization=(&md.dg.metr.atom)=&md.dg.metr.num");
	print2file("dgsa.inp", "   	       end if");
	print2file("dgsa.inp", "   	    end if");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "      end loop embed");
	print2file("dgsa.inp", "      do (x = x * \$dgscale) (known)");
	print2file("dgsa.inp", "      do (y = y * \$dgscale) (known)");
	print2file("dgsa.inp", "      do (z = z * \$dgscale) (known)");
	print2file("dgsa.inp", "      minimize lbfgs");
	print2file("dgsa.inp", "   	 nstep=&md.dg.step drop=1. nprint=25");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", "      \@CNS_NMRMODULE:printdg ( md=&md;");
	print2file("dgsa.inp", "                               output=\$nmr;");
	print2file("dgsa.inp", "                               pdb=&pdb; )");
	print2file("dgsa.inp", "   end loop dg");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "{- initialize and set scaling factors for simulated annealing -}");
	print2file("dgsa.inp", "set seed=&md.seed end");
	print2file("dgsa.inp", "evaluate (\$nmr.trial.count = 0)    {- Initialize current structure number   -}");
	print2file("dgsa.inp", "evaluate (\$nmr.dg.count = 0)");
	print2file("dgsa.inp", "evaluate (\$nmr.accept.count = 0)   {- Initialize number accepted            -}");
	print2file("dgsa.inp", "evaluate (\$nmr.counter = 0)");
	print2file("dgsa.inp", "evaluate (\$coor_count_init=0.)");
	print2file("dgsa.inp", "evaluate (\$coor_input_count=0.)");
	print2file("dgsa.inp", "if (&flg.dg.flag=true) then");
	print2file("dgsa.inp", "   evaluate (\$coor_input_count=&pdb.end.count)");
	print2file("dgsa.inp", "else");
	print2file("dgsa.inp", "   evaluate (\$coor_input_count=&pdb.dg.count)");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "\@CNS_NMRMODULE:initave  ( flag=&flg;");
	print2file("dgsa.inp", "                          ave=\$ave;");
	print2file("dgsa.inp", "                          ave2=\$ave2;");
	print2file("dgsa.inp", "                          cv=\$cv;");
	print2file("dgsa.inp", "                          ener1=\$ener1;");
	print2file("dgsa.inp", "                          ener2=\$ener2;");
	print2file("dgsa.inp", "                          nmr.prot=&nmr.prot; )");
	print2file("dgsa.inp", "        ");
	print2file("dgsa.inp", "{- scaling of nmr restraint data during regularization -}");
	print2file("dgsa.inp", "\@CNS_CUSTOMMODULE:scalehotedited ( md=&md;");
	print2file("dgsa.inp", "                          nmr=&nmr;");
	print2file("dgsa.inp", "                          input.noe.scale=&md.cool.noe;");
	print2file("dgsa.inp", "                          input.cdih.scale=&md.hot.cdih; )");
	print2file("dgsa.inp", "if (&nmr.dani.axis = \"harm\") then");
	print2file("dgsa.inp", "   do (harmonic=20.0) (resid 500 and name OO)");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name Z )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name X )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name Y )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (not (resid 500))");
	print2file("dgsa.inp", "   restraints harmonic exponent=2 end");
	print2file("dgsa.inp", "elseif (&nmr.sani.axis = \"harm\") then");
	print2file("dgsa.inp", "   do (harmonic=20.0) (resid 500 and name OO)");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name Z )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name X )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (resid 500 and name Y )");
	print2file("dgsa.inp", "   do (harmonic=0.0) (not (resid 500))");
	print2file("dgsa.inp", "   restraints harmonic exponent=2 end");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "{- Increase the disulfide bond force constants to their full strength -}");
	print2file("dgsa.inp", "   parameter");
	print2file("dgsa.inp", "      bonds ( name SG ) ( name SG ) 1000. TOKEN ");
	print2file("dgsa.inp", "   end");
	print2file("dgsa.inp", "{- Regularize structures generated by distance geometry - loop until done -}");
	print2file("dgsa.inp", "if (&flg.dgsa.flag=true) then");
	print2file("dgsa.inp", "   while (&pdb.end.count > \$nmr.counter) loop dgsa");
	print2file("dgsa.inp", "      {- Set parameter values -}");
	print2file("dgsa.inp", "      parameter");
	print2file("dgsa.inp", "         nbonds");
	print2file("dgsa.inp", "            repel=0.5");
	print2file("dgsa.inp", "            rexp=2 irexp=2 rcon=1.");
	print2file("dgsa.inp", "            nbxmod=-2");
	print2file("dgsa.inp", "            wmin=0.01");
	print2file("dgsa.inp", "            cutnb=4.5 ctonnb=2.99 ctofnb=3.");
	print2file("dgsa.inp", "            tolerance=0.5");
	print2file("dgsa.inp", "         end");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", "      evaluate (\$nmr.trial.count = \$nmr.trial.count + 1)");
	print2file("dgsa.inp", "      if (\$nmr.trial.count <= \$coor_input_count) then");
	print2file("dgsa.inp", "         evaluate (\$nmr.dg.count=\$nmr.dg.count+1)");
	print2file("dgsa.inp", "         evaluate (\$coor_count_init=0.)");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "         evaluate (\$coor_count_init=\$coor_count_init+1)");
	print2file("dgsa.inp", "         if (\$coor_count_init > \$coor_input_count ) then");
	print2file("dgsa.inp", "            evaluate (\$coor_count_init=1)");
	print2file("dgsa.inp", "         end if");
	print2file("dgsa.inp", "   	 evaluate (\$nmr.dg.count=\$coor_count_init)");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "   {- \$prefix is generated in the macro printdg -}");
	print2file("dgsa.inp", "      if (&flg.dg.flag=true) then");
	print2file("dgsa.inp", "         evaluate (\$filename=\$nmr.prefix+encode(\$nmr.dg.count)+\".pdb\")");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "         evaluate (\$filename=&pdb.in.name+\"_\"+encode(\$nmr.dg.count)+\".pdb\")");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "         ");
	print2file("dgsa.inp", "      {- Test for correct enantiomer -}");
	print2file("dgsa.inp", "      for \$image in ( 1 -1 ) loop imag");
	print2file("dgsa.inp", "         set remarks=reset end ");
	print2file("dgsa.inp", "   	 coor initialize end");
	print2file("dgsa.inp", "   	 coor \@\@\$filename");
	print2file("dgsa.inp", "   	 do (x=x * \$image) ( known )");
	print2file("dgsa.inp", "   	 identity (store1) (not known)");
	print2file("dgsa.inp", "   	 coor copy end");
	print2file("dgsa.inp", "   	 do (x=refx) ( all )");
	print2file("dgsa.inp", "   	 do (y=refy) ( all )");
	print2file("dgsa.inp", "   	 do (z=refz) ( all )");
	print2file("dgsa.inp", "   	 for \$id in id ( tag ) loop fit");
	print2file("dgsa.inp", "   	    coordinates");
	print2file("dgsa.inp", "   	       fit select = ( byresidue (id \$id) and not store1 )");
	print2file("dgsa.inp", "   	    end");
	print2file("dgsa.inp", "   	   coor copy selection=( byresidue (id \$id) ) end");
	print2file("dgsa.inp", "   	 end loop fit");
	print2file("dgsa.inp", "   	 coor swap end");
	print2file("dgsa.inp", "         if (&nmr.dani.axis = \"fixed\" ) then");
	print2file("dgsa.inp", "            fix");
	print2file("dgsa.inp", "               select=(resname ANI)");
	print2file("dgsa.inp", "            end");
	print2file("dgsa.inp", "         elseif (&nmr.sani.axis = \"fixed\" ) then");
	print2file("dgsa.inp", "            fix");
	print2file("dgsa.inp", "               select=(resname ANI)");
	print2file("dgsa.inp", "            end");
	print2file("dgsa.inp", "         end if");
	print2file("dgsa.inp", "   	 parameter");
	print2file("dgsa.inp", "   	    nbonds");
	print2file("dgsa.inp", "   	       nbxmod=-2");
	print2file("dgsa.inp", "   	       repel=0.5");
	print2file("dgsa.inp", "   	    end");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 flags exclude * include bond vdw noe cdih coup oneb ");
	print2file("dgsa.inp", "   	                         carb ncs dani sani harm end");
	print2file("dgsa.inp", "   	 igroup interaction (all) (all) weights * 1.  vdw 20. end end");
	print2file("dgsa.inp", "   	 minimize lbfgs nstep=100 nprint=10 end");
	print2file("dgsa.inp", "   	 flags include angl end");
	print2file("dgsa.inp", "   	 minimize lbfgs nstep=100 nprint=10 end");
	print2file("dgsa.inp", "   	 flags include impr dihe end");
	print2file("dgsa.inp", "   	 evaluate (\$nstep1 = int(&md.hot.step/8))");
	print2file("dgsa.inp", "   	 evaluate (\$nstep2 = int(&md.hot.step/2))");
	print2file("dgsa.inp", "   	 do ( vx = maxwell(0.5) ) ( all )");
	print2file("dgsa.inp", "   	 do ( vy = maxwell(0.5) ) ( all )");
	print2file("dgsa.inp", "   	 do ( vz = maxwell(0.5) ) ( all )");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 0.1 impr 0.05 vdw 20. end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep1");
	print2file("dgsa.inp", "   	    nprint=\$nstep1");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 0.2 impr 0.1  vdw 20. end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep1");
	print2file("dgsa.inp", "   	    nprint=\$nstep1");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 parameter  nbonds repel=0.9   end  end");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 0.2 impr 0.2 vdw 0.01 end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep1");
	print2file("dgsa.inp", "   	    nprint=\$nstep1");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 parameter nbonds nbxmod=-3  end  end");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 0.4 impr 0.4 vdw 0.003 end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep2");
	print2file("dgsa.inp", "   	    nprint=\$nstep2");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 igroup inter (all) (all) weights * 1.0 impr 1.0 vdw 0.003 end end");
	print2file("dgsa.inp", "   	 dynamics cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling=false");
	print2file("dgsa.inp", "   	    tcoupling=true");
	print2file("dgsa.inp", "   	    timestep=&md.hot.ss");
	print2file("dgsa.inp", "   	    nstep=\$nstep1");
	print2file("dgsa.inp", "   	    nprint=\$nstep1");
	print2file("dgsa.inp", "   	    temperature=&md.hot.temp");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 if (\$image = 1) then");
	print2file("dgsa.inp", "   	    do (store7=x) ( all )");
	print2file("dgsa.inp", "   	    do (store8=y) ( all )");
	print2file("dgsa.inp", "   	    do (store9=z) ( all )");
	print2file("dgsa.inp", "   	    do (store4=vx) ( all )");
	print2file("dgsa.inp", "   	    do (store5=vy) ( all )");
	print2file("dgsa.inp", "   	    do (store6=vz) ( all )");
	print2file("dgsa.inp", "   	 end if");
	print2file("dgsa.inp", "      end loop imag");
	print2file("dgsa.inp", "      {- Establish the correct handedness of the structure -}");
	print2file("dgsa.inp", "      energy end");
	print2file("dgsa.inp", "      evaluate (\$e_minus=\$ener)");
	print2file("dgsa.inp", "      coor copy end");
	print2file("dgsa.inp", "      do (x=store7) ( all )");
	print2file("dgsa.inp", "      do (y=store8) ( all )");
	print2file("dgsa.inp", "      do (z=store9) ( all )");
	print2file("dgsa.inp", "      energy end");
	print2file("dgsa.inp", "      evaluate (\$e_plus=\$ener)");
	print2file("dgsa.inp", "      if ( \$e_plus > \$e_minus ) then");
	print2file("dgsa.inp", "   	 evaluate (\$hand=-1 )");
	print2file("dgsa.inp", "   	 coor swap end");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "   	 evaluate (\$hand= 1 )");
	print2file("dgsa.inp", "   	 do (vx=store4) ( all )");
	print2file("dgsa.inp", "   	 do (vy=store5) ( all )");
	print2file("dgsa.inp", "   	 do (vz=store6) ( all )");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "   {- Slow-cooling with cartesian dynamics -}");
	print2file("dgsa.inp", "      parameter");
	print2file("dgsa.inp", "   	 nbonds");
	print2file("dgsa.inp", "   	    repel=0.80");
	print2file("dgsa.inp", "   	    rexp=2 irexp=2 rcon=1.");
	print2file("dgsa.inp", "   	    nbxmod=3");
	print2file("dgsa.inp", "   	    wmin=0.01");
	print2file("dgsa.inp", "   	    cutnb=6.0 ctonnb=2.99 ctofnb=3.");
	print2file("dgsa.inp", "   	    tolerance=0.5");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", "      flags include plan end");
	print2file("dgsa.inp", "      evaluate (\$final_t = 0)");
	print2file("dgsa.inp", "      evaluate (\$ncycle = int((&md.hot.temp-\$final_t)/&md.cool.tmpstp))");
	print2file("dgsa.inp", "      evaluate (\$nstep = int(&md.cool.step/\$ncycle))");
	print2file("dgsa.inp", "      evaluate (\$vdw_step=(&md.cool.vdw.finl/&md.cool.vdw.init)^(1/\$ncycle))");
	print2file("dgsa.inp", "      evaluate (\$rad_step=(&md.cool.init.rad-&md.cool.fina.rad)/\$ncycle)");
	print2file("dgsa.inp", "      evaluate (\$radius=&&md.cool.init.rad)");
	print2file("dgsa.inp", "      {- set up nmr restraint scaling -}");
	print2file("dgsa.inp", "      evaluate (\$kdani.inter.flag=false)");
	print2file("dgsa.inp", "      evaluate (\$ksani.inter.flag=false)");
	print2file("dgsa.inp", "      evaluate (\$kdani.cart.flag=false)");
	print2file("dgsa.inp", "      evaluate (\$ksani.cart.flag=false)");
	print2file("dgsa.inp", "      \@CNS_CUSTOMMODULE:scalecoolsetupedited ( kdani=\$kdani;");
	print2file("dgsa.inp", "                                      ksani=\$ksani;");
	print2file("dgsa.inp", "                                      nmr=&nmr;");
	print2file("dgsa.inp", "                                      input.noe.scale=&md.cool.noe;");
	print2file("dgsa.inp", "                                      input.cdih.scale=&md.cool.cdih;");
	print2file("dgsa.inp", "                                      input.ncycle=\$ncycle; )");
	print2file("dgsa.inp", "      evaluate (\$bath=&md.hot.temp)");
	print2file("dgsa.inp", "      evaluate (\$k_vdw=&md.cool.vdw.init)");
	print2file("dgsa.inp", "      evaluate (\$i_cool = 0)");
	print2file("dgsa.inp", "      while (\$i_cool <= \$ncycle) loop cool");
	print2file("dgsa.inp", "   	 evaluate (\$i_cool = \$i_cool + 1)");
	print2file("dgsa.inp", "   	 igroup");
	print2file("dgsa.inp", "   	    interaction (chemical h*) (all) weights * 1 vdw 0. elec 0. end");
	print2file("dgsa.inp", "   	    interaction (not chemical h*) (not chemical h*) weights * 1 vdw \$k_vdw end");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 dynamics  cartesian");
	print2file("dgsa.inp", "   	    cmremove=true");
	print2file("dgsa.inp", "   	    vscaling = true");
	print2file("dgsa.inp", "   	    tcoup = false");
	print2file("dgsa.inp", "   	    timestep = &md.cool.ss");
	print2file("dgsa.inp", "   	    nstep = \$nstep");
	print2file("dgsa.inp", "   	    nprint = \$nstep");
	print2file("dgsa.inp", "   	    temperature = \$bath");
	print2file("dgsa.inp", "   	 end");
	print2file("dgsa.inp", "   	 evaluate (\$radius=max(&md.cool.fina.rad,\$radius-\$rad_step))");
	print2file("dgsa.inp", "   	 parameter  nbonds repel=\$radius   end end");
	print2file("dgsa.inp", "   	 evaluate (\$k_vdw=min(&md.cool.vdw.finl,\$k_vdw*\$vdw_step))");
	print2file("dgsa.inp", "   	 evaluate (\$bath=\$bath-&md.cool.tmpstp)");
	print2file("dgsa.inp", "         \@CNS_NMRMODULE:scalecool ( kdani=\$kdani;");
	print2file("dgsa.inp", "                                    ksani=\$ksani;");
	print2file("dgsa.inp", "                                    nmr=&nmr; )");
	print2file("dgsa.inp", "      end loop cool");
	print2file("dgsa.inp", "   {- Final minimization -}");
	print2file("dgsa.inp", "      {- turn on proton chemical shifts -}");
	print2file("dgsa.inp", "      flags include prot end");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "      if (\$nmr.nucl.num > 0) then");
	print2file("dgsa.inp", "         flags include elec end");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "      noe             ");
	print2file("dgsa.inp", "         scale * &md.pow.noe ");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", "        ");
	print2file("dgsa.inp", "      restraints dihedral  ");
	print2file("dgsa.inp", "         scale = &md.pow.cdih  ");
	print2file("dgsa.inp", "      end");
	print2file("dgsa.inp", " 													");
	print2file("dgsa.inp", "      igroup interaction ( all ) ( all ) weights * 1 end end");
	print2file("dgsa.inp", "      evaluate (\$count=0 )");
	print2file("dgsa.inp", "      while (&md.pow.cycl > \$count) loop pmini");
	print2file("dgsa.inp", "         evaluate (\$count=\$count + 1)");
	print2file("dgsa.inp", "         minimize lbfgs nstep=&md.pow.step drop=10.0 nprint=25 end");
	print2file("dgsa.inp", "      end loop pmini");
	print2file("dgsa.inp", "      evaluate (\$nmr.min.num = \$count * &md.pow.step)");
	print2file("dgsa.inp", "      {- translate the geometric center of the structure to the origin -}");
	print2file("dgsa.inp", "      if (\$num.dani > 0. ) then");
	print2file("dgsa.inp", "      elseif (\$num.sani > 0. ) then");
	print2file("dgsa.inp", "      else");
	print2file("dgsa.inp", "         show ave ( x ) ( all )");
	print2file("dgsa.inp", "         evaluate (\$geom_x=-\$result)");
	print2file("dgsa.inp", "         show ave ( y ) ( all )");
	print2file("dgsa.inp", "         evaluate (\$geom_y=-\$result)");
	print2file("dgsa.inp", "         show ave ( z ) ( all )");
	print2file("dgsa.inp", "         evaluate (\$geom_z=-\$result)");
	print2file("dgsa.inp", "         coor translate vector=( \$geom_x \$geom_y \$geom_z ) selection=( all ) end");
	print2file("dgsa.inp", "      end if");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "      \@CNS_NMRMODULE:printaccept ( ave=\$ave;");
	print2file("dgsa.inp", "                                   ave2=\$ave2;");
	print2file("dgsa.inp", "                                   cv=\$cv;");
	print2file("dgsa.inp", "                                   ener1=\$ener1;");
	print2file("dgsa.inp", "                                   ener2=\$ener2;");
	print2file("dgsa.inp", "                                   flag=&flg;");
	print2file("dgsa.inp", "                                   md=&md;");
	print2file("dgsa.inp", "                                   nmr=&nmr;");
	print2file("dgsa.inp", "                                   num=\$num;");
	print2file("dgsa.inp", "                                   output=\$nmr;");
	print2file("dgsa.inp", "                                   pdb=&pdb;  )");
	print2file("dgsa.inp", "   end loop dgsa");
	print2file("dgsa.inp", "   \@CNS_NMRMODULE:calcave ( ave=\$ave;                 ");
	print2file("dgsa.inp", "                            ave2=\$ave2;               ");
	print2file("dgsa.inp", "                            cv=\$cv;                   ");
	print2file("dgsa.inp", "                            ener1=\$ener1;               ");
	print2file("dgsa.inp", "                            ener2=\$ener2;             ");
	print2file("dgsa.inp", "                            flag=&flg;               ");
	print2file("dgsa.inp", "                            md=&md;");
	print2file("dgsa.inp", "                            nmr=&nmr;");
	print2file("dgsa.inp", "                            num=\$num;                 ");
	print2file("dgsa.inp", "                            output=\$nmr;           ");
	print2file("dgsa.inp", "                            pdb=&pdb;  )");
	print2file("dgsa.inp", "	");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "      ");
	print2file("dgsa.inp", "end if");
	print2file("dgsa.inp", "stop");
}

sub write_cns_generate_seq_file{
	print2file("gseq.inp", "{+ file: generate_seq.inp +}");
	print2file("gseq.inp", "{+ directory: general +}");
	print2file("gseq.inp", "{+ description: Generate structure file for protein, dna/rna, water, ");
	print2file("gseq.inp", "                ligands and/or carbohydrate from sequence information only +}");
	print2file("gseq.inp", "{+ comment: modified by Brian Smith (Edinburgh University) to allow protein");
	print2file("gseq.inp", "            residue renumbering +}");
	print2file("gseq.inp", "{+ authors: Paul Adams, and Axel Brunger +}");
	print2file("gseq.inp", "{+ copyright: Yale University +}");
	print2file("gseq.inp", "{- Guidelines for using this file:");
	print2file("gseq.inp", "   - all strings must be quoted by double-quotes");
	print2file("gseq.inp", "   - logical variables (true/false) are not quoted");
	print2file("gseq.inp", "   - do not remove any evaluate statements from the file -}");
	print2file("gseq.inp", "{- Special patches will have to be entered manually at the relevant points");
	print2file("gseq.inp", "   in the file - see comments throughout the file -}");
	print2file("gseq.inp", "{- begin block parameter definition -} define(");
	print2file("gseq.inp", "{============ protein topology, linkage, and parameter files =============}");
	print2file("gseq.inp", "{* topology files *}");
	print2file("gseq.inp", "{===>} topology_infile_1=\"CNS_TOPPAR:protein.top\";");
	print2file("gseq.inp", "{===>} topology_infile_2=\"CNS_TOPPAR:dna-rna.top\";");
	print2file("gseq.inp", "{===>} topology_infile_3=\"CNS_TOPPAR:water.top\";");
	print2file("gseq.inp", "{===>} topology_infile_4=\"CNS_TOPPAR:ion.top\";");
	print2file("gseq.inp", "{===>} topology_infile_5=\"CNS_TOPPAR:carbohydrate.top\";");
	print2file("gseq.inp", "{===>} topology_infile_6=\"\";");
	print2file("gseq.inp", "{===>} topology_infile_7=\"\";");
	print2file("gseq.inp", "{===>} topology_infile_8=\"\";");
	print2file("gseq.inp", "{* linkage files for linear, continuous polymers (protein, DNA, RNA) *}");
	print2file("gseq.inp", "{===>} link_infile_1=\"CNS_TOPPAR:protein.link\";");
	print2file("gseq.inp", "{===>} link_infile_2=\"CNS_TOPPAR:dna-rna-pho.link\";");
	print2file("gseq.inp", "{===>} link_infile_3=\"\";");
	print2file("gseq.inp", "{* parameter files *}");
	print2file("gseq.inp", "{===>} parameter_infile_1=\"CNS_TOPPAR:protein.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_2=\"CNS_TOPPAR:dna-rna_rep.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_3=\"CNS_TOPPAR:water_rep.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_4=\"CNS_TOPPAR:ion.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_5=\"CNS_TOPPAR:carbohydrate.param\";");
	print2file("gseq.inp", "{===>} parameter_infile_6=\"\";");
	print2file("gseq.inp", "{===>} parameter_infile_7=\"\";");
	print2file("gseq.inp", "{===>} parameter_infile_8=\"\";");
	print2file("gseq.inp", "{====================== other linkages and modifications  ==================}");
	print2file("gseq.inp", "{* extra linkages and modifications by custom patches *}");
	print2file("gseq.inp", "{===>} patch_infile=\"\";");
	print2file("gseq.inp", "{============================= sequence files ==============================}");
	print2file("gseq.inp", "{* multiple sequence files of the same type can be defined by duplicating");
	print2file("gseq.inp", "   the entries below and incrementing the file number *}");
	print2file("gseq.inp", "{* protein sequence file 1 *}");
	print2file("gseq.inp", "{===>} prot_sequence_infile_1=\"input.seq\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prot_segid_1=\"\";");
	print2file("gseq.inp", "{* start residue numbering at *}");
	print2file("gseq.inp", "{===>} renumber_1=1;");
	print2file("gseq.inp", "{* protein sequence file 2 *}");
	print2file("gseq.inp", "{===>} prot_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prot_segid_2=\"\";");
	print2file("gseq.inp", "{* start residue numbering at *}");
	print2file("gseq.inp", "{===>} renumber_2=1;");
	print2file("gseq.inp", "{* protein sequence file 3 *}");
	print2file("gseq.inp", "{===>} prot_sequence_infile_3=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prot_segid_3=\"\";");
	print2file("gseq.inp", "{* start residue numbering at *}");
	print2file("gseq.inp", "{===>} renumber_3=1;");
	print2file("gseq.inp", "{* protein sequence file 4 *}");
	print2file("gseq.inp", "{===>} prot_sequence_infile_4=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prot_segid_4=\"\";");
	print2file("gseq.inp", "{* start residue numbering at *}");
	print2file("gseq.inp", "{===>} renumber_4=1;");
	print2file("gseq.inp", "{* nucleic acid sequence file 1 *}");
	print2file("gseq.inp", "{===>} nucl_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} nucl_segid_1=\"\";");
	print2file("gseq.inp", "{* nucleic acid sequence file 2 *}");
	print2file("gseq.inp", "{===>} nucl_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} nucl_segid_2=\"\";");
	print2file("gseq.inp", "{* nucleic acid sequence file 3 *}");
	print2file("gseq.inp", "{===>} nucl_sequence_infile_3=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} nucl_segid_3=\"\";");
	print2file("gseq.inp", "{* nucleic acid sequence file 4 *}");
	print2file("gseq.inp", "{===>} nucl_sequence_infile_4=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} nucl_segid_4=\"\";");
	print2file("gseq.inp", "{* water sequence file 1 *}");
	print2file("gseq.inp", "{===>} water_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} water_segid_1=\"\";");
	print2file("gseq.inp", "{* water sequence file 2 *}");
	print2file("gseq.inp", "{===>} water_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} water_segid_2=\"\";");
	print2file("gseq.inp", "{* carbohydrate sequence file 1 *}");
	print2file("gseq.inp", "{===>} carbo_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} carbo_segid_1=\"\";");
	print2file("gseq.inp", "{* carbohydrate sequence file 2 *}");
	print2file("gseq.inp", "{===>} carbo_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} carbo_segid_2=\"\";");
	print2file("gseq.inp", "{* prosthetic group sequence file 1 *}");
	print2file("gseq.inp", "{===>} prost_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prost_segid_1=\"\";");
	print2file("gseq.inp", "{* prosthetic group sequence file 2 *}");
	print2file("gseq.inp", "{===>} prost_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} prost_segid_2=\"\";");
	print2file("gseq.inp", "{* ligand sequence file 1 *}");
	print2file("gseq.inp", "{===>} lig_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} lig_segid_1=\"\";");
	print2file("gseq.inp", "{* ligand sequence file 2 *}");
	print2file("gseq.inp", "{===>} lig_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} lig_segid_2=\"\";");
	print2file("gseq.inp", "{* ion sequence file 1 *}");
	print2file("gseq.inp", "{===>} ion_sequence_infile_1=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} ion_segid_1=\"\";");
	print2file("gseq.inp", "{* ion sequence file 2 *}");
	print2file("gseq.inp", "{===>} ion_sequence_infile_2=\"\";");
	print2file("gseq.inp", "{* segid *}");
	print2file("gseq.inp", "{===>} ion_segid_2=\"\";");
	print2file("gseq.inp", "{============================= output files ================================}");
	print2file("gseq.inp", "{* output structure file *}");
	print2file("gseq.inp", "{===>} structure_outfile=\"extended.mtf\";");
	print2file("gseq.inp", "{=========================== disulphide bonds ==============================}");
	print2file("gseq.inp", "{* Select pairs of cysteine residues that form disulphide bonds *}");
	print2file("gseq.inp", "{* First 2 entries are the segid and resid of the first cysteine (CYS A). *}");
	print2file("gseq.inp", "{* Second 2 entries are the segid and resid of the second cysteine (CYS B). *}");
	print2file("gseq.inp", "{+ table: rows=8 numbered");
	print2file("gseq.inp", "   cols=5 \"use\" \"segid CYS A\" \"resid CYS A\" \"segid CYS B\" \"resid CYS B\" +}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_1=true;");
	print2file("gseq.inp", "{===>} ss_i_segid_1=\"\"; ss_i_resid_1=11;");
	print2file("gseq.inp", "{===>} ss_j_segid_1=\"\"; ss_j_resid_1=27;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_2=true;");
	print2file("gseq.inp", "{===>} ss_i_segid_2=\"\"; ss_i_resid_2=45;");
	print2file("gseq.inp", "{===>} ss_j_segid_2=\"\"; ss_j_resid_2=73;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_3=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_3=\"\"; ss_i_resid_3=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_3=\"\"; ss_j_resid_3=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_4=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_4=\"\"; ss_i_resid_4=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_4=\"\"; ss_j_resid_4=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_5=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_5=\"\"; ss_i_resid_5=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_5=\"\"; ss_j_resid_5=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_6=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_6=\"\"; ss_i_resid_6=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_6=\"\"; ss_j_resid_6=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_7=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_7=\"\"; ss_i_resid_7=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_7=\"\"; ss_j_resid_7=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} ss_use_8=false;");
	print2file("gseq.inp", "{===>} ss_i_segid_8=\"\"; ss_i_resid_8=0;");
	print2file("gseq.inp", "{===>} ss_j_segid_8=\"\"; ss_j_resid_8=0;");
	print2file("gseq.inp", "{=========================== carbohydrate links  ===========================}");
	print2file("gseq.inp", "{* Select pairs of residues that are linked *}");
	print2file("gseq.inp", "{* First entry is the name of the patch residue. *}");
	print2file("gseq.inp", "{* Second and third entries are the resid and segid for the atoms");
	print2file("gseq.inp", "   referenced by \"-\" in the patch. *}");
	print2file("gseq.inp", "{* Fourth and fifth entries are the resid and segid for the atoms");
	print2file("gseq.inp", "   referenced by \"+\" in the patch *}");
	print2file("gseq.inp", "{+ table: rows=6 numbered");
	print2file("gseq.inp", "          cols=6 \"use\" \"patch name\" \"segid -\" \"resid -\" \"segid +\" \"resid +\" +}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_1=false;");
	print2file("gseq.inp", "{===>} carbo_patch_1=\"B1N\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_1=\"BBBB\"; carbo_i_resid_1=401;");
	print2file("gseq.inp", "{===>} carbo_j_segid_1=\"AAAA\"; carbo_j_resid_1=56;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_2=false;");
	print2file("gseq.inp", "{===>} carbo_patch_2=\"B1N\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_2=\"BBBB\"; carbo_i_resid_2=402;");
	print2file("gseq.inp", "{===>} carbo_j_segid_2=\"AAAA\"; carbo_j_resid_2=182;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_3=false;");
	print2file("gseq.inp", "{===>} carbo_patch_3=\"\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_3=\"\"; carbo_i_resid_3=0;");
	print2file("gseq.inp", "{===>} carbo_j_segid_3=\"\"; carbo_j_resid_3=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_4=false;");
	print2file("gseq.inp", "{===>} carbo_patch_4=\"\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_4=\"\"; carbo_i_resid_4=0;");
	print2file("gseq.inp", "{===>} carbo_j_segid_4=\"\"; carbo_j_resid_4=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_5=false;");
	print2file("gseq.inp", "{===>} carbo_patch_5=\"\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_5=\"\"; carbo_i_resid_5=0;");
	print2file("gseq.inp", "{===>} carbo_j_segid_5=\"\"; carbo_j_resid_5=0;");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} carbo_use_6=false;");
	print2file("gseq.inp", "{===>} carbo_patch_6=\"\";");
	print2file("gseq.inp", "{===>} carbo_i_segid_6=\"\"; carbo_i_resid_6=0;");
	print2file("gseq.inp", "{===>} carbo_j_segid_6=\"\"; carbo_j_resid_6=0;");
	print2file("gseq.inp", "{========================= generate parameters =============================}");
	print2file("gseq.inp", "{* hydrogen flag - determines whether hydrogens will be retained *}");
	print2file("gseq.inp", "{* must be true for NMR, atomic resolution X-ray crystallography ");
	print2file("gseq.inp", "   or modelling.  Set to false for most X-ray crystallographic ");
	print2file("gseq.inp", "   applications at resolution > 1A *}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} hydrogen_flag=false;");
	print2file("gseq.inp", "{* set bfactor flag *}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} set_bfactor=true;");
	print2file("gseq.inp", "{* set bfactor value *}");
	print2file("gseq.inp", "{===>} bfactor=15.0;");
	print2file("gseq.inp", "{* set occupancy flag *}");
	print2file("gseq.inp", "{+ choice: true false +}");
	print2file("gseq.inp", "{===>} set_occupancy=true;");
	print2file("gseq.inp", "{* set occupancy value *}");
	print2file("gseq.inp", "{===>} occupancy=1.0;");
	print2file("gseq.inp", "{===========================================================================}");
	print2file("gseq.inp", "{         things below this line do not need to be changed unless           }");
	print2file("gseq.inp", "{         you need to apply patches - at the appropriate places marked      }");
	print2file("gseq.inp", "{===========================================================================}");
	print2file("gseq.inp", " ) {- end block parameter definition -}");
	print2file("gseq.inp", " checkversion 1.3");
	print2file("gseq.inp", " evaluate (\$log_level=quiet)");
	print2file("gseq.inp", " {- read parameter files -}");
	print2file("gseq.inp", " parameter");
	print2file("gseq.inp", "  evaluate (\$counter=1)");
	print2file("gseq.inp", "  evaluate (\$done=false)");
	print2file("gseq.inp", "  while ( \$done = false ) loop read");
	print2file("gseq.inp", "   if ( &exist_parameter_infile_\$counter = true ) then");
	print2file("gseq.inp", "      if ( &BLANK%parameter_infile_\$counter = false ) then");
	print2file("gseq.inp", "         \@\@&parameter_infile_\$counter");
	print2file("gseq.inp", "      end if");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "    evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", "   evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "  end loop read");
	print2file("gseq.inp", " end");
	print2file("gseq.inp", " {- read topology files -}");
	print2file("gseq.inp", " topology");
	print2file("gseq.inp", "  evaluate (\$counter=1)");
	print2file("gseq.inp", "  evaluate (\$done=false)");
	print2file("gseq.inp", "  while ( \$done = false ) loop read");
	print2file("gseq.inp", "   if ( &exist_topology_infile_\$counter = true ) then");
	print2file("gseq.inp", "      if ( &BLANK%topology_infile_\$counter = false ) then");
	print2file("gseq.inp", "         \@\@&topology_infile_\$counter");
	print2file("gseq.inp", "      end if");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", "   evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "  end loop read");
	print2file("gseq.inp", " end");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop prot");
	print2file("gseq.inp", "   if ( &exist_prot_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%prot_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           evaluate (\$count=1)");
	print2file("gseq.inp", "           evaluate (\$done2=false)");
	print2file("gseq.inp", "           while ( \$done2 = false ) loop read");
	print2file("gseq.inp", "             if ( &exist_link_infile_\$count = true ) then");
	print2file("gseq.inp", "               if ( &BLANK%link_infile_\$count = false ) then");
	print2file("gseq.inp", "                  \@\@&link_infile_\$count");
	print2file("gseq.inp", "               end if");
	print2file("gseq.inp", "             else");
	print2file("gseq.inp", "               evaluate (\$done2=true)");
	print2file("gseq.inp", "             end if");
	print2file("gseq.inp", "             evaluate (\$count=\$count+1)");
	print2file("gseq.inp", "           end loop read");
	print2file("gseq.inp", "           sequence \@\@&prot_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=\"T^\" + encode(\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     if ( &exist_renumber_\$counter = true ) then");
	print2file("gseq.inp", "         if ( &BLANK%renumber_\$counter = false ) then");
	print2file("gseq.inp", "           evaluate (\$segid=\"T^\" + encode(\$counter))");
	print2file("gseq.inp", "           do ( resid = adjustl(format(\"I4\",decode(resid) + &renumber_\$counter - 1))) ");
	print2file("gseq.inp", "              ( (attr refx=9999) and segid \$segid )");
	print2file("gseq.inp", "         end if");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop prot");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop nseg");
	print2file("gseq.inp", "   if ( &exist_prot_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%prot_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       evaluate (\$segtmp=\"T^\" + encode(\$counter))");
	print2file("gseq.inp", "       do (segid=capitalize(&prot_segid_\$counter)) (segid \$segtmp)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop nseg");
	print2file("gseq.inp", " evaluate (\$ssc=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop ssbr");
	print2file("gseq.inp", "   if ( &exist_ss_use_\$ssc = true ) then");
	print2file("gseq.inp", "     if ( &ss_use_\$ssc = true ) then");
	print2file("gseq.inp", "       evaluate (\$segidtmp1=capitalize(&ss_i_segid_\$ssc))");
	print2file("gseq.inp", "       evaluate (\$segidtmp2=capitalize(&ss_j_segid_\$ssc))");
	print2file("gseq.inp", "       patch disu");
	print2file("gseq.inp", "         reference=1=(segid \$QUOTE%segidtmp1 and resid &ss_i_resid_\$ssc)");
	print2file("gseq.inp", "         reference=2=(segid \$QUOTE%segidtmp2 and resid &ss_j_resid_\$ssc)");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$ssc=\$ssc+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop ssbr");
	print2file("gseq.inp", " {* any special protein patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop nucl");
	print2file("gseq.inp", "   if ( &exist_nucl_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%nucl_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           evaluate (\$count=1)");
	print2file("gseq.inp", "           evaluate (\$done2=false)");
	print2file("gseq.inp", "           while ( \$done2 = false ) loop read");
	print2file("gseq.inp", "             if ( &exist_link_infile_\$count = true ) then");
	print2file("gseq.inp", "               if ( &BLANK%link_infile_\$count = false ) then");
	print2file("gseq.inp", "                  \@\@&link_infile_\$count");
	print2file("gseq.inp", "               end if");
	print2file("gseq.inp", "             else");
	print2file("gseq.inp", "               evaluate (\$done2=true)");
	print2file("gseq.inp", "             end if");
	print2file("gseq.inp", "             evaluate (\$count=\$count+1)");
	print2file("gseq.inp", "           end loop read");
	print2file("gseq.inp", "           sequence \@\@&nucl_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&nucl_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop nucl");
	print2file("gseq.inp", " {* patch rna sugars to dna here if needed - select the residues *}");
	print2file("gseq.inp", " {===>} ");
	print2file("gseq.inp", " for \$resid in () loop dna");
	print2file("gseq.inp", "   patch deox reference=nil=(resid \$resid) end");
	print2file("gseq.inp", " end loop dna");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " {* any special nucleic acid patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop carbo");
	print2file("gseq.inp", "   if ( &exist_carbo_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%carbo_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&carbo_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&carbo_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop carbo");
	print2file("gseq.inp", " evaluate (\$carc=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop cabr");
	print2file("gseq.inp", "   if ( &exist_carbo_use_\$carc = true ) then");
	print2file("gseq.inp", "     if ( &carbo_use_\$carc = true ) then");
	print2file("gseq.inp", "       evaluate (\$segidtmp1=capitalize(&carbo_i_segid_\$carc))");
	print2file("gseq.inp", "       evaluate (\$segidtmp2=capitalize(&carbo_j_segid_\$carc))");
	print2file("gseq.inp", "       patch &carbo_patch_\$carc");
	print2file("gseq.inp", "         reference=-=(segid \$QUOTE%segidtmp1 and");
	print2file("gseq.inp", "                      resid &carbo_i_resid_\$carc)");
	print2file("gseq.inp", "         reference=+=(segid \$QUOTE%segidtmp2 and");
	print2file("gseq.inp", "                      resid &carbo_j_resid_\$carc)");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$carc=\$carc+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop cabr");
	print2file("gseq.inp", " {* any special carbohydrate patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop prost");
	print2file("gseq.inp", "   if ( &exist_prost_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%prost_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&prost_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&prost_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop prost");
	print2file("gseq.inp", " {* any special prosthetic group patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop liga");
	print2file("gseq.inp", "   if ( &exist_lig_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%lig_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&lig_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&lig_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop liga");
	print2file("gseq.inp", " {* any special ligand patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop ion");
	print2file("gseq.inp", "   if ( &exist_ion_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%ion_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&ion_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&ion_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop ion");
	print2file("gseq.inp", " {* any special ion patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " evaluate (\$counter=1)");
	print2file("gseq.inp", " evaluate (\$done=false)");
	print2file("gseq.inp", " while ( \$done = false ) loop water");
	print2file("gseq.inp", "   if ( &exist_water_sequence_infile_\$counter = true ) then");
	print2file("gseq.inp", "     if ( &BLANK%water_sequence_infile_\$counter = false ) then");
	print2file("gseq.inp", "       do (refx=0) (all)");
	print2file("gseq.inp", "       segment");
	print2file("gseq.inp", "         chain");
	print2file("gseq.inp", "           sequence \@\@&water_sequence_infile_\$counter end");
	print2file("gseq.inp", "         end");
	print2file("gseq.inp", "       end");
	print2file("gseq.inp", "       do (segid=capitalize(&water_segid_\$counter)) (attr refx=9999)");
	print2file("gseq.inp", "     end if");
	print2file("gseq.inp", "     evaluate (\$counter=\$counter+1)");
	print2file("gseq.inp", "   else");
	print2file("gseq.inp", "     evaluate (\$done=true)");
	print2file("gseq.inp", "   end if");
	print2file("gseq.inp", " end loop water");
	print2file("gseq.inp", " {* any special water patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " {* any final patches can be applied here *}");
	print2file("gseq.inp", " {===>}");
	print2file("gseq.inp", " {<===}");
	print2file("gseq.inp", " if (&hydrogen_flag=false) then");
	#	print2file("gseq.inp", "   delete selection=( hydrogen ) end");
	print2file("gseq.inp", " end if");
	print2file("gseq.inp", " if (&set_bfactor=true) then");
	print2file("gseq.inp", "   do (b=&bfactor) ( all )");
	print2file("gseq.inp", " end if");
	print2file("gseq.inp", " if (&set_occupancy=true) then");
	print2file("gseq.inp", "   do (q=&occupancy) ( all )");
	print2file("gseq.inp", " end if");
	print2file("gseq.inp", " write structure output=&structure_outfile end");
	print2file("gseq.inp", " stop");
}

sub write_cns_generate_extended_file{
	print2file("extn.inp", "{+ file: generate_extended.inp +}");
	print2file("extn.inp", "{+ directory: nmr_calc +}");
	print2file("extn.inp", "{+ description: Generates an extended strand with ideal geometry ");
	print2file("extn.inp", "                for each connected polymer.  ");
	print2file("extn.inp", "                The molecular structure file must not contain any ");
	print2file("extn.inp", "                closed loops except disulfide bonds which are automatically");
	print2file("extn.inp", "                excluded from the generation of the strand conformation.  +}");
	print2file("extn.inp", "{+ authors: Axel T. Brunger +}");
	print2file("extn.inp", "{+ copyright: Yale University +}");
	print2file("extn.inp", "{- begin block parameter definition -} define(");
	print2file("extn.inp", "{======================= molecular structure =========================}");
	print2file("extn.inp", "{* structure file(s) *}");
	print2file("extn.inp", "{===>} structure_file=\"extended.mtf\";");
	print2file("extn.inp", "{* parameter file(s) *}");
	# CNS_TOPPAR:protein-allhdg5-4.param
	print2file("extn.inp", "{===>} par_1=\"CNS_TOPPAR:protein.param\";");
	print2file("extn.inp", "{===>} par_2=\"\";");
	print2file("extn.inp", "{===>} par_3=\"\";");
	print2file("extn.inp", "{===>} par_4=\"\";");
	print2file("extn.inp", "{===>} par_5=\"\";");
	print2file("extn.inp", "{======================= input parameters ============================}");
	print2file("extn.inp", "{* maximum number of trials to generate an acceptable structure *}");
	print2file("extn.inp", "{===>} max_trial=10;");
	print2file("extn.inp", "{=========================== output files ============================}");
	print2file("extn.inp", "{* output coordinates *}");
	print2file("extn.inp", "{===>} output_coor=\"extended.pdb\";");
	print2file("extn.inp", "                                  ");
	print2file("extn.inp", "{===========================================================================}");
	print2file("extn.inp", "{        things below this line do not normally need to be changed          }");
	print2file("extn.inp", "{===========================================================================}");
	print2file("extn.inp", " ) {- end block parameter definition -}");
	print2file("extn.inp", " checkversion 1.3");
	print2file("extn.inp", " evaluate (\$log_level=quiet)");
	print2file("extn.inp", " structure \@&structure_file end");
	print2file("extn.inp", " parameter");
	print2file("extn.inp", "   if (&par_1 # \" \") then");
	print2file("extn.inp", "      \@\@&par_1");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   if (&par_2 # \" \") then");
	print2file("extn.inp", "      \@\@&par_2");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   if (&par_3 # \" \") then");
	print2file("extn.inp", "      \@\@&par_3");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   if (&par_4 # \" \") then");
	print2file("extn.inp", "      \@\@&par_4");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   if (&par_5 # \" \") then");
	print2file("extn.inp", "      \@\@&par_5");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", " end");
	print2file("extn.inp", "{ Set force constants for S-S bond lengths and angles to zero  }");
	print2file("extn.inp", "parameter");
	print2file("extn.inp", "   bonds ( name SG ) ( name SG ) 0. 1. ");
	print2file("extn.inp", "end");
	print2file("extn.inp", "igroup interaction=(all) (all) end");
	print2file("extn.inp", "ident (x) ( all )");
	print2file("extn.inp", "do (x=x/5.) ( all )");
	print2file("extn.inp", "do (y=random(0.5) ) ( all )");
	print2file("extn.inp", "do (z=random(0.5) ) ( all )");
	print2file("extn.inp", "flags exclude * include bond angle impr dihe vdw end");
	print2file("extn.inp", "parameter");
	print2file("extn.inp", "   nbonds");
	print2file("extn.inp", "      rcon=50. nbxmod=-3 repel=0.8 cutnb=6. ");
	print2file("extn.inp", "      rexp=2 irexp=2 inhibit=0.0 wmin=0.1 tolerance=0.5");
	print2file("extn.inp", "   end");
	print2file("extn.inp", "end");
	print2file("extn.inp", "evaluate (\$count=1) ");
	print2file("extn.inp", "while (\$count < 10 ) loop l1");
	print2file("extn.inp", "   do (x=x+gauss(0.1)) ( all ) ");
	print2file("extn.inp", "   do (y=y+gauss(0.1)) ( all ) ");
	print2file("extn.inp", "   do (z=z+gauss(0.1)) ( all ) ");
	print2file("extn.inp", "   minimize lbfgs nstep=200 nprint=10 end");
	print2file("extn.inp", "   evaluate (\$count=\$count+1)");
	print2file("extn.inp", "end loop l1");
	print2file("extn.inp", "evaluate (\$accept=false) ");
	print2file("extn.inp", "evaluate (\$trial=1) ");
	print2file("extn.inp", "while (\$accept=false) loop accp");
	print2file("extn.inp", "   for \$1 in id ( tag ) loop resi");
	print2file("extn.inp", "      igroup ");
	print2file("extn.inp", "         interaction=( byresidue (id \$1 ) and not name SG ) ");
	print2file("extn.inp", "                     ( not name SG ) ");
	print2file("extn.inp", "      end");
	print2file("extn.inp", "      evaluate (\$accept=true) ");
	print2file("extn.inp", "      print thres=0.1 bonds");
	print2file("extn.inp", "      if (\$violations > 0) then");
	print2file("extn.inp", "         evaluate (\$accept=false) ");
	print2file("extn.inp", "      end if");
	print2file("extn.inp", "      print thres=10. angles ");
	print2file("extn.inp", "      evaluate (\$angles=\$result)");
	print2file("extn.inp", "      if (\$violations > 0) then");
	print2file("extn.inp", "         evaluate (\$accept=false) ");
	print2file("extn.inp", "      end if");
	print2file("extn.inp", "      print thres=10. improper");
	print2file("extn.inp", "      if (\$violations > 0) then");
	print2file("extn.inp", "         evaluate (\$accept=false) ");
	print2file("extn.inp", "      end if");
	print2file("extn.inp", "      if (\$accept=false) then");
	print2file("extn.inp", "         do (x=x+gauss(0.3)) ( byresidue (id \$1 ) ) ");
	print2file("extn.inp", "         do (y=y+gauss(0.3)) ( byresidue (id \$1 ) ) ");
	print2file("extn.inp", "         do (z=z+gauss(0.3)) ( byresidue (id \$1 ) ) ");
	print2file("extn.inp", "      end if");
	print2file("extn.inp", "   end loop resi");
	print2file("extn.inp", "   igroup interaction=( all ) ( all ) end");
	print2file("extn.inp", "   parameter");
	print2file("extn.inp", "      nbonds");
	print2file("extn.inp", "         rcon=50. nbxmod=-3 repel=3. cutnb=10. ");
	print2file("extn.inp", "      end");
	print2file("extn.inp", "   end");
	print2file("extn.inp", "   flags exclude angle improper end");
	print2file("extn.inp", "   ");
	print2file("extn.inp", "   minimize lbfgs nstep=200 nprint=10 end");
	print2file("extn.inp", "   parameter");
	print2file("extn.inp", "      nbonds");
	print2file("extn.inp", "         rcon=50. nbxmod=-3 repel=0.8 cutnb=6. ");
	print2file("extn.inp", "      end");
	print2file("extn.inp", "   end");
	print2file("extn.inp", "   flags include angle improper end");
	print2file("extn.inp", "   ");
	print2file("extn.inp", "   evaluate (\$count=1) ");
	print2file("extn.inp", "   while (\$count < 5 ) loop l2");
	print2file("extn.inp", "      do (x=x+gauss(0.05)) ( all ) ");
	print2file("extn.inp", "      do (y=y+gauss(0.05)) ( all ) ");
	print2file("extn.inp", "      do (z=z+gauss(0.05)) ( all ) ");
	print2file("extn.inp", "      minimize lbfgs nstep=200 nprint=10 end");
	print2file("extn.inp", "      evaluate (\$count=\$count+1)");
	print2file("extn.inp", "   end loop l2");
	print2file("extn.inp", "   ");
	print2file("extn.inp", "   parameter");
	print2file("extn.inp", "      nbonds");
	print2file("extn.inp", "         rcon=50. nbxmod=3 repel=0.8 cutnb=6. ");
	print2file("extn.inp", "      end");
	print2file("extn.inp", "   end");
	print2file("extn.inp", "   ");
	print2file("extn.inp", "   minimize lbfgs nstep=300 nprint=10 end   ");
	print2file("extn.inp", "   minimize lbfgs nstep=300 nprint=10 end");
	print2file("extn.inp", "   igroup interaction=( not name SG ) ( not name SG ) end");
	print2file("extn.inp", "   energy end");
	print2file("extn.inp", "   evaluate (\$accept=true) ");
	print2file("extn.inp", "   print thres=0.05 bonds");
	print2file("extn.inp", "   evaluate (\$bonds=\$result)");
	print2file("extn.inp", "   if (\$violations > 0) then");
	print2file("extn.inp", "      evaluate (\$accept=false) ");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   print thres=10. angles ");
	print2file("extn.inp", "   evaluate (\$angles=\$result)");
	print2file("extn.inp", "   if (\$violations > 0) then");
	print2file("extn.inp", "      evaluate (\$accept=false) ");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   print thres=10. improper");
	print2file("extn.inp", "   evaluate (\$impr=\$result)");
	print2file("extn.inp", "   if (\$violations > 0) then");
	print2file("extn.inp", "      evaluate (\$accept=false) ");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "   print thres=180. dihedral ");
	print2file("extn.inp", "   evaluate (\$dihe=\$result)");
	print2file("extn.inp", "   evaluate (\$trial=\$trial + 1) ");
	print2file("extn.inp", "   if (\$trial > &max_trial ) then");
	print2file("extn.inp", "      exit loop accp");
	print2file("extn.inp", "   end if");
	print2file("extn.inp", "end loop accp");
	print2file("extn.inp", "remarks extended strand(s) generation");
	print2file("extn.inp", "remarks input molecular structure file=&structure_file ");
	print2file("extn.inp", "remarks final rms deviations (excluding disulfide bonds): ");
	print2file("extn.inp", "remarks    bonds=	 \$bonds[F8.4] A  ");
	print2file("extn.inp", "remarks    angles=	 \$angles[F8.4] degrees");
	print2file("extn.inp", "remarks    impropers= \$impr[F8.4] degrees");
	print2file("extn.inp", "remarks    dihedrals= \$dihe[F8.4] degrees (not used in some parameter sets!)");
	print2file("extn.inp", "remarks final van der Waals (repel) energy=\$vdw kcal/mole");
	print2file("extn.inp", "write coordinates output=&output_coor format=PDBO end  ");
	print2file("extn.inp", "stop");
}

sub print_usage{
	my $param_info = <<EOF;
--------------------------------------------------------------------------------
CONFOLD v1.0 - build protein 3D models using contacts and secondary structures
--------------------------------------------------------------------------------

PARAMETER          DEFAULT  DESCRIPTION 
o        : string : *     : Output directory (mandatory)
seq      : string : *     : Input FASTA file (mandatory)
pair     : string : *     : Beta-pairing file
rr       : string : *     : Input contacts RR file
ss       : string : *     : Input FASTA format secondary structure file;
                            Either 'rr' or 'ss' parameter must be provided
selectrr : string : all   : Number of contacts to use from the input rr file;
                            Options : 'all', '0.1L', '0.2L', ..., '2.0L', etc.
rrtype   : string : cb    : 'ca' for Cα-Cα contacts and 'cb' for Cβ-Cβ contacts
lambda   : float  : 0.4   : Range for secondary structure restraints to define; 
                            '0.4' for predicted contacts
                            '1.0' for reconstructions 
stage2   : int    : 0     : Run stage 2 for contact filtering and sheet-detection;
                            '0' - don't run stage2
                            '1' - stage2 with both contact filtering and sheet-detection
                            '2' - stage2 with contact filtering only
                            '3' - stage2 with sheet-detection only
contwt   : float  : 10    : Weight for contact restraints
sswt     : float  : 5     : Weight of secondary structure (SS) restraints 
                            relative to contact restraints; 
                            '1' implies equal weight for contacts and SS restraints
mcount   : int    : 20    : Number of models to generate in each stage;
                            Leave it to 20 unless CPU resource is very limited;
                            Must be between 5 and 50;
rep2     : float  : 0.85  : Final van der Waals repel radius in dgsa.inp script;
                            '0.8'  - for reconstructions using true contacts
                            '0.85' - for predicted contacts
pthres   : float  : 7.0   : Distance threshold for detecting beta pairs

REFERENCE:
"CONFOLD: Residue-Residue Contact-guided ab initio Protein Folding",
Proteins: Structure, Function, and Bioinformatics, 2015.
B. Adhikari, D. Bhattacharya, R. Cao, J. Cheng. 
EOF
	print $param_info;
	print "\n";
	exit 1;
}
