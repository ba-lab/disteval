# How to upgrade CONFOLD to accept distances?

## Install CNS Suite
#### 1. Download CNS suite   
Provide your academic profile related information at http://cns-online.org/cns_request/. An email with (a) link to download, (b) login, and (c) password will be sent to you. Follow the link, possibly http://cns-online.org/download/, and download CNS suite "cns_solve_1.3_all_intel-mac_linux.tar.gz".
#### 2. Unzip  
```bash
tar xzvf cns_solve_1.3_all_intel-mac_linux.tar.gz
```
#### 3. Change directory to cns_solve  
```bash
cd cns_solve_1.3
```
#### 4. Unhide the file '.cns_solve_env_sh'  
```bash
mv .cns_solve_env_sh cns_solve_env.sh
```
#### 5. Edit path in 'cns_solve_env.sh'  
Replace '_CNSsolve_location_' with CNS installation directory. For instance, if your CNS installation path is '/home/user/programs/cns_solve_1.3' replace '_CNSsolve_location_' with this path
```bash
vim cns_solve_env.sh
```

#### 6. Increase the value for ‘nrestraints’ (maximum number of restraints it can take)
Change the code at line 60 of the module `cns_solve_1.3/modules/nmr/readdata`
```bash
vim cns_solve_1.3/modules/nmr/readdata # change 20000 to 200000 (by adding a zero)
```

#### 7. Test CNS installation  
```bash
source cns_solve_env.sh
cd test 
../bin/run_tests -tidy noe*.inp
```
#### 8. Change directory  
```bash
cd ../../
```

## Configure and test 'distfold.pl'  
#### 1. Update paths for the variables '$cns_suite' and '$program_dssp'  
```bash
vim distfold.pl
```
```perl
my $program_dssp   = "$DIR_BASE/dssp-2.0.4-linux-amd64";
my $cns_suite      = "/home/badri/DISTFOLD/cns_solve_1.3";
```
#### 2. Test  
##### Example with distances as input (1a3aA)  
a. Build models  
```bash
cd test
rm -r output-1a3aA
perl ../distfold.pl -rr ./1a3aA.dist.rr -ss ./1a3aA.ss_sa -o ./output-1a3aA -mcount 20 -selectrr 1.0L
```
b. Evaluate the models
```bash
perl ./eval-using-tmscore.pl 1guu.pdb output-1guu/ all header
```
##### Example with contacts as input (1guu)  

```bash
perl ../confold.pl -seq ./1guu.fasta -rr ./1guu.rr -o ./output-1guu -mcount 20 -selectrr all
```

```perl
# Line #282 throws errors sometimes
system_cmd("./job.sh", "job.log");
# Replace with:
`./job.sh > job.log`;
```

Comment the following line in line 109
```
#system_cmd("cp $dir_out/stage1/$id.ss     ./");
```

```perl
sub rr2tbl{
	my $file_rr  = shift;
	my $file_tbl = shift;
	my $rrtype   = shift;
	confess ":(" if not -f $file_rr;
	confess ":(" if not $file_tbl;
	confess ":(" if not ($rrtype eq "ca" or $rrtype eq "cb");
	my %r1a1r2a2 = rr2r1a1r2a2($file_rr, $rrtype);
	my %lowerbound = rr2contacts_hash($file_rr, 100000, "lowerbound");
	my %upperbound = rr2contacts_hash($file_rr, 100000);
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
```

```perl
sub rr2contacts_hash{
	my $file_rr = shift;
	my $count = shift;
	my $hashvalue = shift;
	my $seq_sep = 2;
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
```

a. Build models
```bash
rm -r output-1guu
perl ../confold.pl -rr ./1guu.rr -o ./output-1guu -mcount 20 -selectrr all
```

With secondary structures:
```bash
perl ../distfold.pl -rr ./1guu.rr -ss ./1guu.ss -o ./output-1guu -mcount 20 -selectrr 1.0L
```

b. Evaluate the models
```bash
perl ./eval-using-tmscore.pl 1a3aA.pdb output-1a3aA all header
```
