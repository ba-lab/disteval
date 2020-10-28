# DISTFOLD
Updated version of the CONFOLD program to accept distance constraints.

## Step 1. Install CNS Suite

### 1. Download CNS suite   
Provide your academic profile related information at http://cns-online.org/cns_request/. An email with (a) link to download, (b) login, and (c) password will be sent to you. Follow the link, possibly http://cns-online.org/download/, and download CNS suite "cns_solve_1.3_all_intel-mac_linux.tar.gz".
### 2. Unzip  
```bash
tar xzvf cns_solve_1.3_all_intel-mac_linux.tar.gz
```
### 3. Change directory to cns_solve  
```bash
cd cns_solve_1.3
```
### 4. Unhide the file '.cns_solve_env_sh'  
```bash
mv .cns_solve_env_sh cns_solve_env.sh
```
### 5. Edit path in 'cns_solve_env.sh'  
Replace '_CNSsolve_location_' with CNS installation directory. For instance, if your CNS installation path is '/home/user/programs/cns_solve_1.3' replace '_CNSsolve_location_' with this path
```bash
pwd
vim cns_solve_env.sh
```
### 6. Increase the value for ‘nrestraints’ (maximum number of restraints it can take)
Change the code at line 60 of the module `cns_solve_1.3/modules/nmr/readdata`
```bash
vim modules/nmr/readdata # change 20000 to 200000 (by adding a zero)
```
### 7. Test CNS installation  
```bash
source cns_solve_env.sh
cd test 
../bin/run_tests -tidy noe*.inp
```
### 8. Change directory  
```bash
cd ../../
```

## Step 2. Download `distfold.pl`

```
wget ???????
```

## Step 3. Test `distfold.pl`

```bash
wget https://raw.githubusercontent.com/multicom-toolbox/CONFOLD/master/test/input/1guu.fasta
wget https://raw.githubusercontent.com/multicom-toolbox/CONFOLD/master/test/input/1guu.rr
wget https://raw.githubusercontent.com/multicom-toolbox/CONFOLD/master/test/input/1guu.ss
# This will take one or two minutes
perl ./confold.pl -seq ./1guu.fasta -rr ./1guu.rr -ss ./1guu.ss -o ./output-1guu -mcount 20 -selectrr 1.0L
```
