'''
Author: Badri Adhikari, University of Missouri-St. Louis, 10-25-2020
File: Reconstruct 3D models with CONFOLD using a predicted distance map and evaluate using TM-score
Options for building models:
    a) a plain 2D numpy distance map
    b) trRosetta predicted distance map
    c) CASP RR file
    d) use distances from the PDB itself
'''

import argparse
import sys
import numpy as np
import re
import os
from math import sqrt
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics import precision_score

if sys.version_info < (3,0,0):
    print('Python 3 required!!!')
    sys.exit(1)

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', type=str, required = False, dest = 'native',    help="true PDB file")
    parser.add_argument('-f', type=str, required = False, dest = 'fasta',     help="FASTA file of the input distance map (for building models)")
    parser.add_argument('-d', type=str, required = False, dest = 'dmap',      help="Predicted distance map as a 2D numpy array")
    parser.add_argument('-t', type=int, required = True,  dest = 'threshold', help="Distance cutoff threshold")
    parser.add_argument('-c', type=str, required = False, dest = 'inputrr',   help="CASP RR file as input (all input rows are used)") 
    parser.add_argument('-r', type=str, required = False, dest = 'trrosetta', help="trRosetta prediction")
    parser.add_argument('-o', type=str, required = True,  dest = 'jobdir',    help="Output directory")
    parser.add_argument('-s', type=str, required = False, dest = 'ss',        help="3-class (H/C/E) Secondary structure file in FASTA format (for building models)")
    parser.add_argument('-m', type=int, required = True,  dest = 'minsep',    help="Minimum sequence separation (24 for long-range & 12 for medium+long-range)")
    parser.add_argument('-p',           required = False, dest = 'truedmap',   action = 'store_true',  help='Use true distances from the PDB as input')
    parser.add_argument('-b',           required = False, dest = 'modeling3d', action = 'store_true', help='Build 3D models using CNS')
    args = parser.parse_args()
    return args

args = get_args()
print(args)

TM = '/home/badri/disteval-using-confold/TMscore'
CONFOLD = '/home/badri/disteval-using-confold/confold.pl'

native     = None
file_fasta = None
dmap       = None
trrosetta  = None
inputrr    = None
ss         = None
truedmap   = False
modeling3d = False

threshold  = args.threshold
job_dir    = args.jobdir
minsep     = args.minsep

if args.native is not None: native = os.path.abspath(args.native)
if args.fasta is not None: file_fasta = os.path.abspath(args.fasta)
if args.dmap is not None: dmap = os.path.abspath(args.dmap)
if args.trrosetta is not None: trrosetta  = os.path.abspath(args.trrosetta)
if args.inputrr is not None: inputrr = os.path.abspath(args.inputrr)
if args.ss is not None: ss = os.path.abspath(args.ss)
if args.truedmap is True: truedmap = True
if args.modeling3d is True: modeling3d = True

valid_amino_acids = {
    'LLP': 'K', 'TPO': 'T', 'CSS': 'C', 'OCS': 'C', 'CSO': 'C', 'PCA': 'E', 'KCX': 'K', \
    'CME': 'C', 'MLY': 'K', 'SEP': 'S', 'CSX': 'C', 'CSD': 'C', 'MSE': 'M', \
    'ALA': 'A', 'ASN': 'N', 'CYS': 'C', 'GLN': 'Q', 'HIS': 'H', 'LEU': 'L', \
    'MET': 'M', 'MHO': 'M', 'PRO': 'P', 'THR': 'T', 'TYR': 'Y', 'ARG': 'R', 'ASP': 'D', \
    'GLU': 'E', 'GLY': 'G', 'ILE': 'I', 'LYS': 'K', 'PHE': 'F', 'SER': 'S', \
    'TRP': 'W', 'VAL': 'V', 'SEC': 'U'
    }

def dmin_dmax_for_d(d):
    dev = 0 # 0.125 * d
    return (d - dev / 2.0, d + dev / 2.0)

def check_pdb_valid_row(valid_amino_acids, l):
    if (get_pdb_rname(l) in valid_amino_acids.keys()) and (l.startswith('ATOM') or l.startswith('HETA')):
        return True
    return False

def get_pdb_atom_name(l):
    return l[12: 16].strip()

def get_pdb_rnum(l):
    return int(l[22: 27].strip())

def get_pdb_rname(l):
    return l[17: 20].strip()

def get_pdb_xyz_cb(lines):
    xyz = {}
    for l in lines:
        if get_pdb_atom_name(l) == 'CB':
            xyz[get_pdb_rnum(l)] = (float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip()))
    for l in lines:
        if (get_pdb_rnum(l) not in xyz) and get_pdb_atom_name(l) == 'CA':
            xyz[get_pdb_rnum(l)] = (float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip()))
    return xyz

def get_pdb_xyz_ca(lines):
    xyz = {}
    for l in lines:
        if get_pdb_atom_name(l) == 'CA':
            xyz[get_pdb_rnum(l)] = (float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip()))
    return xyz

def pdb2dmap(pdbfile):
    f = open(pdbfile, mode = 'r')
    flines = f.read()
    f.close()
    lines = flines.splitlines()
    templines = flines.splitlines()
    for l in templines:
        if not l.startswith('ATOM'): lines.remove(l)
    # We have filtered out all non ATOMs at this point
    rnum_rnames = {}
    for l in lines:
        atom = get_pdb_atom_name(l)
        if atom != 'CA': continue
        if not get_pdb_rname(l) in valid_amino_acids.keys():
            print ('' + get_pdb_rname(l) + ' is unknown amino acid in ' + l)
            return
        rnum_rnames[int(get_pdb_rnum(l))] = valid_amino_acids[get_pdb_rname(l)]
    seq = ""
    for i in range(max(rnum_rnames.keys())):
        if i+1 not in rnum_rnames:
            #print (rnum_rnames)
            #print ('Warning! residue not defined for rnum = ' + str(i+1))
            seq += '-'
        else:
            seq += rnum_rnames[i+1]
    L = len(seq)
    xyz_cb = get_pdb_xyz_cb(lines)
    total_valid_residues = len(xyz_cb)
    if len(xyz_cb) != L:
        print(rnum_rnames)
        for i in range(L):
            if i+1 not in xyz_cb: print('XYZ not defined for ' + str(i+1))
        print ('Warning! Something went wrong - len of cbxyz != seqlen!! ' + str(len(xyz_cb)) + ' ' +  str(L))
    cb_map = np.full((L, L), np.nan)
    for r1 in sorted(xyz_cb):
        (a, b, c) = xyz_cb[r1]
        for r2 in sorted(xyz_cb):
            (p, q, r) = xyz_cb[r2]
            cb_map[r1 - 1, r2 - 1] = sqrt((a-p)**2+(b-q)**2+(c-r)**2)
    return (total_valid_residues, cb_map, rnum_rnames)

def seqfasta(fasta):
    f = open(fasta)
    lines = f.readlines()
    f.close()
    seq = ''
    for l in lines:
        if l.startswith('>'): continue
        seq += l.strip()
    return seq

def trrosetta_probindex2dist(index):
    d = 1.75
    for k in range(1, 37):
        d += 0.5
        if index == k:
            return d
    return d

def trrosetta2maps(trrosetta):
    x = np.load(trrosetta)
    a = x['dist']
    if len(a[0, 0, :]) != 37:
        print('ERROR! This does not look like a trRosetta prediction')
        return
    D = np.full((len(a), len(a)), 21.0)
    for i in range(len(a)):
        for j in range(len(a)):
            maxprob_value = 0.0
            for k in range(37):
                if maxprob_value < a[i, j, k]:
                    maxprob_value = a[i, j, k]
                    D[i, j] = trrosetta_probindex2dist(k)
    C = np.full((len(a), len(a)), 0.0)
    for i in range(len(a)):
        for j in range(i, len(a)):
            for k in range(1, 13):
                C[i, j] += a[i, j, k]
    return (D, C)

def calc_dist_errors(P, Y, L, dist_thres = 8.0, min_sep = 24, top_l_by_x = 5, pred_limit = 20):
    if Y is None:
        print('ERROR! Y is None!')
        return
    if P is None:
        print('ERROR! P is None!')
        return
    if np.isnan(Y).all():
        print('ERROR! Y is all NaNs!')
        return
    if np.isnan(P).all():
        print('ERROR! P is all NaNs!')
        return
    errors = {}
    errors['mae'] = np.nan
    errors['mse'] = np.nan
    errors['rmse'] = np.nan
    pred_dict = {}
    true_dict = {}
    for p in range(len(Y)):
        for q in range(len(Y)):
            if q - p < min_sep: continue
            if np.isnan(P[p, q]): continue
            if np.isnan(Y[p, q]): continue
            if Y[p, q] >= dist_thres: continue
            if P[p, q] >= pred_limit: continue
            if np.isnan(Y[p, q]): continue
            pred_dict[(p, q)] = P[p, q]
            true_dict[(p, q)] = Y[p, q]
    xl = int (L / top_l_by_x)
    pred_list = []
    true_list = []
    for pair in sorted(pred_dict.items(), key=lambda x: x[1]):
        if pair[0] not in true_dict: continue
        pred_list.append(pred_dict[pair[0]])
        true_list.append(true_dict[pair[0]])
        xl -= 1
        if xl == 0: break
    if len(pred_list) > 0:
        errors['mae'] = round(mean_absolute_error(true_list, pred_list), 4)
        errors['mse'] = round(mean_squared_error(true_list, pred_list), 4)
        errors['rmse'] = round(sqrt(errors['mse']), 4)
        errors['count'] = len(pred_list)
    return errors

def calc_dist_errors_various_xl(P, Y, L, separation = [12, 24]):
    all_metrics = {}
    dist_thres = ['1000'] # ['08', '12']
    topxl = {5:'l5', 1:'l1', 0.001:'nc'}
    for dt in dist_thres:
        for sep in separation:
            for xl in topxl.keys():
                results = calc_dist_errors(P = P, Y = Y, L = L, dist_thres = int(dt), min_sep = int(sep), top_l_by_x = xl)
                if len(dist_thres) > 1:
                    all_metrics["dthres: " + dt + " sep: " + str(sep) + " xL: " + topxl[xl]] = results
                else:
                    all_metrics["sep: " + str(sep) + " xL: " + topxl[xl]] = results
    return all_metrics

def calc_contact_errors_various_xl(CPRED, CTRUE, separation = [12, 24]):
    all_metrics = {}
    topxl = {0.2:'l5', 1:'l1', 10:'nc'}
    for sep in separation:
        for xl in topxl.keys():
            results = calculate_contact_precision(CPRED = CPRED, CTRUE = CTRUE, minsep = sep, topxl = xl)
            all_metrics["sep: " + str(sep) + " xL: " + topxl[xl]] = results
    return all_metrics

def rr2dmap(filerr):
    f = open(filerr)
    lines = f.readlines()
    f.close()
    # Detect RMODE (if present)
    mode = 0
    for l in lines:
        if l.startswith('RMODE 1') or l.startswith('RMODE  1'):
            mode = 1
            break
        if l.startswith('RMODE 2') or l.startswith('RMODE  2'):
            mode = 2
            break
    pindex = 4
    if mode == 1 or mode == 2:
        pindex = 2
    # Detect target ID (if present)
    tgt = ''
    for l in lines:
        if l.startswith('TARGET'):
            cols = l.strip().split()
            tgt = cols[1]
    # Extract sequence
    seq = ''
    for l in lines:
        if l.startswith('PFRMAT '): continue
        if l.startswith('TARGET '): continue
        if l.startswith('AUTHOR '): continue
        if l.startswith('REMARK '): continue
        if l.startswith('METHOD '): continue
        if l.startswith('RMODE '): continue
        if l.startswith('MODEL '): continue
        if l.startswith('END'): continue
        if l[0].isalpha():
            seq += l.strip()
    # Try to download if it is a CASP14 target
    if mode !=0 and len(seq) < 1:
        os.system('wget -O ' + tgt  + '.fasta --content-disposition "http://predictioncenter.org/casp14/target.cgi?target=' + tgt + '&view=sequence\"')
        f = open(tgt + '.fasta')
        f.readline()
        seq = f.readline().strip()
        f.close()
    L = len(seq)
    print(L)
    print(seq)
    # Absent values are NaNs
    C = np.full((L, L), np.nan)
    D = None  
    if mode != 2:
        for l in lines:
            if not l[0].isdigit(): continue
            c = l.split()
            C[int(c[0]) - 1, int(c[1]) - 1] = float(c[pindex])
            #D[int(c[0]) - 1, int(c[1]) - 1] = 4.0 / (float(c[pindex]) + 0.01)
    else:
        # Absent values are NaNs
        D = np.full((L, L), np.nan)
        for l in lines:
            if not l[0].isdigit():
                continue
            c = l.strip().split()
            i = int(c[0]) - 1
            j = int(c[1]) - 1
            max_prob_value = 0.0
            max_prob_index = -1
            if len(c) != 13:
                print('ERROR! Unexpected number of columns in line:', l)
                sys.exit()
            # identify the maximum probability
            for position in range(2, len(c)):
                if float(c[position]) > max_prob_value:
                    max_prob_value = float(c[position])
                    max_prob_index = position
            d = 3.0
            for position in range(2, len(c)):
                d += 2.0
                if max_prob_index == position:
                    break
            D[i,j] = d
        # To Do: The case of CASP14
        for i in range(L):
            for j in range(i, L):
                for k in range(1, 13):
                    if not np.isnan(C[i, j]): C[i, j] += 0
                    if np.isnan(C[i, j]): C[i, j] = 0
        D = np.clip(D, 0.0, 100.0)
    return (D, C)

def calculate_contact_precision(CPRED, CTRUE, minsep, topxl, LPDB = None):
    L = len(CPRED)
    if LPDB is None: LPDB = len(np.where(~np.isnan(np.diagonal(CTRUE)))[0])
    # The number of valid true values must be <= predicted
    num_true = 0
    for j in range(0, L):
        for k in range(j, L):
            if np.isnan(CTRUE[j, k]): continue
            if abs(j - k) < minsep: continue
            if CTRUE[j, k] > 1.0 or CTRUE[j, k] < 0.0: print("WARNING!! True contact at "+str(j)+" "+str(k)+" is "+str(CTRUE[j, k]))
            num_true += 1
    num_pred = 0
    for j in range(0, L):
        for k in range(j, L):
            if np.isnan(CPRED[j, k]): continue
            if abs(j - k) < minsep: continue
            if CPRED[j, k] > 1.0 or CPRED[j, k] < 0.0: print("WARNING!! Predicted probability at "+str(j)+" "+str(k)+" is "+str(CPRED[j, k]))
            num_pred += 1
    if num_true < 1:
        return float('nan')
    # Put predictions in a dictionary so they can be sorted
    p_dict = {}
    for j in range(0, L):
        for k in range(j, L):
            if np.isnan(CTRUE[j, k]): continue
            if np.isnan(CPRED[j, k]): continue
            if abs(j - k) < minsep: continue
            p_dict[(j, k)] = CPRED[j, k]
    # Obtain nc, the total number of contacts in the PDB
    nc_count = 0
    for j in range(0, L):
        for k in range(j, L):
            if np.isnan(CTRUE[j, k]): continue
            if abs(j - k) < minsep: continue
            if CTRUE[j, k] != 1: continue
            nc_count += 1
    if nc_count < 1:
        return float('nan')
    # Obtain top xL predictions
    xl = int(LPDB * topxl)
    if xl > nc_count: xl = nc_count
    pred_list = []
    true_list = []
    for pair in reversed(sorted(p_dict.items(), key=lambda x: x[1])):
        if np.isnan(CTRUE[pair[0][0], pair[0][0]]): continue
        pred_list.append(1) # This is assumed to be a +ve prediction
        true_list.append(CTRUE[pair[0][0], pair[0][1]])
        xl -= 1
        if xl == 0: break
        ####### PUT NEW CONTACT METRICS HERE ######
    return round(precision_score(true_list, pred_list), 5)

def dmap2rr(P, seq, file_rr):
    f = open(file_rr, 'w')
    f.write(seq + '\n')
    for j in range(0, len(P)):
        for k in range(j, len(P)):
            if abs(j - k) < minsep: continue
            if P[j][k] > threshold: continue
            (dmin, dmax) = dmin_dmax_for_d(P[j][k]) 
            f.write("%d %d %0.2f %.2f 1.0\n" %(j+1, k+1, dmin, dmax) )
    f.close()

def evaltm(pred, native):
    os.system(TM + " " + pred + " " + native + " | grep -e RMSD\\ of -e TM-score\\ \\ \\  -e MaxSub-score -e GDT-TS-score -e GDT-HA-score > x.tmp")
    f = open('x.tmp')
    lines = f.readlines()
    f.close()
    os.system('rm x.tmp')
    rmsd = -1
    tmsc = -1
    gdts = -1
    for l in lines:
        if l.startswith('RMSD'): rmsd = float(re.sub("[^\d\.]", "", l)[:5])
    for l in lines:
        if l.startswith('TM-score'): tmsc = float(re.sub("[^\d\.]", "", l)[:5])
    for l in lines:
        if l.startswith('GDT-TS'): gdts = float(re.sub("[^\d\.]", "", l)[:5])
    return (rmsd, tmsc, gdts)

# Prepare the native
l_for_xL = None
ND = None
NC = None
rnum_rnames = None
if native:
    print('')
    print('Load PDB..')
    (l_for_xL, ND, rnum_rnames) = pdb2dmap(native)
    print(ND.shape, l_for_xL)
    NC = np.copy(ND)
    NC [NC < 8.0] = 1
    NC [NC >= 8.0] = 0

# Prepare the predicted distances and contacts
D = None
C = None 
if trrosetta is not None:
    print('Load the input trRosetta prediction..')
    (D, C) = trrosetta2maps(trrosetta)
elif dmap is not None:
    print('Load the input 2D distance map..')
    D = np.load(dmap)
    L = len(D)
    if D.ndim == 3:
        D = D.reshape((L, L))
    C = 4.0 / (D + 0.001)
elif truedmap:
    print('Obtaining a true distance map from the input PDB..')
    D = np.copy(ND)
    C = 4.0 / (D + 0.001)
else:
    (D, C) = rr2dmap(inputrr)
    if np.isnan(C).all(): sys.exit('ERROR! C is all NaNs!')
    if D is not None and np.isnan(D).all(): sys.exit('ERROR! D is all NaNs!')

# Check for NaNs
print('')
if C is not None:
    print('C.shape', C.shape)
    print('Contact nans:', np.count_nonzero(np.isnan(C)), 'of', str(len(C) * len(C)))
if D is not None:
    print('D.shape', D.shape)
    print('Distance nans:', np.count_nonzero(np.isnan(D)), 'of', str(len(D) * len(D)))
if C is None and D is None:
    sys.exit('ERROR!! Could not load contact or distance!')

# Evaluate distances
if D is not None:
    print('')
    print('Evaluating distances..')
    all_metrics = calc_dist_errors_various_xl(P = D, Y = ND, L = l_for_xL) #, separation = [minsep])
    for k in all_metrics:
        print(os.path.basename(native), k, all_metrics[k])

if C is not None:
    print('')
    print('Evaluating contacts..')
    all_metrics = calc_contact_errors_various_xl(CPRED = C, CTRUE = NC) #, separation = [minsep])
    for k in all_metrics:
        print(os.path.basename(native), k, all_metrics[k])

if not modeling3d:
    sys.exit()

file_rr = job_dir + '/x.rr'
if inputrr is None:
    seq = seqfasta(file_fasta)
    dmap2rr(C, seq, file_rr)
else:
    os.system('cp ' + inputrr + ' ' + job_dir + '/x.rr')

f = open(job_dir + '/x.rr')
lines = f.readlines()
f.close()
restraint_count = 0
for l in lines:
    if l[0].isdigit(): restraint_count += 1

print('Restraints:')
os.system('head ' + file_rr)

if restraint_count < 1:
    print('ERROR!! No restraints to pass on to CONFOLD! Exiting..')
    sys.exit(1)

print('')
print('Run CONFOLD')
ssparam = ''
if ss is not None: ssparam = ' -ss ' + ss
os.system(f"perl {CONFOLD} -seq {file_fasta} -rr {file_rr} -o {job_dir} -mcount 20 -selectrr all" + ssparam)

if not os.path.exists(native):
    sys.exit()

print('')
print('Run TM-scores')

os.chdir(job_dir + '/stage1/')
tmscores = {}
for pdb in os.listdir('./'):
    if not pdb.endswith('pdb'): continue
    tmscores[pdb] = evaltm(native, pdb)

for pdb in sorted(tmscores.items(), key=lambda kv: kv[1][1]):
    p = pdb[0]
    (r,t,g) = tmscores[pdb[0]]
    print(f'{t:0.3f} {r:0.2f} {g:0.3f} {p}')

