# DISTEVAL

## Prerequisites (for distance and contact evaluation)
1. Python3
1. Numpy

## Prerequisites (for building 3D models)
1. Install csh
   ```
   sudo apt install csh
   ```
1. Download 'dssp-2.0.4-linux-amd64' from https://osf.io/qydjv/
   ```
   chmod +x dssp-2.0.4-linux-amd64
   ```
1. Download TM-score from https://zhanglab.ccmb.med.umich.edu/TM-score/
    ```
    gunzip TM-score.gz
    chmod +x TM-score
    ```
1. CONFOLD
    - Follow instructions [here](CONFOLD-CHANGES.md) to download, install, and modify CONFOLD

## Test cases (without model building)
Download all the files in the test folder, for example:
   ```
   wget https://raw.githubusercontent.com/ba-lab/disteval/main/test/1guuA.contact.rr
   ```

1. Evaluate a predicted RR contacts file
   ```
   python3 disteval.py 
   ```
   Expected output:
   ```
   
   ```
1. Evaluate a predicted distance map
   ```
   python3 ../disteval-using-confold.py -n 1guuA.pdb -d 1guuA.predicted.npy
   ```
   Expected output:
   ```
   ...
   sep: 24 xL: Top-L/5 {'mae': 1.8154, 'mse': 4.6469, 'rmse': 2.1557, 'count': 10}
   sep: 24 xL: Top---L {'mae': 2.1541, 'mse': 8.1816, 'rmse': 2.8603, 'count': 50}
   ...
   sep: 24 xL: Top-L/5 0.5
   sep: 24 xL: Top---L 0.38462
   ...
   ```
1. Evaluate trRosetta prediction

## Test cases (with model building)
1. predicted contacts + SS
1. predicted distances with high seq sep + SS
1. native dmap
1. very high accuracy native dmap reconstruction
1. predicted 2D numpy
1. trRosetta
1. CASP distances RR


## Contact  
Badri Adhikari  
adhikarib@umsl.edu  
University of Missouri-St. Louis  
