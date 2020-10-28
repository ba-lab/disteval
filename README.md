# DISTEVAL

## Prerequisites (for distance and contact evaluation)
1. Python3
1. Numpy

## Prerequisites (building 3D models)
1. csh
   ```
   sudo apt install csh
   ```
1. DSSP 2.0.4
   - Download 'dssp-2.0.4-linux-amd64' from https://osf.io/qydjv/
   ```
   chmod +x dssp-2.0.4-linux-amd64
   ```
1. TM-score 
    - Download from https://zhanglab.ccmb.med.umich.edu/TM-score/
    ```
    gunzip TM-score.gz
    chmod +x TM-score
    ```
1. CONFOLD
    - Follow instructions [here](CONFOLD-CHANGES.md) to download, install, and modify CONFOLD

## Test cases
To Do: include the best TM-score for each case:
1. predicted contacts only
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
