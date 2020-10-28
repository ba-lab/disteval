
# DISTEVAL: For evaluating protein inter-residue distances
Protein inter-residue contact and distance prediction is a key intermediate step towards accurate protein structure prediction. Distance prediction comes in two forms: real-valued distances and 'binned' distograms, a more finely grained variant of the binary contact prediction problem. Importantly, the later has been introduced as a new challenge in the 14th Critical Assessment of Techniques for Protein Structure Prediction (2020) experiment. Despite the recent proliferation of methods for predicting distances, few methods exist for evaluating these predictions.  Currently only numerical metrics, which evaluate the entire prediction at once, are used.  These give no insight into the structural details of a prediction and, as such, new methods and tools are needed. We have developed a web-server for evaluating predicted inter-residue distances. Our server, DISTEVAL, accepts predicted contacts, distances, and true structure as optional inputs to generate informative heatmaps, chord diagrams, and 3D models which facilitate visual and qualitative assessment. The server also evaluates predictions using mean absolute error (MAE) and the standard 'contact precision' metric. DISTEVAL will be useful for researchers in the field of protein structure prediction. The visualizations generated complement each other and collectively serve as a powerful tool for both quantitative and qualitative assessments of predicted contacts and distances, even in the absence of a true 3D structure.

**Webserver:** [http://deep.cs.umsl.edu/disteval/](http://deep.cs.umsl.edu/disteval/)

<p align="center">
<img src="disteval.png" alt="DISTEVAL BANNER" width=500/>
</p>

# Distance/contact evaluation using the `disteval.py`

## Download
Download from [https://github.com/ba-lab/disteval/releases](https://github.com/ba-lab/disteval/releases)

## Prerequisites
- [x] Python3
- [x] Numpy
- [x] Scikit-learn

## Test

### Example 0. See help
   ```bash
   python3 ./disteval.py -h
   ```

### Example 1. Evaluate a predicted RR contacts file
   ```bash
   python3 ./disteval.py -n ./test/1guuA.pdb -c ./test/1guuA.contact.rr
   ```
   Expected output:
   ```
   Evaluating contacts..
   min-seq-sep: 12 xL: Top-L/5 {'precision': 1.0, 'count': 9}
   min-seq-sep: 12 xL: Top-L   {'precision': 1.0, 'count': 9}
   min-seq-sep: 12 xL: Top-NC  {'precision': 1.0, 'count': 9}
   min-seq-sep: 24 xL: Top-L/5 {'precision': 1.0, 'count': 1}
   min-seq-sep: 24 xL: Top-L   {'precision': 1.0, 'count': 1}
   min-seq-sep: 24 xL: Top-NC  {'precision': 1.0, 'count': 1}
   ```
### Example 2. Evaluate a predicted distance map
   ```bash
   python3 ./disteval.py -n ./test/1guuA.pdb -d ./test/1guuA.predicted.npy
   ```
   Expected output:
   ```
   Evaluating distances..
   min-seq-sep: 12 xL: Top-L/5 {'mae': 0.9403, 'mse': 1.5143, 'rmse': 1.2306, 'count': 10}
   min-seq-sep: 12 xL: Top-L   {'mae': 1.7522, 'mse': 5.6841, 'rmse': 2.3841, 'count': 50}
   min-seq-sep: 12 xL: Top-NC  {'mae': 1.9263, 'mse': 6.6872, 'rmse': 2.586, 'count': 603}
   min-seq-sep: 24 xL: Top-L/5 {'mae': 1.8154, 'mse': 4.6469, 'rmse': 2.1557, 'count': 10}
   min-seq-sep: 24 xL: Top-L   {'mae': 2.1541, 'mse': 8.1816, 'rmse': 2.8603, 'count': 50}
   min-seq-sep: 24 xL: Top-NC  {'mae': 2.4536, 'mse': 9.6231, 'rmse': 3.1021, 'count': 295}
   Evaluating contacts..
   min-seq-sep: 12 xL: Top-L/5 {'precision': 0.9, 'count': 10}
   min-seq-sep: 12 xL: Top-L   {'precision': 0.6, 'count': 30}
   min-seq-sep: 12 xL: Top-NC  {'precision': 0.6, 'count': 30}
   min-seq-sep: 24 xL: Top-L/5 {'precision': 0.5, 'count': 10}
   min-seq-sep: 24 xL: Top-L   {'precision': 0.38462, 'count': 13}
   min-seq-sep: 24 xL: Top-NC  {'precision': 0.38462, 'count': 13}
   ```
### Example 3. Evaluate trRosetta prediction
   ```bash
   python3 ./disteval.py -n ./test/1guuA.pdb -r ./test/1guuA.npz 
   ```
   Expected output:
   ```
   Evaluating distances..
   min-seq-sep: 12 xL: Top-L/5 {'mae': 0.5485, 'mse': 0.5375, 'rmse': 0.7331, 'count': 10}
   min-seq-sep: 12 xL: Top-L   {'mae': 0.6789, 'mse': 0.7678, 'rmse': 0.8762, 'count': 50}
   min-seq-sep: 12 xL: Top-NC  {'mae': 1.2951, 'mse': 3.8733, 'rmse': 1.9681, 'count': 741}
   min-seq-sep: 24 xL: Top-L/5 {'mae': 0.537, 'mse': 0.4237, 'rmse': 0.6509, 'count': 10}
   min-seq-sep: 24 xL: Top-L   {'mae': 0.6691, 'mse': 0.6725, 'rmse': 0.8201, 'count': 50}
   min-seq-sep: 24 xL: Top-NC  {'mae': 1.2281, 'mse': 3.2863, 'rmse': 1.8128, 'count': 351}

   Evaluating contacts..
   min-seq-sep: 12 xL: Top-L/5 {'precision': 1.0, 'count': 10}
   min-seq-sep: 12 xL: Top-L   {'precision': 0.8, 'count': 30}
   min-seq-sep: 12 xL: Top-NC  {'precision': 0.8, 'count': 30}
   min-seq-sep: 24 xL: Top-L/5 {'precision': 1.0, 'count': 10}
   min-seq-sep: 24 xL: Top-L   {'precision': 0.84615, 'count': 13}
   min-seq-sep: 24 xL: Top-NC  {'precision': 0.84615, 'count': 13}
   ```

### Example 4. Evaluate a CASP14 RR file
   ```bash
   wget http://deep.cs.umsl.edu/disteval/static/data/casp14/T1024/RaptorX_RR1
   wget http://deep.cs.umsl.edu/disteval/static/data/casp14/casp14_pdbs/T1024.pdb
   python3 ./disteval.py -n ./T1024.pdb -c ./RaptorX_RR1
   ```
   Expected output:
   ```
   Evaluating distances..
   min-seq-sep: 12 xL: Top-L/5 {'mae': 1.7837, 'mse': 4.9053, 'rmse': 2.2148, 'count': 78}
   min-seq-sep: 12 xL: Top-L   {'mae': 2.4797, 'mse': 13.0069, 'rmse': 3.6065, 'count': 392}
   min-seq-sep: 12 xL: Top-NC  {'mae': 3.6061, 'mse': 16.4059, 'rmse': 4.0504, 'count': 5459}
   min-seq-sep: 24 xL: Top-L/5 {'mae': 1.7837, 'mse': 4.9053, 'rmse': 2.2148, 'count': 78}
   min-seq-sep: 24 xL: Top-L   {'mae': 2.4398, 'mse': 12.8404, 'rmse': 3.5834, 'count': 392}
   min-seq-sep: 24 xL: Top-NC  {'mae': 3.6114, 'mse': 16.4634, 'rmse': 4.0575, 'count': 4906}
   Evaluating contacts..
   min-seq-sep: 12 xL: Top-L/5 {'precision': 0.9359, 'count': 78}
   min-seq-sep: 12 xL: Top-L   {'precision': 0.82143, 'count': 392}
   min-seq-sep: 12 xL: Top-NC  {'precision': 0.68562, 'count': 633}
   min-seq-sep: 24 xL: Top-L/5 {'precision': 0.9359, 'count': 78}
   min-seq-sep: 24 xL: Top-L   {'precision': 0.80357, 'count': 392}
   min-seq-sep: 24 xL: Top-NC  {'precision': 0.68631, 'count': 577}
   ```

# Evaluation through 3D modeling using `disteval.py`

## Prerequisites
- [x] Install csh
   ```bash
   sudo apt install csh
   ```
- [x] Download 'dssp-2.0.4-linux-amd64' from https://osf.io/qydjv/
   ```bash
   chmod +x dssp-2.0.4-linux-amd64
   ```
- [x] Download TM-score from https://zhanglab.ccmb.med.umich.edu/TM-score/TMscore.gz
    ```bash
    wget https://zhanglab.ccmb.med.umich.edu/TM-score/TMscore.gz
    gunzip TMscore.gz
    chmod +x TMscore
    ```
- [x] DISTFOLD
    - Follow instructions [here](DISTFOLD.md) to download DISTFOLD, an updated version of CONFOLD.

## Test
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
