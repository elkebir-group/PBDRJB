# PBDRJD -- Phylogenetic  Bulk  DNA-RNA  Joint  Deconvolutio
PBDRJD is an algorithm for phylogenetic inference of clone-specific expression and mutationprofiles using matched RNA and DNA bulk sequencing samples.
![Overview of PBDRJD](doc/overview.png)

### Dependencies    
PBDRJD is implemented in Python and compatible with Python 2.7+ and Python 3. In addition, PBDRJD has the following dependencies:
- [Numpy](http://www.numpy.org/)
- [CPLEX](https://www.ibm.com/analytics/cplex-optimizer)

### Example
An example execution
```
python src/pbdrjb.py --input_dir example --output_dir example
```

### Usage 
```
usage: pbdrjb.py [-h] [--input_dir INPUT_DIR] [--output_dir OUTPUT_DIR]
                 [--alpha ALPHA] [--eps EPS]

optional arguments:
  -h, --help            show this help message and exit
  --input_dir INPUT_DIR
                        directory of input directory (default: None)
  --output_dir OUTPUT_DIR
                        directory of output directory (default: None)
  --alpha ALPHA         edges with neutral probability greater than this
                        threshold will be removed (default: 0.99)
  --eps EPS             neglectable constant used to avoid numerical issues
                        (default: 1e-05)

```

### Input format
The `input_dir` directory must have the following input files:
 - `B.txt`: a `n x m` matrix describing the clone tree where `n` is the number of clones and `m` is the number of genes. The file should have `n` lines and each line contains `m` numbers.
 - `u.txt`: a `1 x n` vector describing the clone composition. The file has one line with `n` numbers.
 - `e.txt`: a `1 x m` vector describing the RNA abundance in the matched normal sample. The file has one line with `m` numbers.
 - `d.txt`: a `1 x m` vector describing the RNA abundance in the tumor sample. The file has one line with `m` numbers.
 - `Z_minus`: a `m x m` matrix describing the down-regulation probability given by xseq.  The file has `m` lines and each line has `m` numbers. The entry at the `q`-th row and `p`-th column represents the probability of down-regulation of gene `p` due to the mutation in gene `q`.
 - `Z_zero`: a `m x m` matrix describing the neutral probability given by xseq.  The file has `m` lines and each line has `m` numbers. The entry at the `q`-th row and `p`-th column represents the probability of being neutral of gene `p` due to the mutation in gene `q`.
 - `Z_plus`: a `m x m` matrix describing the up-regulation probability given by xseq.  The file has `m` lines and each line has `m` numbers. The entry at the `q`-th row and `p`-th column represents the probability of up-regulation of gene `p` due to the mutation in gene `q`.

### Output format
The following output files will be saved in the `output_dir`:
- `output_summary.txt`: a human-readable text file describing the input and output of the PDBRJB problem.
- `output_C.txt`: a `n x m` matrix describing the clone-specific expression profile. The j-th value in the i-th row is the relative abundance of gene j in clone i.
- `output_model.lp`: the CPLEX model file.
