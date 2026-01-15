# bammap2
![](https://img.shields.io/github/v/tag/wangyibin/bammap2) ![](https://img.shields.io/github/commit-activity/m/wangyibin/bammap2) ![](https://img.shields.io/github/last-commit/wangyibin/bammap2) ![platforms](https://img.shields.io/badge/platforms-aarch64%20|%20x86-blue) [](https://img.shields.io/github/downloads/wangyibin/bammap2/total?style=flat)

## Introduction 
A lightweight wrapper to run minimap2 (v2.30) alignments directly on BAM file, powered by a robust Rust wrapper for minimap2 ([minimap2-rs](https://github.com/jguhlin/minimap2-rs)).

### Motivation
- **Efficiency**: The `dorado` binary is excessively large for users who only require read alignment.

- **Flexibility**: `dorado` provides limited access to `minimap2` parameters, restricting fine-grained control over the alignment process.



## Current supported parameters
```shell
Usage: bammap2 [OPTIONS] <reference> <query>...

Arguments:
  <reference>  target.fa or target.idx
  <query>...   query sequences with BAM Format, only support BAM files

Options:
  -h, --help     Print help
  -V, --version  Print version

Indexing:
  -k <INT>      k-mer size (no larger than 28) [15]
  -w <INT>      minimizer window size [5]
  -I <NUM>      split index for every ~NUM input bases [16G] [default: 16G]

Mapping:
  -f <FLOAT>          filter out top FLOAT fraction of repetitive minimizers [0.0002]
  -g <INT>            stop chain enlongation if there are no minimizers in INT-bp [5000]
  -G <INT>            max intron length (effective with -xsplice; changing -r) [200k]
  -F <INT>            max fragment length (effective with -xsr or in the fragment mode) [800]
  -r <INT,[INT]>      chaining/alignment bandwidth and long-join bandwidth [500,20000]
  -n <INT>            minimal number of minimizers on a chain [3]
  -m <INT>            minimal chaining score (matching bases minus log gap penalty) [40]
  -p <FLOAT>          Min secondary-to-primary score ratio [0.8]
  -N <INT>            Retain at most N secondary alignments [5]

Alignments:
  -A <INT>            matching score [2]
  -B <INT>            mismatch penalty (larger value for lower divergence) [4]
  -O <INT,[INT]>      gap open penalty [4,24]
  -E <INT,[INT]>      gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]
  -z <INT,[INT]>      Z-drop score and inversion Z-drop score [400,200]
  -s <INT>            minimal peak DP alignment score [80]

Input/Output:
  -o <FILE>                    output file path, currently only support output in BAM formst [stdout] [default: -]
  -Y                           use soft clipping for supplementary alignments
      --seed <INT>             Integer seed for randomizing equally best hits. Minimap2 hashes INT and read name when choosing between equally best hits. [11]
      --max-qlen <INT>         skip reads longer than INT [0 (disabled)]
      --secondary <secondary>  Whether to output secondary alignments [default: yes] [possible values: yes, no]
  -t <INT>                     number of threads [default: 8]
  -K <STR>                     minibatch size logging when mapping [default: 10k]

Presets:
  -x <STR>      - lr:hq - accurate long reads (error rate <1%) against a reference genome
                - splice/splice:hq - spliced alignment for long reads/accurate long reads
                - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
                - sr - short reads against a reference
                - map-pb/map-hifi/map-ont - CLR/HiFi/Nanopore vs reference mapping
                - ava-pb/ava-ont - PacBio CLR/Nanopore read overlap
                 [default: lr:hq] [possible values: lr:hq, map-hifi, map-pb, map-ont, asm5, asm10, asm20, sr, splice, splice:hq]
```

## Citation
Please cite the minimap2 papers.
### Minimap2 Papers

> Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
> *Bioinformatics*, **34**:3094-3100. [doi:10.1093/bioinformatics/bty191][doi]

and/or:

> Li, H. (2021). New strategies to improve minimap2 alignment accuracy.
> *Bioinformatics*, **37**:4572-4574. [doi:10.1093/bioinformatics/btab705][doi2]
