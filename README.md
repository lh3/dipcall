## Getting Started

```sh
wget https://github.com/lh3/dipcall/releases/download/v0.1/dipcall-0.1_x64-linux.tar.bz2
tar -jxf dipcall-0.1_x64-linux.tar.bz2
# for female
dipcall.kit/run-dipcall prefix hs38.fa pat.fa.gz mat.fa.gz > prefix.mak
# or for male, requiring PAR regions in BED
# dipcall.kit/run-dipcall -x dipcall.kit/hs38.PAR.bed prefix hs38.fa pat.fa.gz mat.fa.gz > prefix.mak
make -j2 -f prefix.mak
```

## Introduction

Dipcall is a reference-based variant calling pipeline for a pair of phased
haplotype assemblies. It was originally developed for constructing the
[syndip][syndip] benchmark dataset and has been applied to other phased
assemblies, too. Dipcall can call small variants and long INDELs as long as
they are contained in minimap2 alignment.

If you use dipcall, please cite

> Li H, Bloom JM, Farjoun Y, Fleharty M, Gauthier L, Neale B, MacArthur D
> (2018) A synthetic-diploid benchmark for accurate variant-calling evaluation.
> Nat Methods, 15:595-597. [PMID:30013044]

[syndip]: https://github.com/lh3/CHM-eval
