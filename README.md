## Getting Started

```sh
wget https://github.com/lh3/dipcall/releases/download/v0.3/dipcall-0.3_x64-linux.tar.bz2
tar -jxf dipcall-0.3_x64-linux.tar.bz2
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

## Output

The final output of dipcall includes two files: *prefix*.dip.vcf.gz and
*prefix*.dip.bed. A raw variant call is made in the VCF if a non-reference
allele is observed in any alignments >=50kb and mapped with mapping quality >=5.
In the VCF, dipcall sets a `HET1` filter if parent1 is a heterozygous at the
raw variant; dipcall sets a `GAP1` filter if no >=50kb alignment from parent1
covers the variant. Note that if a call is covered by multiple >=50kb
alignments in the same parent but the alignments all have the same allele, the
call is not filtered in the VCF.

The BED file gives the confident regions. A base is included in the BED if 1)
it is covered by one >=50kb alignment with mapQ>=5 from each parent and 2)
it is not covered by other >=10kb alignments in each parent. Nearly all calls
filtered in the VCF are excluded in the BED, except very rare edge cases.
However, a fraction of unfiltered calls in the VCF may also be excluded by the
BED. The BED file applies more stringent filter.

The above is applied to autosomes and female chrX. For a male sample, parent1
is assumed to be the father and parent2 the mother. Dipcall treats PARs the
same way as autosomes. However, outside PARs, dipcall filters out chrX regions
covered by father contigs and filters out chrY regions covered by mother
contigs. To make proper calls on sex chromosomes, users should hard mask PARs
on the reference chrY.

[syndip]: https://github.com/lh3/CHM-eval
