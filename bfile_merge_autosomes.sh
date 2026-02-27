#!/bin/bash
# merges plink binaries of different chromosomes if they belong to the same subject
# inspired by https://apol1.blogspot.com/2013/10/merging-plink-binary-files.html
cp chr1.fam autosomes.fam
for chrom in {1..22}; do cat chr${chrom}.bim; done > autosomes.bim
(echo -en "\x6C\x1B\x01"; for chrom in {1..22}; do tail -c chr${chrom}.bed; done) > autosomes.bed
