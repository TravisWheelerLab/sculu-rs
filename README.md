# Rust implementation of SCULU (Subfamily Clustering Using Label Uncertainty)

Original implementation starts with alignment using `cross_match`:

```
cross_match test_set.fa consensi.fa \
-gap_init -25 \
-gap_ext -5 \
-minmatch 7 \
-bandwidth 14 \
-masklevel 101 \
-matrix 25p41g.matrix \
-minscore 200 \
-alignments \
> test_set.fa.cm
```

New method will use `rmblastn` (https://www.repeatmasker.org/rmblast/).

NOTE: Transpose the matrices because `cross_match` and NCBI BLAST have different ideas about query/database matrix order.

Download 

* `RepeatModeler`: https://github.com/Dfam-consortium/RepeatModeler
* `RepeatMasker`: https://www.repeatmasker.org/RepeatMasker/

I placed these into parallel directories, e.g., in _$HOME/work_.
Run `RepeatModeler/util/align.pl` as follows:

```
export PERL5LIB=$HOME/work/RepeatMasker 
perl $HOME/work/RepeatModeler/util/align.pl \
-rmblast \
-gap_init -25 \
-extension -5 \
-minmatch 7 \
-bandwidth 14 \
-masklevel 101 \
-matrix $HOME/work/RepeatMasker/Matrices/ncbi/nt/25p41g.matrix \
-minscore 200 \
-alignments \
tests/inputs/test_set.fa tests/inputs/consensi.fa > test_set.fa.ali
```

Where:

* _test_set.fa_: 50 representative sequences from AluY, AluYb8, and AluYm1 scattered across Hg38
* _consensi.fa_: 10 consensus sequences for AluY, AluYa5, AluYb8, AluYb9, and AluYm1

The expectation is that each target (test_set) sequence will map to exactly one consensus sequence.

The beginning of the output file _test_set.fa.ali_ shows how to run `rmblastn`:

```
export BLASTMAT=$HOME/work/RepeatMasker/Matrices/ncbi/nt
rmblastn \
-num_alignments 9999999 \
-db tests/inputs/consensi.fa \
-query tests/inputs/test_set.fa \
-gapopen 20 \
-gapextend 5 \
-mask_level 101 \
-complexity_adjust \
-word_size 7 \
-xdrop_ungap 400 \
-xdrop_gap_final 200 \
-xdrop_gap 100 \
-min_raw_gapped_score 200 \
-dust no \
-outfmt="6 score perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen sstrand sseqid sstart send slen kdiv cpg_kdiv transi transv cpg_sites qseq sseq" \
-num_threads 4 \
-matrix 25p41g.matrix \
> test_set.fa.ali
```

Each alignment in _test_set.fa.ali_ starts with a left-justified score, here `2170`:

```
2170 6.61 0.30 7.74 AluY__hg38_chr10:116711986-116711654 1 333 (0) AluY 1 310 (1)

  AluY__hg38_ch          1 GGCCGGGCACAGTGGCTCACGCCTGTAATCCCAGCAGTTTGGGAGGCCGA 50
                                   i i                         v
  AluY                   1 GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGA 50

  AluY__hg38_ch         51 GGCCAGCGGATCACAAGGTCAGGAGATCGAGACCATCCTGGCTAACACAG 100
                              vi         i                                 i
  AluY                  51 GGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGG 100

  AluY__hg38_ch        101 TGAAATCCTGTCTCTACTAAAAATAAAAAAATTTAAAAATTTAAAAAAAA 150
                                i  i                v     -------------------
  AluY                 101 TGAAACCCCGTCTCTACTAAAAATACAAAAA------------------- 131

  AluY__hg38_ch        151 AAAACATTAGCTGGGCGCGTTGGCAGGCGCCTGTAGTCCCAGCTACTCGG 200
                           -----      i     i v    i
  AluY                 132 -----ATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGG 176

  AluY__hg38_ch        201 GAGGCTGAGGCAGGAGAATGTCGTGAACCTGGGAGGTGGAGCTTGCAGTG 250
                                               v        i      i
  AluY                 177 GAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTG 226

  AluY__hg38_ch        251 AGCCGAGATCGCACCACTGCACT-CAGCCTGGGCGACAGAGCGAGACTCA 299
                                       i          -                         v
  AluY                 227 AGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCC 276

  AluY__hg38_ch        300 GTCTCGAAAAAAAAAAAAAAGAAAAAAAAAGAAA 333
                                i              i         i
  AluY                 277 GTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAA 310

Matrix = 25p41g.matrix
Kimura (with divCpGMod) = 4.34
CpG sites = 16, Kimura (unadjusted) = 7.55
Transitions / transversions = 2.67 (16/6)
Gap_init rate = 0.08 (25 / 332), avg. gap size = 1.00 (25 / 25)
```

These are all we care about right now, so:

```
$ grep -E '^\d+\s+' test_set.fa.ali | awk 'BEGIN {OFS="\t"} {print $1, $5, $9}' > ali.scores
$ wc -l ali.scores
    2108 ali.scores
$ head ali.scores | csvlook -H
|     a | b                                    | c      |
| ----- | ------------------------------------ | ------ |
| 2,170 | AluY__hg38_chr10:116711986-116711654 | AluY   |
| 2,147 | AluY__hg38_chr10:116711986-116711654 | AluYm1 |
| 2,082 | AluY__hg38_chr10:116711986-116711654 | AluYa5 |
| 2,031 | AluY__hg38_chr10:116711986-116711654 | AluYb8 |
| 2,008 | AluY__hg38_chr10:116711986-116711654 | AluYb9 |
|   336 | AluY__hg38_chr10:116711986-116711654 | AluY   |
|   336 | AluY__hg38_chr10:116711986-116711654 | AluYm1 |
|   355 | AluY__hg38_chr10:116711986-116711654 | AluYa5 |
|   311 | AluY__hg38_chr10:116711986-116711654 | AluYb8 |
|   300 | AluY__hg38_chr10:116711986-116711654 | AluYb9 |
```

I wrote a Python script to find the best score for each target:

```
$ ./scripts/best-score.py -h
usage: best-score.py [-h] FILE

Find best score

positional arguments:
  FILE        Input file

options:

$ ./scripts/best-score.py ali.scores
{'AluY': defaultdict(<class 'int'>,
                     {'AluY': 2370,
                      'AluYa5': 2294,
                      'AluYb8': 2214,
                      'AluYb9': 2191,
                      'AluYm1': 2347}),
 'AluYb8': defaultdict(<class 'int'>,
                       {'AluY': 2449,
                        'AluYa5': 2361,
                        'AluYb8': 2663,
                        'AluYb9': 2640,
                        'AluYm1': 2426}),
 'AluYm1': defaultdict(<class 'int'>,
                       {'AluY': 2437,
                        'AluYa5': 2349,
                        'AluYb8': 2307,
                        'AluYb9': 2303,
                        'AluYm1': 2444})}
target AluY       query AluY       score 2370
target AluYb8     query AluYb8     score 2663
target AluYm1     query AluYm1     score 2444
```

## Authors

Based on prior work by Audrey Shingleton <audrey.shingleton@umontana.edu>

* Ken Youens-Clark <kyclark@arizona.edu>
* Travis Wheeler <twheeler@arizona.edu>
