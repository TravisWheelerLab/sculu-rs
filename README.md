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

Download Perl source code:

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

**Using Alu from Travis**:

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
data/alu/subfams/*.fa
data/alu/alu_consensi.fa
```

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

NOTE: The Python implementation next filters all the alignments for those "with scores that are 10 bits less than the maximum score of the region." Is this necessary?

These alignment scores are all we care about right now.
I can use Python module `csvkit` to look at the results:

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

I wrote a Python script to find families that need to be merged:

```
$ ./scripts/best-score.py ali.scores
...
Clean Winners
defaultdict(<class 'int'>, {'AluY': 37, 'AluYb8': 50, 'AluYm1': 36})

Winning Sets
defaultdict(<class 'int'>,
            {('AluY', 'AluYa5'): 12,
             ('AluY', 'AluYb8'): 6,
             ('AluY', 'AluYm1'): 46,
             ('AluYa5', 'AluYb8'): 6,
             ('AluYa5', 'AluYm1'): 4})
Independence AluYa5/AluYm1: 0.00
Independence AluYa5/AluYb8: 0.00
Independence AluYa5/AluY  : 0.00
Independence AluYm1/AluY  : 0.44
Independence AluY  /AluYm1: 0.45
Independence AluY  /AluYa5: 0.76
Independence AluY  /AluYb8: 0.86
Independence AluYb8/AluY  : 0.89
Independence AluYb8/AluYa5: 0.89
Independence AluYm1/AluYa5: 0.90
Independence AluY  /AluYb9: 1.00
Independence AluYa5/AluYb9: 1.00
Independence AluYb8/AluYb9: 1.00
Independence AluYb8/AluYm1: 1.00
Independence AluYb9/AluY  : 1.00
Independence AluYb9/AluYa5: 1.00
Independence AluYb9/AluYb8: 1.00
Independence AluYb9/AluYm1: 1.00
Independence AluYm1/AluYb8: 1.00
Independence AluYm1/AluYb9: 1.00
Done
```

Downsample to keep only the longest 100 instances for each family:

```
mkdir data/alu/longest_100
for file in data/alu/subfams/*.fa; do 
cargo run --bin sculu-downsample -- \
-o data/alu/longest_100/$(basename $file) $file
done
```

Randomly sample 50 from each set to select 100 for the input to the MSA:

```
mkdir data/alu/random_50
for file in data/alu/longest_100/*; do 
seqkit sample -n 50 -o data/alu/random_50/$(basename $file) $file
done
```

Create a MSA (multiple sequence alignment)/consensus using the sequences from the least independent pair, e.g., AluYa5/AluYm1:

```
export PERL5LIB=$HOME/work/RepeatMasker 
mkdir refiner
cat data/alu/random_50/{AluYa5,AluYm1}.fa > refiner/tmp.fa
$HOME/work/RepeatModeler/Refiner -threads 4 \
--rmblast_dir $HOME/.local/bin refiner/tmp.fa
```

This generates:

```
$ ls -1 refiner
tmp.fa
tmp.fa.njs
tmp.fa.refiner.stk  # MSA in Stockholm format
tmp.fa.refiner_cons # Consensus sequence in FASTA format
```

## Authors

Based on prior work by Audrey Shingleton <audrey.shingleton@umontana.edu>

* Ken Youens-Clark <kyclark@arizona.edu>
* Travis Wheeler <twheeler@arizona.edu>
