# Rust implementation of SCULU (Subfamily Clustering Using Label Uncertainty)

From Audrey's MS thesis:

> Biological sequence annotation is typically performed by aligning a sequence to a database of known sequence elements.
> For transposable elements, these known sequences represent subfamily consensus sequences.
> When many of the subfamily models in the database are highly similar to each other, a sequence belonging to one subfamily can easily be mistaken as belonging to another, causing non-reproducible subfamily annotation.
> Because annotation with subfamilies is expected to give some amount of insight into a sequence’s evolutionary history, it is important that such annotation be reproducible.
> Here, we present our software tool, SCULU, which builds upon our previously- described methods for computing annotation confidence, and uses those confidence estimates to find and collapse pairs of subfamilies that have a high risk of annotation collision.
> The result is a reduced set of subfamilies, with increased expected subfamily annotation reliability.

## Setup

Download and install:

* `RepeatModeler`: https://github.com/Dfam-consortium/RepeatModeler
* `RepeatMasker`: https://www.repeatmasker.org/RepeatMasker/
* `rmblastn`: https://www.repeatmasker.org/rmblast/

Alignment calls `RepeatModeler/util/align.pl`, and multiple sequence alignment (MSA) calls `RepeatModeler/Refiner`.
For instance, I placed these into parallel directories, e.g., in _$HOME/work_.

## Discussion

To run from source, install Rust (https://www.rust-lang.org/tools/install) and use **`cargo run`**.
Following is the usage detailing the arguments:

```
$ cargo run -- -h
SCULU subfamily clustering tool

Usage: sculu [OPTIONS] --consensus <CONSENSUS> --instances <INSTANCES>... 
       --aligner <ALIGNER> --refiner <REFINER>

Options:
      --consensi <CONSENSI>           FASTA file of subfamily consensi
      --instances <INSTANCES>...      Directory of instance files for each subfamily
  -o, --outdir <OUTDIR>               Output directory [default: sculu-out]
      --lambda <LAMBDA>               Lambda value [default: 0.1227]
  -i, --independence-threshold <IND>  Independence threshold [default: 0.5]
  -c, --confidence-margin <CONF>      Confidence margin [default: 3]
      --aligner <ALIGNER>             Path to RepeatModeler/util/align.pl
      --refiner <REFINER>             Path to RepeatModeler/Refiner
      --threads <THREADS>             Number of threads for rmblastn/Refiner [default: 4]
      --perl5lib <PERL5LIB>           PERL5LIB, e.g., to find RepeatMasker/RepeatModeler
      --rmblast-dir <RMBLAST_DIR>     Path to rmblastn
      --alignment-matrix <MATRIX>     Alignment matrix
      --align-gap-init <GAPINIT>      Alignment gap init [default: -25]
      --align-extension <EXT>         Alignment extension [default: -5]
      --align-min-match <MINMATCH>    Alignment minimum match [default: 7]
      --align-bandwidth <BANDWIDTH>   Alignment bandwidth [default: 14]
      --align-mask-level <MASKLEVEL>  Alignment mask level [default: 101]
      --align-min-score <MINSCORE>    Alignment minimum score [default: 200]
  -l, --log <LOG>                     Log level [possible values: info, debug]
  -h, --help                          Print help
  -V, --version                       Print version
```

The 

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

Randomly sample 50 from each set for the input to the MSA:

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
