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

Alignment calls `rmblastn`, and multiple sequence alignment (MSA) calls `RepeatModeler/Refiner`.
For instance, I placed these into parallel directories, e.g., in _$HOME/work_.

## Discussion

To run from source, install Rust (https://www.rust-lang.org/tools/install) and use **`cargo run`**.
Following is the usage detailing the arguments:

```
$ cargo run -- -h
SCULU subfamily clustering tool

Usage: sculu [OPTIONS] --consensi <CONSENSI> --instances <INSTANCES>...

Options:
  -c, --consensi <CONSENSI>
          FASTA file of subfamily consensi
  -i, --instances <INSTANCES>...
          Instance files for each subfamily
  -o, --outfile <OUTFILE>
          Output file [default: families.fa]
      --outdir <OUTDIR>
          Output directory (if you want working files preserved)
      --lambda <LAMBDA>
          Lambda value [default: 0.1227]
      --independence-threshold <IND>
          Independence threshold [default: 0.5]
      --confidence-margin <CONF>
          Confidence margin [default: 3]
      --aligner <ALIGNER>
          Path to rmblastn
      --refiner <REFINER>
          Path to RepeatModeler/Refiner
      --num-threads <THREADS>
          Number of threads for rmblastn/Refiner
      --perl5lib <PERL5LIB>
          PERL5LIB, e.g., to find RepeatMasker/RepeatModeler
      --rmblast-dir <RMBLAST_DIR>
          Path to rmblastn
      --align-matrix <MATRIX>
          Alignment matrix
      --align-gap-open <ALIGN_GAP_OPEN>
          Alignment gap open penalty [default: 20]
      --align-gap-extension <ALIGN_GAP_EXT>
          Alignment gap extension penalty [default: 5]
      --align-word-size <ALIGN_WORD_SIZE>
          Alignment word size [default: 7]
      --align-xdrop-gap <ALIGN_XDROP_GAP>
          Alignment x-drop gap [default: 100]
      --align-xdrop-ungap <ALIGN_XDROP_UNGAP>
          Alignment x-drop ungap [default: 400]
      --align-xdrop-final <ALIGN_XDROP_FINAL>
          Alignment x-drop final [default: 200]
      --align-dust
          Alignment dust option
      --align-complexity-adjust
          Alignment complexity adjust
      --align-mask-level <MASKLEVEL>
          Alignment mask level [default: 101]
      --align-min-score <MINSCORE>
          Alignment minimum score [default: 200]
  -h, --help
          Print help
  -V, --version
          Print version
```

The `--consensi` file should be a FASTA file with sequence IDs of the families, e.g.:

```
$ grep -E '^>' tests/inputs/consensi.fa | sort
>AluY
>AluYa5
>AluYb8
>AluYb9
>AluYm1
```

The file names of the `--instances` must match a consensi IDs, e.g.:

```
$ ls -1 tests/inputs/instances
AluY.fa
AluYa5.fa
AluYb8.fa
AluYb9.fa
AluYm1.fa
```

Note:

* The `--outfile` argument is optional and defaults to _families.fa_. The program will halt if the outfile already exists and has a nonzero size.
* The `--outdir` argument is optional and defaults to a hidden temporary directory that will be removed after the program runs. Specify an output directory if you wish to retain all the build artefacts include the _debug.log_ file with verbose output of the program's choices.
* If you wish to watch the progress, use `tail -f <outdir>/debug.log`.

## Testing

Run `make run` in the root directory to verify that SCULU correctly clusters a few Alu families.

```
$ make run
cargo run -- \
		--outfile                test.fa \
		--outdir                 out-test \
		--consensi               tests/inputs/consensi.fa \
		--instances              tests/inputs/instances/*.fa \
		--independence-threshold .8 \
		--confidence-margin      3 \
		--rmblast-dir            /Users/kyclark/.local/bin \
		--align-matrix /Users/kyclark/work/RepeatMasker/Matrices/ncbi/nt/25p41g.matrix
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.04s
	 Running `target/debug/sculu --outfile test.fa --outdir out-test --consensi
		tests/inputs/consensi.fa --instances tests/inputs/instances/AluY.fa
		tests/inputs/instances/AluYa5.fa tests/inputs/instances/AluYb8.fa
		tests/inputs/instances/AluYb9.fa tests/inputs/instances/AluYm1.fa
		--independence-threshold .8 --confidence-margin 3 --rmblast-dir
		/Users/kyclark/.local/bin --align-matrix
		/Users/kyclark/work/RepeatMasker/Matrices/ncbi/nt/25p41g.matrix`
Done, see final families in "test.fa"
```

The resulting _test.fa_ file should have the following:

```
$ cat test.fa
>(AluYb9,AluYb8):0.74;
GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCATGAGGTCAGGAGATCGA
GACCATCCTGGCTAACAAGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGCGGTGGCGGGCGCCTGT
AGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAAGCGGAGCTTGCAGTGAGCCGAGATTGCGC
CACTGCAGTCCGCAGTCCGGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAA
>AluYa5
GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGA
GACCATCCCGGCTAAAACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTAGTGGCGGGCGCCTGT
AGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGC
CACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>(AluYm1,AluY):0.74;
GGCCGGGCGCGGTGGCTCANGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGA
GACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGT
AGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGC
CACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

## Method

First, SCULU will do the following setup:

- Optionally create the `--outdir` or use a temporary directory for all working files.
- Copy the `--instances` files into an _instances_ directory inside the output directory.
- Check that the consensi names and instances match and there are no duplicated family names.
- Run `rmblastn` on all the consensi to determine if there are clusters of families to run in batches.

For each batch of similar families:
- Concatenate their instance sequences a file _batchXXX/all_seqs.fa_ for use in alignments.
- Copy the family consensi sequences to _batchXXX/consensi.fa_.

Inside each batch, consensi are merged when `rmblastn` is unable to use the family instances to discriminate the family assignments.
In `round000`, the consensi sequence names are replaced with integer values and the ID/family names are moved to the description field to prevent errors when family names are concatenated into IDs that are too long for `makeblastdb` to handle.

Next:

- Create a _roundXXX_ directory for the current round's output.
- Run alignment of _all_seqs.fa_ to the consensi to create _roundXXX/blast.out_. This file can be large but is retained for record keeping.
- Extract the scores from the _blast.out_ into _roundXXX/alignment-scores.tsv_, a tab-separated file that lists the query, target, and alignment score. 
- For each target in the scores file, note the highest score for each query (consensi). Then iterate over the targets to figure out which consensi are "clear winners" (no contention with other sequences) and which are involved in "winning sets" with other sequences, meaning they lack discriminatory power.
- Use the total number of "clear winners" and "winning sets" to determine the independence of each pair of consensi, sorted from least independent to most with a cutoff of the `--independence-threshold`, e.g., 0.5 means pairs found to be less than 50% independent are marked for merging.
- If all consensi are determined to be independent, exit the loop.
- Merge each pair of non-independent consensi. Sample the instance sequences for input to a multiple sequence alignment. The new consenus sequence ID will be a Newick-formatted string noting the two families and their independence value.
- On the next round, use the new consensi file that contains only the newly merged sequences. (There would be no new information from re-aligning the other sequences.) When extracting the alignment scores, include the previous round's score file to add in all the previous targets that weren't included in this new file.

Finally, when the loop exits after all consensi merges, write a new consensi file.

## Testing

- Run **`tar xvf tests/inputs/alignment.tgz`**
- Run **`cargo test`**

## Authors

* Audrey Shingleton <audrey.shingleton@umontana.edu> (prior Python implementation for 2022 MS thesis)
* Ken Youens-Clark <kyclark@arizona.edu> (Rust implementation)
* Travis Wheeler <twheeler@arizona.edu> (Advisor)
