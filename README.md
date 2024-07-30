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
...
SCULU subfamily clustering tool

Usage: sculu [OPTIONS] --consensi <CONSENSI> --instances <INSTANCES>...

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
      --log-file <LOG_FILE>           Log file
  -h, --help                          Print help
  -V, --version                       Print version
```

The `--consensi` file should be a FASTA file with sequence IDs of the families, e.g.:

```
$ grep -E '^>' tests/inputs/consensi.fa
>AluYm1
>AluYb8
>AluYb9
>AluYa5
>AluY
```

The `--instances` should each have a name matching the consensi IDs, e.g.:

```
$ ls tests/inputs/instances
AluY.fa   AluYa5.fa AluYb8.fa AluYb9.fa AluYm1.fa
```

## Method

First, SCULU will do the following setup:

- Create the `--outdir`, which defaults to _sculu-out_, for all program output.
- Take the 100 longest sequences from each of the instances and place them into the same filename in the _outdir/instances_100_ directory.
- Check that the consensi names and instances match and there are no duplicated family names.
- Concatenate all the 100 longest sequences into the file _outdir/all_seqs.fa_ for use in alignments.
- Write a new version of the consensi sequences to _outdir/consensi.fa_ where the IDs are replaced with integer values and the ID/family names are moved to the description field. This prevents errors when family names are concatenated into IDs that are too long for `makeblastdb` to handle.

Next, SCULU enters a loop with the following actions:

- Create a _roundXXX_ directory for the current round's output.
- Run alignment of _all_seqs.fa_ to the consensi to create _roundXXX/alignment.txt_. This file can be large but is retained for record keeping.
- Extract the scores from the _alignment.txt_ into _roundXXX/alignment-scores.tsv_, a tab-separated file that lists the query, target, and alignment score. 
- For each target in the scores file, note the highest score for each query (consensi). Then iterate over the targets figure out which consensi are "clear winners" (no contention with other sequences) or are involved in "winning sets" with other sequences, meaning they lack discriminatory power.
- Use the total number of "clear winners" and "winning sets" to determine the independence of each pair of consensi, sorted from least independent to most with a cutoff of the `--independence-threshold`, e.g., 0.5 means pairs found to be less than 50% independent are marked for merging.
- If all consensi are determined to be independent, exit the loop.
- Merge each pair of non-independent consensi. Sample the instance sequences for input to a multiple sequence alignment. The new consenus sequence ID will be a Newick-formatted string noting the two families and their independence value.
- On the next round, use the new consensi file that contains only the newly merged sequences. (There would be no new information from re-aligning the other sequences.) When extracting the alignment scores, include the previous round's score file to add in all the previous targets that weren't included in this new file.

Finally, when the loop exits after all consensi merges, write a new consensi file.

NOTE: The Python implementation filters all the alignments for those "with scores that are 10 bits less than the maximum score of the region." Is this necessary?

## Testing

- Run **`tar xvf tests/inputs/alignment.tgz`**
- Run **`cargo test`**

## Authors

* Audrey Shingleton <audrey.shingleton@umontana.edu> (prior Python implementation for 2022 MS thesis)
* Ken Youens-Clark <kyclark@arizona.edu> (Rust implementation)
* Travis Wheeler <twheeler@arizona.edu> (Advisor)
