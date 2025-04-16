# Rust implementation of SCULU (Subfamily Clustering Using Label Uncertainty)

The original SCULU tool was written by Audrey Shingleton in Python and handled DNA.
This is a reimplementation in Rust that extends to handling proteins, too.

## Background

From Audrey's MS thesis:

> Biological sequence annotation is typically performed by aligning a sequence
> to a database of known sequence elements.  For transposable elements, these
> known sequences represent subfamily consensus sequences.  When many of the
> subfamily models in the database are highly similar to each other, a
> sequence belonging to one subfamily can easily be mistaken as belonging to
> another, causing non-reproducible subfamily annotation.  Because annotation
> with subfamilies is expected to give some amount of insight into a
> sequence’s evolutionary history, it is important that such annotation be
> reproducible.  Here, we present our software tool, SCULU, which builds upon
> our previously- described methods for computing annotation confidence, and
> uses those confidence estimates to find and collapse pairs of subfamilies
> that have a high risk of annotation collision.  The result is a reduced set
> of subfamilies, with increased expected subfamily annotation reliability.

## Setup

SCULU requires several external tools to compare DNA/protein sequences and
generate new consensi sequences:

* Rust: https://www.rust-lang.org/tools/install
* `RepeatModeler`: https://github.com/Dfam-consortium/RepeatModeler
* `RepeatMasker`: https://www.repeatmasker.org/RepeatMasker/
* `rmblastn`: https://www.repeatmasker.org/rmblast/
* `blastp`: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
* `mafft`: https://mafft.cbrc.jp/alignment/software/
* `hmmer`: http://hmmer.org/documentation.html

## Overview

SCULU has several distinct steps:

1. Optionally generate a configuration file
2. Group the families into components by aligning the consensi sequences to
   themselves. The product will be an optional _singletons_ file containing 
   the sequences that don't group with any others and several _component-N_
   files listing the families that should be considered together.
3. Cluster all the members of each component until the resulting consensi are
   reliably differentiated.
4. Concat the singleton sequences and the newly merged component consensi into
   a final output file.

These steps can be view when running **`sculu help`**:

```
SCULU subfamily clustering tool

Usage: sculu [OPTIONS] [COMMAND]

Commands:
  config      Generate config TOML
  components  Align consensi to self and create "components"
  cluster     Cluster output from "components"
  concat      Concatenate singletons from "components" and output from "cluster"
  run         Run all steps (components, cluster, concat)
  help        Print this message or the help of the given subcommand(s)

Options:
      --logfile <LOGFILE>      Log output
      --num-threads <THREADS>  Number of threads for rmblastn/Refiner
  -h, --help                   Print help
  -V, --version                Print version
```

NOTE: The "run" action will perform the whole pipeline.

Very large datasets may take a long time to process, and the distinct steps
allow for parallelization of the "cluster" step using HPC resources. 
(Given the need for multiple cores for this step, you will likely need
separate physical machines to effectively parallelize. If running on a single
maching, you will likely need to cluster each component serially.)

## Config

You will likely want to begin by running **`sculu config`** to create a configuration file to change the runtime settings.
Run with the `-h|--help` flag to see the arguments:

```
$ sculu config -h
Generate config TOML

Usage: sculu config [OUTFILE]

Arguments:
  [OUTFILE]  Output file [default: sculu.toml]

Options:
  -h, --help  Print help
```

The file should be in TOML (Tom's Obvious Markup Language, https://toml.io/en/) format, and **`sculu config`** will generate a sample file with the default values.

```
$ sculu config example.toml
Wrote "example.toml"

$ cat example.toml
[general]
confidence_margin = 3
independence_threshold = 0.5
lambda = 0.1227
max_num_instances = 100
min_align_cover = 0.9
min_consensus_coverage = 5
min_instance_sequence_length_dna = 30
min_instance_sequence_length_prot = 12
min_len_similarity = 0.9
min_num_instances_dna = 10
min_num_instances_prot = 1

[rmblastn]
matrix = "/path/to/matrix.txt"
gap_open = 20
gap_extend = 5
word_size = 7
mask_level = 101
dust = false
complexity_adjust = false
min_raw_gapped_score = 400
xdrop_gap = 800
xdrop_ungap = 200
xdrop_gap_final = 400

[blastp]
matrix = "/path/to/matrix.txt"
gap_open = 20
gap_extend = 5
word_size = 7
mask_level = 101
dust = false
complexity_adjust = false
min_raw_gapped_score = 400
xdrop_gap = 800
xdrop_ungap = 200
xdrop_gap_final = 400
```

NOTE: TOML requires string values such as the _matrix_ to be in quotes. All
the other values are not quoted because they are interpreted as numeric
values, whether integer/whole numbers or floats.

## Running in One Step

The simplest way to execute SCULU is using the `run` command:

```
$ sculu run -h
Run all steps (components, cluster, concat)

Usage: sculu run [OPTIONS] --alphabet <ALPHABET> --consensi <CONSENSI> 
  --instances <INSTANCES> --outdir <OUTDIR>

Options:
  -a, --alphabet <ALPHABET>    Sequence alphabet [possible values: dna, protein]
      --consensi <CONSENSI>    FASTA file of subfamily consensi
      --instances <INSTANCES>  Directory of instance files for each subfamily
      --component <COMPONENT>  One file from the "components" action
      --outdir <OUTDIR>        Output directory
      --outfile <OUTFILE>      Output file [default: families.fa]
      --config <CONFIG>        Config file
  -h, --help                   Print help
```

The `--outfile` will contain the results of clustering and merging all the
families into a new set of consensi sequences with the power to differentiate
the provided instances.

If you find it takes longer than you like to use `run`, then you may wish to
execute the discrete steps that follow to parallelize the clustering step over
additional hardware.
A [Nextflow](https://github.com/TravisWheelerLab/sculu-nf) wrapper has been provided to assist in distributing this step on an HPC.

## Generate Components/Singletons

If you intend to run the steps manually, start with _components_:

```
$ sculu components -h
Align consensi to self and create "components"

Usage: sculu components [OPTIONS] --alphabet <ALPHABET> --consensi <CONSENSI> 
  --instances <INSTANCES>

Options:
  -a, --alphabet <ALPHABET>    Sequence alphabet [possible values: dna, protein]
      --consensi <CONSENSI>    FASTA file of subfamily consensi
      --instances <INSTANCES>  Directory of instance files for each subfamily
      --outdir <OUTDIR>        Output directory [default: sculu-out]
      --config <CONFIG>        Config file
  -h, --help                   Print help
```

NOTE: The consensi sequence names *must* match the filenames of the _instances_ files.

This step aligns the consensi sequences to themselves using `rmblastn` for DNA
and `blastp` for protein.

The directory structure of this step will be:

```
$ tree sculu-out -L 1
sculu-out
├── components
├── consensi.fa
├── consensi_cluster
├── instances
└── select
```

* _components_: This directory contains the consensi found to be sufficiently similar according to the settings in the configuration. Consensi sequences that are not similar to *any* others are placed into a _singletons_ file that is needed for the final _concat_ step.  The files are named _component-N_ (for N from 0 to the number of components) and will contain the consensi sequence/family names. *These files are the inputs to the next step of clustering.*
* _consensi.fa_: This file contains only consensi supported by sufficient instances.
* _consensi_cluster_: This directory contains the result of clustering the consensi.
* _instances_: This directory contains downsampled instances filtered to those of the best quality.
* _select_: This directory contains the artefacts of selecting the best instances.

NOTE: You should use the _consensi.fa_ and _instances_ from this step in following clustering step.

## Cluster and Merge Components

The next step is clustering and merging the members of each component:

```
$ sculu cluster -h
Cluster output from "components"

Usage: sculu cluster [OPTIONS] --alphabet <ALPHABET> --consensi <CONSENSI> 
  --instances <INSTANCES> --component <COMPONENT> --outdir <OUTDIR>

Options:
  -a, --alphabet <ALPHABET>    Sequence alphabet [possible values: dna, protein]
      --consensi <CONSENSI>    The filtered FASTA consensi from "components"
      --instances <INSTANCES>  Directory of filtered instance files from "components"
      --component <COMPONENT>  A file from "components" action
      --outdir <OUTDIR>        Output directory
      --config <CONFIG>        Config file
  -h, --help                   Print help
```

NOTE: When running manually, you should use the _filtered_ consensi and
instances produced by the _components_ step.

All the families in a given component are compared to find sequences that are
not reliably differentiated and merged.
The process is repeated until no more sequences need to be merged.

This step can be parallelized over multiple machines, but it would probably be
a very bad idea to parallelize on a single machine because each process will
likely require all available cores, which would lead to resource contention
and degraded performance.
The Nextflow wrapper will handle distribution of the component clustering to
separate nodes; however, very small components (e.g., having 10 or fewer
families) will likely finish very quickly while extremely large components
(having hundreds or thousands of families) may take many hours.

The result of this step will be new _component-N_ directories in the output
directory. 
Each of these will contain various _roundXXX_ directories containing the
artefacts of clustering and merging the consensi until no more merges are
needed.
Each directory should contain a _final.fa_ file that has the end result of
merging the families.

## Concatenation

The final step is the concatenation of the original consensi sequences that
landed in the _singletons_ file with the final merged products of the
components from cluster step:

```
$ sculu concat -h
Concatenate singletons from "components" and output from "cluster"

Usage: sculu concat [OPTIONS] --outfile <OUTFILE> --consensi <CONSENSI>

Options:
  -o, --outfile <OUTFILE>             Output file
      --consensi <CONSENSI>           The filtered FASTA consensi from "components"
      --singletons <SINGLETONS>       Singletons file from "components"
      --components [<COMPONENTS>...]  Merged components from "cluster"
  -h, --help                          Print help
```

The _singletons_ file is optional because the _components_ step may not
produce this file.
The `--components` argument should be all the _component-N/final.fa_ files
from the _cluster_ step.
The final output file will be in FASTA format.
The sequence headers will contain the families that were merged in Newick
format, e.g.:

```
>(AluYh7,(AluYh9,AluYh3):0.11):0.25;
```

In the preceding example, AluYh9 and AluYh3 were found to be only 11%
independent and were merged first.
The resulting family was found to be only 25% independent from AluYh7 and
subsequently merged.

## Discussion

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

## See Also

* https://github.com/TravisWheelerLab/sculu-nf

## Authors

* Audrey Shingleton <audrey.shingleton@umontana.edu> (prior Python implementation for 2022 MS thesis)
* Ken Youens-Clark <kyclark@arizona.edu> (Rust implementation)
* Travis Wheeler <twheeler@arizona.edu> (Advisor)
