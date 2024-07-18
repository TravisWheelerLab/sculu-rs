MATRIX = /Users/kyclark/work/RepeatMasker/Matrices/ncbi/nt/25p41g.matrix

run:
	cargo run --bin sculu -- \
		--debug \
		--perl5lib /Users/kyclark/work/RepeatMasker \
		--aligner /Users/kyclark/work/RepeatModeler/util/align.pl \
		--alignment-matrix $(MATRIX) \
		--consensus tests/inputs/consensi.fa \
		--instances tests/inputs/instances/*.fa

run2:
	cargo run --bin sculu -- \
		--debug \
		--perl5lib /Users/kyclark/work/RepeatMasker \
		--aligner /Users/kyclark/work/RepeatModeler/util/align.pl \
		--alignment-matrix $(MATRIX) \
		--consensus data/alu/alu_consensi.fa \
		--instances data/alu/subfams/*.fa

ds:
	cargo run --bin sculu-downsample -- data/alu/subfams/AluJb.fa

filter:
	cargo run --bin filter-alignments -- test_set.fa.ali

scores:
	./scripts/best-score.py ali.scores
