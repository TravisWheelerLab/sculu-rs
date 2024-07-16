run:
	cargo run --bin sculu -- \
		--perl5lib $(HOME)/work/RepeatMasker \
		--aligner $(HOME)/work/RepeatModeler/util/align.pl \
		--consensus tests/inputs/consensi.fa \
		--instances tests/inputs/test_set.fa

run2:
	cargo run --bin sculu -- \
		--perl5lib $(HOME)/work/RepeatMasker \
		--aligner $(HOME)/work/RepeatModeler/util/align.pl \
		--consensus data/alu/alu_consensi.fa \
		--instances data/alu/subfams/*.fa


ds:
	cargo run --bin sculu-downsample -- data/alu/subfams/AluJb.fa

filter:
	cargo run --bin filter-alignments -- test_set.fa.ali

scores:
	./scripts/best-score.py ali.scores
