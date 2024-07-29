MATRIX = /Users/kyclark/work/RepeatMasker/Matrices/ncbi/nt/25p41g.matrix

run:
	cargo run --bin sculu -- \
		--log debug \
		--independence-threshold .8 \
		--confidence-margin 3 \
		--perl5lib /Users/kyclark/work/RepeatMasker \
		--refiner /Users/kyclark/work/RepeatModeler/Refiner \
		--threads 8 \
		--rmblast-dir /Users/kyclark/.local/bin \
		--aligner /Users/kyclark/work/RepeatModeler/util/align.pl \
		--alignment-matrix $(MATRIX) \
		--consensus tests/inputs/consensi.fa \
		--instances tests/inputs/instances/*.fa

alu:
	cargo run --bin sculu -- \
		--log debug \
		--outdir ./sculu-alu \
		--independence-threshold .4 \
		--confidence-margin 3 \
		--perl5lib /Users/kyclark/work/RepeatMasker \
		--refiner /Users/kyclark/work/RepeatModeler/Refiner \
		--threads 8 \
		--rmblast-dir /Users/kyclark/.local/bin \
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

clean:
	rm -rf RM_* makedb.log
