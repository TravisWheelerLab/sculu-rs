MATRIX = /Users/kyclark/work/RepeatMasker/Matrices/ncbi/nt/25p41g.matrix

run:
	cargo run -- \
		--consensi  tests/inputs/consensi.fa \
		--instances tests/inputs/instances/*.fa \
		--log debug \
		--independence-threshold .8 \
		--confidence-margin 3 \
        --rmblast-dir /Users/kyclark/.local/bin \
		--alignment-matrix $(MATRIX) 

#		--log-file log.out \
#        --perl5lib /Users/kyclark/work/RepeatMasker \
#        --refiner /Users/kyclark/work/RepeatModeler/Refiner \
#        --aligner /Users/kyclark/work/RepeatModeler/util/align.pl \

alu:
	cargo run -- \
		--consensi  data/alu/alu_consensi.fa \
		--instances data/alu/subfams/*.fa \
		--log debug \
		--outdir ./sculu-alu \
		--rmblast-dir /Users/kyclark/.local/bin \
		--alignment-matrix $(MATRIX) 

#--independence-threshold .5 \
#--confidence-margin 3 \
#--perl5lib /Users/kyclark/work/RepeatMasker \
#--refiner /Users/kyclark/work/RepeatModeler/Refiner \
#--aligner /Users/kyclark/work/RepeatModeler/util/align.pl \

clean:
	rm -rf RM_* makedb.log tests/inputs/consensi.fa.*
