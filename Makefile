MATRIX = $(shell pwd)/tests/inputs/matrices/25p41g.matrix

# --rmblast-dir            /Users/kyclark/.local/bin \
# --rmblast-dir            /usr/local/bin \

run:
	cargo run -- \
		--outfile                test.fa \
		--outdir                 out-test \
		--consensi               tests/inputs/consensi.fa \
		--instances              tests/inputs/instances/*.fa \
		--independence-threshold .8 \
		--confidence-margin      3 \
		--align-matrix $(MATRIX) 

alu:
	cargo run -- \
		--consensi     data/alu/alu_consensi.fa \
		--instances    data/alu/subfams/*.fa \
		--outfile      final-alu.fa \
		--outdir       ./alu-out \
		--align-matrix $(MATRIX) 

mixed:
	cargo run -- \
		--force \
		--outfile      final-mixed.fa \
		--outdir       ./mixed-out \
		--consensi     data/mixed/consensi.fa \
		--instances    data/mixed/instances/*.fa \
		--align-matrix $(MATRIX) 

clean:
	rm -rf RM_* makedb.log tests/inputs/consensi.fa.*
