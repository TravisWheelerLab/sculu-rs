MATRIX = $(shell pwd)/tests/inputs/matrices/25p41g.matrix

cluster:
	cargo run -- --outdir out-test --align-matrix $(MATRIX) \
		cluster tests/inputs/consensi.fa

mixclust:
	cargo run -- --outdir out-mixed --align-matrix $(MATRIX) \
		cluster data/mixed/consensi.fa

mixmerge:
	cargo run -- --outdir out-mixed --align-matrix $(MATRIX) \
		merge \
		--components   out-mixed/components.json \
		--consensi     data/mixed/consensi.fa \
		--instances    data/mixed/instances \
		--outfile      final-mixed.fa 

run:
	cargo run -- \
		--outfile                test.fa \
		--outdir                 out-test \
		--consensi               tests/inputs/consensi.fa \
		--instances              tests/inputs/instances/*.fa \
		--independence-threshold .8 \
		--confidence-margin      3 \
		--align-matrix           $(MATRIX) 

alu:
	cargo run -- \
		--consensi     data/alu/alu_consensi.fa \
		--instances    data/alu/subfams/*.fa \
		--outfile      final-alu.fa \
		--outdir       ./alu-out \
		--align-matrix $(MATRIX) 


clean:
	rm -rf RM_* makedb.log tests/inputs/consensi.fa.*
