MATRIX = $(shell pwd)/tests/inputs/matrices/25p41g.matrix

run:
	cargo run -- \
		--outdir                 out-test \
		--outfile                out-test/test.fa \
		--consensi               tests/inputs/consensi.fa \
		--instances              tests/inputs/instances \
		--independence-threshold .8 \
		--confidence-margin      3 \
		--align-matrix           $(MATRIX) 

prot:
	cargo run -- \
		--outdir                 out-protein \
		--outfile                out-protein/test.fa \
		--consensi               data/proteins/kazel-10.fa \
		--instances              tests/inputs/instances \
		--independence-threshold .8 \
		--confidence-margin      3 \
		--align-matrix           $(MATRIX) 

dup:
	cargo run -- \
		--instances tests/inputs/instances \
		--consensi tests/inputs/dup_consensi.fa 

alu:
	cargo run -- \
		--consensi     data/alu/alu_consensi.fa \
		--instances    data/alu/subfams \
		--outfile      ./alu-out/final-alu.fa \
		--outdir       ./alu-out \
		--align-matrix $(MATRIX) 

mix:
	cargo run -- \
		--outdir     out-mixed \
		--consensi   data/mixed/consensi.fa \
		--instances  data/mixed/instances \
		--outfile    final-mixed.fa \
		--align-matrix $(MATRIX)


tua1:
	cargo run -- \
        --align-matrix $(MATRIX) \
		--build-components-only \
		--logfile   - \
        --outdir    tuatara-out \
        --outfile   tuatara-out/new-consensi.fa \
        --instances ~/wheelerlab/tuatara/instances \
        --consensi  ~/wheelerlab/tuatara/consensi.fa

tua2:
	cargo run -- \
        --align-matrix $(MATRIX) \
        --outdir       tuatara-out \
        --logfile      - \
        --component    tuatara-out/components/component-007 \
        --instances    tuatara-out/instances \
        --consensi     tuatara-out/consensi.fa

clean:
	rm -rf RM_* makedb.log
