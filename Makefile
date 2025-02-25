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
        --outdir    tuatara-out \
        --outfile   tuatara-out/new-consensi.fa \
        --instances data/tuatara/instances \
        --consensi  data/tuatara/consensi.fa

tua2:
	cargo run -- \
        --align-matrix $(MATRIX) \
        --outdir       tuatara-out \
		--logfile      - \
        --outfile      tuatara-out/new-consensi.fa \
        --component    tuatara-out/components/component-0 \
        --instances    data/tuatara/instances \
        --consensi     data/tuatara/consensi.fa

clean:
	rm -rf RM_* makedb.log
