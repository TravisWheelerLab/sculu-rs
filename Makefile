filter:
	cargo run --bin filter-alignments -- test_set.fa.ali

scores:
	./scripts/best-score.py ali.scores
