python -m scripts.datapipe_save_batches \
  --num-samples 1000000 \
  --cs-path ./data/enamine_rxn115_chemspace.bin \
  --out-dir ./results/samples/

python -m scripts.process_batches \
  --cs-path data/enamine_rxn115_chemspace.bin \
  --smarts-path ./data/rxn115.txt \
  --batch-dir results/samples/ \
  --output results/processed_1M.jsonl
