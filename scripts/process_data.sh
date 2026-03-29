python -m scripts.datapipe_save_batches \
  --num-samples 100000 \
  --cs-path ./data/enamine_rxn115_chemspace.bin \
  --out-dir ./results/samples/

python -m scripts.process_batches \
  --cs-path data/enamine_rxn115_chemspace.bin \
  --smarts-path ./data/rxn115.txt \
  --batch-dir results/samples/ \
  --output results/processed_100k.jsonl


for smi in ../SynLlama/synllama_test_sets/*.smi; do
    testset_name=$(basename "$smi" .smi)
    python -m scripts.diversity_ceiling \
        --jsonl ./results/processed_100k.jsonl \
        --testset "$smi" \
        --out-dir "results/diversity_ceiling/${testset_name}/" \
        --n-generated 100000
done
