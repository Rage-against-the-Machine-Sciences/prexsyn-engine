pixi run python -m scripts.datapipe_save_batches \
    --num-samples 10000000 \
    --cs-path ./data/enamine_rxn115_chemspace.bin \
    --out-dir ./results/samples/

pixi run python -m scripts.process_batches \
    --cs-path data/enamine_rxn115_chemspace.bin \
    --smarts-path ./data/rxn115.txt \
    --batch-dir results/samples/ \
    --output results/processed_10M.jsonl

# outside pixi
source /workspace/miniconda3/bin/activate
conda activate synllama_env

for smi in ../SynLlama/synllama_test_sets/*.smi; do
    testset_name=$(basename "$smi" .smi)
    python -m scripts.diversity_ceiling \
        --jsonl ./results/processed_10M.jsonl \
        --testset "$smi" \
        --out-dir "results/diversity_ceiling/${testset_name}/" \
        --n-generated 10000000
done
