pixi run python -m scripts.datapipe_save_batches \
    --num-samples 1000000 \
    --cs-path ./data/enamine_rxn115_chemspace.bin \
    --out-dir ./results/samples_guided/ \
    --bb-weights ./results/bb_dist_1M.npy

pixi run python -m scripts.process_batches \
    --cs-path data/enamine_rxn115_chemspace.bin \
    --smarts-path ./data/rxn115.txt \
    --batch-dir results/samples_guided/ \
    --output results/processed_guided_1M.jsonl

# outside pixi
source /workspace/miniconda3/bin/activate
conda activate synllama_env

for smi in ../SynLlama/synllama_test_sets/*.smi; do
    testset_name=$(basename "$smi" .smi)
    python -m scripts.diversity_ceiling \
        --jsonl ./results/processed_guided_1M.jsonl \
        --testset "$smi" \
        --out-dir "results/diversity_ceiling/${testset_name}/" \
        --n-generated 1000000
done
