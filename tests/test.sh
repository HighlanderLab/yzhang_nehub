conda run -n ts_tools python ../tests/generate_test_data.py \
                --outdir ${{ env.TEST_DATA }} \
                --n_per_pop 3 \
                --seq_length 1_000_000

# ── Run all methods ───────────────────────────────────────────────────
nextflow run main.nf \
            -profile github,github_small \
            -params-file ./ci_params.json \
            --input_dir ${{ env.TEST_DATA }} \
            --output_dir ../tests/output \
            --chromosomes 1 \
            --population_file popLabel.txt \
            -resume