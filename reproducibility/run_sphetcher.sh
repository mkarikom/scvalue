i=10000 # sketch size=10000, for example
for dataset in "pbmc" "tma" "human_gut" "Spleen" "mouse_T" "human_PAC" "gorilla" "monkey" "tumour" "SEA_AD" "mouse_CNS" "human_fetal"
do
    csv_file="experiments/${dataset}/Sphetcher.${i}.time.csv"
    echo "sketch_size,sketch_time" > $csv_file
    start_time=$(date +%s)
    echo "Start sketching $i at $(date)"
    # Run Sphetcher
    baselines/sphetcher/src/sphetcher data/processed/${dataset}.pca.csv $i experiments/${dataset}/Sphetcher_indicator.${i}.csv
    end_time=$(date +%s)
    echo "End sketching $i at $(date)"
    duration=$((end_time - start_time))
    echo "${i},${duration}" >> $csv_file
    echo "Duration written to $csv_file"
done
