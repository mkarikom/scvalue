#!/bin/bash
dataset="tma"
lr=0.0001
ne=50
ms=128
for method in "Uniform" "GeoSketch" "Sphetcher" "Hopper" "KH" "scSampler" "scValue"
do
    for size in 863 1727 2591 3455 4319
    do
        for i in $(seq 0 9)
        do    
            seed=$((42 + i))
            echo "Learning on $method with size=$size seed=$seed lr=$lr ms=$ms ne=$ne"
            /mnt/alamo01/users/liasktao/anaconda3/envs/actinn/bin/python actinn_predict_sketch.py -trs experiments/${dataset}/${method}.${size}.h5ad -ts actinn_data/test.age.h5ad -lr $lr -ne $ne -ms $ms -rs $seed -sm $method -ss $size
            /mnt/alamo01/users/liasktao/anaconda3/envs/actinn/bin/python eval_actinn_pred_sketch.py -rs $seed -sm $method -ss $size
        done
    done
done         