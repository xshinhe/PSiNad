#!/bin/bash


start=2
end=10
main_dir="cas"

for i in $(seq $start $end)
do
    echo "Running job $i"
    sbatch task_lhc.sh $i $main_dir
done