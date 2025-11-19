#!/bin/bash
echo $FAIRCHEM_PATH
EXECUTABLE_PATH="$FAIRCHEM_PATH/src/fairchem/core/scripts/create_uma_finetune_dataset.py"
TRAIN_DIR="$PWD/fine_tuning_data/train/"
VAL_DIR="$PWD/fine_tuning_data/test/"
OUTPUT_DIR="$PWD/model_output"

cd $FAIRCHEM_PATH
python3 "$EXECUTABLE_PATH" \
  --train-dir $TRAIN_DIR \
  --val-dir $VAL_DIR \
  --output-dir $OUTPUT_DIR \
  --uma-task=odac \
  --regression-task ef
