#!/usr/bin/env python3
import sys
import os
import glob
import shutil
import random
import argparse
from ase.io import read, write, Trajectory

parser = argparse.ArgumentParser()

parser.add_argument(
        "-d", "--data",
        type=str,
        default="OUTCAR",
        help="Specifies the path to the data you are finetuning with.",
        )

parser.add_argument(
        "-s", "--split",
        type=int,
        default=80,
        help="Specifies the percentage of your data to be considered for training (the rest used for testing).",
        )

args = parser.parse_args()

# Reading data and making the directory where all our data will be
outcar = read(args.data, index=":")
os.makedirs("fine_tuning_data", exist_ok=True)

# generating individual *.traj files from OUTCAR
for i, image in enumerate(outcar):
        traj_writer = Trajectory(f"fine_tuning_data/structure_{i:04}.traj", mode='w')
        traj_writer.write(image)

traj_writer.close()

#Making train and test trajectories of *.traj files
train_dir = 'fine_tuning_data/train'
test_dir = 'fine_tuning_data/test'
os.makedirs(train_dir, exist_ok=True)
os.makedirs(test_dir, exist_ok=True)
files = glob.glob(os.path.join("fine_tuning_data", "*.traj"))
print("files: ",files) 
random.shuffle(files)

# 80-20 traj split
split_index = int(len(files) * 0.80) 
train_files = files[:split_index]
test_files = files[split_index:]
print("train_files: ", train_files)
# Moving traj files to respective directories
for file_path in train_files:
    print(file_path)
    shutil.move(file_path, os.path.join(train_dir, os.path.basename(file_path)))

for file_path in test_files:
    shutil.move(file_path, os.path.join(test_dir, os.path.basename(file_path)))
