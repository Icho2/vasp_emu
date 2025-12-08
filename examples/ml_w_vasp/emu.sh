#!/bin/bash
#SBATCH -J sf0_24
#SBATCH --gres=shard:6
#SBATCH --mail-type=all 
#SBATCH --mail-user=ma62425@my.utexas.edu

vasp_emu
