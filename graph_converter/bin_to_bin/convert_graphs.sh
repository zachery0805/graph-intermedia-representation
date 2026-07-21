#!/bin/bash
# convert_grid.sh
# Usage: ./convert_grid.sh
# Converts grid_1000_100000_sym.bin to forward, backward, and undirected binary files

# 输入文件路径
INPUT_BIN="/home/cxl6029/projects/GIR/data/grid_1000_100000/grid_1000_100000_sym.bin"

echo "=== Step 1: Forward edge parser ==="
./bin_to_bin_fw_int $INPUT_BIN

echo "=== Step 2: Backward edge parser ==="
./bin_to_bin_bw_int $INPUT_BIN

echo "=== Step 3: Convert to undirected ==="
# 输入文件顺序：begin.bin adjacent.bin head.bin out_degree.bin
./bin_to_bin_undirected fw_begin.bin fw_adjacent.bin head.bin out_degree.bin

echo "=== All steps completed successfully ==="
