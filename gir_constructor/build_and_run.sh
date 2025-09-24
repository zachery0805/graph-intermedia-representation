#!/bin/bash

# 默认不开启日志
DEBUG=0
OUTPUT_TIME=0

# 支持命令行参数开启日志
# 使用方法: ./build_and_run.sh --debug --time
for arg in "$@"; do
    case $arg in
        --debug)
            DEBUG=1
            shift
            ;;
        --time)
            OUTPUT_TIME=1
            shift
            ;;
        *)
            ;;
    esac
done

echo "Compiling GIR project..."
make clean
make all DEBUG=$DEBUG OUTPUT_TIME=$OUTPUT_TIME

if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

echo "Running GIR..."
./gir
