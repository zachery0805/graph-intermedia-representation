#!/bin/bash

# 编译
echo "Compiling GIR project..."
make clean
make all

# 检查是否编译成功
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

# 运行
echo "Running GIR..."
./gir

# 如果需要，可以添加运行后清理
# make clean
