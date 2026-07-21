# 编译器和编译选项
CXX = g++
CXXFLAGS = -O3

# 所有的源文件（每个都含 main）
TOOLS = \
    graph_clean \
    tuple_text_to_bin \
    txt_to_bin_bw \
    txt_to_bin_bw_int \
    txt_to_bin_fw \
    txt_to_bin_fw_int \
    txt_to_bin_undirected

# 对应的可执行文件
BINARIES = $(TOOLS)

# 默认目标
all: $(BINARIES)

# 每个工具独立编译成一个可执行文件
%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# 清理目标
clean:
	rm -f $(BINARIES) *.o
