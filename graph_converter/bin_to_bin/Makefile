CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall

# 所有的源文件和目标文件
SRCS = bin_to_bin_bw_int.cpp bin_to_bin_fw_int.cpp bin_to_bin_undirected.cpp bin_parser.cpp
BINS = $(SRCS:.cpp=)

all: $(BINS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(BINS) *.o

.PHONY: all clean
