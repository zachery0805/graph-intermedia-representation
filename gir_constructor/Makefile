CXX = g++
CXXFLAGS = -std=c++17 -O3 -fopenmp -I. -I./includes -I./includes/scc_cpu

# make DEBUG=1 / make OUTPUT_TIME=1
ifeq ($(DEBUG),1)
    CXXFLAGS += -DDEBUG
endif
ifeq ($(OUTPUT_TIME),1)
    CXXFLAGS += -DOUTPUT_TIME
endif

CORE = gir.o gir_query.o graph.o graph_undirected.o

# jemalloc when installed (make JEMALLOC= to build without): the batched store
# updates allocate many small per-vertex lists from all threads at once, and
# glibc malloc's locking roughly halves their throughput.
JEMALLOC ?= $(firstword $(wildcard /lib/x86_64-linux-gnu/libjemalloc.so.2 /usr/lib/x86_64-linux-gnu/libjemalloc.so.2))
LDLIBS = $(JEMALLOC)

all: gir bench_gir gir_selftest el_to_csr

%.o: %.cpp gir.h gir_iter.h graph.h graph_undirected.h block_cut_tree.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# gir_test.cpp is #included by the drivers, not compiled on its own.
gir: main.o $(CORE)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

bench_gir: bench_main.o $(CORE)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

gir_selftest: gir_selftest.o $(CORE)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

main.o bench_main.o gir_selftest.o: gir_test.cpp graph_tests.h dynamic_repair.h

# Differential correctness tests on generated graphs + data/example_graph.
check: gir_selftest
	./gir_selftest /tmp/gir_selftest_graphs ../data/example_graph

clean:
	rm -f *.o gir bench_gir gir_selftest

.PHONY: all check clean

# SNAP edge list -> binary CSR (data/<name>/fw_*.bin, bw_*.bin)
el_to_csr: el_to_csr.cpp
	$(CXX) -std=c++17 -O3 -fopenmp $< -o $@
